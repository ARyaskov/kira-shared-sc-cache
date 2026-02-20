use std::fs::File;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use thiserror::Error;

const HEADER_SIZE: usize = 256;
const ALIGNMENT: usize = 64;
const MAGIC: &[u8; 4] = b"KORG";
const VERSION_MAJOR: u16 = 1;
const VERSION_MINOR: u16 = 0;
const ENDIAN_TAG: u32 = 0x1234_5678;
const CRC64_ECMA_POLY: u64 = 0x42F0_E1EB_A9EA_3693;

#[derive(Debug, Error)]
pub enum SharedCacheError {
    #[error("I/O error in {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("cache format error in {path:?}: {message}")]
    Format { path: PathBuf, message: String },
}

#[derive(Debug)]
pub struct SharedCacheMmap {
    mmap: Mmap,
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    col_ptr_offset: usize,
    row_idx_offset: usize,
    values_u32_offset: usize,
}

impl SharedCacheMmap {
    pub fn col_ptr(&self) -> &[u64] {
        let len = self.n_cells + 1;
        let bytes = &self.mmap[self.col_ptr_offset..self.col_ptr_offset + len * 8];
        // SAFETY: validated section bounds in mmap_shared_cache.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u64, len) }
    }

    pub fn row_idx(&self) -> &[u32] {
        let bytes = &self.mmap[self.row_idx_offset..self.row_idx_offset + self.nnz * 4];
        // SAFETY: validated section bounds in mmap_shared_cache.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u32, self.nnz) }
    }

    pub fn values_u32(&self) -> &[u32] {
        let bytes = &self.mmap[self.values_u32_offset..self.values_u32_offset + self.nnz * 4];
        // SAFETY: validated section bounds in mmap_shared_cache.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u32, self.nnz) }
    }
}

#[derive(Debug, Clone)]
pub struct SharedCacheOwned {
    pub path: PathBuf,
    pub n_genes: u64,
    pub n_cells: u64,
    pub nnz: u64,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    pub col_ptr: Vec<u64>,
    pub row_idx: Vec<u32>,
    pub values_u32: Vec<u32>,
}

pub struct SharedCacheWriteInput<'a> {
    pub genes: &'a [String],
    pub barcodes: &'a [String],
    pub col_ptr: &'a [u64],
    pub row_idx: &'a [u32],
    pub values_u32: &'a [u32],
}

pub fn validate_dimensions(path: &Path) -> Result<(u64, u64), SharedCacheError> {
    let mapped = mmap_shared_cache(path)?;
    Ok((mapped.n_genes as u64, mapped.n_cells as u64))
}

pub fn read_shared_cache_owned(path: &Path) -> Result<SharedCacheOwned, SharedCacheError> {
    let mapped = mmap_shared_cache(path)?;
    let col_ptr = mapped.col_ptr().to_vec();
    let row_idx = mapped.row_idx().to_vec();
    let values_u32 = mapped.values_u32().to_vec();
    Ok(SharedCacheOwned {
        path: path.to_path_buf(),
        n_genes: mapped.n_genes as u64,
        n_cells: mapped.n_cells as u64,
        nnz: mapped.nnz as u64,
        genes: mapped.genes,
        barcodes: mapped.barcodes,
        col_ptr,
        row_idx,
        values_u32,
    })
}

pub fn mmap_shared_cache(path: &Path) -> Result<SharedCacheMmap, SharedCacheError> {
    let file = File::open(path).map_err(|source| SharedCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let metadata = file.metadata().map_err(|source| SharedCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() < HEADER_SIZE as u64 {
        return Err(format_error(path, "file too small"));
    }

    let mmap = unsafe {
        Mmap::map(&file).map_err(|source| SharedCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?
    };

    let header = &mmap[..HEADER_SIZE];
    validate_header(path, header, metadata.len() as usize)?;

    let n_genes = read_u64(header, 16)? as usize;
    let n_cells = read_u64(header, 24)? as usize;
    let nnz = read_u64(header, 32)? as usize;

    let genes_table_offset = read_u64(header, 40)? as usize;
    let genes_table_bytes = read_u64(header, 48)? as usize;
    let barcodes_table_offset = read_u64(header, 56)? as usize;
    let barcodes_table_bytes = read_u64(header, 64)? as usize;
    let col_ptr_offset = read_u64(header, 72)? as usize;
    let row_idx_offset = read_u64(header, 80)? as usize;
    let values_u32_offset = read_u64(header, 88)? as usize;

    let genes = parse_string_table(
        path,
        &mmap,
        genes_table_offset,
        genes_table_bytes,
        n_genes,
        "genes",
    )?;
    let barcodes = parse_string_table(
        path,
        &mmap,
        barcodes_table_offset,
        barcodes_table_bytes,
        n_cells,
        "barcodes",
    )?;

    let col_ptr_bytes = (n_cells + 1)
        .checked_mul(8)
        .ok_or_else(|| format_error(path, "col_ptr length overflow"))?;
    let row_idx_bytes = nnz
        .checked_mul(4)
        .ok_or_else(|| format_error(path, "row_idx length overflow"))?;
    let values_bytes = nnz
        .checked_mul(4)
        .ok_or_else(|| format_error(path, "values_u32 length overflow"))?;

    check_range(path, &mmap, col_ptr_offset, col_ptr_bytes, "col_ptr")?;
    check_range(path, &mmap, row_idx_offset, row_idx_bytes, "row_idx")?;
    check_range(path, &mmap, values_u32_offset, values_bytes, "values_u32")?;

    let col_ptr = {
        let bytes = &mmap[col_ptr_offset..col_ptr_offset + col_ptr_bytes];
        // SAFETY: checked range and alignment invariants in validate_header/check_range.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u64, n_cells + 1) }
    };
    let row_idx = {
        let bytes = &mmap[row_idx_offset..row_idx_offset + row_idx_bytes];
        // SAFETY: checked range and alignment invariants in validate_header/check_range.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u32, nnz) }
    };
    validate_csc(path, col_ptr, row_idx, n_genes, nnz)?;

    Ok(SharedCacheMmap {
        mmap,
        n_genes,
        n_cells,
        nnz,
        genes,
        barcodes,
        col_ptr_offset,
        row_idx_offset,
        values_u32_offset,
    })
}

pub fn write_shared_cache(
    path: &Path,
    input: &SharedCacheWriteInput<'_>,
) -> Result<(), SharedCacheError> {
    let n_genes = input.genes.len();
    let n_cells = input.barcodes.len();
    if input.col_ptr.len() != n_cells + 1 {
        return Err(format_error(path, "col_ptr length must be n_cells + 1"));
    }
    if input.row_idx.len() != input.values_u32.len() {
        return Err(format_error(path, "row_idx and values_u32 length mismatch"));
    }

    let nnz = input.values_u32.len();
    validate_csc(path, input.col_ptr, input.row_idx, n_genes, nnz)?;

    let genes_table = encode_string_table(path, input.genes)?;
    let barcodes_table = encode_string_table(path, input.barcodes)?;

    let genes_table_offset = align_to(HEADER_SIZE, ALIGNMENT);
    let barcodes_table_offset = align_to(genes_table_offset + genes_table.len(), ALIGNMENT);
    let col_ptr_offset = align_to(barcodes_table_offset + barcodes_table.len(), ALIGNMENT);
    let col_ptr_bytes = input.col_ptr.len() * 8;
    let row_idx_offset = align_to(col_ptr_offset + col_ptr_bytes, ALIGNMENT);
    let row_idx_bytes = input.row_idx.len() * 4;
    let values_u32_offset = align_to(row_idx_offset + row_idx_bytes, ALIGNMENT);
    let values_u32_bytes = input.values_u32.len() * 4;
    let file_bytes = values_u32_offset + values_u32_bytes;

    let mut header = [0u8; HEADER_SIZE];
    header[0..4].copy_from_slice(MAGIC);
    write_u16(&mut header, 4, VERSION_MAJOR);
    write_u16(&mut header, 6, VERSION_MINOR);
    write_u32(&mut header, 8, ENDIAN_TAG);
    write_u32(&mut header, 12, HEADER_SIZE as u32);
    write_u64(&mut header, 16, n_genes as u64);
    write_u64(&mut header, 24, n_cells as u64);
    write_u64(&mut header, 32, nnz as u64);
    write_u64(&mut header, 40, genes_table_offset as u64);
    write_u64(&mut header, 48, genes_table.len() as u64);
    write_u64(&mut header, 56, barcodes_table_offset as u64);
    write_u64(&mut header, 64, barcodes_table.len() as u64);
    write_u64(&mut header, 72, col_ptr_offset as u64);
    write_u64(&mut header, 80, row_idx_offset as u64);
    write_u64(&mut header, 88, values_u32_offset as u64);
    write_u64(&mut header, 96, 0);
    write_u64(&mut header, 104, 0);
    write_u64(&mut header, 112, file_bytes as u64);
    write_u64(&mut header, 128, 0);
    let crc = crc64_ecma(&header);
    write_u64(&mut header, 120, crc);

    let mut file = File::create(path).map_err(|source| SharedCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    file.set_len(file_bytes as u64)
        .map_err(|source| SharedCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    write_at(&mut file, 0, &header, path)?;
    write_at(&mut file, genes_table_offset, &genes_table, path)?;
    write_at(&mut file, barcodes_table_offset, &barcodes_table, path)?;

    write_u64_slice(&mut file, col_ptr_offset, input.col_ptr, path)?;
    write_u32_slice(&mut file, row_idx_offset, input.row_idx, path)?;
    write_u32_slice(&mut file, values_u32_offset, input.values_u32, path)?;

    file.sync_all().map_err(|source| SharedCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    Ok(())
}

fn validate_header(path: &Path, header: &[u8], file_len: usize) -> Result<(), SharedCacheError> {
    if &header[0..4] != MAGIC {
        return Err(format_error(path, "invalid magic"));
    }

    if read_u16(header, 4)? != VERSION_MAJOR {
        return Err(format_error(path, "unsupported major version"));
    }
    if read_u16(header, 6)? != VERSION_MINOR {
        return Err(format_error(path, "unsupported minor version"));
    }
    if read_u32(header, 8)? != ENDIAN_TAG {
        return Err(format_error(path, "invalid endian tag"));
    }
    if read_u32(header, 12)? as usize != HEADER_SIZE {
        return Err(format_error(path, "unsupported header size"));
    }

    let file_bytes = read_u64(header, 112)? as usize;
    if file_bytes != file_len {
        return Err(format_error(path, "file_bytes mismatch"));
    }

    let n_blocks = read_u64(header, 96)?;
    let blocks_offset = read_u64(header, 104)?;
    let data_crc64 = read_u64(header, 128)?;
    if n_blocks != 0 || blocks_offset != 0 || data_crc64 != 0 {
        return Err(format_error(
            path,
            "v1 requires n_blocks=0, blocks_offset=0, data_crc64=0",
        ));
    }

    let mut crc_header = header.to_vec();
    crc_header[120..128].fill(0);
    let expected_crc = read_u64(header, 120)?;
    let actual_crc = crc64_ecma(&crc_header);
    if actual_crc != expected_crc {
        return Err(format_error(path, "header_crc64 mismatch"));
    }

    for offset in [
        read_u64(header, 40)? as usize,
        read_u64(header, 56)? as usize,
        read_u64(header, 72)? as usize,
        read_u64(header, 80)? as usize,
        read_u64(header, 88)? as usize,
    ] {
        if offset % ALIGNMENT != 0 {
            return Err(format_error(path, "section offset is not 64-byte aligned"));
        }
        if offset >= file_len {
            return Err(format_error(path, "section offset out of bounds"));
        }
    }

    Ok(())
}

fn validate_csc(
    path: &Path,
    col_ptr: &[u64],
    row_idx: &[u32],
    n_genes: usize,
    nnz: usize,
) -> Result<(), SharedCacheError> {
    if col_ptr.is_empty() {
        return Err(format_error(path, "col_ptr must not be empty"));
    }
    if col_ptr[0] != 0 {
        return Err(format_error(path, "col_ptr[0] must be 0"));
    }
    if *col_ptr.last().unwrap_or(&0) as usize != nnz {
        return Err(format_error(path, "col_ptr[n_cells] must equal nnz"));
    }

    for obs in 0..(col_ptr.len() - 1) {
        if col_ptr[obs] > col_ptr[obs + 1] {
            return Err(format_error(path, "col_ptr must be monotonic"));
        }
        let start = col_ptr[obs] as usize;
        let end = col_ptr[obs + 1] as usize;
        if end > nnz {
            return Err(format_error(path, "col_ptr entry exceeds nnz"));
        }

        let mut prev: Option<u32> = None;
        for &row in &row_idx[start..end] {
            if row as usize >= n_genes {
                return Err(format_error(path, "row_idx out of bounds"));
            }
            if let Some(last) = prev
                && row <= last
            {
                return Err(format_error(
                    path,
                    "row_idx must be strictly increasing inside each column",
                ));
            }
            prev = Some(row);
        }
    }

    Ok(())
}

fn parse_string_table(
    path: &Path,
    bytes: &[u8],
    offset: usize,
    table_bytes: usize,
    expected_count: usize,
    field_name: &str,
) -> Result<Vec<String>, SharedCacheError> {
    check_range(path, bytes, offset, table_bytes, field_name)?;

    let table = &bytes[offset..offset + table_bytes];
    if table.len() < 8 {
        return Err(format_error(path, &format!("{field_name} table too small")));
    }

    let count = read_u32(table, 0)? as usize;
    if count != expected_count {
        return Err(format_error(
            path,
            &format!("{field_name} count mismatch: expected {expected_count}, got {count}"),
        ));
    }

    let offsets_len = (count + 1)
        .checked_mul(4)
        .ok_or_else(|| format_error(path, &format!("{field_name} offsets overflow")))?;
    let blob_offset = 4 + offsets_len;
    if table.len() < blob_offset {
        return Err(format_error(
            path,
            &format!("{field_name} offsets truncated"),
        ));
    }

    let blob_len = read_u32(table, 4 + count * 4)? as usize;
    if blob_offset + blob_len != table.len() {
        return Err(format_error(
            path,
            &format!("{field_name} blob length mismatch"),
        ));
    }

    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let start = read_u32(table, 4 + i * 4)? as usize;
        let end = read_u32(table, 4 + (i + 1) * 4)? as usize;
        if start > end || end > blob_len {
            return Err(format_error(
                path,
                &format!("{field_name} offsets out of bounds"),
            ));
        }
        let s = std::str::from_utf8(&table[blob_offset + start..blob_offset + end])
            .map_err(|_| format_error(path, &format!("{field_name} utf8 decode error")))?;
        out.push(s.to_string());
    }

    Ok(out)
}

fn encode_string_table(path: &Path, strings: &[String]) -> Result<Vec<u8>, SharedCacheError> {
    let mut offsets = Vec::with_capacity(strings.len() + 1);
    offsets.push(0u32);

    let mut blob = Vec::new();
    for s in strings {
        blob.extend_from_slice(s.as_bytes());
        let next = u32::try_from(blob.len())
            .map_err(|_| format_error(path, "string table exceeds u32 range"))?;
        offsets.push(next);
    }

    let mut out = Vec::with_capacity(4 + offsets.len() * 4 + blob.len());
    out.extend_from_slice(&(strings.len() as u32).to_le_bytes());
    for off in offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    out.extend_from_slice(&blob);

    Ok(out)
}

fn check_range(
    path: &Path,
    bytes: &[u8],
    offset: usize,
    len: usize,
    field_name: &str,
) -> Result<(), SharedCacheError> {
    let end = offset
        .checked_add(len)
        .ok_or_else(|| format_error(path, &format!("{field_name} range overflow")))?;
    if end > bytes.len() {
        return Err(format_error(
            path,
            &format!("{field_name} section out of bounds"),
        ));
    }
    Ok(())
}

fn write_at(
    file: &mut File,
    offset: usize,
    data: &[u8],
    path: &Path,
) -> Result<(), SharedCacheError> {
    file.seek(SeekFrom::Start(offset as u64))
        .map_err(|source| SharedCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(data).map_err(|source| SharedCacheError::Io {
        path: path.to_path_buf(),
        source,
    })
}

fn write_u64_slice(
    file: &mut File,
    offset: usize,
    data: &[u64],
    path: &Path,
) -> Result<(), SharedCacheError> {
    let mut buf = Vec::with_capacity(data.len() * 8);
    for &v in data {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    write_at(file, offset, &buf, path)
}

fn write_u32_slice(
    file: &mut File,
    offset: usize,
    data: &[u32],
    path: &Path,
) -> Result<(), SharedCacheError> {
    let mut buf = Vec::with_capacity(data.len() * 4);
    for &v in data {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    write_at(file, offset, &buf, path)
}

fn align_to(value: usize, align: usize) -> usize {
    ((value + align - 1) / align) * align
}

fn write_u16(header: &mut [u8], offset: usize, value: u16) {
    header[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn write_u32(header: &mut [u8], offset: usize, value: u32) {
    header[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn write_u64(header: &mut [u8], offset: usize, value: u64) {
    header[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn read_u16(buf: &[u8], offset: usize) -> Result<u16, SharedCacheError> {
    let mut arr = [0u8; 2];
    arr.copy_from_slice(&buf[offset..offset + 2]);
    Ok(u16::from_le_bytes(arr))
}

fn read_u32(buf: &[u8], offset: usize) -> Result<u32, SharedCacheError> {
    let mut arr = [0u8; 4];
    arr.copy_from_slice(&buf[offset..offset + 4]);
    Ok(u32::from_le_bytes(arr))
}

fn read_u64(buf: &[u8], offset: usize) -> Result<u64, SharedCacheError> {
    let mut arr = [0u8; 8];
    arr.copy_from_slice(&buf[offset..offset + 8]);
    Ok(u64::from_le_bytes(arr))
}

fn format_error(path: &Path, message: &str) -> SharedCacheError {
    SharedCacheError::Format {
        path: path.to_path_buf(),
        message: message.to_string(),
    }
}

pub fn crc64_ecma(bytes: &[u8]) -> u64 {
    let mut crc = 0u64;
    for &byte in bytes {
        crc ^= (byte as u64) << 56;
        for _ in 0..8 {
            if (crc & 0x8000_0000_0000_0000) != 0 {
                crc = (crc << 1) ^ CRC64_ECMA_POLY;
            } else {
                crc <<= 1;
            }
        }
    }
    crc
}

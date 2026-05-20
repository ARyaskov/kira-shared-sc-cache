use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use thiserror::Error;

const MAGIC: &[u8; 8] = b"KIRAMTX\0";
const VERSION: u32 = 1;
const HEADER_LEN: usize = 24;
const MODE_MASK: u32 = 0xFF;

/// Buffer size for the writer. 1 MiB amortises syscalls when writing the
/// dense `values` payload (which dominates the file).
const WRITE_BUF: usize = 1 << 20;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ExprCacheMode {
    Unknown,
    Sample,
    Cluster,
    Cell,
}

impl ExprCacheMode {
    fn to_flags(self) -> u32 {
        match self {
            Self::Unknown => 0,
            Self::Sample => 1,
            Self::Cluster => 2,
            Self::Cell => 3,
        }
    }

    fn from_flags(flags: u32) -> Self {
        match flags & MODE_MASK {
            1 => Self::Sample,
            2 => Self::Cluster,
            3 => Self::Cell,
            _ => Self::Unknown,
        }
    }
}

#[derive(Debug, Error)]
pub enum ExprBinError {
    #[error("I/O error reading {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("invalid expr magic in {path:?}")]
    InvalidMagic { path: PathBuf },
    #[error("unsupported expr version {version} in {path:?}")]
    UnsupportedVersion { path: PathBuf, version: u32 },
    #[error("truncated expr cache {path:?}")]
    Truncated { path: PathBuf },
    #[error("size mismatch in {path:?}: expected {expected} bytes, got {actual}")]
    SizeMismatch {
        path: PathBuf,
        expected: usize,
        actual: usize,
    },
}

#[derive(Debug)]
pub struct ExprBinMmap {
    mmap: Mmap,
    pub genes: usize,
    pub samples: usize,
    pub mode: ExprCacheMode,
}

impl ExprBinMmap {
    /// View the dense expression buffer as `&[f32]` without copying.
    ///
    /// SAFETY contract: `HEADER_LEN = 24` is a multiple of `align_of::<f32>()`
    /// (4), and the mmap base is page-aligned by the kernel, so the value
    /// section is correctly aligned for `f32`. `mmap_expr_bin` validates that
    /// the byte length exactly matches `genes * samples * 4`.
    pub fn values(&self) -> &[f32] {
        let values_len = self.genes * self.samples;
        let values_bytes = &self.mmap[HEADER_LEN..HEADER_LEN + values_len * 4];
        debug_assert_eq!(
            (values_bytes.as_ptr() as usize) % std::mem::align_of::<f32>(),
            0,
            "values section must be f32-aligned"
        );
        // SAFETY: alignment is asserted, length matches the file-format
        // invariant validated in `mmap_expr_bin`.
        unsafe { std::slice::from_raw_parts(values_bytes.as_ptr() as *const f32, values_len) }
    }

    #[inline]
    pub fn get(&self, gene: usize, sample: usize) -> f32 {
        self.values()[gene * self.samples + sample]
    }
}

pub fn write_expr_bin(
    path: &Path,
    genes: usize,
    samples: usize,
    values: &[f32],
) -> Result<(), ExprBinError> {
    write_expr_bin_with_mode(path, genes, samples, values, ExprCacheMode::Unknown)
}

pub fn write_expr_bin_with_mode(
    path: &Path,
    genes: usize,
    samples: usize,
    values: &[f32],
    mode: ExprCacheMode,
) -> Result<(), ExprBinError> {
    let expected_values_len =
        genes
            .checked_mul(samples)
            .ok_or_else(|| ExprBinError::SizeMismatch {
                path: path.to_path_buf(),
                expected: usize::MAX,
                actual: values.len(),
            })?;
    if values.len() != expected_values_len {
        return Err(ExprBinError::SizeMismatch {
            path: path.to_path_buf(),
            expected: expected_values_len,
            actual: values.len(),
        });
    }

    let file = File::create(path).map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::with_capacity(WRITE_BUF, file);

    // Header: 24 bytes built on the stack so the BufWriter sees one chunk.
    let mut header = [0u8; HEADER_LEN];
    header[..8].copy_from_slice(MAGIC);
    header[8..12].copy_from_slice(&VERSION.to_le_bytes());
    header[12..16].copy_from_slice(&(genes as u32).to_le_bytes());
    header[16..20].copy_from_slice(&(samples as u32).to_le_bytes());
    header[20..24].copy_from_slice(&mode.to_flags().to_le_bytes());
    write_all(&mut writer, &header, path)?;

    let byte_len = values.len() * std::mem::size_of::<f32>();
    // SAFETY: f32 is plain old data; we only reinterpret contiguous bytes.
    let value_bytes = unsafe { std::slice::from_raw_parts(values.as_ptr() as *const u8, byte_len) };
    write_all(&mut writer, value_bytes, path)?;

    let file = writer.into_inner().map_err(|err| ExprBinError::Io {
        path: path.to_path_buf(),
        source: err.into_error(),
    })?;
    file.sync_all().map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    Ok(())
}

#[inline]
fn write_all<W: Write>(w: &mut W, data: &[u8], path: &Path) -> Result<(), ExprBinError> {
    w.write_all(data).map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })
}

pub fn mmap_expr_bin(path: &Path) -> Result<ExprBinMmap, ExprBinError> {
    let file = File::open(path).map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let metadata = file.metadata().map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() < HEADER_LEN as u64 {
        return Err(ExprBinError::Truncated {
            path: path.to_path_buf(),
        });
    }

    let mmap = unsafe {
        Mmap::map(&file).map_err(|source| ExprBinError::Io {
            path: path.to_path_buf(),
            source,
        })?
    };
    // Hint the kernel to read ahead; we always scan the dense matrix
    // sequentially in downstream tools. memmap2 only exposes `advise` on
    // Unix — on other platforms the kernel uses its default heuristic.
    #[cfg(unix)]
    {
        let _ = mmap.advise(memmap2::Advice::Sequential);
    }

    let header = &mmap[..HEADER_LEN];
    if &header[..8] != MAGIC {
        return Err(ExprBinError::InvalidMagic {
            path: path.to_path_buf(),
        });
    }

    let version = read_u32(header, 8);
    if version != VERSION {
        return Err(ExprBinError::UnsupportedVersion {
            path: path.to_path_buf(),
            version,
        });
    }

    let genes = read_u32(header, 12) as usize;
    let samples = read_u32(header, 16) as usize;
    let mode = ExprCacheMode::from_flags(read_u32(header, 20));

    let values_len = genes
        .checked_mul(samples)
        .ok_or_else(|| ExprBinError::SizeMismatch {
            path: path.to_path_buf(),
            expected: usize::MAX,
            actual: metadata.len() as usize,
        })?;
    let expected_bytes = HEADER_LEN + values_len * 4;
    let actual_bytes = metadata.len() as usize;
    if expected_bytes != actual_bytes {
        return Err(ExprBinError::SizeMismatch {
            path: path.to_path_buf(),
            expected: expected_bytes,
            actual: actual_bytes,
        });
    }

    Ok(ExprBinMmap {
        mmap,
        genes,
        samples,
        mode,
    })
}

#[inline]
fn read_u32(buf: &[u8], offset: usize) -> u32 {
    let mut arr = [0u8; 4];
    arr.copy_from_slice(&buf[offset..offset + 4]);
    u32::from_le_bytes(arr)
}

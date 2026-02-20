use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use thiserror::Error;

const MAGIC: &[u8; 8] = b"KIRAMTX\0";
const VERSION: u32 = 1;
const HEADER_LEN: usize = 24;
const MODE_MASK: u32 = 0xFF;

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
    pub fn values(&self) -> &[f32] {
        let values_len = self.genes * self.samples;
        let values_bytes = &self.mmap[HEADER_LEN..HEADER_LEN + values_len * 4];
        // SAFETY: validated layout in mmap_expr_bin.
        unsafe { std::slice::from_raw_parts(values_bytes.as_ptr() as *const f32, values_len) }
    }

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

    let mut file = File::create(path).map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    file.write_all(MAGIC).map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    file.write_all(&VERSION.to_le_bytes())
        .map_err(|source| ExprBinError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(&(genes as u32).to_le_bytes())
        .map_err(|source| ExprBinError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(&(samples as u32).to_le_bytes())
        .map_err(|source| ExprBinError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(&mode.to_flags().to_le_bytes())
        .map_err(|source| ExprBinError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    // SAFETY: f32 is plain old data; we only reinterpret contiguous bytes.
    let byte_len = values.len() * std::mem::size_of::<f32>();
    let value_bytes = unsafe { std::slice::from_raw_parts(values.as_ptr() as *const u8, byte_len) };
    file.write_all(value_bytes)
        .map_err(|source| ExprBinError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    file.sync_all().map_err(|source| ExprBinError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    Ok(())
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

    let header = &mmap[..HEADER_LEN];
    if &header[..8] != MAGIC {
        return Err(ExprBinError::InvalidMagic {
            path: path.to_path_buf(),
        });
    }

    let version = read_u32(header, 8)?;
    if version != VERSION {
        return Err(ExprBinError::UnsupportedVersion {
            path: path.to_path_buf(),
            version,
        });
    }

    let genes = read_u32(header, 12)? as usize;
    let samples = read_u32(header, 16)? as usize;
    let mode = ExprCacheMode::from_flags(read_u32(header, 20)?);

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
fn read_u32(buf: &[u8], offset: usize) -> Result<u32, ExprBinError> {
    let mut arr = [0u8; 4];
    arr.copy_from_slice(&buf[offset..offset + 4]);
    Ok(u32::from_le_bytes(arr))
}

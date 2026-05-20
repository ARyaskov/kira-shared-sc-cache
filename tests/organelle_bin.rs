//! Round-trip + corruption + validation tests for the shared organelle cache.

use std::fs::{File, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_shared_sc_cache::{
    SHARED_CACHE_BASENAME, SharedCacheError, SharedCacheWriteInput, crc64_ecma, mmap_shared_cache,
    read_shared_cache_owned, resolve_shared_cache_filename, validate_dimensions,
    write_shared_cache,
};

fn temp_file(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("kira_shared_org_{label}_{ts}.bin"))
}

fn tiny_input() -> (Vec<String>, Vec<String>, Vec<u64>, Vec<u32>, Vec<u32>) {
    let genes = vec!["G1".to_string(), "G2".to_string(), "G3".to_string()];
    let barcodes = vec!["C1".to_string(), "C2".to_string()];
    // CSC for a 3×2 matrix:
    //   [[1, 0],
    //    [2, 4],
    //    [0, 5]]
    let col_ptr: Vec<u64> = vec![0, 2, 4];
    let row_idx: Vec<u32> = vec![0, 1, 1, 2];
    let values: Vec<u32> = vec![1, 2, 4, 5];
    (genes, barcodes, col_ptr, row_idx, values)
}

#[test]
fn roundtrip_mmap() {
    let path = temp_file("rt_mmap");
    let (genes, barcodes, col_ptr, row_idx, values) = tiny_input();
    write_shared_cache(
        &path,
        &SharedCacheWriteInput {
            genes: &genes,
            barcodes: &barcodes,
            col_ptr: &col_ptr,
            row_idx: &row_idx,
            values_u32: &values,
        },
    )
    .expect("write");

    let view = mmap_shared_cache(&path).expect("mmap");
    assert_eq!(view.n_genes, 3);
    assert_eq!(view.n_cells, 2);
    assert_eq!(view.nnz, 4);
    assert_eq!(view.genes, genes);
    assert_eq!(view.barcodes, barcodes);
    assert_eq!(view.col_ptr(), col_ptr.as_slice());
    assert_eq!(view.row_idx(), row_idx.as_slice());
    assert_eq!(view.values_u32(), values.as_slice());
}

#[test]
fn roundtrip_owned() {
    let path = temp_file("rt_owned");
    let (genes, barcodes, col_ptr, row_idx, values) = tiny_input();
    write_shared_cache(
        &path,
        &SharedCacheWriteInput {
            genes: &genes,
            barcodes: &barcodes,
            col_ptr: &col_ptr,
            row_idx: &row_idx,
            values_u32: &values,
        },
    )
    .expect("write");

    let owned = read_shared_cache_owned(&path).expect("read_owned");
    assert_eq!(owned.n_genes, 3);
    assert_eq!(owned.n_cells, 2);
    assert_eq!(owned.nnz, 4);
    assert_eq!(owned.col_ptr, col_ptr);
    assert_eq!(owned.row_idx, row_idx);
    assert_eq!(owned.values_u32, values);
    assert_eq!(owned.genes, genes);
    assert_eq!(owned.barcodes, barcodes);
}

#[test]
fn validate_dimensions_returns_expected() {
    let path = temp_file("dims");
    let (genes, barcodes, col_ptr, row_idx, values) = tiny_input();
    write_shared_cache(
        &path,
        &SharedCacheWriteInput {
            genes: &genes,
            barcodes: &barcodes,
            col_ptr: &col_ptr,
            row_idx: &row_idx,
            values_u32: &values,
        },
    )
    .expect("write");

    let (n_genes, n_cells) = validate_dimensions(&path).expect("dims");
    assert_eq!((n_genes, n_cells), (3, 2));
}

#[test]
fn rejects_col_ptr_length_mismatch_on_write() {
    let path = temp_file("bad_colptr");
    let res = write_shared_cache(
        &path,
        &SharedCacheWriteInput {
            genes: &[],
            barcodes: &["C1".to_string(), "C2".to_string()],
            // Should be n_cells + 1 = 3, providing 2 → error.
            col_ptr: &[0, 0],
            row_idx: &[],
            values_u32: &[],
        },
    );
    assert!(matches!(res, Err(SharedCacheError::Format { .. })));
}

#[test]
fn rejects_unsorted_rows_in_column() {
    let path = temp_file("bad_rows");
    let res = write_shared_cache(
        &path,
        &SharedCacheWriteInput {
            genes: &vec!["g".to_string(); 3],
            barcodes: &["c".to_string()],
            col_ptr: &[0, 3],
            // Rows must be strictly increasing inside a column.
            row_idx: &[1, 0, 2],
            values_u32: &[1, 2, 3],
        },
    );
    assert!(matches!(res, Err(SharedCacheError::Format { .. })));
}

#[test]
fn rejects_truncated_file() {
    let path = temp_file("short");
    File::create(&path).unwrap().write_all(&[0u8; 32]).unwrap();
    assert!(mmap_shared_cache(&path).is_err());
}

#[test]
fn rejects_bad_magic() {
    let path = temp_file("magic");
    // Build a header-shaped buffer with the wrong magic.
    let mut buf = vec![0u8; 256];
    buf[0..4].copy_from_slice(b"BADX");
    File::create(&path).unwrap().write_all(&buf).unwrap();
    assert!(mmap_shared_cache(&path).is_err());
}

#[test]
fn rejects_corrupted_header_crc() {
    let path = temp_file("crc");
    let (genes, barcodes, col_ptr, row_idx, values) = tiny_input();
    write_shared_cache(
        &path,
        &SharedCacheWriteInput {
            genes: &genes,
            barcodes: &barcodes,
            col_ptr: &col_ptr,
            row_idx: &row_idx,
            values_u32: &values,
        },
    )
    .expect("write");

    // Flip a single header byte (within the n_genes field) — CRC should fail.
    let mut f = OpenOptions::new().write(true).open(&path).unwrap();
    f.seek(SeekFrom::Start(16)).unwrap();
    f.write_all(&[0xFF]).unwrap();
    drop(f);

    let err = mmap_shared_cache(&path).unwrap_err();
    match err {
        SharedCacheError::Format { message, .. } => {
            assert!(message.contains("crc"), "unexpected: {message}");
        }
        other => panic!("expected Format error, got {other:?}"),
    }
}

#[test]
fn shared_cache_filename_resolution() {
    assert_eq!(resolve_shared_cache_filename(None), SHARED_CACHE_BASENAME);
    assert_eq!(
        resolve_shared_cache_filename(Some("")),
        SHARED_CACHE_BASENAME
    );
    assert_eq!(
        resolve_shared_cache_filename(Some("sample")),
        format!("sample.{SHARED_CACHE_BASENAME}")
    );
}

#[test]
fn crc64_ecma_known_value() {
    // CRC64-ECMA of "123456789" is 0x6c40df5f0b497347
    assert_eq!(crc64_ecma(b"123456789"), 0x6c40df5f0b497347);
}

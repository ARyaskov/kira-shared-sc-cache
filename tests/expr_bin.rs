//! Round-trip and corruption tests for the expression cache.

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_shared_sc_cache::{
    ExprBinError, ExprCacheMode, mmap_expr_bin, write_expr_bin, write_expr_bin_with_mode,
};

fn temp_file(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("kira_shared_expr_{label}_{ts}.bin"))
}

#[test]
fn roundtrip_legacy_writer_keeps_unknown_mode() {
    let path = temp_file("legacy");
    let values: Vec<f32> = (0..6).map(|i| i as f32 * 0.5).collect();
    write_expr_bin(&path, 2, 3, &values).expect("write");

    let view = mmap_expr_bin(&path).expect("mmap");
    assert_eq!(view.genes, 2);
    assert_eq!(view.samples, 3);
    assert_eq!(view.mode, ExprCacheMode::Unknown);
    assert_eq!(view.values(), values.as_slice());
    assert_eq!(view.get(0, 0), 0.0);
    assert_eq!(view.get(1, 2), 2.5);
}

#[test]
fn roundtrip_each_mode() {
    for mode in [
        ExprCacheMode::Sample,
        ExprCacheMode::Cluster,
        ExprCacheMode::Cell,
    ] {
        let path = temp_file(&format!("mode_{mode:?}"));
        let values = vec![1.0, 2.0, 3.0, 4.0];
        write_expr_bin_with_mode(&path, 2, 2, &values, mode).expect("write");
        let view = mmap_expr_bin(&path).expect("mmap");
        assert_eq!(view.mode, mode);
        assert_eq!(view.values(), values.as_slice());
    }
}

#[test]
fn size_mismatch_at_write() {
    let path = temp_file("size_w");
    // Declared dims 3x3 but values has 4 entries.
    let err = write_expr_bin(&path, 3, 3, &[1.0, 2.0, 3.0, 4.0]).unwrap_err();
    matches!(err, ExprBinError::SizeMismatch { .. })
        .then_some(())
        .expect("expected SizeMismatch");
}

#[test]
fn truncated_file_is_detected() {
    let path = temp_file("trunc");
    let mut file = File::create(&path).unwrap();
    file.write_all(&[0u8; 8]).unwrap();
    drop(file);

    let err = mmap_expr_bin(&path).unwrap_err();
    matches!(err, ExprBinError::Truncated { .. })
        .then_some(())
        .expect("expected Truncated");
}

#[test]
fn invalid_magic_is_rejected() {
    let path = temp_file("magic");
    let mut file = File::create(&path).unwrap();
    file.write_all(b"BADMAGIC").unwrap();
    file.write_all(&1u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    drop(file);

    let err = mmap_expr_bin(&path).unwrap_err();
    matches!(err, ExprBinError::InvalidMagic { .. })
        .then_some(())
        .expect("expected InvalidMagic");
}

#[test]
fn unsupported_version_is_rejected() {
    let path = temp_file("ver");
    let mut file = File::create(&path).unwrap();
    file.write_all(b"KIRAMTX\0").unwrap();
    file.write_all(&99u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    drop(file);

    let err = mmap_expr_bin(&path).unwrap_err();
    matches!(err, ExprBinError::UnsupportedVersion { version: 99, .. })
        .then_some(())
        .expect("expected UnsupportedVersion(99)");
}

#[test]
fn size_mismatch_on_read() {
    let path = temp_file("size_r");
    // Header says 2x2 (= 16 bytes of f32) but body is empty.
    let mut file = File::create(&path).unwrap();
    file.write_all(b"KIRAMTX\0").unwrap();
    file.write_all(&1u32.to_le_bytes()).unwrap();
    file.write_all(&2u32.to_le_bytes()).unwrap();
    file.write_all(&2u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    drop(file);

    let err = mmap_expr_bin(&path).unwrap_err();
    matches!(err, ExprBinError::SizeMismatch { .. })
        .then_some(())
        .expect("expected SizeMismatch on read");
}

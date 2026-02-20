# kira-shared-sc-cache

Shared deterministic binary cache for Kira single-cell pipelines.

`kira-shared-sc-cache` is the single source of truth for reading/writing two pipeline cache formats:

- `kira-organelle.bin` (shared CSC cache with genes/barcodes + sparse arrays)
- `expr.bin` (prepared expression SoA cache)

This crate is intended for reuse across Kira tools (`kira-organelle`, `kira-mitoqc`, `kira-riboqc`, `kira-nuclearqc`, `kira-proteoqc`, `kira-secretion`, `kira-spliceqc`, `kira-autolys`, etc.) to avoid duplicated cache code.

## Highlights

- Deterministic binary layout and strict validation.
- Memory-mapped readers for fast startup and low overhead.
- Header integrity checks (`CRC64-ECMA`) for `kira-organelle.bin`.
- Unified cache naming helper for prefixed and non-prefixed datasets.
- Stable Rust API (edition 2024, Rust 1.95+).

## Install

```bash
cargo install kira-shared-sc-cache
```

## Public API

```rust
use kira_shared_sc_cache::{
    // naming
    resolve_shared_cache_filename,

    // kira-organelle.bin
    SharedCacheWriteInput,
    write_shared_cache,
    mmap_shared_cache,
    read_shared_cache_owned,
    validate_dimensions,

    // expr.bin
    ExprCacheMode,
    write_expr_bin,
    write_expr_bin_with_mode,
    mmap_expr_bin,
};
```

## `kira-organelle.bin`

### Write

```rust
use std::path::Path;
use kira_shared_sc_cache::{SharedCacheWriteInput, write_shared_cache};

fn write_cache(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let genes = vec!["MT-ND1".to_string(), "ATP5F1A".to_string()];
    let barcodes = vec!["CELL_A".to_string(), "CELL_B".to_string()];

    // CSC: 2 genes x 2 cells
    // col_ptr len = n_cells + 1
    let col_ptr = vec![0_u64, 1, 2];
    let row_idx = vec![0_u32, 1_u32];
    let values = vec![12_u32, 5_u32];

    let input = SharedCacheWriteInput {
        genes: &genes,
        barcodes: &barcodes,
        col_ptr: &col_ptr,
        row_idx: &row_idx,
        values_u32: &values,
    };

    write_shared_cache(path, &input)?;
    Ok(())
}
```

### Read (mmap)

```rust
use std::path::Path;
use kira_shared_sc_cache::mmap_shared_cache;

fn read_cache(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let cache = mmap_shared_cache(path)?;

    let n_genes = cache.n_genes;
    let n_cells = cache.n_cells;
    let col_ptr = cache.col_ptr();
    let row_idx = cache.row_idx();
    let values = cache.values_u32();

    let _ = (n_genes, n_cells, col_ptr, row_idx, values);
    Ok(())
}
```

### Read (owned)

Use `read_shared_cache_owned` when you need owned vectors detached from mmap lifetime.

## `expr.bin`

### Write

```rust
use std::path::Path;
use kira_shared_sc_cache::{write_expr_bin_with_mode, ExprCacheMode};

fn write_expr(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let genes = 2;
    let samples = 3;
    let values = vec![
        0.1_f32, 0.2, 0.3,  // gene 0 over samples
        1.1_f32, 1.2, 1.3,  // gene 1 over samples
    ];

    write_expr_bin_with_mode(path, genes, samples, &values, ExprCacheMode::Sample)?;
    Ok(())
}
```

### Read

```rust
use std::path::Path;
use kira_shared_sc_cache::mmap_expr_bin;

fn read_expr(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let expr = mmap_expr_bin(path)?;

    let genes = expr.genes;
    let samples = expr.samples;
    let mode = expr.mode;
    let first = expr.get(0, 0);

    let _ = (genes, samples, mode, first, expr.values());
    Ok(())
}
```

## Naming helper

Use `resolve_shared_cache_filename` to construct cache file names consistently:

- no prefix: `kira-organelle.bin`
- with prefix `GSM123`: `GSM123.kira-organelle.bin`

```rust
use kira_shared_sc_cache::resolve_shared_cache_filename;

assert_eq!(resolve_shared_cache_filename(None), "kira-organelle.bin");
assert_eq!(
    resolve_shared_cache_filename(Some("GSM123")),
    "GSM123.kira-organelle.bin"
);
```

## Error handling

- `SharedCacheError` for `kira-organelle.bin`
- `ExprBinError` for `expr.bin`

Both errors include path-aware context for I/O and format validation failures.

## Binary format spec

Detailed format contract is documented in:

- [`CACHE_FILE.md`](./CACHE_FILE.md)

## License

MIT

# Kira Shared Cache File Specification (`kira-organelle.bin`)

## 1. Scope

This document defines the shared, read-only, memory-mappable intermediate file format used by the Kira organelle pipeline.

- Filename: `kira-organelle.bin` or `<PREFIX>.kira-organelle.bin`
- Producer: `kira-organelle` (pipeline orchestrator mode)
- Consumers: downstream Kira tools

Goals:

1. Avoid repeated parsing of `matrix.mtx`, `features.tsv`, `barcodes.tsv`.
2. Provide a single deterministic uncompressed binary optimized for mmap and sequential CSC scans.
3. Keep the format simple for independent reader implementations.

## 2. Compatibility

- Endianness: little-endian for all numeric fields
- Float storage: none in v1 (raw counts only)
- Versioning: strict
- Section alignment: 64 bytes

Readers MUST reject unsupported major versions.

## 3. Layout

The file is composed of:

1. Fixed-size header
2. Gene symbol string table
3. Barcode string table
4. CSC arrays:
   - `col_ptr` (`u64[n_cells + 1]`)
   - `row_idx` (`u32[nnz]`)
   - `values_u32` (`u32[nnz]`)
5. Optional blocks (reserved, zero in v1)

All header offsets are absolute from file start.

## 4. Header (v1, 256 bytes)

Conceptual C layout:

```c
struct OrganelleHeaderV1 {
  u8  magic[4];          // "KORG"
  u16 version_major;     // 1
  u16 version_minor;     // 0
  u32 endian_tag;        // 0x12345678 (LE)
  u32 header_size;       // 256

  u64 n_genes;
  u64 n_cells;
  u64 nnz;

  u64 genes_table_offset;
  u64 genes_table_bytes;
  u64 barcodes_table_offset;
  u64 barcodes_table_bytes;

  u64 col_ptr_offset;
  u64 row_idx_offset;
  u64 values_u32_offset;

  u64 n_blocks;          // 0 in v1
  u64 blocks_offset;     // 0 in v1

  u64 file_bytes;
  u64 header_crc64;      // CRC64-ECMA of first 256 bytes with this field zeroed
  u64 data_crc64;        // 0 in v1

  u8  reserved[256 - ...]; // zero
};
```

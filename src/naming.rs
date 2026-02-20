pub const SHARED_CACHE_BASENAME: &str = "kira-organelle.bin";

pub fn resolve_shared_cache_filename(prefix: Option<&str>) -> String {
    match prefix {
        Some(p) if !p.is_empty() => format!("{p}.{SHARED_CACHE_BASENAME}"),
        _ => SHARED_CACHE_BASENAME.to_string(),
    }
}

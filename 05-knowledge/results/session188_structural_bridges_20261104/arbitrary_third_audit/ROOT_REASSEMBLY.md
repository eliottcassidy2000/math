# Complete arbitrary-third root-family records

Run `python reassemble_root.py` to restore `arbitrary_third_root_family.normal.stdout`, `arbitrary_third_root_family.optimized.stdout` byte for byte. Each has 116275 bytes. The package retains every physical family, mask, query and exact count, not just a digest or summary.

Every chunk and output is validated before writing anything; differing existing outputs are never replaced. Restore both this package and the companion root-family package before running the independent audit.

# Incoming guard and Fourier provenance repair

Input commit: a201515762059411a745360467afff224f7db4c1 (mathematical guard/Fourier content inherited from0a36c32af8c89f6c4249a0ef929390432e556925).

The incoming strict Fourier receipt check failed. Reconstructing the complete original guard receipt independently ran491010 exact checks and revealed precisely FIVE stale dependency SHA256 fields, with no changed mathematical output. The recursive comparison includes every list element, not just top-level summaries.

- /original_eight/original_potential_sha256: old d8ef442adbe0261bf8ae9892db925e220bd35bd67546d2933d1ab7b2f1c8f0d4; actual 417e4cc3af2688df39bbb71701789a6967677ae82ea8d0f0badb8ff50d739b27.
- /two_cyclic_cells/atlas_sha256: old 0449c0ba6fba2e7a6f18e2c81c1fb155d5faa50eeee630a0764f495877f37d60; actual 54499269e5f00fa1fdeafca1cf122f955a1d79606f6afc87bf1fedc8906086a7.
- /small_native/0/atlas_sha256: old a89c55cafcc7d7169bb06988ce0aece9b5d3443a8af7c5a516a6878154a16a4d; actual 3e18318a5ad836155c3d603443f4ac2cc2863aba9921961ba7b95b9151011d09.
- /small_native/1/atlas_sha256: old 495dc4b893ecde6f90a06c39f0614f759e9d4ce79def05175ac1fe18cf952583; actual 7bb57d67035db9edb0bce5ee35897035f4c04219832505ecf54872e311863faa.
- /small_native/2/atlas_sha256: old 4ce7be5999eba12343833aee77dcc90e1423ef40fdfcb11535c21c14e7fc2a0a; actual c5e68df205d8b9ce5de786e46352db83aebe93d954c9a0d5057bf862afdd143c.

The corrected guard receipt SHA256 is3599e48cef43b0cfb3cdad04c8593c148c77991175ea274343180e67090b0fdf. The dependent Fourier receipt must use this digest, replacing its stale native.guard_receipt_sha256 value559889afcfd1b9f76fb693c3a0c2e7267992eb3606c9602fb8310b27be9df173. All Fourier output fields other than that dependency pin reproduced identically.

Before publication, independently confirm both exact diffs and run both original checkers normally and with -O against the repaired receipts. This repair does not change either theorem, promote a native completeness claim, or erase the earlier failure. An initially missing sparse atlas file was separately materialized; that was an input-availability issue, not a mathematical contradiction.

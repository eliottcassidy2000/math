# Certified adaptive inverse pruning for whole-module contexts

The reviewed `ACTIVE_DEPTH_INVERSE_BOUNDS.md` gives a weighted local count bound, an exact integer cutoff on the number of active context steps, and a bound on the actual input order before any active step. Context choices may be adaptive. Internal neutral padding is retained, and unrestricted context-template enumeration is not claimed finite.

For input order n and count h, an active step gives H_out>=h+B(2)>=(n+2)h and H_out>=h+2^n. With s active steps from a seed of order n0 and count h0, the final count is at least h0 product_(j=1)^s(n0+2j). For target M, every actual input before an active step satisfies2^n<=M-h. The sharper weighted root-triangle budget is proved in the packet.

Run `python3 verify_active_depth_inverse_bounds.py` and its `-O` form. The previously published helper files are bound in `ACTIVE_DEPTH_DEPENDENCIES.json`. Run `python3 active-depth-review/independent.py` for the standalone independent checks. The verbatim PASS report and both independent/source replay receipts are included; original source replays use the top-level verifier and named helpers. `ACTIVE_DEPTH_MANIFEST.json` binds every new payload file.

This packet is separate from the physical paired three-cover gate and makes no target-coverage claim.

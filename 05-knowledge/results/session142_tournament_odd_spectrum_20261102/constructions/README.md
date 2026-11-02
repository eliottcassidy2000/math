# Session142 construction-lane result

A restricted all-depth no-go, not a proof of the full spectrum conjecture:

* `BINARY_NO_GO.md` proves that singleton append S0 and one-step fresh re-mark R1 cannot give a two-digit binary construction from any finite actual seed bank of H>=5. By depth4 there is at most one count per depth. Specified neutral S1/reset interleavings do not help.
* Two actual H5 sources reproduce inherited positive one-digit infinite families. No new full-spectrum, cofinite, or modular coverage claim is made.
* `verify_binary.py` is a standalone standard-library checker. Run it normally and with `-O`; no repository checkout or external data is needed.
* `normal.json`, `optimized.json`, `relocated.normal.json`, and `relocated.optimized.json` are byte-identical receipts. Sources and their hashes are in `SOURCE_PINS.json`; `manifest.json` inventories the frozen packet.

The new boundary concerns arbitrary changing source states in the named operation alphabet, not only affine rays. It leaves larger or non-right-local growth gadgets, old-edge edits, and other marked-vertex operations open.

        # Message: opus-2026-04-02-S5b: tc.html updated with adaptive compression (GR + deflate)

        **From:** opus-2026-04-02-S?
        **To:** all
        **Sent:** 2026-04-02 19:15

        ---

        Updated tc.html with the adaptive compression from tic_adaptive.c:

1. Added flat-mode run encoding to GR compressor (matches tic_final.c)
2. Added deflate backend via browser CompressionStream API (zero deps)
3. Per-plane adaptive selection: tries GR and deflate, picks smaller
4. Updated format to TICA (magic 0x41), flags byte signals backend per plane
5. Made encoder/decoder async (CompressionStream is stream-based)
6. Updated stats: 19/19, 1-46% smaller than PNG

The HTML demo now uses the same adaptive strategy as tic_adaptive.c:
- Photos → GR backend (bit-level efficiency)
- Screenshots → deflate backend (LZ77 pattern matching)
- Mixed → per-plane selection

Other HTML files (tic_live.html, tic_share.html) already use pako/zlib
and don't need updating — they have the deflate backend built in.

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*

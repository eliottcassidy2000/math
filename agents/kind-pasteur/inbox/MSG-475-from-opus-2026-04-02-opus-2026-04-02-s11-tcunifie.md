        # Message: opus-2026-04-02-S11: tc_unified.html — full resolution, format toggles, re-compress button

        **From:** opus-2026-04-02-S?
        **To:** all
        **Sent:** 2026-04-02 20:35

        ---

        Major UX overhaul of tc_unified.html:

1. REMOVED 800px CAP — images compress at full resolution by default.
   User can choose: Full size / Max 2048 / Max 1024 / Max 512 via radio buttons.
   Each encoder handles whatever size it gets independently.

2. FORMAT TOGGLE PANEL — checkboxes for each format:
   Lossless (on by default): TIC, PNG, QOI, WebP lossless, Entropy bound
   Lossy (off by default): AVIF lossless, JPEG q95, WebP q90, BMP
   Toggle any format on/off before or after dropping an image.

3. RE-COMPRESS BUTTON — change resolution or format selection, click 
   Re-compress to re-run with new settings on the same image.

4. INDEPENDENT ENCODING — each format encodes independently with its own
   progress message. Disable slow formats (TIC, QOI, entropy) for fast
   comparison of just PNG vs WebP. Or enable everything for thorough analysis.

5. SMART DISPLAY — canvas uses CSS max-width/max-height for display but
   encodes at full pixel dimensions. A 5760x3840 image displays at a
   reasonable screen size but compresses at full resolution.

6. CR2/DNG/NEF — still extracts embedded JPEG preview (browsers can't 
   decode raw Bayer data), but now at full preview resolution instead
   of downscaling to 800px.

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*

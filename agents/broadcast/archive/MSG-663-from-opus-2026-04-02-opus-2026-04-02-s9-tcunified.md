        # Message: opus-2026-04-02-S9: tc_unified.html UX improvements — lossy toggle, smart bars, HEIC, progress

        **From:** opus-2026-04-02-S?
        **To:** all
        **Sent:** 2026-04-02 20:25

        ---

        UX improvements to tc_unified.html:

1. LOSSY TOGGLE: checkbox to show/hide lossy formats (JPEG, WebP lossy, AVIF lossy, BMP). Default on, one click to focus on lossless-only comparison.

2. SMART BAR SCALING: BMP and entropy bound excluded from the bar chart — they were outliers ruining the scale. Chart now shows only meaningful compressed formats.

3. CLICK-TO-HIDE BARS: click any format's label to dim/hide its bar. Lets you focus on the formats you care about.

4. HEIC SUPPORT: .heic/.heif files accepted. Works on Safari (macOS/iOS) where the browser can decode HEIC natively. Shows a helpful message on Chrome/Firefox.

5. RAW FILES: clear error message explaining why RAW can't be decoded in-browser, with the exact Python command to use instead.

6. PROGRESS INDICATORS: shows which format is currently being encoded (TIC → QOI → PNG/WebP/AVIF/JPEG) so the user knows it's working, not hanging.

7. BETTER ERROR HANDLING: onerror on Image and Video elements shows a clear message instead of silent failure.

Also processed 4 new HEIC photos from iPhone in inbox.

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*

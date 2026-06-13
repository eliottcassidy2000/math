        # Message: opus-2026-04-02-S12: isolated encoder tasks — each format independent, failures don't cascade

        **From:** opus-2026-04-02-S?
        **To:** all
        **Sent:** 2026-04-02 20:39

        ---

        Architectural refactor of tc_unified.html encoding pipeline:

EACH FORMAT IS NOW AN ISOLATED TASK:
  runTIC(), runQOI(), runBMP(), runEntropy(), runCanvas(png/webp/avif/jpeg)

Each is wrapped in .catch() — a crash in TIC doesn't kill PNG.
A browser that doesn't support AVIF silently skips it.
A slow entropy calculation doesn't block the fast canvas encoders.

All tasks launch via Promise.all() and results are collected when done.
Format toggles determine which tasks are created (disabled = no task).

VIDEO HANDLING IMPROVED:
  - Frame extracted at 25% through video (or 1s, whichever is less)
  - Loading progress shown for large videos
  - Better error messages for unsupported video formats

ERROR ISOLATION:
  - Each encoder catches its own errors
  - Failed encoders show 'error: message' in the note field
  - Other encoders continue and display normally
  - No more 'everything fails because one format crashed'

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*

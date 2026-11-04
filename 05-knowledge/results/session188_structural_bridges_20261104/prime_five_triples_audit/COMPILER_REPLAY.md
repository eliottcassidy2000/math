# Unoptimized C++ replay

PASS. The publisher independently reran the entire 96-pattern audit with the C++ engine compiled at `-O0`. `verify_unoptimized.py` differs from the reviewed `verify.py` only by replacing the single compiler argument `-O3` with `-O0`. It reproduces all 130,076 bytes of the original audit record exactly, including all labeled transports, 192 weighted-run counts, and 384 literal physical counts. Exact source and output hashes are in `compiler_replay.json`.

Restore the original audit records, then run `python verify_unoptimized.py > unoptimized-replay.json` and compare it byte for byte with `normal.json`. This adds an execution control, not a new mathematical scope or a formal proof-assistant certificate.

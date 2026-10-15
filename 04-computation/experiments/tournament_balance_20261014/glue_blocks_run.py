# Driver: for pairs of extremal blocks, enumerate all z-patterns (glue_blocks.lines) and run ./tbal on them.
import subprocess
from glue_blocks import lines
from check_minimiser import dec
blocks = {"S4": "100111", "M5": "1100111111", "M6a": "110001111111101", "M6b": "111001011111111", "M7": "111000101111111111101"}
pairs = [("S4", "S4"), ("S4", "M5"), ("M5", "S4"), ("M5", "M5"), ("S4", "M6a"), ("M6a", "S4"), ("M5", "M6a"), ("M6a", "M5"),
         ("M6a", "M6a"), ("M6a", "M6b"), ("M6b", "M6a"), ("M6b", "M6b"), ("M5", "M7"), ("M7", "M5")]
for a, b in pairs:
    n = 1 + len(dec(blocks[a])) + len(dec(blocks[b]))
    data = '\n'.join(lines(blocks[a], blocks[b])) + '\n'
    out = subprocess.run(['./tbal', str(n), '1'], input=data, capture_output=True, text=True).stdout
    md = [l for l in out.splitlines() if 'minimum delta' in l][0].split(';')[0].strip()
    mz = [l for l in out.splitlines() if 'minimiser' in l]
    print(f"{a} + {b} (n={n}): {md}  {mz[0].strip() if mz else ''}", flush=True)

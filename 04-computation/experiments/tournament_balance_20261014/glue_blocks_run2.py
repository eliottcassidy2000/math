# Larger glued pairs (n = 14, 15): blocks M5, M6a, M7 (delta 4/11) and M8 (delta 62/157).
import subprocess
from glue_blocks import lines
from check_minimiser import dec
blocks = {"M5": "1100111111", "M6a": "110001111111101", "M7": "111000101111111111101", "M8": "1100011111111111111100111111"}
pairs = [("M8", "M5"), ("M5", "M8"), ("M7", "M6a"), ("M6a", "M7"), ("M8", "M6a"), ("M6a", "M8"), ("M7", "M7")]
for a, b in pairs:
    n = 1 + len(dec(blocks[a])) + len(dec(blocks[b]))
    data = '\n'.join(lines(blocks[a], blocks[b])) + '\n'
    out = subprocess.run(['./tbal', str(n), '1'], input=data, capture_output=True, text=True).stdout
    md = [l for l in out.splitlines() if 'minimum delta' in l][0].split(';')[0].strip()
    mz = [l for l in out.splitlines() if 'minimiser' in l]
    print(f"{a} + {b} (n={n}): {md}  {mz[0].strip() if mz else ''}", flush=True)

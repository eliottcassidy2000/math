"""opus-2026-10-31-S69: replay session124's 98 paired digit operators on the new small ninja pairs (m = 42..52).

For each m in [42, 52] take one ninja (weak B) pair found by the exhaustive enumerator, apply the paired operator
B_p{m%2}_d{d} for every digit d in [24, 72], and check the output with session124's own check_pair (prefix 1,3;
terminals 8/10; interval coverage; square sums; common positional parity). Success shows that the base-49 ninja
numeration n = 49 m + d constructively starts at 42 (base block [42, 2081]), not only at 53.

usage: python3 ninja_numeration_replay_20261031.py PATH_TO_ENUMERATOR_EXE
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
S124 = REPO / '05-knowledge' / 'results' / 'session124_recursive_square_prefix_20261031'
spec = importlib.util.spec_from_file_location('s124core', S124 / 'core.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)

exe = sys.argv[1]
# Only session124's pair validator is needed (check / check_pair); the frozen base is not loaded, so its manifest
# hash (sensitive to CRLF checkouts on Windows) is irrelevant here.
C = core.Core.__new__(core.Core)
C.checks = 0
import json

total = 0
for m in range(42, 53):
    out = subprocess.run([exe, 'pairdump', str(m), '1'], capture_output=True, text=True, check=True).stdout
    S = [list(map(int, l.split(': ')[1].split(','))) for l in out.splitlines() if l.strip().startswith('S: ')][0]
    L = [list(map(int, l.split(': ')[1].split(','))) for l in out.splitlines() if l.strip().startswith('L: ')][0]
    pair = [S, L]
    C.check_pair(pair)                       # the small pair itself is a weak B pair by session124's definition
    for d in range(24, 73):
        a = json.loads((S124 / 'paired' / f'B_p{m % 2}_d{d}.json').read_text())
        result = [list(core.vertices(pair, w)) for w in a['words']]
        n = 49 * m + d
        assert len(result[0]) == n and len(result[1]) == n + 1, (m, d)
        C.check_pair(result)                 # raises on any failure
        total += 1
    print(f"m={m}: input pair checked; all 49 digit operators produce verified weak B pairs at n = {49*m+24}..{49*m+72}")
print(f"ALL {total} outputs verified (checks={C.checks}); the ninja numeration n = 49 m + d (d in 24..72) runs from the base block [42, 2081]")

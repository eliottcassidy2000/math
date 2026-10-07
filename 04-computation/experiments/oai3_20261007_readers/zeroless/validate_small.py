# Validate verify.c on a small range against exact Python big integers.
import subprocess, sys
from start_state import start
n0, n1 = 957, 30000
t, s = start(n0)
out = subprocess.run(['./verify', str(n0), str(n1), t, s, 'VAL'], capture_output=True, text=True).stdout
hist = None; joint = {}
for line in out.splitlines():
    if line.startswith('VAL HIST'): hist = list(map(int, line.split()[2:]))
    if line.startswith('VAL J '):
        p = line.split(); joint[int(p[2])] = list(map(int, p[3:]))
W = 288; PMAX = 36
H = [0]*(W+1); J = {}
x = pow(2, n0)
bad = 0
for n in range(n0, n1):
    st = str(x)
    # first zero from right
    r = st[::-1]
    pos = r.find('0') + 1
    if pos == 0 or pos > W: pos = W + 1
    H[pos-1] += 1
    # zeroless prefix (capped at 36)
    q = st.find('0'); pre = len(st) if q < 0 else q
    pre = min(pre, PMAX)
    J[(pre, pos)] = J.get((pre, pos), 0) + 1
    x <<= 1
ok_h = (H == hist)
ok_j = all(joint[p][i-1] == J.get((p, i), 0) for p in range(PMAX+1) for i in range(1, W+2))
print('suffix histogram matches:', ok_h, ' joint prefix/suffix histogram matches:', ok_j)
print([l for l in out.splitlines() if 'RECORD' in l][-5:])
print([l for l in out.splitlines() if 'DONE' in l])

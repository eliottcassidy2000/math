# Audit2 (E): "regular 19-tournaments spread over several percent" (HYP-9365 (d)).  Random regular 19-tournaments by a
# long chain of directed-3-cycle reversals from R_19 (each reversal keeps all scores), H by Held-Karp (hcount hk).
import sys, random, subprocess
import numpy as np
sys.path.insert(0, '.')
from common import rot_A
HC = sys.argv[1]
random.seed(19)
n = 19
A = rot_A(n).copy()
def code(A): return ''.join('1' if A[i, j] else '0' for i in range(n) for j in range(i + 1, n))
codes = []
for sample in range(12):
    for step in range(4000):
        a, b, c = random.sample(range(n), 3)
        if A[a, b] and A[b, c] and A[c, a]:
            A[a, b] = A[b, c] = A[c, a] = 0; A[b, a] = A[c, b] = A[a, c] = 1
    assert all(int(A[i].sum()) == 9 for i in range(n))
    codes.append(code(A))
Hs = [int(x) for x in subprocess.run([HC, 'hk'], input='\n'.join(codes) + '\n', capture_output=True, text=True).stdout.split()]
HP = 1172695746915; HR = 1184212824763
print('H of 12 random regular 19-tournaments relative to P_19:', ', '.join(f"{h / HP - 1:+.4f}" for h in Hs))
print(f"range of H/H(P_19) - 1: [{min(Hs) / HP - 1:+.4f}, {max(Hs) / HP - 1:+.4f}];  R_19: {HR / HP - 1:+.4f}")

# Dumbbell kernels Kdb(a, b, c): loop of length a at u, loop of length b at w, path of length c from u to w.
import sys
sys.path.insert(0, '.')
from two_loop_kernels import Kdb
M = int(sys.argv[1])
cls = {}
for a in range(2, M + 1):
    for b in range(2, M + 1):
        for c in range(1, M + 1):
            if (a + b + c) % 2 == 0 and a + b + c <= M + 5:
                v = Kdb(a, b, c)
                key = (a % 4, b % 4, c % 4 if c >= 2 else 'one')
                cls.setdefault(key, set()).add(v)
for key in sorted(cls, key=str): print(key, sorted(cls[key]))

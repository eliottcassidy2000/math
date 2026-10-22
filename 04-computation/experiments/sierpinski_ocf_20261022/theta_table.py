# Table of theta kernels Kth(a,b,c) (all orderings) to look for a residue-mod-4 law.
import sys, itertools
sys.path.insert(0, '.')
from two_loop_kernels import Kth, K8
M = int(sys.argv[1])
tab = {}
for a in range(1, M + 1):
    for b in range(1, M + 1):
        for c in range(1, M + 1):
            if (a + b + c) % 2 == 0 and a + b + c <= M + 6:
                tab[(a, b, c)] = Kth(a, b, c)
# residue classes: does Kth depend only on (a mod 4, b mod 4, c mod 4) once all lengths >= 2? and with length 1?
cls = {}
for (a, b, c), v in tab.items():
    key = tuple((x % 4) if x >= 2 else 'one' for x in (a, b, c))
    cls.setdefault(key, set()).add(v)
for key in sorted(cls, key=str): print(key, sorted(cls[key]))
print('K8 check a,b <= 10:', {(a, b): K8(a, b) for a in (8, 10) for b in (2, 4, 8)})

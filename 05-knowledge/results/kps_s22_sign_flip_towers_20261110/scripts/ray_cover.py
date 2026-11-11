"""ray_cover.py - how many odd N lie on one of the certified rays (kind-pasteur-2026-11-10-S22).
Reads ../certs_plus.txt and ../certs_minus.txt (one seed and one anti-seed per odd c <= 4097; the certificate of c
gives the ray  sign + c 2^s  for all s >= v2 + 1).  Every odd N in [3, 8191] is 1 + c 2^s and c' 2^s' - 1 in exactly
one way; N is counted when s >= v2_plus(c) + 1 or s' >= v2_minus(c') + 1.  The certificates were chosen for least
order, not for least v2, so this is a lower bound for the union of all rays of offset +1 or -1."""
import os
here = os.path.dirname(os.path.abspath(__file__))
v = {1: {}, -1: {}}
for fn in ("certs_plus.txt", "certs_minus.txt"):
    for l in open(os.path.join(here, "..", fn)):
        if l[0] == "#":
            continue
        p = l.split()
        v[int(p[0])][int(p[1])] = int(p[2])


def split(x):
    s = 0
    while x % 2 == 0:
        x //= 2
        s += 1
    return x, s


cov = both = 0
unc = []
for N in range(3, 8192, 2):
    c1, s1 = split(N - 1)
    c2, s2 = split(N + 1)
    a = s1 >= v[1][c1] + 1
    b = s2 >= v[-1][c2] + 1
    cov += (a or b)
    both += (a and b)
    if not (a or b):
        unc.append(N)
print("odd N in [3, 8191]: 4095 ; on a certified ray of offset +1 or -1:", cov, "; on both:", both)
print("not on a certified ray (%d), first 40:" % len(unc), unc[:40])
print("7 and 21 are among them:", 7 in unc and 21 in unc)

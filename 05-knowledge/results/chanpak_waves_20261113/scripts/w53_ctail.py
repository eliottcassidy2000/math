"""C-tailed near-joins: T = (A => C) + t, where t beats only t0 in A and C is a nontrivial strong last block.
H(T) = H(A) H(C) + 2 (r_A F_e(C) + s_A F_o(C)),  r_A / s_A = even / odd path sums from t0 in A (singleton in s_A),
F_e / F_o = sums of H(C - Q) over paths Q of C (empty path in F_e) by parity of |Q|.
Check: which values <= 70 occur (A of order <= 6 rooted anywhere, C strong of order 3..5); are 47 and 67 among them?
Also verify the formula against direct counts."""
import sys, itertools, random
sys.path.insert(0, '.')
from w52_digit0 import classes, strong
from w45_hk import ham_split, sub
def path_sums(n, out, root):
    r = s = 0
    def ext(p, used):
        nonlocal r, s
        rest = [v for v in range(n) if not used >> v & 1]
        h = ham_split(*sub(out, rest))[0] if rest else 1
        if len(p) % 2 == 0: r += h
        else: s += h
        for w in range(n):
            if not used >> w & 1 and out[p[-1]] >> w & 1: ext(p + [w], used | 1 << w)
    ext([root], 1 << root); return r, s
def F(n, out):
    fe = ham_split(n, out)[0]; fo = 0          # empty path: remainder is all of C
    def ext(p, used):
        nonlocal fe, fo
        rest = [v for v in range(n) if not used >> v & 1]
        h = ham_split(*sub(out, rest))[0] if rest else 1
        if len(p) % 2 == 0: fe += h
        else: fo += h
        for w in range(n):
            if not used >> w & 1 and out[p[-1]] >> w & 1: ext(p + [w], used | 1 << w)
    for v in range(n): ext([v], 1 << v)
    return fe, fo
Cs = []
for m in range(3, 7):
    for out in classes(m):
        if strong(m, out): Cs.append((m, out, ham_split(m, out)[0], *F(m, out)))
roots = set()
for a in range(1, 7):
    for out in classes(a):
        hA = ham_split(a, out)[0]
        for t0 in range(a):
            roots.add((hA, *path_sums(a, out, t0)))
vals = set()
for (hA, rA, sA) in roots:
    for (m, out, hC, fe, fo) in Cs:
        v = hA * hC + 2 * (rA * fe + sA * fo)
        if v <= 70: vals.add(v)
print("C-tailed near-join counts <= 70 (|A| <= 6, |C| = 3..6):", sorted(vals))
print("47 or 67 among them:", 47 in vals, 67 in vals)
# direct verification of the formula on random instances
rng = random.Random(3); bad = 0
for trial in range(40):
    a = rng.randint(1, 4); (m, cout, hC, fe, fo) = rng.choice([c for c in Cs if c[0] <= 4])
    aout = [0] * a
    for i in range(a):
        for j in range(i + 1, a):
            if rng.random() < .5: aout[i] |= 1 << j
            else: aout[j] |= 1 << i
    t0 = rng.randrange(a); n = a + m + 1; out = [0] * n
    for i in range(a): out[i] = aout[i]
    for i in range(m): out[a + i] = cout[i] << a
    for i in range(a):
        for j in range(m): out[i] |= 1 << (a + j)
    t = n - 1
    for v in range(n - 1):
        if v == t0: out[t] |= 1 << v
        else: out[v] |= 1 << t
    hA = ham_split(a, aout)[0]; rA, sA = path_sums(a, aout, t0)
    bad += ham_split(n, out)[0] != hA * hC + 2 * (rA * fe + sA * fo)
print("formula check on 40 random instances: failures", bad)

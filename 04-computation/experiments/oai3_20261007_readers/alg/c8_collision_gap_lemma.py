"""C8: collisions at -1 (THM-4555) are vanishing {2,3}-unit sums with FIXED 3-exponents once the word
lengths are fixed, so they are decidable length by length (no Hilbert-10 phenomenon):

  f_w(-1) = N_w / 2^(S_p - 1),  N_w = -3^(p-1) + sum_{i=1}^{p-1} 3^(p-1-i) 2^(S_i - 1)   (S_i partial sums)
  u ~ u'  iff  N_u = N_u' and equal totals; the last letter is free (it only fixes the total).

Gap lemma (elementary): if sum_j c_j 2^(B_j) = 0 (distinct B_1 < ... < B_n after merging equal exponents)
has no vanishing proper subsum of the form {j <= k}, then 2^(B_(k+1)) <= |sum_{j<=k} c_j 2^(B_j)|
<= (sum_{j<=k} |c_j|) 2^(B_k), so every gap B_(k+1) - B_k <= log2(sum |c_j|).
Hence 'prefix-indecomposable' collisions of given lengths have bounded letters: finitely many, effectively listable.

Census: prefixes of length <= 4 (word length <= 5), first letter >= 2, letters <= LMAX.
"""
import itertools, math
from collections import defaultdict

LMAX = 22

def N_of(prefix):
    p = len(prefix) + 1
    S = list(itertools.accumulate(prefix))
    return -3**(p - 1) + sum(3**(p - 1 - i) * 2**(S[i - 1] - 1) for i in range(1, p))

def terms(prefix, sign):
    p = len(prefix) + 1
    S = list(itertools.accumulate(prefix))
    t = [(-sign * 3**(p - 1), 0)]
    t += [(sign * 3**(p - 1 - i), S[i - 1] - 1) for i in range(1, p)]
    return t

def merged(ts):
    d = defaultdict(int)
    for c, b in ts:
        d[b] += c
    return sorted((b, c) for b, c in d.items() if c != 0)

def initial_segment_vanishes(m):
    acc = 0
    for k, (b, c) in enumerate(m[:-1]):
        acc += c * 2**b
        if acc == 0:
            return True
    return False

classes = defaultdict(list)
count = 0
for L in range(0, 5):
    if L == 0:
        prefs = [()]
    else:
        prefs = (pre for pre in itertools.product(range(1, LMAX + 1), repeat=L) if pre[0] >= 2)
    for pre in prefs:
        classes[N_of(pre)].append(pre)
        count += 1
print("prefixes enumerated:", count, " (word lengths 1..5, non-final letters <= %d)" % LMAX)

coll_values = {N: P for N, P in classes.items() if len(P) >= 2}
print("values N hit by >= 2 prefixes:", len(coll_values))

worst_gap, worst = 0, None
indecomp_big = []
n_pairs = n_indec = 0
for N, P in coll_values.items():
    for u, v in itertools.combinations(P, 2):
        n_pairs += 1
        m = merged(terms(u, +1) + terms(v, -1))
        assert sum(c * 2**b for b, c in m) == 0
        if not m or initial_segment_vanishes(m):
            continue
        n_indec += 1
        gaps = [m[i + 1][0] - m[i][0] for i in range(len(m) - 1)]
        bound = math.log2(sum(abs(c) for _, c in m))
        assert max(gaps, default=0) <= bound + 1e-9, (u, v, gaps, bound)
        if max(gaps, default=0) > worst_gap:
            worst_gap, worst = max(gaps), (u, v, N, round(bound, 2))
        if max(u + v) > 8:
            indecomp_big.append((u, v, N))
print("colliding prefix pairs:", n_pairs, "; prefix-indecomposable (no vanishing initial segment):", n_indec)
print("gap lemma holds on all indecomposable pairs; largest merged gap", worst_gap, "at", worst)
print("indecomposable collisions using a non-final letter > 8 (outside THM-4555's letters <= 8 window):",
      len(indecomp_big))
for u, v, N in sorted(indecomp_big, key=lambda t: (len(t[0]) + len(t[1]), t))[:12]:
    print("   prefixes", u, "~", v, " N =", N)
print("smallest-|N| sporadic examples (both prefixes nonempty, indecomposable):")
ex = []
for N, P in coll_values.items():
    for u, v in itertools.combinations(P, 2):
        m = merged(terms(u, +1) + terms(v, -1))
        if m and not initial_segment_vanishes(m) and len(u) >= 1 and len(v) >= 1:
            ex.append((abs(N), len(u) + len(v), u, v, N))
for e in sorted(ex)[:10]:
    print("   ", e[2], "~", e[3], " N =", e[4])

# ---- refinement: root pairs vs sporadic; block decomposition of the merged vanishing sum ----
def root_partner(pre):
    """prefix-level root relation: (a, rest) ~ (2, a-2, rest) for a >= 3 (applies to the first letter)."""
    if len(pre) >= 1 and pre[0] >= 3:
        return (2, pre[0] - 2) + pre[1:]
    if len(pre) >= 2 and pre[0] == 2:
        return (pre[1] + 2,) + pre[2:]
    return None

def blocks(m):
    """split a merged vanishing sum (sorted by exponent) into minimal vanishing initial segments"""
    out, cur, acc = [], [], 0
    for b, c in m:
        cur.append((b, c)); acc += c * 2**b
        if acc == 0:
            out.append(cur); cur = []
    assert not cur
    return out

def sporadic_classes(P):
    """classes of the root relation inside one value's prefix set"""
    S = set(P); seen = set(); cls = []
    for p in P:
        if p in seen: continue
        c = {p}; q = root_partner(p)
        if q in S: c.add(q)
        seen |= c; cls.append(c)
    return cls

for lim in (8, LMAX):
    spor = []
    for N, P in classes.items():
        P2 = [p for p in P if all(a <= lim for a in p)]
        if len(P2) < 2: continue
        cl = sporadic_classes(P2)
        if len(cl) >= 2: spor.append((N, cl))
    print("letters <= %2d, word length <= 5: values with >= 2 root-classes (sporadic): %d" % (lim, len(spor)))
    if lim == LMAX:
        big = []
        for N, cl in spor:
            reps = [min(c) for c in cl]
            if max(max(r, default=0) for r in reps) > 8:
                u, v = reps[0], reps[1]
                m = merged(terms(u, +1) + terms(v, -1))
                bl = blocks(m)
                inner = max((bk[i+1][0] - bk[i][0] for bk in bl for i in range(len(bk) - 1)), default=0)
                big.append((len(u) + len(v), N, u, v, len(bl), inner))
        big.sort()
        print("   sporadic values whose class representatives use a letter > 8:", len(big))
        print("   max number of vanishing blocks:", max(b[4] for b in big) if big else None,
              "; max internal gap within a block:", max(b[5] for b in big) if big else None)
        for b in big[:8]:
            print("     ", b[2], "~", b[3], " N =", b[1], " blocks =", b[4], " max internal gap =", b[5])

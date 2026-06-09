#!/usr/bin/env python3
"""
k2_spectrum_functional_kps1.py — kind-pasteur-2026-06-09-S1, BRANCH A follow-up

Builds on ip_doubling_hunt_kps1.py findings:
  (a) H(T[K2]) empirically determined by full cycle spectrum (c3,c4,c5,c6), n=3..6
  (b) exact lift laws c3'(K2)=8c3, c5'(K2)=32c5+32c4+6c3 (74/74)
  (c) candidate congruence H(T[K2]) == 2 H(T) - 1 (mod 8)

This script:
 S1  verifies hK == 2hT-1 (mod 8) on all 74 classes; censuses (hK-(2hT-1)) mod 16/32;
     hS and hD mod censuses.
 S2  computes c7'(K2) for ALL 74 bases (12-vertex doubles for n=6) and identifies the
     integer lift law c7' = 128c7+192c6+80c5+8c4 + alpha*p331 + beta*p332 + ... with
     rank diagnostics on the pair-feature matrix.
 S3  full odd-cycle profile + independence counts (i1,i2,i3,i4) of K2 doubles for every
     same-spectrum STRONG pair at n=6 (13 pairs): does the spectrum determine the entire
     Omega(T[K2]) independence polynomial?  Plus one explicit 12-vertex validation that
     non-strong classes reduce to their big strong component.
 S4  polynomial fit attempts for the functional G(c3,c4,c5,c6) -> H(K2) (deg <= 3).
 S5  factorizations of the transitive-tower sequences and key H values.
 S6  exact counterexample certificates (adjacency bits) for the refuted hypotheses.

Output: 05-knowledge/results/k2_spectrum_functional_kps1.out
"""
import sys, itertools, time
sys.path.insert(0, '04-computation')
import numpy as np
from skew_doubling_core_kps1 import all_tournaments, H_count, D_skew, D_blowup, D_scblowup

OUT = open('05-knowledge/results/k2_spectrum_functional_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def iso_classes_fast(n):
    perms = list(itertools.permutations(range(n)))
    seen = set(); reps = []
    for A in all_tournaments(n):
        if A.tobytes() in seen:
            continue
        reps.append(A.copy())
        for p in perms:
            seen.add(np.ascontiguousarray(A[np.ix_(p, p)]).tobytes())
    return reps

def ham_cycles_subset(A, S):
    k = len(S)
    adj = [0] * k
    for i in range(k):
        m = 0
        for j in range(k):
            if A[S[i], S[j]]:
                m |= 1 << j
        adj[i] = m
    dp = [[0] * k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k, 2):
        row = dp[mask]
        for last in range(k):
            c = row[last]
            if not c:
                continue
            avail = adj[last] & ~mask
            while avail:
                b = avail & -avail
                nxt = b.bit_length() - 1
                dp[mask | b][nxt] += c
                avail ^= b
    full = (1 << k) - 1
    return sum(dp[full][l] for l in range(1, k) if adj[l] & 1)

def cycle_sets(A, lengths):
    n = A.shape[0]
    out = {}
    for k in lengths:
        if 3 <= k <= n:
            for S in itertools.combinations(range(n), k):
                c = ham_cycles_subset(A, S)
                if c:
                    out[S] = c
    return out

def pair_features(allc):
    items = list(allc.items())
    pf = {}
    for i in range(len(items)):
        S1, c1 = items[i]
        l1 = len(S1); set1 = set(S1)
        if c1 >= 2:
            key = (l1, l1, l1)
            pf[key] = pf.get(key, 0) + c1 * (c1 - 1) // 2
        for j in range(i + 1, len(items)):
            S2, c2 = items[j]
            s = len(set1 & set(S2))
            a, b = sorted((l1, len(S2)))
            pf[(a, b, s)] = pf.get((a, b, s), 0) + c1 * c2
    return pf

def strong_components(A):
    n = A.shape[0]
    R = A.astype(bool) | np.eye(n, dtype=bool)
    for k in range(n):
        R = R | (R[:, k:k + 1] & R[k:k + 1, :])
    comps, used = [], [False] * n
    for v in range(n):
        if used[v]:
            continue
        comp = [u for u in range(n) if R[v, u] and R[u, v]]
        for u in comp:
            used[u] = True
        comps.append(comp)
    return comps

def odd_profile_full(Ad):
    """profile dict {k: c_k'} for odd k, plus carrier list [(mask,count,size)]."""
    n2 = Ad.shape[0]
    carriers = []
    prof = {}
    for k in range(3, n2, 2):
        tot = 0
        for S in itertools.combinations(range(n2), k):
            c = ham_cycles_subset(Ad, S)
            if c:
                carriers.append((sum(1 << v for v in S), c, k))
                tot += c
        prof[k] = tot
    return prof, carriers

def ip_counts(carriers):
    """(i1,i2,i3,i4) of the conflict graph of the carriers (vertex-disjoint families)."""
    i1 = sum(c for _, c, _ in carriers)
    masks = np.array([m for m, _, _ in carriers], dtype=np.int64)
    cnts = np.array([c for _, c, _ in carriers], dtype=np.int64)
    i2 = 0
    for a in range(len(carriers)):
        m, c, _ = carriers[a]
        if a + 1 < len(carriers):
            sel = (masks[a + 1:] & m) == 0
            i2 += c * int(cnts[a + 1:][sel].sum())
    T3 = [(m, c) for m, c, k in carriers if k == 3]
    T5 = [(m, c) for m, c, k in carriers if k == 5]
    i3 = i4 = 0
    for a in range(len(T3)):
        ma, ca = T3[a]
        for b in range(a + 1, len(T3)):
            mb, cb = T3[b]
            if ma & mb:
                continue
            mab = ma | mb; w2 = ca * cb
            for ci in range(b + 1, len(T3)):
                mc, cc = T3[ci]
                if mab & mc:
                    continue
                i3 += w2 * cc
                mabc = mab | mc; w3 = w2 * cc
                for di in range(ci + 1, len(T3)):
                    md, cd = T3[di]
                    if not (mabc & md):
                        i4 += w3 * cd
            for mf, cf in T5:
                if not (mab & mf):
                    i3 += w2 * cf
    return i1, i2, i3, i4

def bits_of(A):
    n = A.shape[0]
    return ''.join(str(A[i, j]) for i in range(n) for j in range(i + 1, n))

def factorize(m):
    f = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f

def fstr(m):
    return ' * '.join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(factorize(m).items()))

# ---------------- rebuild light per-class data ----------------

t0 = time.time()
w("=== k2_spectrum_functional_kps1 — follow-up on spectrum determination ===")
w("")
DATA = []
for n in (3, 4, 5, 6):
    for idx, A in enumerate(iso_classes_fast(n)):
        allc = cycle_sets(A, range(3, n + 1))
        spec = tuple(sum(c for S, c in allc.items() if len(S) == k) for k in (3, 4, 5, 6, 7))
        hT = H_count(A)
        hK = H_count(D_blowup(A)[0])
        hD = H_count(D_skew(A)[0])
        hS = H_count(D_scblowup(A)[0])
        DATA.append(dict(n=n, idx=idx, A=A, spec=spec, hT=hT, hK=hK, hD=hD, hS=hS,
                         pf=pair_features(allc),
                         csz=tuple(sorted(len(c) for c in strong_components(A)))))
w(f"[setup {time.time()-t0:.1f}s] {len(DATA)} classes n=3..6")
w("")

# ---------------- S1: congruences ----------------

w("--- S1: congruence laws ---")
bad = [d for d in DATA if (d['hK'] - (2 * d['hT'] - 1)) % 8 != 0]
w(f"  H(K2) == 2 H(T) - 1 (mod 8): {len(DATA)-len(bad)}/{len(DATA)} hold"
  + ("" if not bad else "  COUNTEREXAMPLES: " + str([(d['n'], d['idx']) for d in bad[:5]])))
cen16 = {}
for d in DATA:
    cen16.setdefault((d['hK'] - (2 * d['hT'] - 1)) % 16, 0)
    cen16[(d['hK'] - (2 * d['hT'] - 1)) % 16] += 1
w(f"  (H(K2) - (2H(T)-1)) mod 16 census: {dict(sorted(cen16.items()))}")
cen32 = {}
for d in DATA:
    cen32.setdefault((d['hK'] - (2 * d['hT'] - 1)) % 32, 0)
    cen32[(d['hK'] - (2 * d['hT'] - 1)) % 32] += 1
w(f"  ... mod 32 census: {dict(sorted(cen32.items()))}")
for m in (4, 8):
    tab = {}
    for d in DATA:
        tab.setdefault((d['hT'] % m, d['hS'] % m), 0)
        tab[(d['hT'] % m, d['hS'] % m)] += 1
    w(f"  (H(T) mod {m}, H(SC) mod {m}) census: {dict(sorted(tab.items()))}")
# D vs 2H-1?
badD = sum(1 for d in DATA if (d['hD'] - (2 * d['hT'] - 1)) % 8 != 0)
w(f"  H(D) == 2H(T)-1 (mod 8)? holds in {len(DATA)-badD}/{len(DATA)} (expect: no law)")
w("")

# ---------------- S2: c7'(K2) law over all 74 bases ----------------

w("--- S2: c7'(K2) exact lift law (74 bases; n=6 doubles live on 12 vertices) ---")
t1 = time.time()
rows = []
for d in DATA:
    Ad = D_blowup(d['A'])[0]
    c7p = sum(ham_cycles_subset(Ad, S) for S in itertools.combinations(range(Ad.shape[0]), 7)) \
        if Ad.shape[0] >= 7 else 0
    d['c7p'] = c7p
    rows.append(d)
w(f"  [c7' computed for all {len(rows)} doubles in {time.time()-t1:.1f}s]")
featk = [(3, 3, 1), (3, 3, 2), (3, 4, 1), (3, 4, 2), (3, 4, 3)]
X = np.array([[d['pf'].get(k, 0) for k in featk] for d in rows], float)
Y = np.array([d['c7p'] - (128 * d['spec'][4] + 192 * d['spec'][3] + 80 * d['spec'][2]
                          + 8 * d['spec'][1]) for d in rows], float)
w(f"  residual = c7' - (128c7 + 192c6 + 80c5 + 8c4); regress on {featk}")
w(f"  pair-feature matrix rank = {np.linalg.matrix_rank(X)} (of {X.shape[1]})")
coef, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
pred = X @ coef
w(f"  lstsq coeffs = {np.round(coef, 6)}   exact={np.allclose(pred, Y, atol=1e-6)} "
  f"max|res|={np.max(np.abs(pred-Y)):.4g}")
ri = np.round(coef).astype(int)
predi = X @ ri
w(f"  integer-rounded {list(ri)}: exact={np.array_equal(predi.astype(int), Y.astype(int))}")
if not np.array_equal(predi.astype(int), Y.astype(int)):
    # try integer least squares by brute small search around lstsq? just report residual table
    off = predi - Y
    w(f"  integer version max|res|={np.max(np.abs(off)):.0f}; nonzero rows:")
    cnt = 0
    for d, o in zip(rows, off):
        if o != 0 and cnt < 10:
            w(f"    n={d['n']} idx={d['idx']} spec={d['spec']} pf33/34="
              f"{[d['pf'].get(k,0) for k in featk]} resid_target={int(Y[rows.index(d)])} off={int(o)}")
            cnt += 1
w("")

# ---------------- S3: full Omega(K2) IP on same-spectrum strong pairs ----------------

w("--- S3: full Omega(T[K2]) profile+IP for same-spectrum STRONG groups at n=6 ---")
groups = {}
for d in DATA:
    groups.setdefault(d['spec'], []).append(d)
strong_groups = [(k, [d for d in ds if len(d['csz']) == 1 and d['n'] == 6])
                 for k, ds in groups.items()]
strong_groups = [(k, ds) for k, ds in strong_groups if len(ds) >= 2]
w(f"  {len(strong_groups)} spectrum groups with >=2 STRONG n=6 members")
mismatch = 0
for k, ds in strong_groups:
    sigs = []
    for d in ds:
        t2 = time.time()
        Ad = D_blowup(d['A'])[0]
        prof, carriers = odd_profile_full(Ad)
        i1, i2, i3, i4 = ip_counts(carriers)
        ocf = 1 + 2 * i1 + 4 * i2 + 8 * i3 + 16 * i4
        sigs.append((tuple(prof.values()), (i1, i2, i3, i4)))
        w(f"    spec={k} idx={d['idx']:>2}: profile'={tuple(prof.values())} "
          f"IP'=(1,{i1},{i2},{i3},{i4}) I'(2)={ocf} {'==' if ocf==d['hK'] else '!='} H(K2)={d['hK']}"
          f"   [{time.time()-t2:.0f}s]")
    if len(set(sigs)) > 1:
        mismatch += 1
        w(f"    *** PROFILE/IP MISMATCH within spectrum group {k} ***")
w(f"  groups with identical full profile+IP across members: "
  f"{len(strong_groups)-mismatch}/{len(strong_groups)}")
w("")

# component-reduction validation on one non-strong n=6 class (idx 4, comps (1,5))
w("  component-reduction validation (n=6 idx4, comps (1,5)):")
d4 = [d for d in DATA if d['n'] == 6 and d['idx'] == 4][0]
Ad = D_blowup(d4['A'])[0]
prof12, car12 = odd_profile_full(Ad)
comp = max(strong_components(d4['A']), key=len)
sub = d4['A'][np.ix_(comp, comp)]
prof10, car10 = odd_profile_full(D_blowup(sub)[0])
w(f"    12-vertex double profile: {dict(prof12)}")
w(f"    10-vertex comp double   : {dict(prof10)}")
w(f"    match on common lengths + zeros beyond: "
  f"{all(prof12[k]==prof10.get(k,0) for k in prof12)}")
w("")

# ---------------- S4: polynomial functional attempts ----------------

w("--- S4: polynomial fit attempts for G(c3,c4,c5,c6) = H(K2) ---")
uniq = {}
for d in DATA:
    uniq[d['spec'][:4]] = d['hK']
sp = sorted(uniq.items())
w(f"  {len(sp)} distinct spectra")
pts = np.array([list(k) for k, v in sp], float)
vals = np.array([v for k, v in sp], float)
def monomials(pts, deg):
    cols, names = [], []
    for ex in itertools.product(range(deg + 1), repeat=4):
        if sum(ex) <= deg:
            cols.append(np.prod([pts[:, i] ** e for i, e in enumerate(ex)], axis=0))
            names.append(ex)
    return np.array(cols).T, names
for deg in (2, 3):
    Xm, names = monomials(pts, deg)
    coef, _, rank, _ = np.linalg.lstsq(Xm, vals, rcond=None)
    pred = Xm @ coef
    w(f"  degree {deg}: {Xm.shape[1]} monomials, rank={rank}, "
      f"exact={np.allclose(pred, vals, atol=1e-6)}, max|res|={np.max(np.abs(pred-vals)):.4g}")
# multiplicativity reminder: G((1,0,0,0))^2 = 45^2 = 2025 = G((2,0,0,0)) (disjoint union) => G not polynomial-bounded
w("  note: disjoint-union multiplicativity G(s1+s2)=G(s1)G(s2) (e.g. 45^2=2025) "
  "forbids any global polynomial law")
w("")

# ---------------- S5: factorizations ----------------

w("--- S5: factorizations ---")
w("  H(D(transitive_n)) n=3..7:")
for v in (13, 95, 1033, 15611, 313285):
    w(f"    {v} = {fstr(v)}")
w("  H(SC(transitive_n)) n=3..7:")
for v in (41, 629, 14937, 513669, 24104937):
    w(f"    {v} = {fstr(v)}")
w("  key H(K2) values (strong classes):")
for v in (45, 393, 2785, 3225, 6069, 8097, 10773, 15565, 27513, 19189, 36709, 64773,
          90117, 59373, 23361, 172881, 85185, 210945, 160449, 131373, 230457, 236433,
          121537, 110049, 188721, 118113, 356977, 386445, 421425, 411513):
    w(f"    {v} = {fstr(v)}")
w("")

# ---------------- S6: counterexample certificates ----------------

w("--- S6: exact counterexample certificates (upper-triangle adjacency bits) ---")
for n, idx in ((5, 4), (5, 6), (4, 1), (4, 3), (5, 1), (5, 3)):
    d = [x for x in DATA if x['n'] == n and x['idx'] == idx][0]
    w(f"  n={n} idx={idx}: bits={bits_of(d['A'])} spec={d['spec']} "
      f"H(T)={d['hT']} H(D)={d['hD']} H(K2)={d['hK']} H(SC)={d['hS']}")
w("  HYP 'IP determines H(K2)' refuted by (5,4) vs (5,6): same I=1+4x, hK 3225 vs 2785")
w("  HYP 'spectrum determines H(D)' refuted by (4,1) vs (4,3): same spec, hD 189 vs 333")
w("  HYP 'spectrum determines H(SC)' refuted by (5,1) vs (5,3): same spec, hS 14961 vs 14973")
w("")
w(f"=== done in {time.time()-t0:.1f}s ===")
OUT.close()

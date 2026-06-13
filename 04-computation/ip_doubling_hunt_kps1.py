#!/usr/bin/env python3
"""
ip_doubling_hunt_kps1.py — kind-pasteur-2026-06-09-S1, BRANCH A
H-formula hunt for the three doublings (HYP-2334, OPEN-Q-045 Q1).

Phases:
 P1  per iso class n=3..6: all-cycle spectrum (c3..c6), odd-cycle conflict-graph
     independence polynomial I(Omega(T),x) (plain + length-typed), H(T) (OCF check),
     H(D(T)), H(T[K2]), H(SC blowup), self-complementarity, strong components.
 P2  determination tests: does {plain IP / typed IP / full cycle spectrum /
     spectrum+pair-intersection features} determine {H(D), H(K2), H(SC)}?
     Exact counterexamples reported.
 P2b STRONG-COMPONENT PRODUCT LAW: H(T[K2]) ?= prod over nontrivial strong
     components C of H(C[K2]).  (Candidate theorem explaining cross-n collisions.)
 P3  Omega structure of the doubles (bases n=3..5, all three doublings):
     odd-cycle profiles of the doubles, OCF check on doubles, exact lift laws
     c3'(K2)=8 c3, c5'(K2)=32 c5+32 c4+6 c3, regression for c7'(K2) incl.
     figure-eight (cycle-pair) features; same regressions for D and SCblow.
 P4  H(K2) functional hunt on strong classes (raw table + linear fits).
 P5  transitive towers n=3..7: H(D), H(K2), H(SC) sequences (for OEIS).
 P6  mod-4/8/16 patterns of H(D), H(K2), H(SC) vs H(T); op-asymmetry census of D.

Output: 05-knowledge/results/ip_doubling_hunt_kps1.out
"""
import sys, itertools, time
sys.path.insert(0, '04-computation')
import numpy as np
from skew_doubling_core_kps1 import (all_tournaments, M_of, A_of, scores, H_count,
                                     D_skew, D_blowup, D_scblowup, is_iso, canon)

OUT = open('05-knowledge/results/ip_doubling_hunt_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s)

# ---------------- helpers ----------------

def iso_classes_fast(n):
    """Orbit-sweep iso classes (much faster than canon-everything at n=6)."""
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
    """# directed Hamiltonian cycles of the sub-tournament induced on S."""
    k = len(S)
    if k == 3:
        a, b, c = S
        return 1 if (A[a, b] and A[b, c] and A[c, a]) or (A[a, c] and A[c, b] and A[b, a]) else 0
    adj = [[int(A[S[i], S[j]]) for j in range(k)] for i in range(k)]
    dp = [[0] * k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        row = dp[mask]
        for last in range(k):
            cnt = row[last]
            if not cnt:
                continue
            al = adj[last]
            for nxt in range(1, k):
                if al[nxt] and not ((mask >> nxt) & 1):
                    dp[mask | (1 << nxt)][nxt] += cnt
    full = (1 << k) - 1
    return sum(dp[full][l] for l in range(1, k) if adj[l][0])

def cycle_sets(A, lengths):
    """dict vertex-tuple -> # distinct directed cycles on exactly that set."""
    n = A.shape[0]
    out = {}
    for k in lengths:
        if k > n or k < 3:
            continue
        for S in itertools.combinations(range(n), k):
            c = ham_cycles_subset(A, S)
            if c:
                out[S] = c
    return out

def omega_ip(cyc):
    """Independence polynomial of the conflict graph of the given cycles.
    cyc: dict vertex-tuple -> multiplicity (distinct cycles on that set).
    Two cycles conflict iff they share a vertex (same-set cycles always conflict).
    Returns (plain coefficient list, typed dict {sorted length tuple: weight})."""
    items = [(sum(1 << v for v in S), c, len(S)) for S, c in cyc.items()]
    plain = {0: 1}
    typed = {(): 1}
    def rec(i, used, weight, lens):
        for j in range(i, len(items)):
            m, c, L = items[j]
            if m & used:
                continue
            w2 = weight * c
            l2 = tuple(sorted(lens + (L,)))
            typed[l2] = typed.get(l2, 0) + w2
            plain[len(l2)] = plain.get(len(l2), 0) + w2
            rec(j + 1, used | m, w2, l2)
    rec(0, 0, 1, ())
    co = [plain.get(k, 0) for k in range(max(plain) + 1)]
    return co, typed

def ip_eval(co, x):
    return sum(c * x ** k for k, c in enumerate(co))

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

def pair_features(allc):
    """dict (l1<=l2, shared) -> weighted count of unordered cycle pairs."""
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
            key = (a, b, s)
            pf[key] = pf.get(key, 0) + c1 * c2
    return pf

def typed_key(typed):
    return tuple(sorted(typed.items()))

def bits_of(A):
    n = A.shape[0]
    return ''.join(str(A[i, j]) for i in range(n) for j in range(i + 1, n))

def transitive(n):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = 1
    return A

# ---------------- P1: master table ----------------

t0 = time.time()
w("=== ip_doubling_hunt_kps1 — BRANCH A: H-formula hunt for the doublings ===")
w("")
w("--- P1: master table, iso classes n=3..6 ---")
w("spec = (c3,c4,c5,c6) all directed cycles; IP = I(Omega(T),x) coeffs (odd cycles only)")
w("typed = length-typed IP monomials; comps = strong component sizes")
w("")
hdr = (f"{'n':>2} {'idx':>3} {'scores':>16} {'comps':>10} {'spec':>16} "
      f"{'IP':>14} {'H(T)':>5} {'H(D)':>7} {'H(K2)':>7} {'H(SC)':>7} {'SC?':>4}")
w(hdr)

DATA = []
K2_cache = {}   # canonical bytes of strong comp -> H(C[K2])

def hK2_of(sub):
    key = canon(sub).tobytes()
    if key not in K2_cache:
        K2_cache[key] = H_count(D_blowup(sub)[0])
    return K2_cache[key]

ocf_fail = 0
for n in (3, 4, 5, 6):
    reps = iso_classes_fast(n)
    for idx, A in enumerate(reps):
        allc = cycle_sets(A, range(3, n + 1))
        oddc = {S: c for S, c in allc.items() if len(S) % 2 == 1}
        spec = tuple(sum(c for S, c in allc.items() if len(S) == k) for k in (3, 4, 5, 6))
        ip, typed = omega_ip(oddc)
        hT = H_count(A)
        if ip_eval(ip, 2) != hT:
            ocf_fail += 1
            w(f"  *** OCF FAILURE n={n} idx={idx}: I(2)={ip_eval(ip,2)} H={hT}")
        hD = H_count(D_skew(A)[0])
        hK = H_count(D_blowup(A)[0])
        hS = H_count(D_scblowup(A)[0])
        comps = strong_components(A)
        csz = tuple(sorted(len(c) for c in comps))
        sc = is_iso(A, A.T)
        prodK = 1
        for c in comps:
            if len(c) >= 2:
                prodK *= hK2_of(A[np.ix_(c, c)])
        pf = pair_features(allc)
        DATA.append(dict(n=n, idx=idx, A=A, bits=bits_of(A), scores=scores(A),
                         spec=spec, ip=tuple(ip), typed=typed_key(typed),
                         hT=hT, hD=hD, hK=hK, hS=hS, sc=sc, csz=csz,
                         prodK=prodK, pf=pf,
                         hDop=H_count(D_skew(A.T)[0])))
        w(f"{n:>2} {idx:>3} {str(scores(A)):>16} {str(csz):>10} {str(spec):>16} "
          f"{str(tuple(ip)):>14} {hT:>5} {hD:>7} {hK:>7} {hS:>7} {('Y' if sc else '-'):>4}")
w("")
w(f"OCF check I(Omega(T),2)==H(T): {'ALL PASS' if ocf_fail==0 else str(ocf_fail)+' FAILURES'}"
  f" over {len(DATA)} classes   [{time.time()-t0:.1f}s]")
w("")

# ---------------- P2: determination tests ----------------

w("--- P2: does invariant X determine H(double)? (exact counterexamples) ---")

def determination_test(name, keyf, include_n_variants=True):
    for valname in ('hD', 'hK', 'hS'):
        groups = {}
        for d in DATA:
            groups.setdefault(keyf(d), []).append(d)
        bad = [(k, ds) for k, ds in groups.items()
               if len({dd[valname] for dd in ds}) > 1]
        multi = sum(1 for k, ds in groups.items() if len(ds) > 1)
        w(f"  [{name}] -> {valname}: {len(groups)} groups, {multi} with >=2 members, "
          f"{len(bad)} BROKEN  => {'DETERMINED' if not bad else 'NOT determined'}")
        for k, ds in bad[:4]:
            vals = sorted({dd[valname] for dd in ds})
            mem = '; '.join(f"n={dd['n']} idx={dd['idx']} sc={dd['scores']} "
                            f"{valname}={dd[valname]}" for dd in ds[:6])
            w(f"      CE key={k}: values {vals} | {mem}")
        if len(bad) > 4:
            w(f"      ... and {len(bad)-4} more broken groups")

determination_test("plain IP, cross-n", lambda d: d['ip'])
w("")
determination_test("typed IP, cross-n", lambda d: d['typed'])
w("")
determination_test("full cycle spectrum (c3,c4,c5,c6), cross-n", lambda d: d['spec'])
w("")
determination_test("(spectrum, n)", lambda d: (d['spec'], d['n']))
w("")
determination_test("(spectrum, pair-features), cross-n",
                   lambda d: (d['spec'], tuple(sorted(d['pf'].items()))))
w("")

# ---------------- P2b: strong-component product law for H(K2) ----------------

w("--- P2b: STRONG-COMPONENT PRODUCT LAW: H(T[K2]) ?= prod_C H(C[K2]) ---")
badP = [d for d in DATA if d['prodK'] != d['hK']]
w(f"  holds in {len(DATA)-len(badP)}/{len(DATA)} classes (n=3..6)")
for d in badP[:10]:
    w(f"  FAIL n={d['n']} idx={d['idx']} scores={d['scores']} comps={d['csz']} "
      f"H(K2)={d['hK']} prod={d['prodK']}")
# sanity: same law for plain H(T) (known) and explicit failure check for D and SC
badH = badD = badS = 0
for d in DATA:
    A = d['A']; comps = strong_components(A)
    pH = pD = pS = 1
    for c in comps:
        if len(c) >= 2:
            sub = A[np.ix_(c, c)]
            pH *= H_count(sub)
            pD *= H_count(D_skew(sub)[0])
            pS *= H_count(D_scblowup(sub)[0])
    badH += (pH != d['hT'])
    badD += (pD != d['hD'])
    badS += (pS != d['hS'])
w(f"  control: product law for H(T) itself: {len(DATA)-badH}/{len(DATA)} hold (known thm)")
w(f"  product law for H(D): {len(DATA)-badD}/{len(DATA)} hold;  for H(SC): {len(DATA)-badS}/{len(DATA)} hold")
w("")

# ---------------- P3: Omega structure of the doubles (bases n=3..5) ----------------

w("--- P3: odd-cycle profiles of the doubles (bases n=3..5; doubles on 6,8,10) ---")
w("profile' = (c3',c5',c7',c9') of the double;  OCF re-checked on every double")
w("")
P3 = []
ocf2_fail = 0
for d in DATA:
    n = d['n']
    if n > 5:
        continue
    A = d['A']
    row = dict(base=d)
    for name, fn, h in (('D ', D_skew, d['hD']), ('K2', D_blowup, d['hK']),
                        ('SC', D_scblowup, d['hS'])):
        Ad = fn(A)[0]
        cyc = cycle_sets(Ad, range(3, 2 * n, 2))   # odd lengths 3..2n-1
        prof = tuple(sum(c for S, c in cyc.items() if len(S) == k) for k in (3, 5, 7, 9))
        ipd, typd = omega_ip(cyc)
        ok = (ip_eval(ipd, 2) == h)
        if not ok:
            ocf2_fail += 1
        row[name.strip()] = prof
        w(f"  n={n} idx={d['idx']:>2} {name}: profile'={prof}  IP'={tuple(ipd)}  "
          f"I'(2)={ip_eval(ipd,2)} {'==' if ok else '!='} H={h}")
    P3.append(row)
w(f"  OCF on doubles: {'ALL PASS' if ocf2_fail==0 else str(ocf2_fail)+' FAILURES'}")
w("")

# exact lift laws for K2
w("  K2 lift-law checks (predictions from twin-decorated closed-walk lifting):")
ok3 = ok5 = tot = 0
for row in P3:
    d = row['base']; c3, c4, c5, c6 = d['spec']
    p = row['K2']
    tot += 1
    ok3 += (p[0] == 8 * c3)
    ok5 += (p[1] == 32 * c5 + 32 * c4 + 6 * c3)
w(f"    c3'(K2) == 8*c3                : {ok3}/{tot}")
w(f"    c5'(K2) == 32*c5 + 32*c4 + 6*c3: {ok5}/{tot}")
# c7'(K2): known single-cycle terms 128c7+192c6+80c5+8c4, regress residual on pair features
w("    c7'(K2) = 128c7+192c6+80c5+8c4 + figure-eight terms; residual regression:")
feat_keys = [(3, 3, 1), (3, 3, 2), (3, 4, 1), (3, 4, 2), (3, 4, 3)]
X, Y = [], []
for row in P3:
    d = row['base']; c3, c4, c5, c6 = d['spec']
    resid = row['K2'][2] - (192 * c6 + 80 * c5 + 8 * c4)
    pf = d['pf']
    X.append([pf.get(k, 0) for k in feat_keys])
    Y.append(resid)
X = np.array(X, dtype=float); Y = np.array(Y, dtype=float)
coef, res, rank, _ = np.linalg.lstsq(X, Y, rcond=None)
pred = X @ coef
exact = np.allclose(pred, Y, atol=1e-6)
w(f"      features {feat_keys}")
w(f"      coeffs   {np.round(coef, 4)}   exact fit: {exact}  "
  f"max|resid|={np.max(np.abs(pred-Y)):.4g}")
if not exact:
    w("      (non-exact => 7-cycle count needs arc-level or finer features)")
    for row, p, y in zip(P3, pred, Y):
        d = row['base']
        w(f"        n={d['n']} idx={d['idx']} resid_target={int(y)} fitted={p:.2f}")
w("")

# c3' laws for D and SC doubles: regress on [1, n, c3, c4, c5]
for nm in ('D', 'SC'):
    w(f"  {nm} double: c3' regression on [1, n, c3, c4, c5]:")
    X, Y = [], []
    for row in P3:
        d = row['base']; c3, c4, c5, c6 = d['spec']
        X.append([1, d['n'], c3, c4, c5])
        Y.append(row[nm][0])
    X = np.array(X, float); Y = np.array(Y, float)
    coef, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    pred = X @ coef
    w(f"      coeffs [1,n,c3,c4,c5] = {np.round(coef,4)}  exact={np.allclose(pred,Y,atol=1e-6)}"
      f"  max|res|={np.max(np.abs(pred-Y)):.4g}")
w("")

# n=6 bases: doubles' c3',c5' only (cheap), to stress-test K2 lift laws at n=12
w("  n=6 bases: K2/D/SC doubles' (c3',c5') and K2 lift-law check at 12 vertices:")
ok3 = ok5 = tot = 0
D6 = []
for d in DATA:
    if d['n'] != 6:
        continue
    A = d['A']; c3, c4, c5, c6 = d['spec']
    rec = {}
    for name, fn in (('K2', D_blowup), ('D', D_skew), ('SC', D_scblowup)):
        Ad = fn(A)[0]
        cy3 = cycle_sets(Ad, [3]); cy5 = cycle_sets(Ad, [5])
        rec[name] = (sum(cy3.values()), sum(cy5.values()))
    D6.append((d, rec))
    tot += 1
    ok3 += (rec['K2'][0] == 8 * c3)
    ok5 += (rec['K2'][1] == 32 * c5 + 32 * c4 + 6 * c3)
w(f"    c3'(K2)==8c3: {ok3}/{tot}    c5'(K2)==32c5+32c4+6c3: {ok5}/{tot}")
# D/SC c3',c5' regressions with n=3..6 data combined
for nm in ('D', 'SC'):
    for ci, lab in ((0, "c3'"), (1, "c5'")):
        X, Y = [], []
        for row in P3:
            d = row['base']; c3, c4, c5, c6 = d['spec']
            X.append([1, d['n'], c3, c4, c5, c6]); Y.append(row[nm][ci])
        for d, rec in D6:
            c3, c4, c5, c6 = d['spec']
            X.append([1, d['n'], c3, c4, c5, c6]); Y.append(rec[nm][ci])
        X = np.array(X, float); Y = np.array(Y, float)
        coef, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
        pred = X @ coef
        w(f"    {nm} {lab} ~ [1,n,c3,c4,c5,c6]: coeffs={np.round(coef,4)} "
          f"exact={np.allclose(pred,Y,atol=1e-6)} max|res|={np.max(np.abs(pred-Y)):.3g}")
w("")

# K2 c5' law re-verified on combined n=3..6 (printed above); also c7'/c9' raw for record
w("  raw double profiles stored in table above (P3) for c7', c9' future laws")
w("")

# ---------------- P4: H(K2) functional hunt on strong classes ----------------

w("--- P4: H(K2) on STRONG classes (component law reduces everything to these) ---")
w(f"{'n':>2} {'idx':>3} {'spec':>16} {'pairs(3,3,*)':>14} {'H(T)':>5} {'H(K2)':>7}")
strong = [d for d in DATA if len(d['csz']) == 1]
for d in strong:
    p33 = tuple(d['pf'].get((3, 3, s), 0) for s in (1, 2, 3))
    w(f"{d['n']:>2} {d['idx']:>3} {str(d['spec']):>16} {str(p33):>14} {d['hT']:>5} {d['hK']:>7}")
w("")
# linear fits of H(K2) over spectra+pair features on strong classes
featf = [
    ("1", lambda d: 1), ("c3", lambda d: d['spec'][0]), ("c4", lambda d: d['spec'][1]),
    ("c5", lambda d: d['spec'][2]), ("c6", lambda d: d['spec'][3]),
    ("c3^2", lambda d: d['spec'][0] ** 2), ("c3c4", lambda d: d['spec'][0] * d['spec'][1]),
    ("c4^2", lambda d: d['spec'][1] ** 2),
    ("p331", lambda d: d['pf'].get((3, 3, 1), 0)), ("p332", lambda d: d['pf'].get((3, 3, 2), 0)),
    ("p341", lambda d: d['pf'].get((3, 4, 1), 0)), ("p342", lambda d: d['pf'].get((3, 4, 2), 0)),
    ("p343", lambda d: d['pf'].get((3, 4, 3), 0)),
    ("HT", lambda d: d['hT']), ("HT^2", lambda d: d['hT'] ** 2),
]
X = np.array([[f(d) for _, f in featf] for d in strong], float)
Y = np.array([d['hK'] for d in strong], float)
coef, _, rank, _ = np.linalg.lstsq(X, Y, rcond=None)
pred = X @ coef
w(f"  lstsq H(K2) ~ [{', '.join(nm for nm,_ in featf)}] on {len(strong)} strong classes:")
w(f"  rank={rank}, coeffs={np.round(coef,3)}")
w(f"  exact={np.allclose(pred,Y,atol=1e-6)}  max|res|={np.max(np.abs(pred-Y)):.4g}")
w("")

# ---------------- P5: transitive towers ----------------

w("--- P5: transitive tournament towers (sequences for OEIS) ---")
seqD, seqK, seqS = [], [], []
for n in range(3, 8):
    A = transitive(n)
    t1 = time.time()
    hD = H_count(D_skew(A)[0]); hK = H_count(D_blowup(A)[0]); hS = H_count(D_scblowup(A)[0])
    seqD.append(hD); seqK.append(hK); seqS.append(hS)
    w(f"  n={n} (double on {2*n}): H(D)={hD}  H(K2)={hK}  H(SC)={hS}   [{time.time()-t1:.1f}s]")
w(f"  H(D(transitive_n)),  n=3..7: {seqD}")
w(f"  H(K2(transitive_n)), n=3..7: {seqK}")
w(f"  H(SC(transitive_n)), n=3..7: {seqS}")
# also the C3-anchored K2 ladder: H(K2) of single-3-cycle classes is 45 at n=3,4,5 — verify at n=6
ladder = sorted({d['hK'] for d in DATA if d['ip'] == (1, 1)})
w(f"  H(K2) over ALL classes with IP=1+x (one odd cycle), n=3..6: {ladder}")
w("")

# ---------------- P6: mod patterns + op asymmetry ----------------

w("--- P6: H(D) mod patterns vs H(T), and op-asymmetry of D ---")
allodd = all(d['hD'] % 2 == 1 and d['hK'] % 2 == 1 and d['hS'] % 2 == 1 for d in DATA)
w(f"  all H(D),H(K2),H(SC) odd (Redei): {allodd}")
for m in (4, 8, 16):
    tab = {}
    for d in DATA:
        tab.setdefault((d['hT'] % m, d['hD'] % m), 0)
        tab[(d['hT'] % m, d['hD'] % m)] += 1
    w(f"  (H(T) mod {m}, H(D) mod {m}) census: {dict(sorted(tab.items()))}")
for m in (4, 8):
    tabK = {}
    for d in DATA:
        tabK.setdefault((d['hT'] % m, d['hK'] % m), 0)
        tabK[(d['hT'] % m, d['hK'] % m)] += 1
    w(f"  (H(T) mod {m}, H(K2) mod {m}) census: {dict(sorted(tabK.items()))}")
asym = [d for d in DATA if d['hD'] != d['hDop']]
w(f"  classes with H(D(T)) != H(D(T^op)): {len(asym)}/{len(DATA)}")
for d in asym[:8]:
    w(f"    n={d['n']} idx={d['idx']} scores={d['scores']}: H(D(T))={d['hD']}  H(D(T^op))={d['hDop']}")
if len(asym) > 8:
    w(f"    ... and {len(asym)-8} more")
w("")
w(f"=== done in {time.time()-t0:.1f}s ===")
OUT.close()

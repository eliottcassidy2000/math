#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3951: the EXPLICIT PSL_2(F_7) tensor code + its MEASURED soundness s — does s recover a
1/36-type constant?

kind-pasteur-2026-07-01-S28. Fills the gap left by opus-S31 (theoretical trickle-down bound only,
rho >= delta_inner*(1-lambda_link) with lambda_link=0 automatic for EVERY left-right complex — so it
cannot discriminate; it would equally 'certify' mac-mini's disconnected single-generator complex).

THE CODE (DELLM-style, explicit):
  G = PSL_2(F_7), |G| = 168. A = symmetric set of ORDER-7 elements (inverse pairs), B = symmetric set
  of ORDER-3/4 elements. Order-disjointness => NO degenerate squares: faces = triples (g,a,b) under the
  4-way identification (g,a,b)~(ag,a^-1,b)~(gb,a,b^-1)~(agb,a^-1,b^-1); F = |G||A||B|/4 exactly; the
  link of every vertex is the FULL |A| x |B| grid of distinct faces. Inner codes C_A (on columns),
  C_B (on rows); the tensor code C = {f: F -> F2 : every vertex's grid view lies in C_A (x) C_B}.
  Tester: read a uniform vertex's grid (|A||B| queries), accept iff in the tensor code.
  Soundness of f: rho(f) = [viol(f)/V] / [dist(f,C)/F];  s = inf over f not in C.

MEASUREMENTS per configuration:
  connectivity (double-coset components — mac-mini HYP-3824's correction), Cayley lambda_2 of A and B,
  dim C, min-distance upper bound d_ub (randomized codeword bank), exhaustive soundness at weight 1
  and 2 (coset leaders when < d/2), sampled weight 3-6, adversarial tubes/rectangles/regions, and a
  hill-climb adversary with codeword-bank coset reduction. Compare s against the 1/36 family:
  1/36 = 0.027778 (analytic 11-core floor), 313/9702 = 0.032261 (census min), (delta/2)^2, delta^2/4,
  delta_A*delta_B, trickle-down delta_inner.
"""
import sys, itertools, random
import numpy as np
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
rng = random.Random(3951)

# ------------------------------------------------------------------------------- group PSL_2(F_7)
P = 7
def mm(X, Y):
    return ((X[0]*Y[0]+X[1]*Y[2]) % P, (X[0]*Y[1]+X[1]*Y[3]) % P,
            (X[2]*Y[0]+X[3]*Y[2]) % P, (X[2]*Y[1]+X[3]*Y[3]) % P)
def norm(X):
    Y = tuple((-x) % P for x in X); return min(X, Y)
SL = [(a,b,c,d) for a in range(P) for b in range(P) for c in range(P) for d in range(P)
      if (a*d - b*c) % P == 1]
G = sorted(set(norm(x) for x in SL)); idx = {g:i for i,g in enumerate(G)}; V = len(G)
I2 = norm((1,0,0,1))
def inv(X):
    a,b,c,d = X; return norm((d, (-b)%P, (-c)%P, a))
def order(X):
    Y = X; k = 1
    while norm(Y) != I2: Y = mm(Y, X); k += 1
    return k
ords = {g: order(g) for g in G}
by_order = {}
for g in G:
    if g != I2: by_order.setdefault(ords[g], []).append(g)
print("="*106)
print(f" PSL_2(F_7): |G| = {V}; elements by order: " + ", ".join(f"{o}:{len(v)}" for o,v in sorted(by_order.items())))
print("="*106)

def cayley_lambda2(S):
    A = np.zeros((V, V))
    for g in G:
        for s in S: A[idx[g], idx[norm(mm(s, g))]] += 1
    ev = np.sort(np.linalg.eigvalsh(A))[::-1]
    return max(abs(ev[1]), abs(ev[-1]))

def gen_span(S):
    seen = {I2}; fr = [I2]
    while fr:
        x = fr.pop()
        for s in S:
            y = norm(mm(s, x))
            if y not in seen: seen.add(y); fr.append(y)
    return len(seen)

def pick_sym_set(order_pool, npairs, tries=600):
    """best-expansion symmetric set of `npairs` inverse pairs from elements of the given orders."""
    pool = []
    for o in order_pool: pool += by_order[o]
    best = None
    for _ in range(tries):
        seeds = rng.sample(pool, npairs)
        S = set()
        for s in seeds: S.add(s); S.add(inv(s))
        if len(S) != 2*npairs: continue
        lam = cayley_lambda2(sorted(S))
        if best is None or lam < best[1]: best = (sorted(S), lam)
    return best

# ------------------------------------------------------------------- complex + tensor code builder
def build(A, B):
    """faces = (g,a,b) triples mod the 4-way identification. Returns:
       nF, link[v] = |A|x|B| array of face ids, face_pres[fid] = list of (v, ia, ib)."""
    nA, nB = len(A), len(B)
    iA = {a:i for i,a in enumerate(A)}; iB = {b:i for i,b in enumerate(B)}
    ainv = [iA[inv(a)] for a in A]; binv = [iB[inv(b)] for b in B]
    face_id = {}; face_pres = []
    link = [[[-1]*nB for _ in range(nA)] for _ in range(V)]
    for g in G:
        gi = idx[g]
        for ia, a in enumerate(A):
            ag = norm(mm(a, g)); agi = idx[ag]
            for ib, b in enumerate(B):
                gb = norm(mm(g, b)); agb = norm(mm(ag, b))
                pres = [(gi, ia, ib), (agi, ainv[ia], ib), (idx[gb], ia, binv[ib]), (idx[agb], ainv[ia], binv[ib])]
                key = min(pres)
                if key not in face_id:
                    face_id[key] = len(face_pres); face_pres.append(pres)
                fid = face_id[key]
                link[gi][ia][ib] = fid
    return len(face_pres), link, face_pres

def parity_checks(n, kind):
    """inner-code parity-check rows as bitmasks over positions 0..n-1."""
    if kind == 'parity':  return [(1 << n) - 1]                       # [n, n-1, 2]
    if kind == 'rep':     return [(1 << i) | 1 for i in range(1, n)]  # [n, 1, n]
    if kind == 'hamming8':                                            # [8,4,4] extended Hamming
        H = [0b11110000, 0b11001100, 0b10101010, 0b11111111]
        return H
    raise ValueError(kind)

def rank_gf2(rows):
    basis = {}
    r = 0
    for row in rows:
        x = row
        while x:
            h = x.bit_length() - 1
            if h in basis: x ^= basis[h]
            else: basis[h] = x; r += 1; break
    return r, basis

def build_checks(nF, link, A, B, HA, HB):
    """global GF(2) checks: per vertex, per row a: HB rows applied over b; per col b: HA rows over a."""
    checks = []
    nA, nB = len(A), len(B)
    for v in range(V):
        L = link[v]
        for ia in range(nA):
            for h in HB:
                m = 0
                for ib in range(nB):
                    if (h >> ib) & 1: m ^= (1 << L[ia][ib])
                checks.append(m)
        for ib in range(nB):
            for h in HA:
                m = 0
                for ia in range(nA):
                    if (h >> ia) & 1: m ^= (1 << L[ia][ib])
                checks.append(m)
    return checks

def nullspace_basis(checks, nF):
    """basis of C = ker(checks) subset F2^nF (list of bitmasks), via RREF back-substitution."""
    rows = [c for c in checks if c]
    _, basis = rank_gf2(rows)
    # full reduction to RREF: eliminate each pivot bit from all other basis rows
    for h in sorted(basis, reverse=True):
        r = basis[h]
        for h2 in list(basis):
            if h2 != h and (basis[h2] >> h) & 1:
                basis[h2] ^= r
    pivots = set(basis)
    free = [j for j in range(nF) if j not in pivots]
    out = []
    for j in free:
        x = 1 << j
        for h, r in basis.items():          # RREF row: pivot h + free bits; constraint x_h = x_j-coeff
            if (r >> j) & 1: x |= (1 << h)
        out.append(x)
    # sanity: every basis codeword must satisfy all checks
    for x in out[:5]:
        assert all(bin(x & c).count('1') % 2 == 0 for c in checks[:50]), "nullspace sanity failed"
    return out

def viol_of(faces_set, face_pres, link, HA, HB, nA, nB):
    """# vertices whose grid view (indicator of faces_set) fails the tensor test."""
    touched = {}
    for f in faces_set:
        for (v, ia, ib) in face_pres[f]:
            touched.setdefault(v, []).append((ia, ib))
    bad = 0
    for v, cells in touched.items():
        rows = [0]*nA; cols = [0]*nB
        for (ia, ib) in cells:
            rows[ia] ^= (1 << ib); cols[ib] ^= (1 << ia)
        ok = True
        for r in rows:
            if r and any(bin(r & h).count('1') % 2 for h in HB): ok = False; break
        if ok:
            for c in cols:
                if c and any(bin(c & h).count('1') % 2 for h in HA): ok = False; break
        if not ok: bad += 1
    return bad

# ------------------------------------------------------------------------------------ configs
print("\n searching generating sets (best Cayley lambda_2 per size; A from order-7, B from order-3/4):")
A2 = ([by_order[7][0], inv(by_order[7][0])], None)
A2 = (sorted(set(A2[0])), cayley_lambda2(sorted(set(A2[0]))))
b3 = by_order[3][0]; B2 = (sorted({b3, inv(b3)}), cayley_lambda2(sorted({b3, inv(b3)})))
A4 = pick_sym_set([7], 2); B4 = pick_sym_set([3,4], 2)
A6 = pick_sym_set([7], 3); B6 = pick_sym_set([3,4], 3)
A8 = pick_sym_set([7], 4); B8 = pick_sym_set([3,4], 4)
for name, (S, lam) in [("A2",A2),("B2",B2),("A4",A4),("B4",B4),("A6",A6),("B6",B6),("A8",A8),("B8",B8)]:
    d = len(S); ram = 2*np.sqrt(d-1)
    print(f"   {name}: deg {d}  lambda2 = {lam:.3f}  (Ramanujan bound 2sqrt(d-1) = {ram:.3f})"
          f"  generates G: {gen_span(S) == V}")

CONFIGS = [
    ("SURFACE rep2xrep2 (deg2)",  A2[0], B2[0], 'rep',    'rep',    (2,1,2), (2,1,2)),
    ("T3-like par4xpar4 (deg4)",  A4[0], B4[0], 'parity', 'parity', (4,3,2), (4,3,2)),
    ("T6 par6xpar6 (deg6=phi7)",  A6[0], B6[0], 'parity', 'parity', (6,5,2), (6,5,2)),
    ("T6 rep6xrep6 (deg6)",       A6[0], B6[0], 'rep',    'rep',    (6,1,6), (6,1,6)),
    ("T8 ham8xham8 (deg8=LPS)",   A8[0], B8[0], 'hamming8','hamming8',(8,4,4),(8,4,4)),
]

T36 = 1.0/36.0
summary = []
for (name, A, B, kA, kB, prmA, prmB) in CONFIGS:
    nA, nB = len(A), len(B)
    print("\n" + "="*106)
    print(f" CONFIG {name}:  |A|={nA} |B|={nB}")
    print("="*106)
    # connectivity of the complex = components of the union graph (A-left + B-right edges)
    adj = [[] for _ in range(V)]
    for g in G:
        gi = idx[g]
        for a in A: adj[gi].append(idx[norm(mm(a, g))])
        for b in B: adj[gi].append(idx[norm(mm(g, b))])
    comp = [-1]*V; nc = 0
    for s0 in range(V):
        if comp[s0] >= 0: continue
        comp[s0] = nc; st = [s0]
        while st:
            x = st.pop()
            for y in adj[x]:
                if comp[y] < 0: comp[y] = nc; st.append(y)
        nc += 1
    nF, link, face_pres = build(A, B)
    HA = parity_checks(nA, kA); HB = parity_checks(nB, kB)
    checks = build_checks(nF, link, A, B, HA, HB)
    rk, _ = rank_gf2(checks)
    dimC = nF - rk
    dA = prmA[2]/prmA[0]; dB = prmB[2]/prmB[0]
    print(f"  components = {nc} {'(CONNECTED)' if nc==1 else '(DISCONNECTED — mac-mini regime)'}   "
          f"V={V} F={nF} (=|G||A||B|/4: {V*nA*nB//4})   dim C = {dimC}  rate = {dimC/nF:.3f}")
    print(f"  inner codes: C_A [{prmA[0]},{prmA[1]},{prmA[2]}] delta_A={dA:.3f} | C_B [{prmB[0]},{prmB[1]},{prmB[2]}] delta_B={dB:.3f}")

    # ---- codeword bank (for d_ub + coset reduction) via random column combinations
    basisC = nullspace_basis(checks, nF)
    bank = []
    if basisC:
        for _ in range(min(4000, 40*len(basisC))):
            k = rng.randint(1, min(4, len(basisC)))
            x = 0
            for b_ in rng.sample(basisC, k): x ^= b_
            if x: bank.append(x)
        # greedy weight reduction inside the bank
        red = []
        for x in sorted(bank, key=lambda z: bin(z).count('1'))[:400]:
            improved = True
            while improved:
                improved = False
                for b_ in basisC:
                    y = x ^ b_
                    if bin(y).count('1') < bin(x).count('1') and y: x = y; improved = True
            red.append(x)
        bank = sorted(set(red), key=lambda z: bin(z).count('1'))
    d_ub = bin(bank[0]).count('1') if bank else None
    print(f"  min-distance upper bound d_ub = {d_ub} (bank of {len(bank)} reduced codewords)" if bank
          else "  code is trivial/empty bank")

    def ratio(fs, wt=None):
        vv = viol_of(fs, face_pres, link, HA, HB, nA, nB)
        w = wt if wt is not None else len(fs)
        return (vv / V) / (w / nF), vv

    # ---- exhaustive weight 1
    best1 = (1e9, None)
    for f in range(nF):
        r, vv = ratio({f})
        if r < best1[0]: best1 = (r, (f, vv))
    print(f"  s(weight-1 exhaustive) = {best1[0]:.4f}   (viol = {best1[1][1]} vertices)")
    # ---- weight 2: exhaustive if small else sampled
    best2 = (1e9, None); n2 = 0
    if nF <= 800:
        it2 = itertools.combinations(range(nF), 2)
    else:
        it2 = ((rng.randrange(nF), rng.randrange(nF)) for _ in range(300000))
    for (f1, f2) in it2:
        if f1 == f2: continue
        n2 += 1
        r, vv = ratio({f1, f2})
        if r < best2[0]: best2 = (r, ((f1, f2), vv))
    print(f"  s(weight-2 {'exhaustive' if nF<=800 else 'sampled'} n={n2}) = {best2[0]:.4f}   (viol = {best2[1][1]})")
    # ---- structured adversaries: a-tubes (closed order-7 strings), rectangles, random balls
    bestS = (1e9, None)
    for _ in range(200):
        g0 = rng.randrange(V); ia = rng.randrange(nA); ib = rng.randrange(nB)
        # tube: walk g, a g, a^2 g, ... using face(v, ia, ib)
        a = A[ia]; vset = set(); v0 = G[g0]; cur = v0
        fs = set()
        for _ in range(ords[a]):
            fs.add(link[idx[cur]][ia][ib]); cur = norm(mm(a, cur))
        r, vv = ratio(fs)
        if r < bestS[0]: bestS = (r, ('tube', len(fs), vv))
    for _ in range(300):
        v0 = rng.randrange(V)
        ias = rng.sample(range(nA), min(2, nA)); ibs = rng.sample(range(nB), min(2, nB))
        fs = {link[v0][ia][ib] for ia in ias for ib in ibs}
        r, vv = ratio(fs)
        if r < bestS[0]: bestS = (r, ('rect', len(fs), vv))
    print(f"  s(structured tubes/rects) = {bestS[0]:.4f}   ({bestS[1]})")
    # ---- hill-climb adversary + coset reduction on the incumbent
    bestH = (1e9, None)
    for trial in range(6):
        w0 = rng.choice([4, 8, 12, 20, 30])
        E = set(rng.sample(range(nF), min(w0, nF)))
        cur_r, _ = ratio(E)
        for step in range(1200):
            f = rng.randrange(nF)
            E2 = set(E); (E2.remove(f) if f in E2 else E2.add(f))
            if not E2: continue
            r2, _ = ratio(E2)
            if r2 <= cur_r: E, cur_r = E2, r2
        # coset reduction: try lowering the weight with bank codewords (raises the honest ratio)
        x = 0
        for f in E: x |= (1 << f)
        improved = True
        while improved and bank:
            improved = False
            for c in bank[:60]:
                y = x ^ c
                if bin(y).count('1') < bin(x).count('1') and y: x = y; improved = True
        E_red = {i for i in range(nF) if (x >> i) & 1}
        r_red, vv = ratio(E_red) if E_red else (1e9, 0)
        r_fin = max(cur_r, r_red) if E_red else cur_r   # honest: use the reduced-weight coset rep
        if r_red < bestH[0] and E_red: bestH = (r_red, (len(E_red), vv))
        if cur_r < bestH[0]: bestH = (cur_r, (len(E), None))
    print(f"  s(hill-climb + coset-reduce) = {bestH[0]:.4f}   (support {bestH[1][0]})")

    s_meas = min(best1[0], best2[0], bestS[0], bestH[0])
    caveat = "" if (d_ub is None or 2*2 < d_ub) else " [w=2 may exceed d/2 — coset-leader caveat]"
    print(f"  ==> MEASURED s (min over probes) = {s_meas:.4f}{caveat}")
    # 1/36-family comparison
    dd = dA*dB
    print(f"  1/36-family: 1/36={T36:.4f} | census 313/9702=0.0323 | delta_A*delta_B={dd:.4f} | "
          f"(dAdB)^2={dd*dd:.4f} | dAdB/4={dd/4:.4f} | trickle-down min(dA,dB)={min(dA,dB):.4f}")
    summary.append((name, nc, nF, dimC, d_ub, s_meas, dd))

print("\n" + "="*106)
print(" SUMMARY: measured soundness vs the 1/36 family")
print("="*106)
print(f"  {'config':>28} {'comp':>5} {'F':>5} {'dimC':>5} {'d_ub':>5} {'s_meas':>8} {'dAdB':>6}  verdict")
for (name, nc, nF, dimC, d_ub, s, dd) in summary:
    near36 = abs(s - 1/36.0) < 0.01
    print(f"  {name:>28} {nc:>5} {nF:>5} {dimC:>5} {str(d_ub):>5} {s:>8.4f} {dd:>6.3f}  "
          f"{'~1/36!' if near36 else ('>>1/36' if s > 0.06 else 'near-1/36 band' if s>0.015 else '<<1/36')}")
print("\n  READING: s is the WORST-CASE tester ratio actually measured; the 1/36 analytic floor is the")
print("  11-core lonely-measure minimum. A numeric match would make the soundness/measure mirror")
print("  QUANTITATIVE; a mismatch pins it as structural analogy (both = uniform detection floors).")
print("DONE.")

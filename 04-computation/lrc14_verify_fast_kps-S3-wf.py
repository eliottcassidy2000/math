"""
FAST adversarial verification of lrc-lever-overlap (LRC14 S3).

Cost insight: Mval on Vmax~400 costs ~2.6s (9k candidates). C(S) needs only Wwidth (cheap).
Strategy:
  PRIMARY screen = C(S) via Wwidth (the PROVED-implication criterion). Cheap.
  Only call expensive Mval on sets where C(S) FAILS (to confirm/deny a real counterexample),
  which the claim says is essentially never.
Plus:
  - large-Vmax persistence test using ONLY Wwidth (cheap even for big speeds? No -- Wwidth ~ sum
    of speeds. Cap Vmax for Wwidth too, but push higher than Mval-bound since Wwidth is ~5x cheaper).
  - WIDTH-INEQ + L_collapse + L_ruler geometry checks (pure arithmetic, instant).

Prints progress so we can see partial results.
kind-pasteur-S3-wf
"""
from fractions import Fraction as F
from math import gcd
import itertools, random, sys, time

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def gcd_all(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g

def case_of(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13*Vmin: return 'S2'
    return 'S3'

def best_C(S):
    best = F(0); arg = None
    for v in S:
        A = [x for x in S if x != v]
        r = Wwidth(A) * 7 * v
        if r > best: best = r; arg = v
    return best, arg

# ---------------- generators ----------------
def gen_S3(rng, vmax_cap=150, spread_lo=14, spread_hi=60):
    for _ in range(300):
        nsmall = rng.randint(8, 12)
        small = sorted(rng.sample(range(1, 14), nsmall))
        nlarge = 13 - nsmall
        if nlarge < 2: continue
        V0 = rng.randint(14, vmax_cap - spread_hi - 1)
        if V0 < 14: continue
        spread = rng.randint(spread_lo, spread_hi)
        large = set(); tries = 0
        while len(large) < nlarge and tries < 100:
            large.add(V0 + rng.randint(0, spread)); tries += 1
        large = sorted(large)
        if len(large) < nlarge: continue
        S = sorted(set(small) | set(large))
        if len(S) != 13: continue
        if max(S) > vmax_cap: continue
        if gcd_all(S) != 1: continue
        if not is_covering(S): continue
        if case_of(S) != 'S3': continue
        return S
    return None

# ============================================================================
def test_geometry():
    print("[GEOM] L_collapse / L_ruler / WIDTH-INEQ exact checks")
    rng = random.Random(7); f1 = f2 = f3 = 0
    # L_ruler single-runner safe-arc formula
    for _ in range(5000):
        r = rng.randint(14, 500); j = rng.randint(0, r-1)
        lo = F(14*j+1, 14*r); hi = F(14*j+13, 14*r)
        for tau in (lo, (lo+hi)/2, hi):
            if nrm(r*tau) < F(1, 14): f1 += 1
    # L_collapse uniform bound
    for _ in range(3000):
        w0 = rng.randint(14, 500); wT = w0 + rng.randint(0, 80); j = rng.randint(0, 5)
        lo = F(14*j+1, 14*w0); hi = F(14*j+13, 14*wT)
        if lo >= hi: continue
        for u in (w0, (w0+wT)//2, wT):
            for tau in (lo, (lo+hi)/2, hi):
                if nrm(u*tau) < F(1, 14): f2 += 1
    # WIDTH-INEQ algebra + sufficiency
    for _ in range(20000):
        w0 = rng.randint(14, 1000); wT = w0 + rng.randint(0, 200); Vmax = wT + rng.randint(0, 500)
        width = F(13, 14)/wT - F(1, 14)/w0
        ineq = (13*w0 - wT)*Vmax > 2*w0*wT
        if (width > F(1, 7*Vmax)) != ineq: f3 += 1
    suff = 0
    for _ in range(20000):
        w0 = rng.randint(14, 1000); wT = rng.randint(14, 11*w0 - 1); Vmax = wT + rng.randint(0, 1000)
        if 11*w0 > wT and Vmax >= wT and not ((13*w0 - wT)*Vmax > 2*w0*wT):
            suff += 1
    print(f"  L_ruler violations={f1}  L_collapse violations={f2}  WIDTH-INEQ algebra-mismatch={f3}  sufficiency-cex={suff}")
    return f1 == 0 and f2 == 0 and f3 == 0 and suff == 0

def test_random_S3(num, vmax_cap, label):
    print(f"[RAND {label}] {num} sets, vmax_cap={vmax_cap}")
    rng = random.Random(hash(label) & 0xffff)
    checked = 0; Cfail = []; minmargin = F(100); worstmarg = None
    t0 = time.time()
    for i in range(num):
        S = gen_S3(rng, vmax_cap=vmax_cap)
        if S is None: continue
        checked += 1
        bc, arg = best_C(S)
        if bc < minmargin: minmargin = bc; worstmarg = (S, arg)
        if bc <= 1:
            m = Mval(S)  # only expensive call when C fails
            Cfail.append((S, bc, m))
        if (i+1) % 500 == 0:
            print(f"    ...{i+1}/{num} checked={checked} Cfails={len(Cfail)} minmarg={float(minmargin):.3f} ({time.time()-t0:.0f}s)")
            sys.stdout.flush()
    print(f"  checked={checked}  C-total-failures={len(Cfail)}  min best-C-margin={float(minmargin):.4f} at {worstmarg}")
    for S, bc, m in Cfail[:10]:
        print(f"    C FAILS: {S} margin={float(bc):.4f} M={m}={float(m):.5f} (>=1/14? {m>=F(1,14)})")
    return Cfail, minmargin, worstmarg

def test_boundary_k2(BCAP):
    print(f"[BOUNDARY k=2] large speeds in [14,{BCAP}], exhaustive over 11-subsets")
    n_total = 0; Cfail = []; minmargin = F(100); worstmarg = None; minM = F(1); worstM = None
    Mbelow = []
    t0 = time.time(); cnt = 0
    for P in itertools.combinations(range(1, 14), 11):
        P = list(P)
        for a in range(14, BCAP+1):
            for b in range(a+1, BCAP+1):
                S = sorted(set(P) | {a, b})
                if len(S) != 13: continue
                if gcd_all(S) != 1: continue
                if not is_covering(S): continue
                if case_of(S) != 'S3': continue
                n_total += 1
                bc, arg = best_C(S)
                if bc < minmargin: minmargin = bc; worstmarg = (S, arg)
                if bc <= 1:
                    m = Mval(S); Cfail.append((S, bc, m))
                    if m < minM: minM = m; worstM = S
                    if m < F(1, 14): Mbelow.append((S, m))
        cnt += 1
        if cnt % 13 == 0:
            print(f"    ...subsets {cnt}/78 total={n_total} Cfails={len(Cfail)} minmarg={float(minmargin):.3f} ({time.time()-t0:.0f}s)")
            sys.stdout.flush()
    print(f"  total covering S3 k=2 sets={n_total}  C-failures={len(Cfail)}  min best-C-margin={float(minmargin):.4f} at {worstmarg}")
    for S, bc, m in Cfail[:20]:
        print(f"    C FAILS: {S} margin={float(bc):.4f} M={m}={float(m):.5f} (>=1/14? {m>=F(1,14)})")
    print(f"  #M<1/14 among C-failures: {len(Mbelow)}")
    for S, m in Mbelow[:20]:
        print("    M<1/14:", S, m)
    return Cfail, minmargin, worstmarg

def test_persistence_largeVmax():
    print("[PERSIST] large-Vmax: does C(S) keep holding as V0 grows? (Wwidth only)")
    # fix offset pattern Delta and small part; grow w0
    small = [1, 2, 3, 5, 6, 7, 9, 11]
    Delta = [0, 3, 9, 12, 15]
    results = []
    for V0 in [200, 500, 1000, 2000, 4000, 8000]:
        large = [V0 + d for d in Delta]
        S = sorted(set(small) | set(large))
        if len(S) != 13:
            print("   skip V0", V0, "len", len(S)); continue
        if not is_covering(S):
            print("   V0", V0, "NOT covering; skip"); continue
        if case_of(S) != 'S3':
            print("   V0", V0, "case", case_of(S), "skip"); continue
        bc, arg = best_C(S)
        results.append((V0, bc, arg))
        print(f"   V0={V0:6d}: best-C-margin={float(bc):.4f} via v={arg}  (C holds? {bc>1})")
        sys.stdout.flush()
    return results

def test_known_exc():
    print("[KNOWN-EXC] [1,2,3,5,7,8,9,10,11,12,13,27,28]")
    S = [1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 27, 28]
    m = Mval(S)
    print(f"  case={case_of(S)} covering={is_covering(S)} M={m}={float(m):.5f} (>=1/14? {m>=F(1,14)}, ==2/19? {m==F(2,19)})")
    for v in S:
        A = [x for x in S if x != v]
        r = Wwidth(A) * 7 * v
        print(f"   remove {v:3d}: ratio={float(r):.4f} {'C-YES' if r > 1 else 'no'}")
    return m

def main():
    print("FAST ADVERSARIAL VERIFICATION — lrc-lever-overlap S3\n")
    g = test_geometry()
    print()
    kx = test_known_exc()
    print()
    p = test_persistence_largeVmax()
    print()
    cf_b, mm_b, wm_b = test_boundary_k2(BCAP=80)
    print()
    cf_r1, mm_r1, wm_r1 = test_random_S3(2000, 150, "tight")
    print()
    cf_r2, mm_r2, wm_r2 = test_random_S3(1500, 120, "wide-ish")
    print()
    print("="*70)
    print("SUMMARY")
    print(f"  geometry (L_ruler/L_collapse/WIDTH-INEQ) all exact-sound: {g}")
    print(f"  known-exc M = {kx} = {float(kx):.5f}  (>=1/14: {kx>=F(1,14)})")
    print(f"  persistence largeVmax: {[(v,float(bc)) for v,bc,a in p]}")
    print(f"  boundary k=2 [14,80]: C-failures={len(cf_b)} min-margin={float(mm_b):.4f}")
    print(f"  random tight: C-failures={len(cf_r1)} min-margin={float(mm_r1):.4f}")
    print(f"  random wide:  C-failures={len(cf_r2)} min-margin={float(mm_r2):.4f}")
    allcf = cf_b + cf_r1 + cf_r2
    real_cex = [(S,m) for S,bc,m in allcf if m < F(1,14)]
    print(f"  TOTAL C-failures across all: {len(allcf)}")
    print(f"  TRUE LRC14 counterexamples (M<1/14): {len(real_cex)}")
    for S,m in real_cex[:20]:
        print("    !!", S, m)

if __name__ == '__main__':
    main()

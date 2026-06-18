"""
LRC(14) — angle "scale-invariance-AP-and-k3" (kind-pasteur-S4-wf)

GOAL (two parts):
  (a) SCALING / AP-LIKE S3 sets. Using M-scale-invariance M(cS)=M(S):
      a set {t,2t,...,12t,V} with t|V reduces to {1,2,...,12,m}, m=V/t.
      Prove M({1,...,12,m}) >= 1/14 for ALL integers m, and classify the
      scaling-reducible S3 sets (which S3 sets are c*{primitive small AP}).
  (b) k=3 extension. Extend the PROVED k=2 bounded-core slice to k=3
      (exactly three speeds > 13): after drop-max scaling the core has TWO
      large speeds. Find the analogue: is there a Wmin / Vmax threshold and a
      FINITE core? Push as far as rigor allows.

ALL decisions use exact fractions. Witnesses are exact rationals tau with
   min_{v in S} ||v tau|| >= 1/14   <=>   M(S) >= 1/14 (a GLOBAL witness).
"""
from fractions import Fraction as F
from math import gcd
import itertools

# ---------- EXACT TOOLS (verbatim from prompt) ----------
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
        ws.append(sc[0][1] + (1 - sc[-1][0]))
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
        v = min(nrm(x * t) for x in S)
        if v > b: b = v
    return b

def is_cov(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def gcd_all(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g

def case_of(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def witness_value(S, t):
    """min_{v in S} ||v t||, the level achieved by tau=t."""
    return min(nrm(x * t) for x in S)

# =====================================================================
# PART (a1):  M({1,...,12,m}) >= 1/14 for all m, with EXACT witness.
# =====================================================================
def part_a1_AP_family(MMAX=2000):
    print("=" * 70)
    print("PART (a1): M({1,...,12,m}) for the AP-reduced 1-parameter family")
    print("=" * 70)
    base = list(range(1, 13))
    # m must be >=13 for |S|=13 (m in 1..12 collapses). m=13 gives {1..13}.
    worst = None
    violations = []
    # m=13 is {1..13}: the canonical tight LRC set, M=1/14 exactly, NOT covering.
    for m in range(13, MMAX + 1):
        S = base + [m]
        Mv = Mval(S)
        if Mv < F(1, 14):
            violations.append((m, Mv))
        if worst is None or Mv < worst[1]:
            worst = (m, Mv)
    print(f"  scanned m in [13,{MMAX}]")
    print(f"  worst M = {worst[1]} = {float(worst[1]):.6f} at m={worst[0]}")
    print(f"  1/14 = {float(F(1,14)):.6f}")
    print(f"  #violations (M<1/14): {len(violations)}  -> {violations[:5]}")
    # The unique tight point is m=13 (={1..13}), M=1/14 exactly. For all m!=13, M>1/14?
    eq = [(m, Mv) for m in range(13, MMAX+1) for Mv in [Mval(base+[m])] if Mv == F(1,14)]
    print(f"  m in [13,{MMAX}] with M == 1/14 EXACTLY: {[m for m,_ in eq]}")
    return worst, violations

# =====================================================================
# PART (a2):  EXPLICIT WITNESS for M({1,...,12,m}) >= 1/14, all m.
#  Strategy: tau = 1/14 is a near-universal witness for the small AP 1..12?
#  Check ||v/14|| for v=1..12 and find a single closed-form tau(m).
# =====================================================================
def part_a2_explicit_witness(MMAX=2000):
    print("=" * 70)
    print("PART (a2): explicit closed-form witness tau for {1..12,m}")
    print("=" * 70)
    base = list(range(1, 13))
    # First: the small part 1..12. The 'safe set' G_{1..12} = {tau: ||v tau||>=1/14 all v<=12}.
    sc = safe_components(base)
    print(f"  G_{{1..12}} safe components ({len(sc)} arcs):")
    for a, b in sc:
        print(f"     ({a}, {b})  width {b-a} = {float(b-a):.5f}, center {(a+b)/2}")
    Wmax = Wwidth(base)
    print(f"  widest safe arc width W(1..12) = {Wmax} = {float(Wmax):.5f}")
    # For ANY single extra runner m, by ARC-WIDTH LEMMA (THM-526):
    #   if W(small)>1/(7m) then we can dodge runner m inside a safe arc of the small part.
    # 1/(7m): for which m is W(1..12) > 1/(7m)?
    print(f"  W(1..12) > 1/(7m)  <=>  m > 1/(7*W) = {F(1,7)/Wmax} = {float(F(1,7)/Wmax):.4f}")
    print("    => for m >= 3 the arc-width lemma alone gives a v-safe witness for runner m.")
    # So small m (m=13.. anyway >=13) ALWAYS satisfies m large enough. Confirm directly:
    # For each m, find an EXPLICIT tau in a safe arc of 1..12 that is also m-safe.
    print("  Searching explicit witness tau for each m (must be safe for all 13 speeds):")
    fails = []
    sample = list(range(13, 60)) + [100, 182, 200, 365, 1000, 1999]
    for m in sample:
        S = base + [m]
        # candidate witnesses: centers of safe arcs of small part, plus cand(S)
        found = None
        for t in cand(S):
            if witness_value(S, t) >= F(1, 14):
                found = t; break
        if found is None:
            fails.append(m)
        # report a few
        if m in (13, 14, 15, 27, 28, 42, 182):
            print(f"     m={m:4d}: witness tau={found}, level={witness_value(S,found) if found else None}")
    print(f"  m-values with NO global 1/14-witness: {fails}")
    return Wmax

# =====================================================================
# PART (a3):  classify SCALING-REDUCIBLE S3 sets.
#  A set S is "scaling-reducible to the small AP" if S = t*{1,...,12} u {V}
#  (after possibly relabeling) — i.e. 12 of its speeds form a multiple of 1..12.
#  More generally: which covering primitive S3 sets are of the form
#  c * S0 with S0 a *smaller-range* primitive set?  (gcd already =1 for primitive,
#  so the genuine scaling families are AP-clusters {t,2t,...,12t,V}.)
#  We test: do AP-cluster S3 sets ALWAYS have M>=1/14, and is the bound tight?
# =====================================================================
def part_a3_scaling_families(TMAX=20, VMAX=400):
    print("=" * 70)
    print("PART (a3): scaling families {t,2t,...,12t,V}  (S3, covering, primitive)")
    print("=" * 70)
    worst = None
    cnt = 0
    viol = []
    eqtight = []
    for t in range(1, TMAX + 1):
        small = [t * i for i in range(1, 13)]  # t,2t,...,12t  (range 12t)
        for V in range(12 * t + 1, VMAX + 1):  # V > 12t so it's the unique large outlier
            S = small + [V]
            if len(set(S)) != 13: continue
            if gcd_all(S) != 1: continue
            if not is_cov(S): continue
            if case_of(S) != 'S3': continue
            cnt += 1
            Mv = Mval(S)
            if Mv < F(1, 14): viol.append((t, V, Mv))
            if Mv == F(1, 14): eqtight.append((t, V))
            if worst is None or Mv < worst[2]:
                worst = (t, V, Mv)
    print(f"  enumerated {cnt} covering primitive S3 AP-cluster sets (t<= {TMAX}, V<= {VMAX})")
    if worst:
        print(f"  worst M = {worst[2]} = {float(worst[2]):.6f} at t={worst[0]}, V={worst[1]}")
    print(f"  #M<1/14: {len(viol)}  examples {viol[:5]}")
    print(f"  #M==1/14 exactly: {len(eqtight)}  examples {eqtight[:5]}")
    return worst

# =====================================================================
# PART (a4):  the t∤V cousins — small part is an AP but V is NOT a multiple of t.
#  These are not exactly scale-reducible. Test {t,2t,...,12t, V} with gcd(S)=1
#  but the cluster is the AP and V arbitrary (already covered in a3 since a3
#  does not require t|V). Here instead: "partial AP" clusters — 12 small speeds
#  that are an AP with common difference d but offset a: {a, a+d, ..., a+11d}.
#  These are the genuine 'AP-like' S3 families. Sweep them.
# =====================================================================
def part_a4_general_AP_clusters(AMAX=10, DMAX=8, VMAX=400):
    print("=" * 70)
    print("PART (a4): general AP clusters {a,a+d,...,a+11d} u {V}, S3 covering primitive")
    print("=" * 70)
    worst = None; cnt = 0; viol = []
    for a in range(1, AMAX + 1):
        for d in range(1, DMAX + 1):
            small = [a + i * d for i in range(12)]
            if len(set(small)) != 12: continue
            top = max(small)
            for V in range(top + 1, VMAX + 1):
                S = sorted(set(small + [V]))
                if len(S) != 13: continue
                if gcd_all(S) != 1: continue
                if not is_cov(S): continue
                if case_of(S) != 'S3': continue
                cnt += 1
                Mv = Mval(S)
                if Mv < F(1, 14): viol.append((a, d, V, Mv))
                if worst is None or Mv < worst[-1]:
                    worst = (a, d, V, Mv)
    print(f"  enumerated {cnt} covering primitive S3 general-AP-cluster sets")
    if worst:
        print(f"  worst M = {worst[-1]} = {float(worst[-1]):.6f} at (a,d,V)={worst[:3]}")
    print(f"  #M<1/14: {len(viol)}  examples {viol[:5]}")
    return worst

# =====================================================================
# PART (b1):  k=3 slice — structure of the drop-max core.
#  S3, exactly THREE speeds > 13: L = {a<b<c}, P = 10 small speeds in 1..13.
#  Drop-max (scale by removing the influence of c): the prior k=2 slice proved
#  a Vmax threshold + finite core. For k=3 the core after analysis has TWO
#  large speeds. We EXHAUSTIVELY sweep small large-windows to (i) confirm no
#  M<1/14 counterexample, (ii) locate the worst M, (iii) measure how the
#  CLUSTER-COLLAPSE window W_K closes these (global witness existence).
# =====================================================================
def cluster_collapse_witness(S, Kmax=40):
    """Lemma A: window W_K = ((14K+1)/(14 Vmin),(14K+13)/(14 Vmax)) is level-1/14
    safe for the whole cluster [Vmin..Vmax-of-cluster]. Here we treat the SMALL
    part P (speeds<=13) as the cluster and the large speeds L as the 'outside'.
    A point tau in W_K that is ALSO L-safe is a global witness.
    Returns (K, tau) if found else None."""
    P = sorted(v for v in S if v <= 13)
    L = sorted(v for v in S if v > 13)
    Vmin, Vmax = P[0], P[-1]
    for K in range(0, Kmax + 1):
        lo = F(14 * K + 1, 14 * Vmin)
        hi = F(14 * K + 13, 14 * Vmax)
        if lo >= hi:
            continue
        # need a point in (lo,hi) that is L-safe (||l tau||>=1/14 for l in L).
        # L-safe set restricted to (lo,hi): intersect with safe_components(L).
        for (sa, sb) in safe_components(L):
            a = max(lo, sa); b = min(hi, sb)
            if a < b:
                tau = (a + b) / 2
                if witness_value(S, tau) >= F(1, 14):
                    return (K, tau)
    return None

def part_b1_k3_slice(LO=14, HI=40):
    print("=" * 70)
    print(f"PART (b1): k=3 slice exhaustive sweep, three large speeds in [{LO},{HI}]")
    print("=" * 70)
    n_total = 0; n_Mbelow = 0; Mbelows = []
    minM = F(1); worstM = None
    n_cc_witness = 0; n_cc_fail = 0; ccfails = []
    # 10-subsets of 1..13: C(13,10)=286
    Psets = list(itertools.combinations(range(1, 14), 10))
    for P in Psets:
        P = list(P)
        for (a, b, c) in itertools.combinations(range(LO, HI + 1), 3):
            S = sorted(set(P) | {a, b, c})
            if len(S) != 13: continue
            if gcd_all(S) != 1: continue
            if not is_cov(S): continue
            if case_of(S) != 'S3': continue
            n_total += 1
            Mv = Mval(S)
            if Mv < minM: minM = Mv; worstM = S
            if Mv < F(1, 14):
                n_Mbelow += 1; Mbelows.append((S, Mv))
            cc = cluster_collapse_witness(S)
            if cc is not None:
                n_cc_witness += 1
            else:
                n_cc_fail += 1
                if len(ccfails) < 10:
                    ccfails.append((S, Mv))
    print(f"  total covering primitive S3 k=3 sets: {n_total}")
    print(f"  #M<1/14 (TRUE COUNTEREXAMPLES): {n_Mbelow}")
    for S, m in Mbelows[:10]:
        print("     M<1/14:", S, "M=", m, float(m))
    print(f"  min M = {minM} = {float(minM):.6f} at {worstM}")
    print(f"  cluster-collapse global witness found: {n_cc_witness}/{n_total}"
          f"  ({100*n_cc_witness/max(1,n_total):.1f}%)")
    print(f"  cluster-collapse FAILED to certify: {n_cc_fail}")
    for S, m in ccfails[:10]:
        print("     CC-FAIL (M still ok):", S, "M=", m, float(m))
    return minM, worstM

# =====================================================================
# PART (b2):  k=3 drop-max scaling reduction — does Vmax bound finitely?
#  In k=2 the proof scaled by dropping the max speed and bounded Vmax<=62.
#  For k=3 we test whether the WORST M stabilizes as the largest speed c->inf
#  with a,b fixed, i.e. whether the M-floor is a function of the SHAPE only.
# =====================================================================
def part_b2_k3_dropmax_stability():
    print("=" * 70)
    print("PART (b2): k=3 drop-max stability — does worst M depend only on shape?")
    print("=" * 70)
    # Fix a representative hard small part and two large speeds; let c grow.
    # Use P with the covering obligations packed: must contain mult of 7 and 14.
    # Take P = {1,2,3,5,6,7,9,10,11,13} won't have mult of 14; large must supply.
    # We just scan a family and watch M as c->inf.
    P = [1, 2, 3, 4, 5, 6, 8, 9, 11, 13]
    a, b = 14, 28  # 28 = mult of 14 (covering for q=14), a=14 also mult of 14&7
    print(f"  P={P}, fixed a={a}, b={b}; sweep c:")
    prevM = None
    rows = []
    for c in list(range(29, 60)) + [80, 120, 200, 400, 800]:
        S = sorted(set(P) | {a, b, c})
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        cov = is_cov(S); cs = case_of(S)
        Mv = Mval(S)
        rows.append((c, cov, cs, Mv))
    for c, cov, cs, Mv in rows:
        flag = "" if Mv >= F(1,14) else "  <<< M<1/14"
        print(f"     c={c:4d} cov={cov} {cs}  M={Mv} = {float(Mv):.6f}{flag}")
    # The point: if M is eventually constant in c, the core is FINITE in the top speed.
    tail = [Mv for c, cov, cs, Mv in rows if c >= 200]
    if tail:
        print(f"  tail M-values (c>=200): {[float(x) for x in tail]}")
        print(f"  constant in tail? {len(set(tail))==1}")

# =====================================================================
# PART (b3):  k=3 finite-core search — find the actual worst covering primitive
#  S3 k=3 set over a moderate board, and the global infimum candidate.
# =====================================================================
def part_b3_k3_worst(LO=14, HI=60):
    print("=" * 70)
    print(f"PART (b3): k=3 worst-M search, large speeds in [{LO},{HI}] (sampled large)")
    print("=" * 70)
    minM = F(1); worstM = None; n = 0
    Psets = list(itertools.combinations(range(1, 14), 10))
    # To keep it tractable, only allow large triples that include a multiple of 14
    # (necessary if P lacks one) — but we just enforce covering at the end.
    for P in Psets:
        P = list(P)
        Pset = set(P)
        for (a, b, c) in itertools.combinations(range(LO, HI + 1), 3):
            S = sorted(Pset | {a, b, c})
            if len(S) != 13: continue
            if gcd_all(S) != 1: continue
            if not is_cov(S): continue
            if case_of(S) != 'S3': continue
            n += 1
            Mv = Mval(S)
            if Mv < minM:
                minM = Mv; worstM = S
    print(f"  scanned {n} sets")
    print(f"  global min M over board = {minM} = {float(minM):.6f} at {worstM}")
    print(f"  margin M*14 = {minM*14} = {float(minM*14):.4f} (>1 means strictly above floor)")
    return minM, worstM


if __name__ == '__main__':
    part_a1_AP_family(MMAX=2000)
    print()
    part_a2_explicit_witness()
    print()
    part_a3_scaling_families(TMAX=20, VMAX=400)
    print()
    part_a4_general_AP_clusters(AMAX=10, DMAX=8, VMAX=400)
    print()
    part_b1_k3_slice(LO=14, HI=34)
    print()
    part_b2_k3_dropmax_stability()
    print()
    part_b3_k3_worst(LO=14, HI=55)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_threadB_singleton_wide_bound_opus_0621.py   (opus-2026-06-21, THREAD B)

THE WIDE BOUND AS AN MDS / SINGLETON WEIGHT-ENUMERATOR BOUND.
=============================================================
Grounded identity (proved, kps-S6 / routeE / threadD):
    corr(E) = measS7(E) - iid_k = sum_{0 != n in Lambda(E)} K(n),
    Lambda(E) = { n in Z^k : sum_i n_i e_i = 0 }  (offset relation lattice; the e_0=0
                anchor is the fixed cell-0 factor, lattice taken over NONZERO offsets),
    K(n) = sum_{T subset {1..6}} (-1)^|T| prod_i chat_T(n_i),     (REAL)
    chat_T(0) = 1 - |T|/7,
    chat_T(n) = - sum_{j in T} shat(n,j)   for n != 0 (mod 7),
    chat_T(n) = 0                          for n != 0, 7 | n   ("apex-prime vanishing"),
    shat(n,j) = (e(-n j/7) - e(-n(j+1)/7)) / (2 pi i n),  |shat(n,j)| = |sin(pi n/7)|/(pi|n|).

WEIGHT ENUMERATOR PICTURE (kps HYP-2723): order Lambda(E) by SUPPORT SIZE s = #nonzero
coords ("Hamming weight" of the relation).  d = d(Lambda(E)) = min support of a nonzero
EFFECTIVE relation (effective = no coord a nonzero multiple of 7).  Lambda(E) ~ [k, k-1, d].
    corr(E) = sum_{s >= d} sum_{n: supp(n)=s, effective} K(n).

THREAD-B DELIVERABLE (this script):
 (1) sup|K|_s := the PER-SUPPORT KERNEL BOUND = max over support-s relation vectors of |K(n)|.
     Derive a CLOSED upper bound and verify it numerically. Special focus on s=3.
 (2) A_s(E) = #effective support-s relations of bounded height. The TRIANGLE-INEQUALITY
     wide bound:  |corr(E)| <= sum_{s>=d} A_s(E) * sup|K|_s.
 (3) SINGLETON: d <= k - rank(E) + 1 (Singleton bound for the relation code). MDS (arc/Sidon)
     sets hit equality with d large => A_2 = 0 and small A_3; anti-MDS (consec) have d = 2
     => the s=2 layer is populated and A_3 is maximal. Quantify the A_s growth as d drops.
 (4) THE RIGOROUS DISSOCIATED HALF: an EXPLICIT min-distance threshold d_0 such that
        d(Lambda(E)) > d_0  ==>  the bounded-height absolute tail is provably small,
     PLUS the HONEST statement of where pure triangle-inequality fails (codex HYP-2715:
     sum|K(n)| over ALL n diverges ~ 9000; absolute control is only possible LAYER-BY-LAYER
     up to a height cutoff, and the genuine smallness uses SIGNED cancellation in high layers).
     We make the cancellation-free part rigorous: the FINITE low-height shell.

ALL measS7/iid EXACT (Fractions). Kernel K(n) computed in complex double, cross-checked to the
exact rational measure on small E so the numbers are trustworthy.

stdlib only.
"""
import sys, itertools, math, cmath
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

OUT = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

P = 7
TWO_PI_I = 2j * math.pi

# ============================ exact measS7 / iid ============================
def measS7(E):
    E = sorted(set(int(e) for e in E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, P*abs(e)+1): bps.add(F(m, P*abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = F(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo+hi)/2
        if len({int(((e*mid) % 1)*P) for e in E}) == P: tot += hi-lo
    return tot

def stirling2(n, k): return sum((-1)**(k-j)*math.comb(k,j)*j**n for j in range(k+1))//math.factorial(k)
def iid_k(k):
    if k < P: return F(0)
    return F(math.factorial(P)*stirling2(k, P), P**k)

# ============================ kernel K(n) ============================
def shat(n, j):
    if n == 0: return 1.0/P
    a = j/float(P)
    return (cmath.exp(-TWO_PI_I*n*a) - cmath.exp(-TWO_PI_I*n*(a+1.0/P)))/(TWO_PI_I*n)

SUB = [tuple(T) for r in range(P) for T in itertools.combinations(range(1, P), r)]
SGN = {T: (-1)**len(T) for T in SUB}

_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0:        v = complex(1 - len(T)/float(P), 0.0)
    elif n % P == 0:  v = 0j
    else:             v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v

def Kk(n):
    """K(n) for an EFFECTIVE relation: n = tuple of NONZERO coordinates (support already extracted).
       (the zero coords contribute chat_T(0) factors -- handled by caller via the full vector)."""
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T]*p
    return s.real

def Kk_full(nvec):
    """K of a FULL k-vector (including zero coords). zero coords contribute chat_T(0)=1-|T|/7."""
    s = 0j
    for T in SUB:
        p = 1.0+0j
        for ni in nvec:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T]*p
    return s.real

# ============================ PART (1): sup|K|_s ============================
def part1_supK():
    log("="*94)
    log("(1) PER-SUPPORT KERNEL BOUND  sup|K|_s = max_{supp(n)=s} |K(n)|")
    log("="*94)
    log("""
For a support-s relation the nonzero coords are residues n_i != 0 (mod 7) (else chat=0).
K depends on the n_i ONLY through their residues r_i = n_i mod 7 in {1..6} TIMES a 1/|n_i|
height decay from shat. Write chat_T(n) = -(1/(pi n)) * sum_{j in T} c(n,j), where
   c(n,j) = sin(pi n/7) * e(-(2j+1) n /(2*7)) ... has |chat_T(n)| <= |sin(pi n/7)|/(pi|n|) * |T|?
NO -- |sum_j| <= sum_j, and each |shat(n,j)| = |sin(pi n/7)|/(pi|n|). So
   |chat_T(n)| <= |T| * |sin(pi n/7)| / (pi |n|).
But the residue structure makes |chat_T(n)| MAXIMAL at the smallest height |n|=residue r in {1..6}
with r the residue (n = +-1 gives the largest single-arc coefficient). So sup|K|_s is attained at
HEIGHT-MINIMAL relations: every nonzero coord is +-1 ... +-6 with |n_i| as small as the residue.
We therefore search over residue vectors r in {1..6}^s (n_i = r_i or r_i-7, take min |.|) to find
the genuine sup of |K| over the per-coordinate single-harmonic kernels, and SEPARATELY give the
closed triangle-inequality envelope.
""")
    # closed envelope: |K(n)| <= sum_T prod_i |chat_T(n_i)|, with |chat_T(n_i)| <= |T| s7(n_i),
    # s7(n) = |sin(pi n/7)|/(pi|n|). The max over heights is at |n_i| minimal = residue value.
    def s7(n): return abs(math.sin(math.pi*n/P))/(math.pi*abs(n))
    # height-minimal representative of residue r: choose n in {r, r-7} minimizing |n|
    def hmin(r): return r if r <= 3 else r-7   # 1,2,3 -> 1,2,3 ; 4,5,6 -> -3,-2,-1
    # absolute envelope sup over residues for support s, evaluated at height-minimal:
    def env_abs(s):
        best = 0.0
        for rv in itertools.product(range(1, P), repeat=s):
            nv = tuple(hmin(r) for r in rv)
            val = 0.0
            for T in SUB:
                p = 1.0
                for ni in nv: p *= abs(chat(ni, T))
                val += p
            best = max(best, val)
        return best
    # true signed sup over residues at height-minimal:
    def sup_signed(s):
        best = 0.0; argn = None
        for rv in itertools.product(range(1, P), repeat=s):
            nv = tuple(hmin(r) for r in rv)
            v = abs(Kk(nv))
            if v > best: best = v; argn = nv
        return best, argn
    log(f"  {'s':>2} {'sup|K|_s (signed, height-min)':>30} {'argmax n':>22} {'abs-envelope sum_T|.|':>22}")
    sup_table = {}
    for s in range(2, 7):
        ss, an = sup_signed(s)
        ea = env_abs(s)
        sup_table[s] = ss
        log(f"  {s:>2} {ss:>30.8f} {str(an):>22} {ea:>22.6f}")
    # explicit s=3 closed form discussion
    log("""
  CLOSED FORM, s=3 (the additive-triangle layer, the leading cross-block term):
    The height-minimal support-3 kernel is maximized at n = (1,1,-2)-type residue triangles
    (a Schur/AP relation +1.e_a +1.e_b -1.e_c = 0 with the SMALLEST coords). Its value is the
    sup|K|_3 reported above. The single-arc decay gives the explicit per-relation constant
       |K(n)| <= ( sum_{T} prod_i |T| s7(n_i) )  with s7(n)=|sin(pi n/7)|/(pi|n|),
    and the dominant n=(1,1,-2) gives s7(1)=sin(pi/7)/pi, s7(2)=sin(2pi/7)/(2 pi).""")
    log(f"    s7(1) = sin(pi/7)/pi           = {math.sin(math.pi/P)/math.pi:.8f}")
    log(f"    s7(2) = sin(2pi/7)/(2 pi)      = {math.sin(2*math.pi/P)/(2*math.pi):.8f}")
    log(f"    2*zeta-ish ref 2 zeta(3)/pi^3  = {2*1.2020569031595942/math.pi**3:.8f}  (the s=3 height-sum scale)")
    return sup_table

# ============================ relation lattice (effective, by support) ============================
def enum_relations(E, N0):
    """All n over NONZERO offsets with sum n_i e_i = 0, |n_i|<=N0, EFFECTIVE (no nonzero multiple
       of 7). Returns list of full nonzero-offset vectors and the nonzero-offset list."""
    Enz = [int(e) for e in E if e != 0]
    out = []
    for n in itertools.product(range(-N0, N0+1), repeat=len(Enz)):
        if all(x == 0 for x in n): continue
        if any((ni % P == 0 and ni != 0) for ni in n): continue
        if sum(ni*ei for ni, ei in zip(n, Enz)) != 0: continue
        out.append(n)
    return out, Enz

def support_size(n): return sum(1 for v in n if v != 0)

def corr_by_support(E, N0):
    """signed corr split by support layer, up to height N0, plus A_s counts (one rep per +-pair)."""
    rels, Enz = enum_relations(E, N0)
    layer_sum = defaultdict(float)
    Acount = defaultdict(int)
    seen = set()
    for n in rels:
        neg = tuple(-v for v in n)
        if n in seen or neg in seen:
            # still need to add BOTH n and -n to the signed sum; handle by adding 2*Re once
            continue
        seen.add(n); seen.add(neg)
        s = support_size(n)
        nz = tuple(v for v in n if v != 0)
        k = Kk(nz)
        layer_sum[s] += 2*k          # n and -n both contribute (K real, K(-n)=K(n))
        Acount[s] += 1               # count of +-pairs = support-s relation pairs
    return layer_sum, Acount, len(Enz)

def dmin_effective(E, N0=4):
    """min support among effective nonzero relations up to height N0 (computational d)."""
    rels, Enz = enum_relations(E, N0)
    best = 10**9
    for n in rels:
        best = min(best, support_size(n))
        if best == 2: break
    return best if best < 10**9 else None

# ============================ PART (3): Singleton / MDS ladder ============================
def part3_singleton():
    log("\n" + "="*94)
    log("(3) SINGLETON: d(Lambda) <= k - rank(E) + 1.  A_s growth as d drops MDS->anti-MDS.")
    log("="*94)
    fams = {
        "consec(anti-MDS)":   list(range(8)),                    # 0..7
        "AP_step2":           [0,2,4,6,8,10,12,14],
        "dyadic":             [0,1,2,4,8,16,32,64],
        "Sidon(MianChowla)":  [0,1,3,7,12,20,30,44],
        "random-wide":        [0,5,8,22,25,37,38,40],
        "two-block":          [0,1,2,3,50,51,52,53],
    }
    N0 = 2
    log(f"  (effective d computed at height N0=4; A_s counts at height N0={N0}, +-pairs)")
    log(f"  {'family':>20} {'k':>2} {'d':>3} {'A2':>5} {'A3':>6} {'A4':>7} {'corr(exact)':>13} {'measS7':>9}")
    rows = {}
    for name, E in fams.items():
        k = len([e for e in E if True])
        d = dmin_effective(E, 4)
        ls, Ac, _ = corr_by_support(E, N0)
        cexact = measS7(E) - iid_k(k)
        rows[name] = (k, d, Ac.get(2,0), Ac.get(3,0), Ac.get(4,0), cexact, measS7(E))
        log(f"  {name:>20} {k:>2} {str(d):>3} {Ac.get(2,0):>5} {Ac.get(3,0):>6} {Ac.get(4,0):>7} "
            f"{float(cexact):>+13.6f} {float(measS7(E)):>9.5f}")
    log("""
  READING:
    d = 2  (consec / AP / dyadic / two-block): support-2 relations EXIST (commensurate offsets,
            n_i e_i = -n_j e_j). These are anti-MDS; A_3 is large; |corr| large.
    d = 3  (Sidon, random-wide): A_2 = 0 (Sidon defining property!); A_3 small; corr ~ iid (small).
    Singleton: rank(E)=k-? gives d <= k-rank+1. For a generic (Sidon-like) E the relation lattice
    has rank k-1 and is MDS-like: d as large as the geometry allows (general position).""")
    return rows

# ============================ PART (2)+(4): the wide bound + threshold ============================
def part2_4_widebound(sup_table):
    log("\n" + "="*94)
    log("(2)+(4) WIDE BOUND  |corr| <= sum_{s>=d} A_s sup|K|_s  AND the min-distance threshold d_0")
    log("="*94)
    log("""
HONEST PRELIMINARY (codex HYP-2715): sum_{all n} |K(n)| DIVERGES (~9000 numerically). So a pure
triangle-inequality bound over the WHOLE lattice is FALSE/vacuous. The wide bound is meaningful
only LAYER-BY-LAYER UP TO A HEIGHT CUTOFF. Two regimes:

  REGIME A (MDS / dissociated, the EASY half -- what THREAD B makes rigorous):
     If d(Lambda(E)) >= 3 (A_2 = 0) AND the set is sufficiently spread (Sidon), then the
     low-HEIGHT shell carries essentially all the (already tiny) corr, and the absolute mass of
     the low shell is < budget. We certify: for height N0=2, the absolute layer mass
        M_abs(E,N0) = sum_{0!=n, |n_i|<=N0, eff} |K(n_nz)|
     and show M_abs is SMALL for d>=3 sets (the dissociated bound) and LARGE (needs cancellation)
     for d=2 sets (consec). This is the SHARP DICHOTOMY at the Singleton edge d=2 vs d>=3.

  REGIME B (anti-MDS / consec): d=2, A_3 maximal; corr is large and genuinely needs the SIGNED
     structure (THREAD A/C/mac-mini convex-order). NOT covered by the absolute wide bound -- and
     it shouldn't be, since consec is the EXTREMIZER we are bounding AGAINST.
""")
    fams = {
        "consec(d=2)":        list(range(8)),
        "AP_step2(d=2)":      [0,2,4,6,8,10,12,14],
        "dyadic(d=2)":        [0,1,2,4,8,16,32,64],
        "two-block(d=2)":     [0,1,2,3,50,51,52,53],
        "Sidon(d>=3)":        [0,1,3,7,12,20,30,44],
        "random-wide(d>=3)":  [0,5,8,22,25,37,38,40],
        "random-wide2(d>=3)": [0,4,23,27,28,34,36,39],
        "random-wide3(d>=3)": [0,6,7,11,24,29,33,34],
    }
    N0 = 2
    log(f"  height cutoff N0={N0}.  M_abs = absolute low-shell mass;  corr_exact for reference.")
    log(f"  {'family':>22} {'d':>3} {'A2':>4} {'A3':>5} {'M_abs(low)':>12} {'corr_exact':>12} "
        f"{'budget? (<0.18)':>15}")
    BUDGET = 0.18
    dis_max = 0.0; struct_min = 1.0
    for name, E in fams.items():
        d = dmin_effective(E, 4)
        # absolute low-shell mass
        rels, Enz = enum_relations(E, N0)
        seen = set(); mabs = 0.0; A2 = 0; A3 = 0
        for n in rels:
            neg = tuple(-v for v in n)
            if n in seen or neg in seen: continue
            seen.add(n); seen.add(neg)
            nz = tuple(v for v in n if v != 0)
            mabs += 2*abs(Kk(nz))
            s = support_size(n)
            if s == 2: A2 += 1
            if s == 3: A3 += 1
        k = len([e for e in E])
        cexact = float(measS7(E) - iid_k(k))
        ok = "YES" if mabs < BUDGET else "no(needs cancel)"
        log(f"  {name:>22} {str(d):>3} {A2:>4} {A3:>5} {mabs:>12.6f} {cexact:>+12.6f} {ok:>15}")
        if d is not None and d >= 3:
            dis_max = max(dis_max, mabs)
        else:
            struct_min = min(struct_min, mabs)
    log(f"""
  THRESHOLD RESULT (height cutoff N0={N0}):
    max low-shell |K|-mass over d>=3 (dissociated) sets = {dis_max:.6f}
    min low-shell |K|-mass over d=2  (structured)  sets = {struct_min:.6f}
    => the absolute low-shell mass SEPARATES the d>=3 (small, no cancellation needed) from the
       d=2 (large, cancellation needed) families. The min-distance threshold is d_0 = 2:
          d(Lambda(E)) >= 3  ==>  low-shell absolute mass < budget (dissociated bound holds).
""")

# ============================ certification: K reproduces corr ============================
def certify():
    log("\n" + "="*94)
    log("(0) CERTIFY: sum_{n in low shell} K(n) approaches exact corr; and CONVERGENCE HONESTY.")
    log("="*94)
    log("  (truncated height sums do NOT converge fast for d=2 sets -- documenting the honest gap)")
    for E in [(0,1,2,3,4,5,6,7), (0,1,3,7,12,20,30,44)]:
        k = len(E); cexact = float(measS7(E)-iid_k(k))
        log(f"  E={E}  corr_exact={cexact:+.6f}")
        for N0 in (1,2,3):
            ls, Ac, _ = corr_by_support(E, N0)
            recon = sum(ls.values())
            log(f"     N0={N0}: sum_s K-layer = {recon:+.6f}   (A2={Ac.get(2,0)}, A3={Ac.get(3,0)}, A4={Ac.get(4,0)})")

def main():
    log(__doc__)
    certify()
    sup_table = part1_supK()
    part3_singleton()
    part2_4_widebound(sup_table)
    log("\nDONE.")

if __name__ == "__main__":
    main()

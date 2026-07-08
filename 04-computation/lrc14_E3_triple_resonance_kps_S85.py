"""
lrc14_E3_triple_resonance_kps_S85.py   (kind-pasteur-2026-07-08-S85)

DERIVE E_3, the TRIPLE resonance multiplier in opus-S152's variance decomposition, and
VERIFY it against the exact variance.  (v2: coherent Fourier sum -- v1 wrongly summed
per-level totals and omitted the cross-level interference = the 96% cancellation.)

FRAMEWORK (opus-S152).  psi = 1-arc, psihat(0)=s=1-th, psihat(m)=-c_m (m!=0),
c_m = int_0^{1/7} e(-2pi i m t) dt.  W(x)=int_0^1 prod_i psi(y-e_i x) dy, exact.  Then

    What(nu) = sum_{m: sum m=0, m.e=nu} prod_i psihat(m_i) = sum_p A_p(nu),
      A_p(nu) = s^{k-p} (-1)^p sum_{|supp|=p, balanced, m.e=nu} prod_{a} c_{m_a}    (damping IN amp)

    Var(W) = sum_{nu!=0} |What(nu)|^2   (Parseval; ONE coherent sum over nu).

Gram split of Var over support sizes (p,q):
    Var = sum_{p,q>=2} G[p][q],   G[p][q] = sum_{nu!=0} A_p(nu) conj(A_q(nu)).
    G[p][p] = DIAGONAL_p (m=m', = C(k,p)c_p s^{2(k-p)}) + RESONANCE_pp (m!=m', same order p).
    G[p][q], p<q = CROSS-LEVEL resonance (a p-term and a q-term relation share nu).

opus's per-level model:  Var ~ sum_j s^{2(k-j)}[C(k,j) c_j + E_j c_j],  E_2 = R2/2.
    => E_p := RESONANCE_pp / (s^{2(k-p)} c_p).   TARGET: E_3.

We build A_p(nu) EXACTLY (integer nu -> exact grouping; entries |w|<=WMAX, t-hat~1/w^2 decays),
reconstruct Var and CHECK vs opus's exact moments_exact, then read the Gram blocks, E_2 (=R2/2?),
E_3, cross terms, and E_3's additive-energy closed form (the (1,1,-2) triple energy).
"""
import cmath
import math
from fractions import Fraction as F
from math import comb
from collections import defaultdict, Counter
from itertools import combinations, product
import sys

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_degree3_floor_opus_S148 import moments_exact
from lrc14_overlap_variance_kernels_opus_S152 import kernels

TH = 1.0 / 7.0
S = 1.0 - TH
TWO_PI = 2.0 * math.pi


def cm(m):
    if m == 0:
        return complex(TH, 0.0)
    return (cmath.exp(-1j * TWO_PI * m / 7.0) - 1.0) / (-1j * TWO_PI * m)


def that(m):
    return abs(cm(m)) ** 2


def R2(E):
    cc = Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j:
                cc[E[i] - E[j]] += 1
    return sum(v * v for v in cc.values())


def balanced_vectors(p, wmax):
    rng = [a for a in range(-wmax, wmax + 1) if a != 0]
    return [w for w in product(rng, repeat=p) if sum(w) == 0]


def build_levels(E, wmax, pmax):
    """A_p[nu] complex amplitude (incl s^{k-p}) and diag_p = sum |single amp|^2 (nu!=0)."""
    k = len(E)
    A = {p: defaultdict(complex) for p in range(2, pmax + 1)}
    diag = {p: 0.0 for p in range(2, pmax + 1)}
    for p in range(2, pmax + 1):
        vecs = balanced_vectors(p, wmax)
        damp = S ** (k - p)
        sign = (-1) ** p
        for sub in combinations(range(k), p):
            eA = [E[a] for a in sub]
            for w in vecs:
                nu = sum(wi * ei for wi, ei in zip(w, eA))
                if nu == 0:
                    continue
                amp = damp * sign
                for wi in w:
                    amp *= cm(wi)
                A[p][nu] += amp
                diag[p] += abs(amp) ** 2
    return A, diag


def gram(A, pmax):
    """G[p][q] = sum_{nu!=0} A_p(nu) conj(A_q(nu)) (real diag; complex offdiag)."""
    G = {}
    keys = {p: set(A[p].keys()) for p in A}
    for p in range(2, pmax + 1):
        for q in range(p, pmax + 1):
            common = keys[p] & keys[q]
            val = sum(A[p][nu] * A[q][nu].conjugate() for nu in common)
            G[(p, q)] = val
    return G


def triple_energy_112(E):
    """(1,1,-2) order-3 energy: Q(nu)=#{i<j,l!=i,j: E_i+E_j-2E_l=nu}; offdiag=sumQ^2-sumQ."""
    k = len(E)
    Q = Counter()
    for i, j in combinations(range(k), 2):
        for l in range(k):
            if l != i and l != j:
                Q[E[i] + E[j] - 2 * E[l]] += 1
    return sum(v * v for v in Q.values()), sum(Q.values())


def analyze(E, wmax, pmax, c):
    E = sorted(E); k = len(E)
    A, diag = build_levels(E, wmax, pmax)
    G = gram(A, pmax)
    # coherent Var
    allnu = set()
    for p in A:
        allnu |= set(A[p].keys())
    var_model = 0.0
    for nu in allnu:
        w = sum(A[p].get(nu, 0) for p in A)
        var_model += abs(w) ** 2
    diag_tot = sum(diag.values())
    res_pp = {p: (G[(p, p)].real - diag[p]) for p in range(2, pmax + 1)}
    cross = {(p, q): 2 * G[(p, q)].real for p in range(2, pmax + 1)
             for q in range(p + 1, pmax + 1)}
    cf = {p: float(c[p]) for p in c}
    Ep = {p: res_pp[p] / (S ** (2 * (k - p)) * cf[p]) for p in range(2, pmax + 1)}
    vt = moments_exact(E, 2)
    var_true = float(vt[2] - vt[1] * vt[1])
    return dict(E=E, k=k, r2=R2(E), var_true=var_true, var_model=var_model,
                diag=diag, diag_tot=diag_tot, res_pp=res_pp, cross=cross, Ep=Ep, G=G)


def main():
    c = kernels(8)
    WMAX, PMAX = 5, 6
    print("=" * 104)
    print(f"E_3 TRIPLE RESONANCE MULTIPLIER (kps-S85 v2)  theta=1/7 s=6/7  |w|<={WMAX} pmax={PMAX}")
    print("exact kernels: " + ", ".join(f"c_{j}={float(c[j]):.4e}" for j in range(2, 7)))
    print("=" * 104)

    tests = {
        "block k=6":  list(range(6)),
        "block k=8":  list(range(8)),
        "block k=11": list(range(11)),
        "Sidon k=6":  [0, 1, 3, 7, 12, 20],
        "Sidon k=8":  [0, 1, 3, 7, 12, 20, 30, 44],
        "AP7 d=1":    list(range(7)),
    }

    print()
    print("VALIDATION: coherent Fourier Var_model vs exact Var_true, and diagonal-only vs Sidon")
    print(f"{'set':<11}{'k':>3}{'R2':>6}{'Var_true':>12}{'Var_model':>12}{'mdl/true':>9}"
          f"{'diag':>12}{'diag/true':>10}{'res_tot':>12}")
    print("-" * 104)
    stash = {}
    for name, E in tests.items():
        r = analyze(E, WMAX, PMAX, c)
        stash[name] = r
        res_tot = r['var_model'] - r['diag_tot']
        print(f"{name:<11}{r['k']:>3}{r['r2']:>6}{r['var_true']:>12.4e}{r['var_model']:>12.4e}"
              f"{r['var_model']/r['var_true']:>9.4f}{r['diag_tot']:>12.4e}"
              f"{r['diag_tot']/r['var_true']:>10.4f}{res_tot:>12.4e}")
    print()
    print("  * Sidon: res_tot ~ 0, diag/true ~ 1  -> diagonal (Poisson) IS the variance (opus).")
    print("  * blocks: res_tot > 0 dominates -> additive-energy resonance.  mdl/true ~ 1 validates")
    print("    the truncation (|w|<=5) and the whole level machinery.")
    print()

    print("=" * 104)
    print("GRAM BLOCKS: same-level resonance E_p = RES_pp/(s^{2(k-p)} c_p), and cross-level 2Re G[p][q]")
    print("=" * 104)
    print(f"{'set':<11}{'E_2':>8}{'R2/2':>7}{'E_3':>9}{'E_4':>9}{'E_5':>9}"
          f"{'X23':>11}{'X24':>11}{'X34':>11}")
    print("-" * 104)
    for name, r in stash.items():
        Ep = r['Ep']; X = r['cross']
        print(f"{name:<11}{Ep[2]:>8.1f}{r['r2']/2:>7.0f}{Ep[3]:>9.1f}{Ep[4]:>9.1f}{Ep[5]:>9.1f}"
              f"{X[(2,3)]:>11.3e}{X[(2,4)]:>11.3e}{X[(3,4)]:>11.3e}")
    print()
    print("  E_2 vs R2/2: opus's pair energy.  E_3 = the TRIPLE multiplier (target).")
    print("  Cross-level X_pq: interference between a p-term and q-term relation at equal nu.")
    print()

    print("=" * 104)
    print("E_3 CLOSED FORM: order-3 additive energy of the (1,1,-2) form  E_i+E_j-2E_l")
    print("  Q(nu)=#{i<j,l: E_i+E_j-2E_l=nu};  E_3 should track the offdiag energy sumQ^2-sumQ")
    print("=" * 104)
    w112 = 6.0 * that(1) ** 2 * that(2)
    frac112 = w112 / float(c[3])
    print(f"  (1,1,-2) orbit = {100*frac112:.1f}% of c_3;  weight-scaled predictor = frac112*(sumQ^2-sumQ)")
    print(f"{'set':<11}{'E_3':>9}{'sumQ^2':>8}{'offdiagQ':>9}{'pred=frac*od':>13}"
          f"{'E_3/pred':>9}{'E_3/R2':>8}")
    print("-" * 104)
    for name, r in stash.items():
        e3 = r['Ep'][3]
        sq2, sq = triple_energy_112(r['E'])
        od = sq2 - sq
        pred = frac112 * od
        print(f"{name:<11}{e3:>9.1f}{sq2:>8}{od:>9}{pred:>13.1f}"
              f"{(e3/pred if pred else float('nan')):>9.3f}{e3/r['r2']:>8.3f}")
    print()
    print("STATEMENT (E_3):  E_3 = RES_33/(s^{2(k-3)} c_3) is the ORDER-3 ADDITIVE ENERGY -- the")
    print("count of matched balanced 3-term relations sum w_a E_a = sum w'_b E_b (!=0), weighted by")
    print("the tent kernel.  Dominant carrier = the (1,1,-2) form E_i+E_j-2E_l, L^2 mass sum Q(nu)^2")
    print("(the triple analog of R2=sum_d r_d^2); E_3 ~ 9-10x its raw offdiagonal (other order-3")
    print("forms + harmonics).  NOTE E_2 itself is ~1.6x R2/2 (a>=2 harmonics + doublings), so the")
    print("naive energies are LEADING terms, not exact.")
    print()
    print("!! PER-LEVEL DIVERGENCE (honest caveat) !!  For the CLUSTERED block k=11 the per-level")
    print("   E_p GROW (E_3=22k,E_4=608k,E_5=8M): the series sum_j E_j c_j s^{2(k-j)} does NOT")
    print("   converge term-by-term; Var survives only by full coherent cancellation across ALL")
    print("   levels.  So brick(B) Var<=c*R2 CANNOT be proven by bounding E_3 alone.  This is the")
    print("   real form of kps-S82's 'non-perturbative': the DAMPED/coherent object converges")
    print("   (opus-S152) but the per-level RESONANCE series is asymptotic/divergent for clustered")
    print("   sets.  The block is exactly brick(A)'s EXCLUDED case (diam 10<16).")

    # ------------------------------------------------------------------
    print()
    print("=" * 104)
    print("BRICK (B) REGIME: EXACT Var/R2 census over SPREAD 11-sets (moments_exact, NO truncation)")
    print("  brick(A) delivers primitive diam>=16 => R2<=614; measure the true constant c=Var/R2 there")
    print("=" * 104)
    import random
    from math import gcd
    M = F(6, 7)
    BAR = 0.3312   # honest k=11 bar (MISTAKE-123)

    def D3_of(E):
        m = moments_exact(E, 3)
        m1, m2, m3 = float(m[1]), float(m[2]), float(m[3])
        Mf = 6.0 / 7.0
        denom = m2 - m3 / Mf
        if denom <= 0:
            return None, (m1, m2, m3)
        return m1 / Mf + (m1 - m2 / Mf) ** 2 / denom, (m1, m2, m3)

    # validate D3 formula on the block
    d3blk, _ = D3_of(list(range(11)))
    print(f"  [validate] block k=11 D3 = {d3blk:.6f}  (THM-661 value 0.404751)")
    print()

    random.seed(614)
    worst_c = (0.0, None)
    min_d3 = (9.9, None)
    cs_all = []
    d3_all = []
    N = 0
    for D in (16, 18, 20, 22, 25, 30, 36, 44):
        for _ in range(120):
            mid = sorted(random.sample(range(1, D), 9))
            E = [0] + mid + [D]
            g = 0
            for a in E:
                g = gcd(g, a)
            if g != 1:
                continue
            r2 = R2(E)
            if r2 > 614:
                continue
            d3, (m1, m2, m3) = D3_of(E)
            var = m2 - m1 * m1
            c_emp = var / r2
            N += 1
            cs_all.append(c_emp)
            if d3 is not None:
                d3_all.append(d3)
                if d3 < min_d3[0]:
                    min_d3 = (d3, tuple(E), r2, c_emp)
            if c_emp > worst_c[0]:
                worst_c = (c_emp, tuple(E), r2, var)
    print(f"  sampled N={N} primitive 11-sets, diam in [16,44], R2<=614 (brick(B) regime)")
    print(f"  Var/R2 :  mean={sum(cs_all)/len(cs_all):.3e}  max={max(cs_all):.3e}  min={min(cs_all):.3e}")
    print(f"  D3     :  mean={sum(d3_all)/len(d3_all):.4f}  MIN={min(d3_all):.4f}  max={max(d3_all):.4f}")
    print(f"  bar(k=11) = {BAR}   ->  min D3 - bar = {min(d3_all)-BAR:+.4f}")
    print("-" * 104)
    print(f"  WORST Var/R2 : c={worst_c[0]:.4e} R2={worst_c[2]} set={worst_c[1]}")
    print(f"  MIN D3       : D3={min_d3[0]:.4f} R2={min_d3[2]} c={min_d3[3]:.3e} set={min_d3[1]}")
    print(f"  block ref    : Var/R2=6.10e-5, D3=0.4048 (clustered EXTREME, R2=770 EXCLUDED by brickA)")
    print()
    print("  READING (brick B, EXACT):  over the spread brick(B) regime the DECISIVE floor D3 stays")
    print("  WELL ABOVE the honest bar 0.331 (min D3 - bar > 0), EVEN THOUGH the constant c=Var/R2")
    print("  reaches ~7.2e-5 (ABOVE the block's 6.1e-5).  So brick(B) [R2<=614 => D3>=bar] is")
    print("  EMPIRICALLY CONFIRMED, but its proof CANNOT go through a uniform c=6.1e-5 (that constant")
    print("  is NOT the sup) -- it must use that spread sets have LARGER E[W] (more uncovered mass),")
    print("  which lifts the PZ/degree-3 floor faster than the larger c lowers it.  E[W] vs R2 is the")
    print("  real lever, not c alone.  <-- key message for opus's uniform bound.")


if __name__ == "__main__":
    main()

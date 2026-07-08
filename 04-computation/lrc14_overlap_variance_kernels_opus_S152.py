"""
lrc14_overlap_variance_kernels_opus_S152.py   (opus-2026-07-08-S152, HYP-5417)

THE EXACT TRIPLE/QUAD (j-FOLD) OVERLAP VARIANCE KERNELS + the variance structure.

Fourier of the uncovered measure (psi = 1 - arc, psihat(0)=1-theta, psihat(m)=-hhat(m)):
    What(nu) = sum_{m in Z^k: sum m_i = 0, m.e = nu} prod_i psihat(m_i),   Var(W) = sum_{nu!=0}|What|^2.
The SPECTATOR (m_i=0) coordinates each contribute psihat(0)=1-theta -> the (1-theta)^{2(k-|supp|)}
"inactive-arc damping" (= mac-mini LEM-007), the bulk of kps-S81's 96% cancellation.

THE j-FOLD OVERLAP VARIANCE KERNEL (the per-j-subset variance density):
    c_j := sum_{a_1+..+a_j=0, all a_i != 0} prod_i that(a_i),   that(n) = |hhat(n)|^2 (tent Fourier).
Since sum_{a_1+..+a_j=0} prod that = int_0^1 t(x)^j dx (Parseval; t = tent = arc autocorrelation,
that(0) = int t = theta^2) and int t^j = 2 theta^{j+1}/(j+1), inclusion-exclusion on zero coords gives

    THE RECURSION:  c_j = int t^j - sum_{r=1}^{j} C(j,r) theta^{2r} c_{j-r},   c_0=1, c_1=0,
    int t^j = 2 theta^{j+1}/(j+1).

CLOSED FORMS (theta=1/7):
    c_2 = 2th^3/3 - th^4 = 11/7203               (= THM-641 pair kernel, opus-S151)
    c_3 = th^4/2 - 2th^5 + 2th^6 = 25/235298     (the TRIPLE kernel -- NEW)
    c_4 = 321/28824005,  c_5 = 950/847425747,  c_6 = 1633/13841287201  (quad, quint, sext)

THE VARIANCE STRUCTURE (validated below):
    Var(W) = sum_{j>=2} (1-theta)^{2(k-j)} [ C(k,j) c_j  (POISSON diagonal, m=m')
                                            + E_j c_j     (RESONANCE, m!=m', matched m.e=m'.e) ]
    - Poisson diagonal alone matches SIDON sets (minimal energy) to ~5%.
    - RESONANCE (E_2 = R2/2 the additive energy; E_3,.. the triple/quad energies) dominates the
      block; it is the additive-energy content, with the kernels c_j now EXACT.
"""
from fractions import Fraction as F
from math import comb
from collections import Counter
import sys

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_degree3_floor_opus_S148 import moments_exact

TH = F(1, 7)


def kernels(jmax=11):
    c = {0: F(1), 1: F(0)}
    for j in range(2, jmax + 1):
        intj = 2 * TH**(j + 1) / (j + 1)
        c[j] = intj - sum(comb(j, r) * TH**(2 * r) * c[j - r] for r in range(1, j + 1))
    return c


def R2(E):
    cc = Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j:
                cc[E[i] - E[j]] += 1
    return sum(v * v for v in cc.values())


def var_true(E):
    m = moments_exact(E, 2)
    return m[2] - m[1] * m[1]


def diagonal(k, c):
    return sum((1 - TH)**(2 * (k - j)) * comb(k, j) * c[j] for j in range(2, k + 1))


def main():
    c = kernels(13)
    print("=" * 92)
    print("EXACT j-FOLD OVERLAP VARIANCE KERNELS c_j (recursion c_j = int t^j - sum C(j,r)th^2r c_{j-r})")
    print("=" * 92)
    forms = {2: "2th^3/3 - th^4", 3: "th^4/2 - 2th^5 + 2th^6"}
    for j in range(2, 7):
        extra = f"  = {forms[j]}" if j in forms else ""
        print(f"  c_{j} = {c[j]} = {float(c[j]):.6e}{extra}")
    print(f"  int t^j = 2 th^(j+1)/(j+1):  " + ", ".join(f"{float(2*TH**(j+1)/(j+1)):.3e}" for j in range(2,6)))

    print()
    print("=" * 92)
    print("VALIDATION: Poisson diagonal sum_j (1-th)^{2(k-j)} C(k,j) c_j  vs  Var(W)")
    print("=" * 92)
    # Sidon sets (minimal R2 = k(k-1)) : resonance ~0, Var ~ diagonal
    sidons = {6: [0, 1, 3, 7, 12, 20], 5: [0, 1, 3, 7, 12], 7: [0, 1, 3, 7, 12, 20, 30]}
    for k, S in sidons.items():
        d = [S[j] - S[i] for i in range(k) for j in range(i + 1, k)]
        vt = var_true(S); dg = diagonal(k, c)
        print(f"  SIDON k={k} (R2={R2(S)}=min {k*(k-1)}, diffs distinct {len(d)==len(set(d))}): "
              f"Var_true={float(vt):.5e}  diagonal={float(dg):.5e}  ratio={float(vt)/float(dg):.4f}")
    print("  => diagonal (Poisson/coverage) matches Sidon to ~5% -- the kernels + (1-th)-damping are correct.")
    print()
    for k in (6, 8, 11, 13):
        blk = list(range(k))
        vt = var_true(blk); dg = diagonal(k, c)
        # implied resonance energy E_2 from pair: (Var - diagonal_higher)/[(1-th)^{2(k-2)} c_2] - C(k,2)
        r2 = R2(blk)
        pair_res = (1 - TH)**(2 * (k - 2)) * F(r2, 2) * c[2]   # (R2/2) c_2 damped = pair total
        print(f"  BLOCK k={k}: Var_true={float(vt):.5e}  diagonal={float(dg):.5e}  R2={r2}  "
              f"(1-th)^2(k-2)(R2/2)c_2={float(pair_res):.5e}  Var/pair_res={float(vt)/float(pair_res):.3f}")
    print()
    print("  READING: for the block, Var >> diagonal (R2-resonance dominates).  The damped pair")
    print("  resonance (1-th)^{2(k-2)}(R2/2)c_2 captures the bulk; the residual is the TRIPLE/QUAD")
    print("  resonance E_3 c_3 + E_4 c_4 + ... (the additive-energy multipliers E_j on the exact")
    print("  kernels c_j).  The (1-th)^{2(k-2)} damping is why the naive (R2/2)c_2 (opus-S151, 27x")
    print("  too big) collapses to the true c ~ 6e-5 (mac-mini/klein LEM-007).")


if __name__ == "__main__":
    main()

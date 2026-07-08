"""
lrc14_fourier_pair_kernel_opus_S151.py   (opus-2026-07-08-S151, HYP-5407)

THE FOURIER PAIR-RESONANCE KERNEL for Var(W) (the resonance lemma Var(W) <= c*R2).

W(x) = uncovered measure; arc indicator a_i(x,y) = chi(e_i x - y), chi = 1_{(-theta,0]},
theta = 1/7.  Fourier: chi(t) = sum_n hhat(n) e(nt), hhat(n) = (1 - e(n theta))/(2 pi i n),
|hhat(n)|^2 = sin^2(pi n theta)/(pi^2 n^2).

THE PAIR TERM (leading, |S|=|T|=2 in kps-S81's W = sum_S (-1)^|S| L_S decomposition; = THM-641
pair overlap mass in Fourier).  The overlap of two arcs at phase-offset delta is the tent
t(delta) = (theta - ||delta||)_+ = autocorrelation of the arc, so its Fourier coeff is
|hhat(n)|^2, and Var(tent) = sum_{n!=0}|hhat(n)|^4 =: c_pair.  Matched-difference pairs
(e_i-e_j = e_k-e_l, the additive energy) share the frequency e(n(e_i-e_j)x), so

    Var_pair(W) = (R2/2) * c_pair,   c_pair = sum_{n!=0}|hhat(n)|^4.

EXACT CLOSED FORM (this session):  c_pair = theta^3 (2/3 - theta) = 11/7203  (theta=1/7),
   = (1/4)[1/30 + (4/3)B4(1/7) - (1/3)B4(2/7)]  (Bernoulli B4(x)=x^4-2x^3+x^2-1/30);
   verified = the k=2 tent variance 2 theta^3/3 - theta^4.

THE SCREENING (kps-S81's ~96% cancellation): the true Var(W) is FAR below Var_pair because the
triple+ overlap terms cancel it.  This computes the exact screening factor s(k) =
Var(block_k)/Var_pair(block_k) and looks for its law (in k, in coverage k*theta).
"""
from fractions import Fraction as F
from math import floor
import sys

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")

TH = F(1, 7)
C_PAIR = TH**3 * (F(2, 3) - TH)   # = 11/7203


def var_block_exact(k):
    """Exact Var(W) for the block {0,..,k-1} via Farey-cell integration (E[W], E[W^2])."""
    from lrc14_pz_degree3_floor_opus_S148 import moments_exact
    m = moments_exact(list(range(k)), 2)
    return m[2] - m[1] * m[1]


def R2_block(k):
    # reduced additive energy of {0..k-1}: sum_{d != 0} r_d^2, r_d = k - |d|
    return (k - 1) * k * (2 * k - 1) // 3


def main():
    print("=" * 92)
    print("THE FOURIER PAIR-RESONANCE KERNEL  c_pair = sum_{n!=0}|hhat(n)|^4")
    print("=" * 92)
    print(f"  c_pair = theta^3(2/3 - theta) = {C_PAIR} = {float(C_PAIR):.8e}  (theta=1/7)")
    print(f"    = k=2 tent variance 2theta^3/3 - theta^4;  Var_pair(W) = (R2/2)*c_pair")
    print()
    print("=" * 92)
    print("THE SCREENING s(k) = Var(block_k) / Var_pair(block_k)  (kps-S81: ~96% cancellation)")
    print("=" * 92)
    print(f"  {'k':>3} {'R2':>6} {'Var_pair=(R2/2)c':>17} {'Var(W) true':>13} "
          f"{'screening s':>12} {'k*theta':>8} {'Var/R2':>10}")
    rows = []
    for k in range(2, 14):
        R2 = R2_block(k)
        vpair = F(R2, 2) * C_PAIR
        vtrue = var_block_exact(k)
        s = float(vtrue) / float(vpair)
        rows.append((k, R2, float(vpair), float(vtrue), s))
        print(f"  {k:>3} {R2:>6} {float(vpair):>17.6f} {float(vtrue):>13.6f} "
              f"{s:>12.4f} {float(k*TH):>8.4f} {float(vtrue)/R2:>10.3e}")
    print()
    print("  READING: s(2) = 1 EXACTLY (pair kernel exact at k=2); s(k) FALLS as coverage k*theta")
    print("  grows (the multiple-overlap screening).  Var/R2 = s(k)*c_pair/2 -- so the effective")
    print("  c(k) = s(k)*c_pair/2; at k=11 (k*theta=1.57) s~0.08 => c ~ 6e-5 (klein's constant).")
    print()
    # look for the screening law: is s(k) ~ some clean function?  test (1-theta)^? , e^{-k theta}, etc.
    print("  screening-law probes (does s(k) match a clean form?):")
    import math
    for (k, R2, vp, vt, s) in rows:
        kt = k * float(TH)
        # candidate laws
        c1 = (1 - float(TH))**(k - 2)      # (6/7)^{k-2}
        c2 = math.exp(-(k - 2) * float(TH))
        c3 = 1.0 / (1 + (k - 2) * float(TH))   # 1/(1+(k-2)theta)
        print(f"    k={k:2d}: s={s:.4f}   (6/7)^(k-2)={c1:.4f}  exp(-(k-2)th)={c2:.4f}  1/(1+(k-2)th)={c3:.4f}")


if __name__ == "__main__":
    main()

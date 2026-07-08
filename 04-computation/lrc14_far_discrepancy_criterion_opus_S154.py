"""
lrc14_far_discrepancy_criterion_opus_S154.py   (opus-2026-07-08-S154)

THE RIGOROUS ABSOLUTE-BOUND CRITERION for  far <= E[W]^2  (the discrepancy route, LEM-005).

From the exact Fourier far (verified in lrc14_far_fourier_discrepancy_opus_S154):
    far = (5/7)^{k+1} + sum_{m in L\0} (5/7)^{k-|S|}(-1)^|S| (prod ahat(m_i)) J(m),
    E[W] = (6/7)^k    + sum_{m in L\0} (6/7)^{k-|S|}(-1)^|S| (prod ahat(m_i)).
Taking ABSOLUTE values (NO cancellation -- this is why it sidesteps the non-perturbative wall
mac-mini/klein hit in the moment expansion) with |J(m)| <= (5/7)2^{|S|} and |ahat(m)| <= 1/(pi|m|):

    |far - (5/7)^{k+1}| <= (5/7)^{k+1} * PHI(E),   PHI(E) := sum_{m in L\0} prod_{i in S} (14/5)|ahat(m_i)|
    |E[W]  - (6/7)^k  | <= (6/7)^k     * PSI(E),   PSI(E) := sum_{m in L\0} prod_{i in S} (7/6) |ahat(m_i)|

Each factor (14/5)|ahat(m)| <= (14/5)/(pi|m|) = 0.891/|m| < 1 and (7/6)|ahat(m)| <= 0.371/|m|, so
BOTH sums are absolutely convergent and dominated by the smallest relations (support>=3: 3-APs
(2,-1,-1) and additive quadruples (1,1,-1,-1)).  Then far <= E[W]^2 is IMPLIED by

    ***  rho_k (1 + PHI(E)) <= (1 - PSI(E))^2,   rho_k := (5/7)^{k+1} / (6/7)^{2k} = (5/7)(35/36)^k  ***

(valid when PSI<1).  rho_11 = 0.5247, so at k=11 the criterion is  PHI <= 0.906  when PSI ~ 0.
This is a SUFFICIENT, EXPLICIT, per-family condition; PHI,PSI -> 0 as the family spreads
(relation lattice thins).  We find the primitive-diameter D0 above which it holds, aiming D0<=24
to meet klein-S184's exhaustive.
"""
import sys, itertools
from math import gcd
from fractions import Fraction as F
import mpmath as mp

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")

mp.mp.dps = 30
TH = mp.mpf(1) / 7
TWO_PI_I = 2 * mp.pi * mp.mpc(0, 1)


def abs_ahat(m):
    if m == 0:
        return TH
    return abs(mp.sinpi(m * TH)) / (mp.pi * abs(m))   # |sin(pi m theta)|/(pi|m|)


def rho(k):
    return (mp.mpf(5) / 7) ** (k + 1) / (mp.mpf(6) / 7) ** (2 * k)


def enumerate_relations(E, Mmax, Smax):
    """All nonzero m in L (support<=Smax, |m_i|<=Mmax), via solve-for-last-two (exact int)."""
    k = len(E)
    rels = []
    for s in range(3, Smax + 1):
        for idx in itertools.combinations(range(k), s):
            ei = [E[i] for i in idx]
            ep, eq = ei[-2], ei[-1]
            den = ep - eq
            free = idx[:-2]
            efree = ei[:-2]
            rng = [v for v in range(-Mmax, Mmax + 1) if v != 0] if s > 2 else [0]
            # free coords take any nonzero value; but a free coord could also be structurally...
            # we allow the s-2 free coords in [-Mmax,Mmax]\{0}; the 2 solved coords may be 0? No:
            # support means all s coords nonzero. We enforce solved coords nonzero below.
            if s == 3:
                free_iter = ([v] for v in rng)
            else:
                free_iter = itertools.product(rng, repeat=s - 2)
            for vals in free_iter:
                A = -sum(vals)                              # m_p + m_q = A
                B = -sum(v * e for v, e in zip(vals, efree))  # m_p ep + m_q eq = B
                num = B - A * eq
                if num % den != 0:
                    continue
                mp_ = num // den
                mq_ = A - mp_
                if mp_ == 0 or mq_ == 0:
                    continue
                if abs(mp_) > Mmax or abs(mq_) > Mmax:
                    continue
                m = list(vals) + [mp_, mq_]
                rels.append((idx, m))
    return rels


def PHI_PSI(E, Mmax, Smax):
    rels = enumerate_relations(E, Mmax, Smax)
    phi = mp.mpf(0)
    psi = mp.mpf(0)
    for idx, m in rels:
        pf = mp.mpf(1)
        pp = mp.mpf(1)
        for v in m:
            a = abs_ahat(v)
            pf *= (mp.mpf(14) / 5) * a
            pp *= (mp.mpf(7) / 6) * a
        phi += pf
        psi += pp
    return phi, psi, len(rels)


def criterion_holds(E, Mmax, Smax):
    k = len(E)
    phi, psi, nr = PHI_PSI(E, Mmax, Smax)
    lhs = rho(k) * (1 + phi)
    rhs = (1 - psi) ** 2 if psi < 1 else mp.mpf(-1)
    return (lhs <= rhs), phi, psi, lhs, rhs, nr


def gcd_all(xs):
    g = 0
    for x in xs:
        g = gcd(g, abs(x))
    return g


def is_primitive(E):
    d = [E[i + 1] - E[i] for i in range(len(E) - 1)]
    return gcd_all(d) == 1


def main():
    K = 11
    print("=" * 100)
    print("RIGOROUS ABSOLUTE-BOUND CRITERION:  rho_k(1+PHI) <= (1-PSI)^2  =>  far <= E[W]^2  (=> Var<=near)")
    print(f"  k={K}: rho_k = {float(rho(K)):.6f};  criterion at PSI~0 is PHI <= {float(1/rho(K)-1):.4f}")
    print("=" * 100)

    # (0) convergence of PHI on a representative spread family
    print("\n(0) CONVERGENCE of PHI,PSI in (Mmax,Smax) on a spread family E=[0,1,3,7,12,20,30,44,60,80,100]:")
    Esp = [0, 1, 3, 7, 12, 20, 30, 44, 60, 80, 100]
    for (M, S) in [(2, 4), (3, 4), (3, 5), (4, 5), (4, 6), (5, 6)]:
        phi, psi, nr = PHI_PSI(Esp, M, S)
        print(f"    Mmax={M} Smax={S}: PHI={float(phi):.5f}  PSI={float(psi):.5f}  #rels={nr}")

    # (1) validate PHI is a true upper bound on the far deviation for small families
    print("\n(1) VALIDATE  |far_exact-(5/7)^{k+1}| <= (5/7)^{k+1} PHI  (PHI is a rigorous upper bound):")
    from lrc14_pz_general_integrator_opus_S148 import pz_exact
    from lrc14_k11_tail_sharp_near_opus_S149 import EW_and_near
    for E in [[0, 1, 3, 7, 12], [0, 2, 3, 4, 5], [0, 1, 4, 9, 11], [0, 1, 3, 7, 12, 20]]:
        k = len(E)
        EW, near = EW_and_near(E)
        _, EW2, _ = pz_exact(E)
        far = EW2 - near
        base = (mp.mpf(5) / 7) ** (k + 1)
        devrel = abs(mp.mpf(float(far)) - base) / base
        phi, psi, nr = PHI_PSI(E, 6, k)
        ok = "OK" if devrel <= phi else "!!!"
        print(f"    E={str(E):22s}: |dev|/base={float(devrel):.4f}  PHI={float(phi):.4f}  [{ok}]")

    # (2) the crossover: scan primitive-diameter bands, worst-case PHI over sampled families
    print("\n(2) CROSSOVER: does the criterion hold across the spread tail? (worst PHI per diam band)")
    print("    per band: max PHI, max PSI, #(criterion holds)/#sampled, worst family")
    import random
    for lo, hi in [(25, 30), (31, 40), (41, 60), (61, 100), (101, 200), (201, 400)]:
        worst_phi = mp.mpf(-1); worst_E = None; nhold = 0; ntot = 0; maxpsi = mp.mpf(0)
        rng = random.Random(lo * 131 + hi)
        for _ in range(300):
            D = rng.randint(lo, hi)
            mid = sorted(rng.sample(range(1, D), K - 2))
            E = [0] + mid + [D]
            if len(set(E)) != K or not is_primitive(E):
                continue
            ntot += 1
            hold, phi, psi, lhs, rhs, nr = criterion_holds(E, 4, 6)
            if hold:
                nhold += 1
            if phi > worst_phi:
                worst_phi = phi; worst_E = E
            if psi > maxpsi:
                maxpsi = psi
        print(f"    prim-diam [{lo:3d},{hi:3d}]: maxPHI={float(worst_phi):.4f} maxPSI={float(maxpsi):.4f}"
              f"  holds {nhold}/{ntot}  worst={worst_E}")

    # (3) the STRUCTURAL worst case: what family maximizes PHI at large diameter?
    print("\n(3) STRUCTURAL worst spread families (near-AP / long-AP embeddings that keep PHI high):")
    cands = {
        "AP_11 (diam 10, COMPACT-ref)": list(range(11)),
        "block10+far {0..9,D=40}": list(range(10)) + [40],
        "2-AP {0,2,..,18, D=41}": list(range(0, 20, 2)) + [41],
        "AP9+2far {0..8,40,80}": list(range(9)) + [40, 80],
        "double-AP {0..5, 40..45+...}": [0, 1, 2, 3, 4, 5, 40, 41, 42, 43, 44],
    }
    for name, E in cands.items():
        E = sorted(set(E))
        if len(E) != K:
            continue
        hold, phi, psi, lhs, rhs, nr = criterion_holds(E, 4, 6)
        pd = "prim" if is_primitive(E) else "NONprim(reduce!)"
        print(f"    {name:34s} {pd:16s}: PHI={float(phi):.4f} PSI={float(psi):.4f}"
              f" lhs={float(lhs):.4f} rhs={float(rhs):.4f} {'HOLDS' if hold else 'FAILS'}")
    print("=" * 100)


if __name__ == "__main__":
    main()

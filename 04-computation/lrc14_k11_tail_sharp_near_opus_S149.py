"""
lrc14_k11_tail_sharp_near_opus_S149.py   (opus-2026-07-08-S149, HYP-5357)

CLOSE THE k=11 TAIL -- the single remaining open density-floor leg.

mac-mini's near/far split (THM-661): E[W^2] = near + far,
  near = 2 int_{1/7}^{2/7} q(L) dL,  q(L) = E_x[sum_i (g_i(x)-L)_+]  (empty-arc prob, arc len L);
  far  = int int_{|y1-y2|>1/7} P(both arcs empty)  ->  <= E[W]^2 for spread families (decorrelation).
Then PZ = E[W]^2/E[W^2] >= E[W]^2/(near + E[W]^2).

mac-mini used the CRUDE near <= (2/7)E[W] (since q(L) <= q(1/7)=E[W]), giving PZ ~ 0.31 at k=11,
JUST under the bar 0.331.  mac-mini: "the last +0.02 needs the near part's STRICT decay
q(L) < q(1/7) for L>1/7, which is real."  THIS computes near EXACTLY (q(L) decays), giving the
sharp PZ tail bound = E[W]^2/(near_exact + E[W]^2), and checks it clears 0.331 over spread k=11
families.

near = 2 E_x[ sum_i h(g_i(x)) ],  h(g) = int_{1/7}^{2/7} (g-L)_+ dL =
   0                         if g <= 1/7
   (g - 1/7)^2 / 2           if 1/7 <= g <= 2/7
   g/7 - 3/98                if g >= 2/7          (= int_{1/7}^{2/7}(g-L)dL)
computed exactly via Farey-cell integration (g_i linear per cell => h(g_i) piecewise-quadratic).
"""
from fractions import Fraction as F
from math import floor, gcd
from itertools import combinations
import sys

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_general_integrator_opus_S148 import BAR

TH = F(1, 7)
TH2 = F(2, 7)


def EW_and_near(E):
    """Exact (E[W], near) for integer family E.  near = 2*E_x[sum_i h(g_i)]."""
    E = sorted(E); k = len(E); e0 = E[0]
    Es = [e - e0 for e in E]
    dens = set(e for e in Es if e > 0)
    for i in range(k):
        for j in range(i + 1, k):
            dens.add(Es[j] - Es[i])
    bps = set([F(0), F(1)])
    for d in dens:
        for p in range(0, d + 1):
            bps.add(F(p, d))
    bps = sorted(bps)
    EW = F(0); NEAR = F(0)
    for c in range(len(bps) - 1):
        x0, x1 = bps[c], bps[c + 1]; xm = (x0 + x1) / 2
        lin = [(F(-floor(e * xm)), F(e)) for e in Es]
        order = sorted(range(k), key=lambda j: lin[j][0] + lin[j][1] * xm)
        sp = [lin[j] for j in order]
        gaps = [(sp[i + 1][0] - sp[i][0], sp[i + 1][1] - sp[i][1]) for i in range(k - 1)]
        gaps.append((F(1) + sp[0][0] - sp[k - 1][0], sp[0][1] - sp[k - 1][1]))
        # sub-breakpoints: each gap crosses 1/7 and 2/7
        subs = set([x0, x1])
        for (a, b) in gaps:
            if b != 0:
                for thr in (TH, TH2):
                    xs = (thr - a) / b
                    if x0 < xs < x1:
                        subs.add(xs)
        subs = sorted(subs)
        for s in range(len(subs) - 1):
            u0, u1 = subs[s], subs[s + 1]; um = (u0 + u1) / 2
            # E[W] contribution: sum over gaps with g(um)>1/7 of int (g-1/7)
            Aw = F(0); Bw = F(0)
            # near contribution: sum over gaps of int h(g); h depends on band of g(um)
            # h(g) linear/quadratic in g=a+bx; integrate over [u0,u1]
            for (a, b) in gaps:
                gm = a + b * um
                if gm > TH:
                    Aw += (a - TH); Bw += b
                # near: h(g)
                if gm >= TH2:
                    # h = g/7 - 3/98 = (a+bx)/7 - 3/98  -> linear
                    A = a / 7 - F(3, 98); B = b / 7
                    NEAR += A * (u1 - u0) + B * (u1 * u1 - u0 * u0) / 2
                elif gm > TH:
                    # h = (g-1/7)^2/2 = (a-1/7 + bx)^2/2 -> quadratic
                    A = a - TH
                    # (A + b x)^2 / 2 ; integrate
                    NEAR += (A * A * (u1 - u0) + A * b * (u1 * u1 - u0 * u0)
                             + b * b * (u1**3 - u0**3) / 3) / 2
                # gm <= 1/7: h=0
            EW += Aw * (u1 - u0) + Bw * (u1 * u1 - u0 * u0) / 2
    return EW, 2 * NEAR


def pz_tail_lb(E):
    """Rigorous PZ lower bound assuming far <= E[W]^2 (spread): E[W]^2/(near + E[W]^2)."""
    EW, near = EW_and_near(E)
    return EW * EW / (near + EW * EW), EW, near


def gcd_all(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g


def main():
    K = 11; bar = BAR[K]
    print("=" * 96)
    print(f"k={K} TAIL: sharp near/far PZ lower bound = E[W]^2/(near + E[W]^2), far<=E[W]^2 (spread)")
    print(f"  bar = {bar} = {float(bar):.6f};  crude near=(2/7)E[W] gave PZ~0.31 (mac-mini, just under)")
    print("=" * 96)
    # (0) iid-limit anchor (fully decorrelated): q(L)=(1-L)^k, near=(1/6)[(6/7)^12-(5/7)^12],
    #     E[W]=(6/7)^11, far=(5/7)^12
    from fractions import Fraction as Fr
    EWiid = Fr(6, 7)**11
    neariid = Fr(1, 6) * (Fr(6, 7)**12 - Fr(5, 7)**12)
    fariid = Fr(5, 7)**12
    pz_iid_sharp = EWiid**2 / (neariid + EWiid**2)   # far<=E[W]^2 used
    pz_iid_exact = EWiid**2 / (neariid + fariid)      # exact iid far
    print(f"  IID limit: E[W]={float(EWiid):.5f} near={float(neariid):.5f} far={float(fariid):.5f}"
          f" (far<=E[W]^2={float(EWiid**2):.5f}: {fariid<=EWiid**2})")
    print(f"    PZ_tail(far<=E[W]^2) = {float(pz_iid_sharp):.5f}   PZ_tail(exact iid far) = {float(pz_iid_exact):.5f}"
          f"   {'CLEARS' if pz_iid_sharp>bar else 'below'} bar {float(bar):.4f}")
    print()
    # (1) exact per-family near/far PZ tail over spread shapes, min per diameter band
    print("  min PZ_tail over spread k=11 families (exhaustive small diam, sampled large):")
    import random
    for D in (12, 14, 16, 18, 20, 25, 30):
        mn = None; arg = None; cnt = 0
        if D <= 16:
            it = ([0] + list(m) + [D] for m in combinations(range(1, D), K - 2))
        else:
            rng = random.Random(D)
            it = ([0] + sorted(rng.sample(range(1, D), K - 2)) + [D] for _ in range(300))
        for E in it:
            E = sorted(set(E))
            if len(E) != K or gcd_all([E[i+1]-E[i] for i in range(len(E)-1)]) != 1:
                continue
            cnt += 1
            pz, EW, near = pz_tail_lb(E)
            if mn is None or pz < mn:
                mn = pz; arg = (E, EW, near)
        if mn is None:
            continue
        E, EW, near = arg
        print(f"    diam {D:3d}: {cnt:6d} shapes, min PZ_tail = {float(mn):.5f}  "
              f"(E[W]={float(EW):.4f} near={float(near):.4f})  "
              f"{'CLEARS' if mn>=bar else '*** BELOW ***'}  argmin {E}")
    print()
    print("  NOTE: PZ_tail assumes far <= E[W]^2 (mac-mini's one decorrelation lemma, holds for")
    print("  spread families).  If min PZ_tail >= bar for all diam >= D*, then [exhaustive")
    print("  compact diam <= D*] + [this tail] CLOSES k=11.  The sharp near (vs crude (2/7)E[W])")
    print("  is the +0.02 mac-mini flagged as needed.")


if __name__ == "__main__":
    main()

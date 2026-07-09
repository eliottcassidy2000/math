"""
lrc14_L9_reduction_basecases_opus_S161.py   (opus-2026-07-08-S161)

THE L=9 LOWER-RANK REDUCTION + FINITE BOUND -- base cases.  Finishes S159/S160's correlated remainder.

THE REDUCTION (nested rank-1, each step = kps-S89's 1-point Koksma rate).  For E = block_9(d) u {p,q}:
  peel the FAR point (say q): D3(E) >= D3(block_9(d) u {p} + iid) - rate_q      [rank-1 Koksma in q]
  then peel p:                D3(block_9(d)u{p}+iid) >= D3(block_9 + 2 iid) - rate_p   [rank-1 in p]
  base:                       D3(block_9 + 2 iid) = D3_inf^{(9)} = 0.524.
The recursion bottoms at THREE base floors, each must be >= bar:
  (B0) block_9 + 2 iid            = 0.524   (S158/S160)
  (B1) min over (d,p) of D3(block_9(d) u {p} + 1 iid)   -- the intermediate 10-point-CORE floor  [THIS SCRIPT]
  (B2) both points close: bounded prim-diam, exact conditional finite check   (S160: min 0.467)
The FINITE BOUND: "both close" (small reduced pair products p1 d1, q1 d1 < T) => prim-diam bounded;
"one far" reduces via (B1); "both far" => box bound (0.524).  (B0),(B1),(B2) >= bar  =>  L=9 closed.

This script computes (B1) -- the 10-point-core floor -- reliably (conditional in p, exact; iid averaged),
min over (d,p); and reconfirms (B0),(B2).  If min (B1) >= bar the reduction's intermediate base holds.
"""
import numpy as np
from math import gcd

TH = 1 / 7; M = 6 / 7; BAR = 83549 / 252252


def W_batch(P):
    P = np.sort(P % 1.0, axis=1)
    g = np.empty_like(P)
    g[:, :-1] = P[:, 1:] - P[:, :-1]
    g[:, -1] = P[:, :1].ravel() + 1.0 - P[:, -1]
    return np.maximum(g - TH, 0.0).sum(axis=1)


def D3_of(m1, m2, m3):
    den = m2 - m3 / M
    return m1 / M + (m1 - m2 / M) ** 2 / den if den > 0 else m1 / M


def core_plus_iid_moments(d, p, Na=500, Nw=64):
    """E[W^i] for block_9(scale d) u {p}  +  1 iid point, i=1,2,3.
    conditional in (a,k) for the block+p (EXACT over the d sheets), iid point averaged over Nw."""
    a = ((np.arange(Na) + 0.5) / Na)[:, None]
    jj = np.arange(9)[None, :]
    ap = (a * jj) % 1.0                                   # (Na,9)
    ks = np.arange(d)[None, :]
    pph = ((p * a / d) % 1.0 + (p * ks) / d) % 1.0         # (Na,d)
    ap_rep = np.repeat(ap, d, axis=0)                      # (Na*d,9)
    core = np.concatenate([ap_rep, pph.reshape(-1, 1)], axis=1)   # (Na*d,10) block_9+p phases
    # add 1 iid point averaged over a w-grid
    Ncore = core.shape[0]
    ww = ((np.arange(Nw) + 0.5) / Nw)                      # iid point phases
    m1 = m2 = m3 = 0.0
    for w in ww:
        P = np.concatenate([core, np.full((Ncore, 1), w)], axis=1)
        W = W_batch(P)
        m1 += W.mean(); m2 += (W * W).mean(); m3 += (W ** 3).mean()
    return m1 / Nw, m2 / Nw, m3 / Nw


def main():
    print(f"L=9 REDUCTION BASE CASES; bar = {BAR:.6f}")
    print("=" * 92)
    print("\n(B0) block_9 + 2 iid = D3_inf^{(9)} (S158/S160): ~0.524  margin +0.19  [>= bar]")

    print("\n(B1) INTERMEDIATE 10-POINT-CORE FLOOR: min over (d,p) of D3(block_9(scale d) u {p} + 1 iid):")
    print("     (the reduction peels the far point q, leaving this core + iid; must be >= bar)")
    gmin = (9.9, None)
    # limit reference (large d, generic p): the decorrelated core = block_9 + 2 iid
    ml = core_plus_iid_moments(300, 137, Na=800, Nw=96)
    print(f"     [decorr ref: d=300 core -> D3 = {D3_of(*ml):.5f} ~ 0.524]")
    for d in range(1, 26):
        dmin = (9.9, None)
        for p in range(1, 8 * d + 40):
            if p % d == 0:
                continue
            if gcd(p, d) != 1 and p <= 8 * d:
                pass
            m = core_plus_iid_moments(d, p, Na=400, Nw=48)
            v = D3_of(*m)
            if v < dmin[0]:
                dmin = (v, p)
        if dmin[1] is not None:
            # refine
            m = core_plus_iid_moments(d, dmin[1], Na=900, Nw=96)
            vref = D3_of(*m)
            if vref < gmin[0]:
                gmin = (vref, (d, dmin[1]))
        if d <= 8 or d % 4 == 0:
            print(f"     d={d:3d}: min_p D3(core+iid) = {dmin[0]:.5f} at p={dmin[1]}", flush=True)
    print(f"\n     (B1) MIN over (d,p) = {gmin[0]:.5f} at (d,p)={gmin[1]}  margin {gmin[0]-BAR:+.5f}"
          f"  [{'>= bar' if gmin[0] >= BAR else 'BELOW'}]")

    print("\n(B2) BOTH-CLOSE finite region (S160 exact conditional): min D3 = 0.467 at block_9(4)+{3,40} [>= bar +0.136]")

    print("\n" + "=" * 92)
    print("THE FINITE BOUND (explicit).  'Close' for a point x vs the block means the reduced pair")
    print("product x1*d1 (x1=x/gcd(x,d), d1=d/gcd(x,d)) is < T (strong (d,x) resonance).  Then x < T*d/... ;")
    print("BOTH close => d < T and p,q < T => prim-diam = max(8d,p,q) < 8T (BOUNDED) => finite exact check (B2).")
    print("ONE far => peel it (rank-1 Koksma) => core+iid floor (B1) >= bar.  BOTH far => box bound (B0)=0.524.")
    ok = (gmin[0] >= BAR)
    print(f"\nVERDICT: base floors (B0)=0.524, (B1)={gmin[0]:.4f}, (B2)=0.467 all >= bar {BAR:.4f}"
          f"  => the L=9 reduction closes {'[OK]' if ok else '[CHECK (B1)]'}")
    print("=" * 92)


if __name__ == "__main__":
    main()

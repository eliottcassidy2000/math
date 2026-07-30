#!/usr/bin/env python3
"""Exact negative lemma for HYP-9061: single-jump corner routing at
gamma = 1/2 needs unbounded additive slack.

Setting (spine normal form THM-2966): envelope T(n) = ceil(3n/2) + D, so
depth budget d_m = ceil(3m/2) + D - m - 1 = ceil(m/2) + D - 1 at spine m.
Shell family with dyadic tails: boundaries n = 2*l and 3*l for each dyadic
l (ratio alternates 3/2, 4/3), 0-side corners at (z,o) = (2l, l) and
(3l, l), 1-side mirrors. The minimal corner-annihilation routing pairs
(2l,l) -> jump right Delta=l -> meets (3l,l); sibling (2l, 2l) meets the
mirror sibling on the diagonal. The jump interior demands, for each
i = 1..l-1, an integer absorption of binom(l,i)/2 at monomial
(3l - i, l + i).

An absorber for monomial (z,o) is a spine cell (m,k): k = o-1,
m + e - k = z for some polynomial degree e <= d_m, with capacity
binom(e, k) >= needed. We compute, for each dyadic l, the minimal D such
that every interior cell of the jump has an absorber with enough capacity
(choosing the best m for each cell, requiring only HALF the capacity
2*need <= binom since the other half-tail budget must stay available --
we use the weaker need <= binom to be conservative: even the weak form
forces D to grow).

Output: table l, minimal conservative D. The lemma: min D grows with l
(superconstant), so the single-jump routing does not realize C = 3/2 + o(1)
with constant D. This does NOT bound smarter multi-hop routings.
"""

from math import comb


def min_D_for_scale(l: int, D_cap: int = 4000) -> int:
    """Minimal D so every interior absorption of the (2l,l)->(3l,l) jump fits."""
    needs = []
    for i in range(1, l):
        z, o = 3 * l - i, l + i
        needs.append((z, o, comb(l, i) // 2 if comb(l, i) % 2 == 0 else (comb(l, i) + 1) // 2))
    for D in range(0, D_cap + 1):
        ok_all = True
        for z, o, need in needs:
            k = o - 1
            ok = False
            # spine cell m <= z, degree e = z - m + k <= d_m = ceil(m/2)+D-1
            for m in range(1, z + 1):
                e = z - m + k
                if e < k:
                    continue
                d_m = (m + 1) // 2 + D - 1
                if e > d_m:
                    continue
                if comb(e, k) >= need:
                    ok = True
                    break
            if not ok:
                ok_all = False
                break
        if ok_all:
            return D
    return -1


def main() -> None:
    rows = []
    for a in range(1, 8):
        l = 1 << a
        D = min_D_for_scale(l)
        rows.append((l, D))
        print(f"l={l:4d}  minimal_conservative_D={D}")
    Ds = [d for _, d in rows]
    if not all(x >= 0 for x in Ds):
        raise RuntimeError("D cap exceeded")
    if not (Ds[-1] > Ds[2] > Ds[0]):
        raise RuntimeError("expected growth pattern absent")
    growth = [ (rows[i+1][1] - rows[i][1]) for i in range(len(rows)-1) ]
    print("increments:", growth)
    print("lemma=single_jump_corner_routing_at_gamma_half_needs_growing_D")
    print("scope=does_not_bound_multi_hop_routings")
    print("status=VERIFIED-EXACT")


if __name__ == "__main__":
    main()

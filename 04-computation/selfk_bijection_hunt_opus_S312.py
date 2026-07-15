# opus-2026-07-15-S312 -- HYP-6895: bijection-structure hunt for 2*selfK = SC.
# Frame: sigma realizes reversal (T(sigma t) = rho(rev(T(t)))), kappa = full flip.
# qf := {t : cls(kappa t) = cls(t)}. Law: #black qf = SC(n), n >= 5.
# Pair-counting gives the Aut-WEIGHTED count Sum_{qf t}|Aut| = M/n!.
# This census decides: weighted vs unweighted; and hunts the natural map
# {black qf tilings} -> {SC classes} by comparing every available invariant.
import sys
from collections import defaultdict
sys.path.insert(0, '04-computation')
from smith_diagram_of_the_metagraph_opus_S307 import build

for n in (5, 6, 7):
    m = n*(n-1)//2 - (n-1)
    B = build(n)
    cls_of, H_of, rcls, x_of = B['cls_of'], B['H_of'], B['rcls'], B['x_of']
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
    tidx = {t: i for i, t in enumerate(tiles)}
    sig = [tidx[(n+1-y, n+1-x)] for (x, y) in tiles]
    FULL = (1 << m) - 1
    def sig_t(t):
        s = 0
        for i in range(m):
            if (t >> i) & 1: s |= 1 << sig[i]
        return s
    fiber = defaultdict(int)
    for t in range(1 << m): fiber[cls_of[t]] += 1
    aut = {c: H_of[c] // fiber[c] for c in fiber}

    qf = [t for t in range(1 << m) if cls_of[t] == cls_of[t ^ FULL]]
    black = [t for t in qf if sig_t(t) != t]
    blue = [t for t in qf if sig_t(t) == t]
    SCs = sorted(c for c in fiber if rcls[c] == c)
    print(f"\n================= n={n}: qf={len(qf)} black={len(black)} "
          f"blue={len(blue)} SC={len(SCs)}  law black==SC: {len(black)==len(SCs)}")
    wtd = sum(aut[cls_of[t]] for t in qf)
    wtd_black = sum(aut[cls_of[t]] for t in black)
    print(f"WEIGHTED: total={wtd}, black={wtd_black} "
          f"(if all-Aut-1 these equal {len(qf)}, {len(black)})")

    # carrier detail (black)
    per = defaultdict(int)
    for t in black: per[cls_of[t]] += 1
    print(f"black carriers: {len(per)} classes")
    for c in sorted(per):
        print(f"   class {c}: qfblack={per[c]}, H={H_of[c]}, Aut={aut[c]}, "
              f"x={x_of[c]}, {'SC' if rcls[c]==c else 'NS(rev=%d)'%rcls[c]}")
    # SC class list with invariants
    print(f"SC classes ({len(SCs)}):")
    for c in SCs:
        print(f"   class {c}: H={H_of[c]}, Aut={aut[c]}, x={x_of[c]}, fiber={fiber[c]}")

    # candidate map: t (black qf) -> cls(t XOR path-flip)?? path arcs aren't tiles;
    # instead try natural derived classes: cls(sig t)=rev cls(t) (known),
    # cls of the HALF-fold?? try: multiset of x-values match?
    xs_black = sorted(x_of[cls_of[t]] for t in black)
    xs_sc = sorted(x_of[c] for c in SCs)
    print(f"x-multiset match (black qf tilings vs SC classes): {xs_black == xs_sc}")
    print(f"   black x: {xs_black}")
    print(f"   SC    x: {xs_sc}")
    hs_black = sorted(H_of[cls_of[t]] for t in black)
    hs_sc = sorted(H_of[c] for c in SCs)
    print(f"H-multiset match: {hs_black == hs_sc}")
    print(f"   black H: {hs_black}")
    print(f"   SC    H: {hs_sc}")

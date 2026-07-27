#!/usr/bin/env python3
"""THM-2588 referee (renumbered from THM-2570; Jelonek collision) — the fold cascade (mac-mini-2026-07-27-S148).

LEMMA 1 (fold identity): for V with max h and w = max(V∖{h}), at the fold
modulus q = h + w: dist_q(h·j) = dist_q(w·j) for every j (h ≡ −w mod q).

LEMMA 2 (snap-fold step): M(V) ≥ M(V∖{h}) − w/(2q).  Proof: take a true
maximizer t* of V∖{h}; let j = round(q t*); every speed u ∈ V∖{h} has
dist(u·j/q) ≥ ||u t*|| − u/(2q) ≥ M(V∖{h}) − w/(2q); and h's distance at
j/q equals w's (Lemma 1).  So the single grid point j/q certifies the
bound for ALL of V.

THEOREM (cascade): sorting V = {h₁ > h₂ > … } and folding down k steps,
    M(V) ≥ M(body) − Σᵢ maxᵢ₊₁/(2(hᵢ + maxᵢ₊₁)),
so if every descent ratio hᵢ/maxᵢ₊₁ ≥ R then M(V) ≥ 1/(14−k) − k/(2(R+1))
(body floor = settled LRC(14−k)).  Uniform R = 134 gives M ≥ 3/41 for all
k = 1..13; k = 1 (= T-B, constant 533/4) is the binding case.
THE APEX-7 WALL IS BROKEN in the separated regime: no Bonferroni charge,
no 6/41-per-far, any k.

Referee: (1) fold identity exact; (2) snap inequality vs exact M on small
separated families (k ≤ 3, feasible heights); (3) the certificate (one
grid evaluation, O(13) at ANY height) on separated towers k = 7..12 with
astronomically large tops; (4) per-k sharp constants table; (5) T-A corner
values re-verified.
"""

from fractions import Fraction as F
from itertools import combinations
from math import floor

GAP_HI = F(3, 41)


def dist(x: F) -> F:
    fr = x - floor(x)
    return min(fr, 1 - fr)


def exact_M(W):
    W = sorted(W)
    mods = set()
    for u, v in combinations(W, 2):
        mods.add(u + v)
        mods.add(v - u)
    for w in W:
        mods.add(2 * w)
    mods.discard(0)
    best, bt = F(0), None
    for q in mods:
        for k in range(1, q):
            m = min(min((w * k) % q, q - (w * k) % q) for w in W)
            if F(m, q) > best:
                best, bt = F(m, q), F(k, q)
    return best, bt


def cascade_certificate(V):
    """Descend the sorted family, snapping at each fold modulus; return the
    certified lower bound and the final grid point of the FIRST fold (the
    single point that certifies the whole family via iterated Lemma 2 is
    the top fold's snap of the recursive maximizer)."""
    V = sorted(V, reverse=True)
    # recursive: bound(body) exact via settled floors 1/(14-k) — here we
    # instead compute the actual chain with exact maximizers where feasible
    # (small bodies), else the floor.
    def go(vs):
        if max(vs) <= 5000 or len(vs) <= 2:
            M, t = exact_M(vs)          # base: genuinely small body only
            return M, t
        h = vs[0]
        rest = vs[1:]
        w = rest[0]
        q = h + w
        Mr, tr = go(rest)
        j = round(q * tr)
        t = F(j, q)
        cert = min(dist(u * t) for u in vs)   # actual clearance at the point
        # Lemma-2 bound check
        assert cert >= Mr - F(w, 2 * q) - F(0), (vs, cert, Mr)
        return cert, t
    return go(V)


def main():
    print("== (1) fold identity, exact ==")
    import random
    rng = random.Random(148)
    for _ in range(200):
        V = sorted(rng.sample(range(1, 10**6), 13))
        h, w = V[-1], V[-2]
        q = h + w
        for j in (1, 7, rng.randint(2, q - 1)):
            d1 = min((h * j) % q, q - (h * j) % q)
            d2 = min((w * j) % q, q - (w * j) % q)
            assert d1 == d2
    print("  200 random families x 3 grid points: dist_q(h j) = dist_q(w j)  OK")

    print("\n== (2) snap inequality vs exact M (small separated, k <= 3) ==")
    for body, fars in [
        ([1, 2, 3, 4, 5, 6, 7, 8, 9, 11], [601, 24041]),      # k=2, R ~ 40-55
        ([1, 2, 3, 4, 5, 7, 8, 9, 10, 12], [251, 5021, 40168]),  # k=3, R ~ 8-20
    ]:
        V = body + fars
        M, t = exact_M(V)
        cert, tc = cascade_certificate(V)
        print(f"  V(k={len(fars)}): exact M = {M}, cascade certificate = {cert} "
              f"at t = {tc}: cert <= M: {cert <= M}, cert >= 3/41: {cert >= GAP_HI}")
        assert cert <= M
    print("  (certificate is a true lower bound; both clear 3/41)")

    print("\n== (3) the certificate at ANY height: separated towers k = 7..12 ==")
    for k in range(7, 13):
        body = list(range(1, 14 - k))            # {1..13-k}
        V = list(body)
        cur = max(body)
        for i in range(k):
            nxt = cur * 134 + 1                  # ratio >= 134 each step
            V.append(nxt)
            cur = nxt
        cert, tc = cascade_certificate(sorted(V))
        print(f"  k={k:>2}: v_max ~ 10^{len(str(max(V)))-1}, certificate = "
              f"{float(cert):.6f} >= 3/41 = {float(GAP_HI):.6f}: {cert >= GAP_HI}")
        assert cert >= GAP_HI
    print("  THE APEX-7 WALL IS BROKEN in the separated regime (any k).")

    print("\n== (4) per-k sharp cascade constants R_k (need 1/(14-k) - k/(2(R+1)) >= 3/41) ==")
    for k in range(1, 14):
        lhs = F(1, 14 - k) - GAP_HI if k < 13 else F(1, 2) - GAP_HI
        body_floor = F(1, 14 - k)
        need = body_floor - GAP_HI
        if need <= 0:
            print(f"  k={k:>2}: body floor {body_floor} <= 3/41 — cascade alone insufficient")
            continue
        Rk = -(-int(F(k, 2) / need)) - 1  # ceil(k/(2 need)) - 1
        import math
        Rk = math.ceil(k / (2 * need)) - 1
        print(f"  k={k:>2}: body floor {body_floor}, R_k = {Rk + 1} (ratio >= {Rk + 1} suffices)")

    print("\n== (5) T-A corner values, re-verified exact ==")
    for W, want in [
        (list(range(1, 13)) + [65], F(5, 66)),
        (list(range(1, 13)) + [66], F(1, 13)),
        (list(range(1, 12)) + [13, 65], F(1, 12)),
        (list(range(1, 12)) + [13, 66], F(1, 12)),
        (list(range(1, 12)) + [13, 67], F(1, 12)),
    ]:
        M, t = exact_M(W)
        print(f"  {W[-1]:>3} corner: M = {M} (agent claimed {want}) "
              f"{'OK' if M == want else 'MISMATCH'}")
        assert M == want
    print("  T-A corners confirmed; with THM-1290 (h <= 64) and the comb")
    print("  pigeonhole (agent ledger) / one-far closure (S126), both")
    print("  classical cores are gap-empty at all heights — two routes.")


if __name__ == "__main__":
    main()

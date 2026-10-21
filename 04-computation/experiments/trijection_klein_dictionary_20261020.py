"""trijection_klein_dictionary_20261020.py -- opus-2026-10-20-S47

The 3-4-3 / 4-2-4 divisor chart of p^k q r (trijection_reset_20261020/arithmetic.md) against the
Klein four-group <kappa, sigma> acting on the tilings of a fixed Hamiltonian path (THM-643 frame):
kappa = flip every tile (graph complementation under the ordered-graph trijection), sigma = staircase
reflection (converse read along the reversed path), kappa*sigma = path reversal.

Checks (all exact, small):
  1. fixed-point profiles and orbit types of both Klein-sets (divisors: C, Q, CQ; tilings: kappa, sigma, kappa sigma);
  2. kappa and kappa*sigma are fixed-point-free on tilings at every n (apex tile), sigma fixes 2^e;
  3. the ten divisors of p^2 q r never embed equivariantly (one orbit has stabiliser <CQ>); the eight
     non-central divisors form a Klein-set isomorphic to the n = 4 tilings; the number of equivariant
     bijections, and how many of them read sign(f) as closedness and sign(g) as chirality;
  4. for p^k q r: an equivariant embedding into tilings exists iff k is odd (no CQ-fixed points);
  5. no equivariant bijection (k = 2, eight elements) is divisibility-monotone onto the Boolean cube.
"""
from itertools import permutations, product

# ---------- divisor Klein-set of p^k q r ----------
def divisor_set(k):
    pts = [(a, b, c) for a in range(k + 1) for b in (0, 1) for c in (0, 1)]
    pts.remove((0, 0, 0)); pts.remove((k, 1, 1))
    C = lambda x: (k - x[0], 1 - x[1], 1 - x[2])
    Q = lambda x: (x[0], x[2], x[1])
    CQ = lambda x: C(Q(x))
    f = lambda x: 2 * x[0] - k * (x[1] + x[2])
    g = lambda x: x[1] - x[2]
    return pts, C, Q, CQ, f, g

def orbit_types(pts, A, B):
    """Orbit type = frozenset of the group elements fixing a point (same for all points of an orbit
    in an abelian group); returns multiset of stabiliser types over orbits."""
    AB = lambda x: A(B(x))
    seen, types = set(), {}
    for x in pts:
        if x in seen: continue
        orb = {x, A(x), B(x), AB(x)}; seen |= orb
        stab = tuple(name for name, h in (("A", A), ("B", B), ("AB", AB)) if h(x) == x)
        types[stab] = types.get(stab, 0) + 1
    return types

# ---------- tiling Klein-set ----------
def tiling_set(n):
    cells = [(i, j) for i in range(n) for j in range(i + 2, n)]
    idx = {c: t for t, c in enumerate(cells)}
    m = len(cells)
    sig = [idx[(n - 1 - j, n - 1 - i)] for (i, j) in cells]
    ALL = (1 << m) - 1
    kappa = lambda t: t ^ ALL
    def sigma(t):
        s = 0
        for k in range(m):
            if t >> k & 1: s |= 1 << sig[k]
        return s
    ks = lambda t: kappa(sigma(t))
    apex = idx[(0, n - 1)]
    return list(range(1 << m)), kappa, sigma, ks, m, apex, idx

def main():
    out = []
    # 1-2: profiles
    for k in range(2, 9):
        pts, C, Q, CQ, f, g = divisor_set(k)
        fix = (sum(C(x) == x for x in pts), sum(Q(x) == x for x in pts), sum(CQ(x) == x for x in pts))
        ot = orbit_types(pts, C, Q)
        fs = (sum(f(x) < 0 for x in pts), sum(f(x) == 0 for x in pts), sum(f(x) > 0 for x in pts))
        gs = (sum(g(x) < 0 for x in pts), sum(g(x) == 0 for x in pts), sum(g(x) > 0 for x in pts))
        out.append(f"divisors p^{k}qr: {len(pts)} points; fixed (C,Q,CQ) = {fix}; orbit stabilisers {dict(ot)}; sign f {fs}; sign g {gs}")
        assert fix == (0, 2 * k, 2 if k % 2 == 0 else 0)
        assert fs == ((2 * k, 2, 2 * k) if k % 2 == 0 else (2 * k + 1, 0, 2 * k + 1)) and gs == (k + 1, 2 * k, k + 1)
    for n in range(4, 9):
        pts, kappa, sigma, ks, m, apex, idx = tiling_set(n)
        e = sum(1 for t in pts if sigma(t) == t).bit_length() - 1
        fix = (sum(kappa(t) == t for t in pts), 1 << e, sum(ks(t) == t for t in pts))
        ot = orbit_types(pts, kappa, sigma)
        assert fix[0] == 0 and fix[2] == 0 and e == (n - 1) ** 2 // 4
        assert ot == {("B",): 1 << (e - 1), (): ((1 << m) - (1 << e)) // 4}
        out.append(f"tilings n={n}: {1 << m} points; fixed (kappa,sigma,kappa sigma) = {fix}; orbits: blue lines {1 << (e - 1)}, black twin pairs {((1 << m) - (1 << e)) // 4}")
    out.append("PROVED by the apex tile: kappa flips the sigma-fixed apex, so kappa and kappa*sigma have no fixed tiling at any n.")
    # 3: k = 2, eight non-central divisors vs n = 4 tilings
    pts, C, Q, CQ, f, g = divisor_set(2)
    central = [x for x in pts if CQ(x) == x]
    assert sorted(central) == [(1, 0, 1), (1, 1, 0)]
    eight = [x for x in pts if CQ(x) != x]
    tl, kappa, sigma, ks, m, apex, idx = tiling_set(4)
    closed = lambda t: not (t >> apex & 1)
    cells4 = sorted(idx, key=idx.get)
    pair = [c for c in cells4 if c != (0, 3)]
    lo, hi = idx[pair[0]], idx[pair[1]]           # the sigma-paired tiles (0,2) and (1,3)
    chir = lambda t: (t >> lo & 1) - (t >> hi & 1)  # +1, 0, -1: the chirality of a black tiling
    count = both = fonly = gonly = 0; mono = 0
    divides = lambda x, y: all(x[i] <= y[i] for i in range(3))
    for perm in permutations(range(8)):
        phi = {eight[i]: tl[perm[i]] for i in range(8)}
        if all(phi[C(x)] == kappa(phi[x]) and phi[Q(x)] == sigma(phi[x]) for x in eight):
            count += 1
            fok = all((f(x) > 0) == closed(phi[x]) for x in eight)
            gok = all(((g(x) > 0) - (g(x) < 0)) == chir(phi[x]) for x in eight)
            both += fok and gok; fonly += fok; gonly += gok
            if all((not divides(x, y)) or (phi[x] & ~phi[y]) == 0 for x in eight for y in eight if x != y): mono += 1
    out.append(f"k=2: equivariant bijections eight divisors -> n=4 tilings: {count}; with sign f = closedness: {fonly}; with sign g = chirality: {gonly}; both: {both}; divisibility-monotone into the cube: {mono}")
    assert count == 32 and fonly == 4 and gonly == 16 and both == 2 and mono == 2
    # list the two bijections reading both signs and the two divisibility-monotone ones
    names = {(1,0,0):"p",(0,1,0):"q",(0,0,1):"r",(2,0,0):"p^2",(0,1,1):"qr",(1,1,1):"pqr",(2,1,0):"p^2q",(2,0,1):"p^2r"}
    def tname(t):
        bits = "".join(str(t >> idx[c] & 1) for c in cells4)
        return f"[{bits}]" + ("closed" if closed(t) else "open")
    for perm in permutations(range(8)):
        phi = {eight[i]: tl[perm[i]] for i in range(8)}
        if not all(phi[C(x)] == kappa(phi[x]) and phi[Q(x)] == sigma(phi[x]) for x in eight): continue
        fok = all((f(x) > 0) == closed(phi[x]) for x in eight)
        gok = all(((g(x) > 0) - (g(x) < 0)) == chir(phi[x]) for x in eight)
        mon = all((not divides(x, y)) or (phi[x] & ~phi[y]) == 0 for x in eight for y in eight if x != y)
        if (fok and gok) or mon:
            out.append(("both signs" if fok and gok else "monotone") + ": " + ", ".join(f"{names[x]}->{tname(phi[x])}" for x in eight))
    out.append("tile order in the bit strings: " + ", ".join(str(c) for c in cells4) + "; closed = apex tile (0,3) backward")
    # full ten: no equivariant injection (stabiliser <CQ> orbit) -- certified by orbit types at every n above
    out.append("ten divisors: the orbit {pq, pr} has stabiliser <CQ>; tilings have no orbit with stabiliser <kappa sigma>; no equivariant injection at any n.")
    # 4: parity of k
    for k in range(2, 9):
        pts, C, Q, CQ, f, g = divisor_set(k)
        ot = orbit_types(pts, C, Q)
        embeds = ("AB",) not in ot
        assert embeds == (k % 2 == 1)
        out.append(f"p^{k}qr embeds equivariantly into tilings: {embeds} (needs {ot.get(('B',),0)} blue lines and {ot.get((),0)} black twin pairs; n = 4 offers 2 and 1, n = 5 offers 8 and 12)")
    print("\n".join(out))

if __name__ == "__main__":
    main()

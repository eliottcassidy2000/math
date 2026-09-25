"""Exact square-sum, signed Pell, and decimal mod17 bridge controls.

Finite universes and hostile controls are printed. Checks survive -O.
The all-index identities are proved in the companion note.
"""
from itertools import combinations
from math import isqrt


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def square_graph(n):
    adj = {i: [] for i in range(1, n+1)}
    for a, b in combinations(adj, 2):
        if isqrt(a+b)**2 == a+b:
            adj[a].append(b)
            adj[b].append(a)
    return adj


def square_controls():
    q14, q15 = square_graph(14), square_graph(15)
    need([v for v in q14 if len(q14[v]) == 1] == [8,9,10], "Q14 leaves")
    need([v for v in q15 if len(q15[v]) == 1] == [8,9], "Q15 leaves")
    found = []

    def search(path, used):
        if len(path) == 15:
            found.append(tuple(path))
            return
        for v in q15[path[-1]]:
            if v not in used:
                search(path + [v], used | {v})

    search([8], {8})
    expected = (8,1,15,10,6,3,13,12,4,5,11,14,2,7,9)
    need(found == [expected], "unique oriented Q15 path")
    need(9 not in q15[8], "endpoints do not form square edge")
    print("Q14 leaves8,9,10; Q15 unique path from8:", expected)
    print("path has14 edges, sums", tuple(a+b for a,b in zip(expected, expected[1:])))
    need(8*9//2 == 36 == 6**2 and 17**2-8*6**2 == 1, "square triangular seed")
    need(8**2+15**2 == 17**2, "8,15,17 triangle")
    # Primitive odd-root chart for the actual accelerated plus edge3->5.
    s,t = 3,5
    need((s*t,(s*s-t*t)//2,(s*s+t*t)//2) == (15,-8,17), "marked edge triangle")
    print("8,9 -> triangular36 -> Pell(17,6); edge3->5 -> marked triple(15,-8,17)")


def half_step(p, modulus=None):
    x,s = p
    result = (x+2*s,x+s)
    return tuple(v % modulus for v in result) if modulus else result


def pell_controls():
    p = (1,0)
    square_rows = []
    for j in range(25):
        x,s = p
        need(x*x-2*s*s == (-1)**j, "signed Pell norm")
        if j % 2 == 0:
            n,q = (x-1)//2,s//2
            need(x % 2 == 1 and s % 2 == 0 and n*(n+1)//2 == q*q, "triangular chart")
            square_rows.append((j//2,n,q,x))
        p = half_step(p)
    print("Pell square branch (index,n,sqrt(Tn),2n+1):", square_rows[:6])
    x,q = 1,0
    xs = [x]
    for ell in range(1,101):
        x,q = 3*x+8*q,x+3*q
        xs.append(x)
        need((x % 17 == 0) == (ell % 4 == 2), "exact17 Pell phase")
        odd_part = ell
        two_part = 1
        while odd_part % 2 == 0:
            odd_part //= 2
            two_part *= 2
        if odd_part > 1:
            need(x % xs[two_part] == 0 and 1 < xs[two_part] < x, "odd-index-factor divisibility")
    need(xs[6] == 19601 == 17*1153, "later composite Pell trace")
    print("100 Pell indices:17 divides trace iff index=2 mod4; odd-index-factor composite controls PASS")


def psi(a, sheet=1):
    z = (3*a+7) % 17
    need(z != 0 and sheet in (-1,1), "Pell chart domain")
    u,v = pow(z,9,17),sheet*pow(z,15,17) % 17
    return ((u+v)*pow(2,-1,17) % 17, (u-v)*pow(12,-1,17) % 17)


def finite_bridge_controls():
    allowed = set(range(17)) - {9}
    norm_union = {(x,s) for x in range(17) for s in range(17) if (x*x-2*s*s) % 17 in (1,16)}
    images = []
    for sheet in (1,-1):
        image = set()
        for a in allowed:
            p = psi(a,sheet)
            z = (3*a+7) % 17
            need(psi((10*a+21) % 17,sheet) == half_step(p,17), "commuting affine/Pell square")
            need((p[0]**2-2*p[1]**2) % 17 == sheet*pow(z,8,17) % 17, "norm character")
            # Inverse uses u=x+6s=z^9, and 9^2=1 modulo16.
            recovered_z = pow((p[0]+6*p[1]) % 17,9,17)
            need((recovered_z-7)*pow(3,-1,17) % 17 == a, "inverse chart")
            image.add(p)
        need(len(image) == 16, "sixteen distinct states")
        images.append(image)
    need(not images[0] & images[1] and images[0] | images[1] == norm_union, "two complete signed-norm orbits")
    for seed in ((1,0),(0,1)):
        p = seed
        for _ in range(8): p = half_step(p,17)
        need(p == tuple(-x % 17 for x in seed), "S8=-I")
        for _ in range(8): p = half_step(p,17)
        need(p == seed, "S16=I")
    a,p = -2,(1,0)
    for k in range(65):
        need(a == (10**k-7)//3, "integer decimal clock")
        need(psi(a) == p, "all sampled clock indices")
        need(((10**(k+8)-7)//3-(1-a)) % 17 == 0, "eight-step affine reflection")
        need((a % 17 == 0) == (k % 16 == 9), "decimal17 factor phase")
        if k in (1,4,8,9): print(f"k={k}: decimal mod17={a%17}, signed Pell state={p}")
        a,p = 10*a+21,half_step(p,17)
    print("full32 signed-norm states = two16-cycles; both decimal charts and inverses PASS")


def triple_reversal_control():
    counts = {0:0,2:0}
    for code in range(8):
        edges = tuple(combinations(range(3),2))
        degrees = [0,0,0]
        for bit,(a,b) in enumerate(edges):
            degrees[a if code >> bit & 1 else b] += 1
        cycles = int(degrees == [1,1,1])
        # Converse keeps a cycle cyclic and a transitive triple transitive.
        contribution = cycles + int([2-d for d in degrees] == [1,1,1])
        counts[contribution] += 1
    need(counts == {0:6,2:2}, "THM060 repaired TypeA control")
    print("all8 triple orientations: full reversal contributes0 in6 cases,2 in2 cases")


if __name__ == '__main__':
    square_controls()
    pell_controls()
    finite_bridge_controls()
    triple_reversal_control()
    print("PASS: finite controls and explicit scoped identities; no global primality or Collatz conclusion")

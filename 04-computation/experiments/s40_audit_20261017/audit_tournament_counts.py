"""S40 audit: T(n) (tournaments) and SC(n) (self-converse tournaments) for n <= 8 by Burnside over S_n acting on ordered pairs.
fix(sigma) = 2^(#orbit pairs) unless some orbit of pi(u,v) = (s u, s v) contains a reversed pair; for sigma(T) = T^op use
pi(u,v) = (s v, s u). Orbits of S_n x Z_2: (sum fix + sum fixconv)/(2 n!) = (T + SC)/2."""
import itertools, math
def count(n, conv):
    tot = 0
    pairs = [(u, v) for u in range(n) for v in range(n) if u != v]
    for s in itertools.permutations(range(n)):
        seen = set(); orbits = []
        for pr in pairs:
            if pr in seen: continue
            orb = []; x = pr
            while x not in seen:
                seen.add(x); orb.append(x)
                x = (s[x[1]], s[x[0]]) if conv else (s[x[0]], s[x[1]])
            orbits.append(frozenset(orb))
        ok = all((v, u) not in o for o in orbits for (u, v) in o)
        if ok:
            tot += 2 ** (len(orbits) // 2)
    return tot // math.factorial(n)
T = [count(n, False) for n in range(1, 9)]
Fc = [count(n, True) for n in range(1, 9)]
SC = [2 * ((t + f) // 1) // 2 - t + 0 for t, f in zip(T, Fc)]
# orbits of S_n x Z2 = (T + Fc)/2 and also = (T + SC)/2, hence SC = Fc
print('T(n), n=1..8:', T)
print('SC(n), n=1..8:', Fc)
print('solo fraction SC/T at n = 8: %.4f; metagraph nodes (T+SC)/2:' % (Fc[-1] / T[-1]), [(t + f) // 2 for t, f in zip(T, Fc)])

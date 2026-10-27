"""window_groups_g.py -- independent generalisation of window_groups.py to THM-4812 (g):
strong S with tc(S) = V (odd V <= 21) and |S| >= NMIN.  Every Aut(S)-orbit is K1/C3/R5/P7 (by (b), (c)),
prod tc(orbits) <= V, |G| >= f(n)/V, index [P:G] <= V|P|/f(n); index-3 subgroups = kernels of homomorphisms to Z_3
(evaluated element by element); index >= 5 is checked to be impossible (must divide |P| and be <= the bound).
Writes jobs to OUTDIR with prefix g_v{V}_."""
import itertools, os, sys
from window_groups import f_busch, build, elements_of_P, check_aut_of_copies, orbitals_between, transitive_projections

def configurations(V, nmin, nmax=80):
    out = []
    for n in range(nmin, nmax + 1):
        for p in range(0, 3):
            for r in range(0, 5):
                if 3 ** r * 9 ** p > V: continue
                for c in range(0, n // 3 + 1):
                    z = n - 3 * c - 5 * r - 7 * p
                    if z < 0: continue
                    Psize = 3 ** c * 5 ** r * 21 ** p
                    if V * Psize < f_busch(n): continue
                    out.append((n, z, c, r, p, Psize))
    return out

def main(outdir, nmin):
    os.makedirs(outdir, exist_ok=True)
    total = 0; refined_dead = 0; per_v = {}
    for V in range(1, 22, 2):
        cfgs = configurations(V, nmin)
        maxn = max((c[0] for c in cfgs), default=None)
        cases = 0
        for (n, z, c, r, p, Psize) in cfgs:
            orbits, arcs = build(n, z, c, r, p)
            els, z3 = elements_of_P(orbits, n)
            assert len(els) == Psize and check_aut_of_copies(orbits, arcs, els)
            bound = V * Psize / f_busch(n)
            prod_tc = 3 ** r * 9 ** p
            # any admissible index d: odd, d | |P|, d <= bound
            for d in range(5, int(bound) + 1, 2):
                assert Psize % d != 0, f"index {d} possible at V={V} n={n} cfg={(z,c,r,p)}: not handled"
            groups = [('P', [g for g, _ in els], 1)]
            if bound >= 3 and z3:
                seen = set()
                for phi in itertools.product(range(3), repeat=len(z3)):
                    if not any(phi): continue
                    lead = next(x for x in phi if x)
                    norm = tuple((x * pow(lead, -1, 3)) % 3 for x in phi)
                    if norm in seen: continue
                    seen.add(norm)
                    K = [g for g, chi in els if sum(a * b for a, b in zip(norm, chi)) % 3 == 0]
                    groups.append(('ker' + ''.join(map(str, norm)), K, 3))
            for name, G, index in groups:
                if len(G) * V < f_busch(n): continue
                if not transitive_projections(orbits, G): continue
                cases += 1
                if index * prod_tc > V: refined_dead += 1
                pairs, module = orbitals_between(orbits, G, n)
                fn = os.path.join(outdir, f"g_v{V}_n{n}_z{z}_c{c}_r{r}_p{p}_{name}.job")
                with open(fn, 'w') as fh:
                    fh.write(f"{n} {len(orbits)} {len(pairs)} {len(G)} {index} {'module' if module else 'nonmodule'}\n")
                    fh.write(' '.join(f"{k}:" + '.'.join(map(str, vs)) for k, vs in orbits) + "\n")
                    fh.write(' '.join(f"{u},{v}" for u, v in arcs) + "\n")
                    for orb in pairs:
                        fh.write(' '.join(f"{u},{v}" for u, v in orb) + "\n")
        per_v[V] = (cases, maxn)
        total += cases
    print("per V: (cases, max order with a configuration):", per_v)
    print(f"total cases {total}; cases excluded a priori by index*prod tc > V: {refined_dead}")

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))

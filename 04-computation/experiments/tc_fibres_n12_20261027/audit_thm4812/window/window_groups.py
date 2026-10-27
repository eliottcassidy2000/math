"""window_groups.py -- independent audit construction of the THM-4812 (d) window cases.

Independent of window_configs.py: configurations are derived from the inequality |P| >= f(n)/21 with my own f;
P is built as the explicit set of its elements (tuples of local affine maps), NOT by generator closure; index-3
subgroups are found as kernels of every homomorphism P -> Z_3, by evaluating the homomorphism on every element;
standard orbit copies differ from the repository's (C3 oriented 0->2->1->0, P7 = non-residue Paley tournament
{3,5,6}, R5 = {1,2}); orbit order on the vertex set is K1's, C3's, R5's, P7 (the repository uses the reverse).

Output: one job file per (configuration, subgroup):
  line 1: n s K |G| index kind   (s = number of orbits, K = number of free orbital pairs)
  line 2: orbit types and vertex lists  e.g. K1:0 C3:1.2.3 P7:4.5.6.7.8.9.10
  line 3: fixed within-orbit arcs u,v
  then K lines: one orbital (list of arcs u,v) per pair of mutually inverse between-orbit orbitals.
Also prints the derived bound index <= 21|P|/f(n), and the refined bound index * prod tc(orbit) <= 21.
"""
import itertools, os, sys

def f_busch(n):
    q, r = divmod(n, 3)
    if r == 0: return 3 * 5 ** (q - 1)
    if r == 1: return 5 ** q
    return 9 * 5 ** (q - 1)

TC = {'K1': 1, 'C3': 1, 'R5': 3, 'P7': 9}
AUT = {'K1': 1, 'C3': 3, 'R5': 5, 'P7': 21}
QR7 = [1, 2, 4]
NQR7 = [3, 5, 6]

def configurations(tcval=21, orders=None):
    out = []
    for n in orders:
        for p in range(0, 3):
            for r in range(0, 5):
                if 3 ** r * 9 ** p > tcval:          # product lemma
                    continue
                for c in range(0, n // 3 + 1):
                    z = n - 3 * c - 5 * r - 7 * p
                    if z < 0: continue
                    Psize = 3 ** c * 5 ** r * 21 ** p
                    if tcval * Psize < f_busch(n):
                        continue
                    out.append((n, z, c, r, p, Psize))
    return out

def local_elements(kind):
    """elements of Aut(standard copy) as (label, map list) on local coordinates 0..m-1, plus the Z_3-character exponent."""
    if kind == 'K1':
        return [((), [0], 0)]
    if kind == 'C3':
        return [((e,), [(x + e) % 3 for x in range(3)], e) for e in range(3)]
    if kind == 'R5':
        return [((e,), [(x + e) % 5 for x in range(5)], 0) for e in range(5)]
    if kind == 'P7':
        els = []
        for ell, a in enumerate([1, 2, 4]):          # a = 2^ell
            for b in range(7):
                els.append(((a, b), [(a * x + b) % 7 for x in range(7)], ell))
        return els
    raise ValueError

def local_arcs(kind):
    if kind == 'K1': return []
    if kind == 'C3': return [(x, (x - 1) % 3) for x in range(3)]            # 0->2->1->0
    if kind == 'R5': return [(x, (x + d) % 5) for x in range(5) for d in (1, 2)]
    if kind == 'P7': return [(x, (x + d) % 7) for x in range(7) for d in NQR7]
    raise ValueError

def build(n, z, c, r, p):
    kinds = ['K1'] * z + ['C3'] * c + ['R5'] * r + ['P7'] * p
    orbits = []; pos = 0
    for k in kinds:
        m = {'K1': 1, 'C3': 3, 'R5': 5, 'P7': 7}[k]
        orbits.append((k, list(range(pos, pos + m)))); pos += m
    assert pos == n
    arcs = []
    for k, vs in orbits:
        for (a, b) in local_arcs(k):
            arcs.append((vs[a], vs[b]))
    return orbits, arcs

def elements_of_P(orbits, n):
    """all elements of P = prod Aut(orbit copy), as (permutation tuple, character vector over Z_3-coordinates)."""
    per_orbit = [local_elements(k) for k, _ in orbits]
    z3coords = [i for i, (k, _) in enumerate(orbits) if k in ('C3', 'P7')]
    els = []
    for combo in itertools.product(*per_orbit):
        perm = [None] * n
        for (k, vs), (_, mp, _) in zip(orbits, combo):
            for x, y in enumerate(mp):
                perm[vs[x]] = vs[y]
        chi = tuple(combo[i][2] for i in z3coords)
        els.append((tuple(perm), chi))
    return els, z3coords

def check_aut_of_copies(orbits, arcs, els):
    arcset = set(arcs)
    for perm, _ in els:
        for (u, v) in arcs:
            if (perm[u], perm[v]) not in arcset:
                return False
    return True

def orbitals_between(orbits, G, n):
    orbit_index = {}
    for i, (_, vs) in enumerate(orbits):
        for v in vs: orbit_index[v] = i
    assigned = {}
    orbitals = []
    for u in range(n):
        for v in range(n):
            if u == v or orbit_index[u] == orbit_index[v] or (u, v) in assigned: continue
            orb = sorted({(g[u], g[v]) for g in G})
            for pr in orb: assigned[pr] = len(orbitals)
            orbitals.append(orb)
    pairs = []; done = set()
    for t, orb in enumerate(orbitals):
        if t in done: continue
        u, v = orb[0]
        ti = assigned[(v, u)]
        assert ti != t
        done.add(t); done.add(ti)
        pairs.append(orb)
    # module test: every orbital pair is a full rectangle
    sizes = {i: len(vs) for i, (_, vs) in enumerate(orbits)}
    module = all(len(orb) == sizes[orbit_index[orb[0][0]]] * sizes[orbit_index[orb[0][1]]] for orb in pairs)
    return pairs, module

def transitive_projections(orbits, G):
    for k, vs in orbits:
        if len({g[vs[0]] for g in G}) != len(vs):
            return False
    return True

def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    orders = list(range(13, 23))
    cfgs = configurations(21, orders)
    # tc = 7 check for all n >= 13 (should be empty)
    cfg7 = configurations(7, list(range(13, 60)))
    print(f"tc = 7 configurations with 7|P| >= f(n), 13 <= n < 60: {cfg7}")
    print(f"tc = 21 configurations at orders 13..22: {len(cfgs)}")
    jobs = 0; ncase = 0; idx3_refined_dead = 0
    summary = []
    for (n, z, c, r, p, Psize) in cfgs:
        orbits, arcs = build(n, z, c, r, p)
        els, z3 = elements_of_P(orbits, n)
        assert len(els) == Psize
        assert check_aut_of_copies(orbits, arcs, els)
        bound = 21 * Psize / f_busch(n)
        prod_tc = 3 ** r * 9 ** p
        cand_groups = [('P', [g for g, _ in els], 1)]
        # every index-d subgroup must have d | |P| and d <= bound; odd d in {3,5} considered explicitly
        if bound >= 3 and len(z3) > 0:
            seen = set()
            for phi in itertools.product(range(3), repeat=len(z3)):
                if not any(phi): continue
                lead = next(x for x in phi if x)
                norm = tuple((x * pow(lead, -1, 3)) % 3 for x in phi)
                if norm in seen: continue
                seen.add(norm)
                K = [g for g, chi in els if sum(a * b for a, b in zip(norm, chi)) % 3 == 0]
                assert len(K) * 3 == Psize
                cand_groups.append(('ker' + ''.join(map(str, norm)), K, 3))
        if bound >= 5:
            assert Psize % 5 != 0, "index-5 subgroups possible: not handled"
        for name, G, index in cand_groups:
            if len(G) * 21 < f_busch(n): continue
            if not transitive_projections(orbits, G): continue
            ncase += 1
            if index * prod_tc > 21: idx3_refined_dead += 1
            pairs, module = orbitals_between(orbits, G, n)
            fn = os.path.join(outdir, f"a_n{n}_z{z}_c{c}_r{r}_p{p}_{name}.job")
            with open(fn, 'w') as fh:
                fh.write(f"{n} {len(orbits)} {len(pairs)} {len(G)} {index} {'module' if module else 'nonmodule'}\n")
                fh.write(' '.join(f"{k}:" + '.'.join(map(str, vs)) for k, vs in orbits) + "\n")
                fh.write(' '.join(f"{u},{v}" for u, v in arcs) + "\n")
                for orb in pairs:
                    fh.write(' '.join(f"{u},{v}" for u, v in orb) + "\n")
            jobs += 1
            summary.append((n, z, c, r, p, Psize, name, len(G), index, len(pairs), module, round(bound, 3), index * prod_tc))
    for s in summary:
        print("n=%d z=%d c=%d r=%d p=%d |P|=%d G=%s |G|=%d index=%d K=%d module=%s 21|P|/f=%.3f index*prod_tc=%d" % s)
    print(f"cases (config, subgroup) = {ncase}; jobs written = {jobs}; module cases = {sum(1 for s in summary if s[10])}; "
          f"non-module = {sum(1 for s in summary if not s[10])}; index-3 cases excluded a priori by index*prod tc > 21: {idx3_refined_dead}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "jobs")

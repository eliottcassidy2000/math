#!/usr/bin/env python3
"""gen_sym11.py -- symmetry-first candidate generator for the order-11 prime cycle-index census.

Every prime tournament of order 11 with a nontrivial automorphism has an automorphism of prime
order whose cycles have length p with at least two p-cycles (a lone p-cycle is a module) or a
single 11-cycle; so its automorphism group contains an element of cycle type 3^2 1^5, 3^3 1^2,
5^2 1 or 11.  Every self-converse tournament of order 11 has an anti-automorphism with exactly
one fixed point and all other cycles even: types 2^5 1, 4 2^3 1, 4^2 2 1, 6 2^2 1, 6 4 1, 8 2 1,
10 1.  We enumerate the invariant tournaments of each type (one orientation bit per orbit of the
cyclic group on pairs), with two reductions: the sub-tournament on a set of five vertices
permuted freely by the centraliser is taken from the 12 isomorphism classes, and 2-cycles (p q)
of an anti-automorphism are oriented p -> q (the swap commutes with psi).  The output is in
gentourng's text format (upper triangle, row by row, '1' at (i,j) means i -> j); duplicates are
removed downstream by canonical form.  Usage: gen_sym11.py TYPE  (or 'all')."""
import sys, subprocess
n = 11
def perm_from_cycles(cycles):
    p = list(range(n))
    for c in cycles:
        for i in range(len(c)): p[c[i]] = c[(i+1) % len(c)]
    return p
def pair_orbits(psi, anti):
    """orbits of <psi> on unordered pairs, each as a list of oriented pairs (x,y) meaning:
    if the representative is oriented rep_x -> rep_y then this member is oriented x -> y."""
    seen = {}; orbits = []
    for u in range(n):
        for v in range(u+1, n):
            if (u, v) in seen: continue
            orb = []; x, y = u, v
            while True:
                k = (min(x, y), max(x, y))
                if k in seen:
                    assert (x, y) == (u, v), "inconsistent orbit (odd stabiliser)"
                    break
                seen[k] = len(orbits); orb.append((x, y))
                x, y = (psi[y], psi[x]) if anti else (psi[x], psi[y])
            orbits.append(orb)
    return orbits, seen
def five_tournaments():
    out = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', '5'], capture_output=True, text=True).stdout.split()
    assert len(out) == 12
    res = []
    for line in out:
        arcs = {}; p = 0
        for i in range(5):
            for j in range(i+1, 5):
                arcs[(i, j)] = (line[p] == '1'); p += 1
        res.append(arcs)
    return res
def emit(orbits, fixed_orient, free_idx, fout):
    """fixed_orient: dict orbit_index -> bool (rep oriented x->y if True); enumerate free orbits."""
    m = len(free_idx); cnt = 0
    adjbase = [[0]*n for _ in range(n)]
    for oi, b in fixed_orient.items():
        for (x, y) in orbits[oi]:
            if b: adjbase[x][y] = 1
            else: adjbase[y][x] = 1
    for mask in range(1 << m):
        adj = [row[:] for row in adjbase]
        for t, oi in enumerate(free_idx):
            b = (mask >> t) & 1
            for (x, y) in orbits[oi]:
                if b: adj[x][y] = 1
                else: adj[y][x] = 1
        s = ''.join('1' if adj[i][j] else '0' for i in range(n) for j in range(i+1, n))
        fout.write(s + '\n'); cnt += 1
    return cnt
def run_type(name, fout):
    T5 = five_tournaments()
    if name == '3^2 1^5':
        psi = perm_from_cycles([[0,1,2],[3,4,5]]); orbits, where = pair_orbits(psi, False)
        F = [6,7,8,9,10]; fixed_orbs = [where[(F[i],F[j])] for i in range(5) for j in range(i+1,5)]
        free = [k for k in range(len(orbits)) if k not in fixed_orbs]; cnt = 0
        for arcs in T5:
            fo = {where[(F[i],F[j])]: arcs[(i,j)] for i in range(5) for j in range(i+1,5)}
            cnt += emit(orbits, fo, free, fout)
        return cnt
    if name == '3^3 1^2':
        psi = perm_from_cycles([[0,1,2],[3,4,5],[6,7,8]]); orbits, where = pair_orbits(psi, False)
        return emit(orbits, {}, list(range(len(orbits))), fout)
    if name == '5^2 1':
        psi = perm_from_cycles([[0,1,2,3,4],[5,6,7,8,9]]); orbits, where = pair_orbits(psi, False)
        return emit(orbits, {}, list(range(len(orbits))), fout)
    if name == '11':
        psi = perm_from_cycles([list(range(11))]); orbits, where = pair_orbits(psi, False)
        return emit(orbits, {}, list(range(len(orbits))), fout)
    if name == '2^5 1':
        psi = perm_from_cycles([[0,1],[2,3],[4,5],[6,7],[8,9]]); orbits, where = pair_orbits(psi, True)
        P = [0,2,4,6,8]; Q = [1,3,5,7,9]
        fixed_orbs = [where[(P[i],Q[i])] for i in range(5)] + [where[(P[i],P[j])] for i in range(5) for j in range(i+1,5)]
        free = [k for k in range(len(orbits)) if k not in fixed_orbs]; cnt = 0
        for arcs in T5:
            fo = {where[(P[i],Q[i])]: True for i in range(5)}
            for i in range(5):
                for j in range(i+1,5):
                    oi = where[(P[i],P[j])]; rep = orbits[oi][0]
                    # rep is (P[i],P[j]) oriented P[i]->P[j] iff arcs[(i,j)]
                    assert rep == (P[i], P[j])
                    fo[oi] = arcs[(i,j)]
            cnt += emit(orbits, fo, free, fout)
        return cnt
    if name in ('6 2^2 1', '10 1'):
        # the even-cycle types with a cycle of length 0 mod 4 (4 2^3 1, 4^2 2 1, 6 4 1, 8 2 1) are impossible:
        # a pair swapped by an even power of psi forces u->v and v->u (see anti_types_check.py)
        cyc = {'6 2^2 1': [[0,1,2,3,4,5],[6,7],[8,9]], '10 1': [[0,1,2,3,4,5,6,7,8,9]]}[name]
        psi = perm_from_cycles(cyc); orbits, where = pair_orbits(psi, True)
        fo = {}
        for c in cyc:
            if len(c) == 2:
                oi = where[(c[0], c[1])]; assert orbits[oi][0] == (c[0], c[1]); fo[oi] = True
        free = [k for k in range(len(orbits)) if k not in fo]
        return emit(orbits, fo, free, fout)
    raise SystemExit('unknown type ' + name)
if __name__ == '__main__':
    name = sys.argv[1]; outp = sys.argv[2]
    with open(outp, 'w') as fout:
        c = run_type(name, fout)
    print(f'{name}: {c} candidates -> {outp}')

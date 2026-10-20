#!/usr/bin/env python3
"""Independent recount of the C15 counterexample (c15_symmetry_counterexample_20261020): H = <(0 1 2)(3 4 5 6 7)> <= S_10
(orbits 3+5+1+1). tau_T(H) = isomorphism types of H-invariant tournaments, tau_E(H) = even types among H-invariant graphs.
Types via nauty labelg; evenness via dreadnaut automorphism generators (sgn is a homomorphism, so generators suffice).
Also: for comparison the same for H' = <(0 1 2)(3 4 5 6 7)> on 8 and 9 vertices and the graph K_{5,5} test."""
import subprocess, itertools, re, sys
def cyc(perm_cycles, n):
    g = list(range(n))
    for c in perm_cycles:
        for i in range(len(c)): g[c[i]] = c[(i + 1) % len(c)]
    return tuple(g)
def mul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def closure(gens, n):
    idt = tuple(range(n)); S = {idt}; fr = [idt]
    while fr:
        x = fr.pop()
        for g in gens:
            y = mul(g, x)
            if y not in S: S.add(y); fr.append(y)
    return S
def pair_orbits(H, n):
    seen = set(); orbs = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in seen: continue
            O = {(min(g[i], g[j]), max(g[i], g[j])) for g in H}; seen |= O; orbs.append(sorted(O))
    return orbs
def to_d6(n, arcs):
    A = [[0] * n for _ in range(n)]
    for u, v in arcs: A[u][v] = 1
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6): s += chr(63 + int(''.join(map(str, bits[k:k + 6])), 2))
    return s
def to_g6(n, edges):
    E = set(edges); bits = []
    for j in range(1, n):
        for i in range(j): bits.append(1 if (i, j) in E else 0)
    while len(bits) % 6: bits.append(0)
    s = chr(63 + n)
    for k in range(0, len(bits), 6): s += chr(63 + int(''.join(map(str, bits[k:k + 6])), 2))
    return s
def labelg(lines):
    out = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
    assert len(out) == len(lines); return out
def aut_generators(n, edges):
    adj = {v: [] for v in range(n)}
    for u, v in edges: adj[u].append(v)
    txt = f"n={n} g\n" + "\n".join(f"{v} : {' '.join(map(str, adj[v]))};" for v in range(n)) + "\n.\nx\nq\n"
    out = subprocess.run(['dreadnaut'], input=txt, capture_output=True, text=True).stdout
    gens = []; chunk = ""
    for line in out.split('\n'):
        if line.startswith('level') or 'grpsize' in line:
            if chunk.strip():
                g = list(range(n))
                for c in re.findall(r'\(([^)]*)\)', chunk):
                    pts = list(map(int, c.split()))
                    for i in range(len(pts)): g[pts[i]] = pts[(i + 1) % len(pts)]
                gens.append(tuple(g))
            chunk = ""
        elif line.startswith('(') or line.startswith(' '):
            chunk += line
    return gens
def is_even(n, edges):
    for g in aut_generators(n, edges):
        s = 0
        for u, v in edges:
            if g[u] > g[v]: s ^= 1
        if s: return False
    return True

for n in (8, 9, 10):
    H = closure([cyc([[0, 1, 2], [3, 4, 5, 6, 7]], n)], n)
    orbs = pair_orbits(H, n)
    reps = [O[0] for O in orbs]
    orb_arcs = [{(g[p[0]], g[p[1]]) for g in H} for p in reps]
    tl, gl, glist = [], [], []
    for choice in itertools.product([0, 1], repeat=len(orbs)):
        arcs = []; edges = []
        for c, O, A in zip(choice, orbs, orb_arcs):
            if c: arcs += list(A); edges += O
            else: arcs += [(v, u) for u, v in A]
        tl.append(to_d6(n, arcs)); gl.append(to_g6(n, edges)); glist.append(edges)
    ttypes = set(labelg(tl)); gcan = labelg(gl)
    rep_of = {}
    for c, e in zip(gcan, glist): rep_of.setdefault(c, e)
    even_types = [c for c, e in rep_of.items() if is_even(n, e)]
    odd_types = [c for c in rep_of if c not in even_types]
    print(f"n = {n}, H = C15 with orbits 3+5+{'1+' * (n - 8)}: pair orbits {len(orbs)}, objects {2 ** len(orbs)}; "
          f"tau_T = {len(ttypes)}, tau_E = {len(even_types)} (graph types {len(rep_of)}, odd {len(odd_types)})")
    if n == 10:
        # identify K_{5,5}-like odd graphs: 25 edges, bipartite complete
        for c in odd_types:
            e = rep_of[c]
            if len(e) == 25:
                print("   odd type with 25 edges (K_{5,5} candidate), generators of Aut:", len(aut_generators(n, e)),
                      "| degrees:", sorted(sum(1 for u, v in e if x in (u, v)) for x in range(n)))

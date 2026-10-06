#!/usr/bin/env python3
"""Orchestrator check of the frame argument: the 2-wide frame (rings 2,3 of 8x8) has no closed knight tour.
Size identity: e(ring3) - e(ring2) = |ring3| - |ring2| = 8 and ring 3 has only 8 internal knight moves, so a closed
tour of the frame uses all 8 ring-3 internal moves and no ring-2 internal move. Then propagate degree forcing."""
N = 8
F = "abcdefgh"
def ring(i, j):
    c = (N - 1) / 2
    return int(max(abs(i - c), abs(j - c)) - 0.5)
nm = lambda s: F[s[0]] + str(s[1] + 1)
V = [(i, j) for i in range(N) for j in range(N) if ring(i, j) >= 2]
Vs = set(V)
E = set()
for (i, j) in V:
    for a, b in ((1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)):
        u = (i + a, j + b)
        if u in Vs:
            E.add(frozenset(((i, j), u)))
in3 = {e for e in E if all(ring(*x) == 3 for x in e)}
in2 = {e for e in E if all(ring(*x) == 2 for x in e)}
print("frame:", len(V), "squares,", len(E), "knight edges;", len(in3), "inside ring 3,", len(in2), "inside ring 2")
allowed = E - in2
forced = set(in3)
changed = True
while changed:
    changed = False
    for v in V:
        inc = [e for e in allowed if v in e]
        fv = [e for e in forced if v in e]
        if len(inc) < 2:
            print("vertex", nm(v), "has <2 allowed edges -> no tour"); raise SystemExit
        if len(fv) > 2:
            print("vertex", nm(v), "has >2 forced edges -> no tour"); raise SystemExit
        if len(inc) == 2 and not set(inc) <= forced:
            forced |= set(inc); changed = True
        if len(fv) == 2:
            extra = [e for e in inc if e not in forced]
            if extra:
                allowed -= set(extra); changed = True
print("forced edges:", len(forced), "; every square has exactly 2 forced edges:", all(sum(v in e for e in forced) == 2 for v in V))
# cycles of the forced 2-factor
adj = {v: [] for v in V}
for e in forced:
    a, b = tuple(e); adj[a].append(b); adj[b].append(a)
seen, cycles = set(), []
for v in V:
    if v in seen: continue
    cyc, prev, cur = [v], None, v
    seen.add(v)
    while True:
        nxt = [w for w in adj[cur] if w != prev][0]
        if nxt == v: break
        cyc.append(nxt); seen.add(nxt); prev, cur = cur, nxt
    cycles.append(cyc)
print("forced 2-factor cycle lengths:", sorted(len(c) for c in cycles), "-> Hamiltonian:", len(cycles) == 1)
for c in cycles:
    print("   ", " ".join(nm(x) for x in c))

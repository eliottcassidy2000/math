#!/usr/bin/env python3
"""Anti-automorphisms of the Cayley-Dickson sign tournaments CD_k (k = 3, 4, 5): search for an involution psi with
psi(a) -> psi(b) iff b -> a, list its fixed vertex and the invariant triples {v0, a, psi(a)} that are cyclic, build the
quotient tournament Z and the set D_1 of Theorem Q, and print the input for halfpaths.c.
Usage: cd_converse.py k [--halfpath-input]"""
import sys
from cayley_dickson import tournament
sys.setrecursionlimit(10000)
def find_anti(out, n, want_involution=True, limit=None):
    arc = [[(out[a] >> b) & 1 for b in range(n)] for a in range(n)]
    found = []
    perm = [-1] * n; used = [False] * n
    # order vertices to maximise early constraints
    order = list(range(n))
    def rec(i):
        if i == n:
            found.append(perm[:]); return limit is not None and len(found) >= limit
        a = order[i]
        for c in range(n):
            if used[c]: continue
            if want_involution and perm[c] != -1 and perm[c] != a: continue
            ok = True
            for j in range(i):
                b = order[j]
                # anti-automorphism: a -> b iff perm[b] -> perm[a]
                if arc[a][b] != arc[perm[b]][c]: ok = False; break
            if not ok: continue
            if want_involution:
                # if c already has an image it must be a (checked); also if a == c fixed point allowed
                pass
            perm[a] = c; used[c] = True
            if want_involution and c != a and perm[c] == -1:
                # force perm[c] = a later: check consistency lazily
                pass
            if rec(i + 1): return True
            perm[a] = -1; used[c] = False
        return False
    rec(0)
    if want_involution:
        found = [p for p in found if all(p[p[a]] == a for a in range(n))]
    return found
if __name__ == "__main__":
    k = int(sys.argv[1]); out, g = tournament(k); n = (1 << k) - 1
    antis = find_anti(out, n, want_involution=False)
    invs = [p for p in antis if all(p[p[a]] == a for a in range(n))]
    if "--halfpath-input" not in sys.argv:
        print("CD_%d (n = %d): anti-automorphisms found: %d; involutive: %d" % (k, n, len(antis), len(invs)))
    if not invs: sys.exit(0)
    psi = invs[0]; fixed = [a for a in range(n) if psi[a] == a]; v0 = fixed[0]
    arc = lambda a, b: (out[a] >> b) & 1
    plus = [a for a in range(n) if a != v0 and arc(v0, a)]        # one from each pair
    D1 = [a for a in plus if arc(a, psi[a])]
    if "--halfpath-input" in sys.argv:
        verts = []
        for a in plus: verts += [a, psi[a]]
        idx = {v: i for i, v in enumerate(verts)}
        print(len(plus))
        for v in verts:
            mask = 0
            for w in verts:
                if w != v and arc(v, w): mask |= 1 << idx[w]
            print(mask, 1 if arc(v0, v) else 0)
        sys.exit(0)
    print("   fixed vertices of the first involution: %s; pairs: %d; cyclic invariant triples (|D_1|): %d" % (fixed, len(plus), len(D1)))
    # all involutions give the same |D_1| (they are conjugate)
    print("   |D_1| over all involutive anti-automorphisms:", sorted(set(sum(1 for a in range(n) if a != [x for x in range(n) if p[x] == x][0] and arc([x for x in range(n) if p[x] == x][0], a) and arc(a, p[a])) for p in invs)))

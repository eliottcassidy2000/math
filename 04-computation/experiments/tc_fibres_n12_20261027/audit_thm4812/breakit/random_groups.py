"""random_groups.py -- generate tournaments invariant under random odd-order permutation groups (orders 13..18).
Generators: random permutations whose cycles have odd lengths in {3,5,7,9} (and possibly products of such);
closure capped at 20000 elements; groups of even order are impossible (odd-order generators can still generate an
even-order group, e.g. two 3-cycles generate A_4 -- such groups have self-paired orbitals or are skipped when an
orbital contains a pair and its reverse).  For every group, up to R random orientations of the orbital pairs.
Writes lines 'n hex(out_0) ... hex(out_{n-1})' to stdout."""
import random, sys

def closure(gens, n, cap=20000):
    ident = tuple(range(n)); G = {ident}; fr = [ident]
    while fr:
        new = []
        for g in fr:
            for h in gens:
                gh = tuple(h[g[x]] for x in range(n))
                if gh not in G:
                    G.add(gh); new.append(gh)
                    if len(G) > cap: return None
        fr = new
    return G

def random_odd_perm(n, rnd):
    pts = list(range(n)); rnd.shuffle(pts)
    perm = list(range(n)); pos = 0
    ncyc = rnd.randint(1, 4)
    for _ in range(ncyc):
        L = rnd.choice([3, 3, 3, 5, 5, 7, 9])
        if pos + L > n: break
        cyc = pts[pos:pos + L]; pos += L
        for i in range(L): perm[cyc[i]] = cyc[(i + 1) % L]
    return tuple(perm)

def main(seed, ngroups, R):
    rnd = random.Random(seed)
    made = 0; tries = 0
    while made < ngroups and tries < 50 * ngroups:
        tries += 1
        n = rnd.randint(13, 18)
        gens = [random_odd_perm(n, rnd) for _ in range(rnd.randint(1, 3))]
        G = closure(gens, n)
        if G is None or len(G) < 3 or len(G) % 2 == 0: continue
        # orbitals
        idx = {}; orbs = []
        bad = False
        for u in range(n):
            for v in range(n):
                if u == v or (u, v) in idx: continue
                o = {(g[u], g[v]) for g in G}
                if (v, u) in o: bad = True; break
                for pr in o: idx[pr] = len(orbs)
                orbs.append(sorted(o))
            if bad: break
        if bad: continue
        pairs = []; done = set()
        for t, o in enumerate(orbs):
            if t in done: continue
            u, v = o[0]; ti = idx[(v, u)]; done.add(t); done.add(ti); pairs.append((t, ti))
        for _ in range(R):
            out = [0] * n
            for (t, ti) in pairs:
                pick = t if rnd.random() < 0.5 else ti
                for (u, v) in orbs[pick]: out[u] |= 1 << v
            print(n, ' '.join(format(x, 'x') for x in out))
        made += 1
    print(f"# groups {made}", file=sys.stderr)

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))

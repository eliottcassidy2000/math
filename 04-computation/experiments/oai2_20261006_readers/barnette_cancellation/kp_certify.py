#!/usr/bin/env python3
"""Certify beta_P >= lambda_P on the knight torus G_n for prescribed paths P with m edges.

For each path representative P (translations x| D4 x reversal), lambda_P = cheapest local
obstruction (knight_paths2.local_bound).  We prove that every (lambda_P - 1)-set D of edges
disjoint from P leaves a Hamiltonian cycle through P: a pool of Hamiltonian cycles through P
(orbit images of stored cycles) is fed to the C enumerator kp_enum, which lists the D's hit by
every pool cycle; each such D is resolved by CP-SAT (new cycle -> pool; none -> a blocking set
smaller than lambda_P, which would refute beta_P = lambda_P).  The upper bound beta_P <= lambda_P
holds by the explicit local obstruction (re-checked here by CP-SAT infeasibility).
Usage: python3 kp_certify.py n m [maxreps]"""
import sys, os, time, subprocess, random
import knight_paths2 as K

HERE = os.path.dirname(os.path.abspath(__file__))


def local_obstruction_set(n, eng, path, lam, why):
    """explicit blocking set of size lam realizing the local bound"""
    adj, eid = eng.adj, eng.eid
    I = set(path[1:-1]); a, b = path[0], path[-1]
    EPs = set(eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:]))
    kind, v = why
    if kind == "starve":
        usable = sorted(w for w in adj[v] if w not in I)
        keep = usable[0]
        return [eid[(min(v, w), max(v, w))] for w in usable if w != keep]
    if kind == "close":
        usable = sorted(w for w in adj[v] if w not in I)
        return [eid[(min(v, w), max(v, w))] for w in usable if w not in (a, b)]
    if kind == "end":
        other = b if v == a else a
        usable = sorted(w for w in adj[v] if w not in I and w != other)
        return [eid[(min(v, w), max(v, w))] for w in usable if eid[(min(v, w), max(v, w))] not in EPs]


def certify(n, eng, path, k, tag, cap=400, log=print):
    E, nV, eid = eng.E, eng.nV, eng.eid
    EP = [eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:])]
    EPs = set(EP)
    cand = [e for e in range(len(E)) if e not in EPs]
    W = (len(E) + 63) // 64
    blocking = []
    rounds = 0
    while True:
        rounds += 1
        pool, seen = eng.pool_for(EP)
        fn = os.path.join(HERE, f"kp_pool_{tag}.txt")
        with open(fn, "w") as f:
            f.write(f"{len(E)} {W} {k} {len(cand)} {len(seen)} {cap}\n")
            f.write(" ".join(map(str, cand)) + "\n")
            for img in seen:
                words = [0] * W
                for e in img:
                    words[e >> 6] |= 1 << (e & 63)
                f.write(" ".join(format(x, "x") for x in words) + "\n")
        t0 = time.time()
        out = subprocess.run([os.path.join(HERE, "kp_enum"), fn], capture_output=True, text=True).stdout
        dt = time.time() - t0
        U = [list(map(int, line.split()[1:])) for line in out.splitlines() if line.startswith("U")]
        done = [line for line in out.splitlines() if line.startswith("DONE")][0]
        if not U:
            log(f"      round {rounds}: pool {len(seen)} cycles; {done} ({dt:.0f}s) => CERTIFIED")
            return True, blocking, len(seen)
        log(f"      round {rounds}: pool {len(seen)}; {len(U)} uncertified listed ({dt:.0f}s); resolving")
        newc = 0
        random.shuffle(U)
        for D in U[:120]:
            c = K.find_hc(E, nV, set(D), EP)
            if c is None:
                blocking.append([E[e] for e in D])
                log(f"      BLOCKING SET of size {k}: {[E[e] for e in D]}")
                return False, blocking, len(seen)
            assert c != "UNKNOWN"
            if eng.add_cycle(c):
                newc += 1
            if newc >= 40:
                break


def main():
    n, m = int(sys.argv[1]), int(sys.argv[2])
    maxreps = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    start = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    T0 = time.time()
    eng = K.Engine(n)
    eng.seed(40)
    tri = sum(1 for (u, v) in eng.E for w in eng.adj[u] & eng.adj[v]) // 3
    print(f"n={n} m={m}: V={eng.nV} E={len(eng.E)} triangles={tri}; seeded {len(eng.cycles)} HCs", flush=True)
    reps = K.path_reps(n, m)[start:start + maxreps]
    print(f"  {len(reps)} path representatives", flush=True)
    summ = {}
    for r, pts in enumerate(reps):
        path = [n * p[0] + p[1] for p in pts]
        lam, why = K.local_bound(n, eng.adj, path)
        B = local_obstruction_set(n, eng, path, lam, why)
        EP = [eng.eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:])]
        assert len(B) == lam and not (set(B) & set(EP))
        blocks = K.find_hc(eng.E, eng.nV, set(B), EP) is None
        t0 = time.time()
        if lam >= 1:
            ok, blk, ps = certify(n, eng, path, lam - 1, f"n{n}m{m}r{r + start}", log=lambda s: print(s, flush=True))
        else:
            ok, blk, ps = True, [], 0
        res = "beta=lambda" if (ok and blocks) else ("UPPER-FAIL" if not blocks else "LOWER-FAIL")
        summ[(lam, res)] = summ.get((lam, res), 0) + 1
        print(f"   P={pts}: lambda={lam} [{why[0]}], local set blocks: {blocks}; all {lam-1}-sets leave an HC "
              f"through P: {ok} -> {res} ({time.time()-t0:.0f}s, stored HCs {len(eng.cycles)})", flush=True)
    print(f"  SUMMARY n={n} m={m}: {summ}  total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Are all minimum (size lambda_P) tour-blocking sets for a prescribed path P local, i.e. all
lambda_P deleted edges at one vertex?  Enumerate every lambda_P-set D disjoint from P that is hit
by all pool cycles and is NOT local (C enumerator with locality filter); resolve each by CP-SAT.
Usage: python3 kp_minsets.py n m rep_index [rep_index ...]"""
import sys, os, time, subprocess, random
import knight_paths2 as K

HERE = os.path.dirname(os.path.abspath(__file__))


def run(n, m, idxs):
    eng = K.Engine(n)
    eng.seed(40)
    reps = K.path_reps(n, m)
    for r in idxs:
        pts = reps[r]
        path = [n * p[0] + p[1] for p in pts]
        lam, why = K.local_bound(n, eng.adj, path)
        E, nV, eid = eng.E, eng.nV, eng.eid
        EP = [eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:])]
        EPs = set(EP)
        cand = [e for e in range(len(E)) if e not in EPs]
        W = (len(E) + 63) // 64
        exotic = set()
        t0 = time.time()
        rounds = 0
        while True:
            rounds += 1
            pool, seen = eng.pool_for(EP)
            fn = os.path.join(HERE, f"kp_min_n{n}m{m}r{r}.txt")
            with open(fn, "w") as f:
                f.write(f"{len(E)} {W} {lam} {len(cand)} {len(seen)} 3000\n")
                f.write(" ".join(map(str, cand)) + "\n")
                for img in seen:
                    words = [0] * W
                    for e in img:
                        words[e >> 6] |= 1 << (e & 63)
                    f.write(" ".join(format(x, "x") for x in words) + "\n")
                f.write("ENDPOINTS\n" + " ".join(f"{u} {v}" for u, v in E) + "\n")
            out = subprocess.run([os.path.join(HERE, "kp_enum2"), fn], capture_output=True, text=True).stdout
            U = [tuple(map(int, line.split()[1:])) for line in out.splitlines() if line.startswith("U")]
            done = [line for line in out.splitlines() if line.startswith("DONE")][0]
            todo = [D for D in U if D not in exotic]
            print(f"    round {rounds}: pool {len(seen)}; {done}; unresolved non-local {len(todo)} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            if not todo:
                break
            random.shuffle(todo)
            newc = 0
            for D in todo[:200]:
                c = K.find_hc(E, nV, set(D), EP)
                if c is None:
                    exotic.add(D)
                    print(f"    EXOTIC minimum blocking set: {[E[e] for e in D]}", flush=True)
                    continue
                assert c != "UNKNOWN"
                if eng.add_cycle(c):
                    newc += 1
                if newc >= 60:
                    break
            if newc == 0 and all(D in exotic for D in U):
                break
        print(f"  P={pts} lambda={lam} [{why[0]}]: non-local minimum blocking sets: {len(exotic)} "
              f"({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    n, m = int(sys.argv[1]), int(sys.argv[2])
    run(n, m, [int(a) for a in sys.argv[3:]])

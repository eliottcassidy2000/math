#!/usr/bin/env python3
"""Path-extendability of the knight torus G_n: is every path with m edges contained in a
Hamiltonian cycle?  (m = 3 is the knight analogue of the preprint's P4-Hamiltonicity corollary.)
For each path representative (translations x| D4 x reversal) we look for a stored cycle whose
orbit contains it, else ask CP-SAT; a CP-SAT INFEASIBLE answer is a non-extendable path.
Usage: python3 kp_extend.py n mmax"""
import sys, time
import knight_paths2 as K


def main():
    n, mmax = int(sys.argv[1]), int(sys.argv[2])
    T0 = time.time()
    eng = K.Engine(n)
    eng.seed(30)
    print(f"n={n}: V={eng.nV}, seeded {len(eng.cycles)} HCs", flush=True)
    for m in range(1, mmax + 1):
        reps = K.path_reps(n, m)
        bad, viacp = [], 0
        for pts in reps:
            path = [n * p[0] + p[1] for p in pts]
            EP = [eng.eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:])]
            # quick: any stored cycle containing an image of P?
            pre = [[inv[e] for e in EP] for (pe, inv) in eng.perms]
            hit = False
            for C in eng.cycles:
                if any(all(k in C for k in pg) for pg in pre):
                    hit = True; break
            if hit:
                continue
            viacp += 1
            c = K.find_hc(eng.E, eng.nV, set(), EP)
            if c is None:
                bad.append(pts)
            else:
                assert c != "UNKNOWN"
                eng.add_cycle(c)
        print(f"  m={m}: {len(reps)} path classes; {viacp} needed a new CP-SAT call; "
              f"NOT contained in any Hamiltonian cycle: {len(bad)} {bad[:5]} ({time.time()-T0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

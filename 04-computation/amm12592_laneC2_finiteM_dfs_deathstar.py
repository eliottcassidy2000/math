#!/usr/bin/env python3
"""Lane C2: exact finite-M feasibility DFS for the AMM 12592 deficit ledger
(THM-2966 frame), death-star coinC2.

Question. For depth law d_m = floor(g1*m/g2) + D0, is there ANY assignment
of doubled deficits delta_c (integers, |delta_c| <= cap_c = binom(d_m,k),
delta_c == cap_c mod 2) for rows 1..M such that the PROVED necessary
envelope

    |D_m(p)| <= (p^{m+1} + q^{m+1}) / 2      for every level m <= M   (E_m)

holds at every bias p in a finite rational set P?  (E_m) at finitely many
biases is a NECESSARY relaxation of extractor feasibility, so an
exhausted-search infeasibility verdict at finite M is a THEOREM: no fair
extractor with this deadline envelope exists.  A feasible assignment is
only a surviving corridor (weaker than an extractor).

Method: depth-first search over cells in row order, exact integers
(denominators cleared at scale den^Lmax per bias), with a sound prune:
a prefix with value S at bias p dies if for some checkpoint m' >= here

    |S| > RestOfCurrentRowWiggle + sum_{rows r in (row, m']} R_r + Env_m',

R_r = max possible |row-r contribution|.  The prune is a valid
relaxation, so exhaustion certifies infeasibility of (E_m) over P.

Because the verdict quantifies over ALL assignments, an infeasible M also
closes lane F (policy hostility) for this (gamma, D0): no anticipatory
policy can satisfy what no assignment satisfies.
"""

import sys
from math import comb

sys.setrecursionlimit(100000)


def run(g1, g2, D0, Mmax, biases, parity=True, node_cap=None, tag=""):
    d = lambda m: (g1 * m) // g2 + D0
    A = lambda m: m + d(m) + 1
    Lmax = A(Mmax)
    P = biases
    nP = len(P)

    rows = []          # rows[m-1] = list of (cap, par, weights)
    for m in range(1, Mmax + 1):
        dm = d(m)
        cells = []
        for k in range(dm + 1):
            cap = comb(dm, k)
            for (z, o) in ((m + dm - k, k + 1), (k + 1, m + dm - k)):
                w = tuple(n**z * (dd - n)**o * dd**(Lmax - z - o)
                          for (n, dd) in P)
                cells.append((cap, cap % 2 if parity else -1, w))
        rows.append(cells)

    R = [tuple(sum(c[0] * c[2][i] for c in rows[m]) for i in range(nP))
         for m in range(Mmax)]
    Env = [tuple((n**(m + 2) + (dd - n)**(m + 2)) * dd**(Lmax - m - 2)
                 for (n, dd) in P) for m in range(Mmax)]

    def choices(cap, par):
        if par == -1:
            vals = list(range(-cap, cap + 1))
        else:
            vals = [v for v in range(-cap, cap + 1) if (v - par) % 2 == 0]
        return sorted(vals, key=abs)

    results = {}
    for Mtry in range(1, Mmax + 1):
        allowed = []
        for m in range(Mtry):
            best = None
            for mp in range(m, Mtry):
                tot = tuple(Env[mp][i] +
                            sum(R[r][i] for r in range(m + 1, mp + 1))
                            for i in range(nP))
                best = tot if best is None else tuple(
                    min(a, b) for a, b in zip(best, tot))
            allowed.append(best)

        # per-row suffix wiggle within the row
        suffix = []
        for m in range(Mtry):
            cells = rows[m]
            sw = [tuple(0 for _ in range(nP))]
            for c in reversed(cells):
                sw.append(tuple(sw[-1][i] + c[0] * c[2][i]
                                for i in range(nP)))
            sw.reverse()
            suffix.append(sw)   # suffix[m][ci] = wiggle of cells ci..end

        state = {"nodes": 0, "cap": node_cap}

        def dfs(m, ci, S):
            state["nodes"] += 1
            if state["cap"] and state["nodes"] > state["cap"]:
                raise TimeoutError
            cells = rows[m]
            if ci == len(cells):
                for i in range(nP):
                    if abs(S[i]) > Env[m][i] or abs(S[i]) > allowed[m][i]:
                        return False
                if m + 1 == Mtry:
                    return True
                return dfs(m + 1, 0, S)
            rest = suffix[m][ci]
            for i in range(nP):
                if abs(S[i]) > rest[i] + allowed[m][i]:
                    return False
            cap, par, w = cells[ci]
            for v in choices(cap, par):
                S2 = tuple(S[i] + v * w[i] for i in range(nP)) if v else S
                if dfs(m, ci + 1, S2):
                    return True
            return False

        try:
            ok = dfs(0, 0, tuple(0 for _ in range(nP)))
            results[Mtry] = ("FEASIBLE" if ok else "INFEASIBLE",
                             state["nodes"])
            print(f"[{tag}] M={Mtry:2d}: "
                  f"{'FEASIBLE  ' if ok else 'INFEASIBLE'} "
                  f"nodes={state['nodes']}", flush=True)
            if not ok:
                print(f"[{tag}] *** INFEASIBLE at M={Mtry} by exhausted "
                      f"search: envelope THEOREM candidate ***", flush=True)
                break
        except TimeoutError:
            print(f"[{tag}] M={Mtry:2d}: NODE-CAP after {state['nodes']} "
                  f"nodes", flush=True)
            results[Mtry] = ("NODE-CAP", state["nodes"])
            break
    return results


if __name__ == "__main__":
    BIASES = [(1, 2), (1, 3), (2, 3), (2, 5), (3, 5), (1, 4), (3, 4),
              (1285, 2181), (8847357, 11821757)]
    print("=== gamma=1/2 D0=0, parity ON ===", flush=True)
    run(1, 2, 0, 14, BIASES, parity=True, node_cap=150_000_000,
        tag="g1/2,par")
    print("=== gamma=1/2 D0=0, parity OFF (control) ===", flush=True)
    run(1, 2, 0, 10, BIASES, parity=False, node_cap=30_000_000,
        tag="g1/2,nopar")
    print("=== gamma=1/3 D0=0, parity ON ===", flush=True)
    run(1, 3, 0, 12, BIASES, parity=True, node_cap=80_000_000,
        tag="g1/3,par")
    print("DONE", flush=True)

#!/usr/bin/env python3
"""Exact symbolic audit of rank-nine 3+3 full-ladder traces.

With paid cost ten and one unit of initial row-rank slack, THM-2850 forces
both rounds to have size three.  After the common-slope gauge alpha=0, an
initial headroom-three first wave F is equivalent to W restricted to the
three survivors vanishing.  These are linear equations in the seven
remaining slope parameters, so generic exact ranks on every resulting
linear family decide whether rank(A)=9 and rank(A_R)=6 can coexist.
"""

from itertools import combinations, combinations_with_replacement
import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

from ladder_all_cost10_kernel_search import make_topology
from ladder_full_f4_symbolic_audit import generator_matrix
from ladder_full_f3_exact_audit import full_cycle_matrix


def main():
    cases = families = possible = 0
    rank_pairs = {}
    witnesses = []
    for positions in combinations_with_replacement(range(6), 3):
        topology = make_topology((0, 1, 2, 3), positions)
        if topology is None:
            continue
        params, A = generator_matrix(topology)
        W = full_cycle_matrix(topology, params)
        for fire in combinations(range(6), 3):
            rest = tuple(v for v in range(6) if v not in fire)
            # Gauge alpha=0 and force every entry of W_R to vanish.
            equations = [params[0]]
            equations.extend(W[i, v] for i in range(4) for v in rest)
            M, rhs = sp.linear_eq_to_matrix(equations, params)
            if any(rhs):
                raise RuntimeError("unexpected inhomogeneous equation")
            basis = M.nullspace()
            cases += 1
            if not basis:
                continue
            families += 1
            z = sp.symbols(f"z0:{len(basis)}")
            subs = {
                params[i]: sum(z[j] * basis[j][i] for j in range(len(basis)))
                for i in range(len(params))
            }
            As = A.subs(subs)
            cols = tuple(2 * v + j for v in rest for j in (0, 1))
            pair = (
                As.to_DM().to_field().rank(),
                As[:, cols].to_DM().to_field().rank(),
            )
            rank_pairs[pair] = rank_pairs.get(pair, 0) + 1
            if pair == (9, 6):
                possible += 1
                witnesses.append((positions, fire, len(basis)))
    print("LADDER_FULL_RANK9_F3_SYMBOLIC_AUDIT")
    print(f"cases={cases};families={families};rank_possible={possible}")
    print(f"rank_pairs={tuple(sorted(rank_pairs.items()))}")
    print(f"possible={tuple(witnesses)}")
    print("status=RANK9_TWO_ROUND_IMPOSSIBLE" if not possible else "status=POTENTIALS_REMAIN")


if __name__ == "__main__":
    main()

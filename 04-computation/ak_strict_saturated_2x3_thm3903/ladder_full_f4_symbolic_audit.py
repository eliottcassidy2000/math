#!/usr/bin/env python3
"""Exact symbolic audit of 4+2 two-round traces in the full 2x3 ladder."""

from itertools import combinations, combinations_with_replacement
import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

from ladder_all_cost10_kernel_search import make_topology, incidence

N = 6


def generator_matrix(topology):
    edges = topology[2]
    params = sp.symbols("r0:8")
    rows = []
    for u, v, p in edges:
        row = [0]*(2*N)
        rho = params[p]
        if u != 6:
            row[2*u], row[2*u+1] = 1, rho
        if v != 6:
            row[2*v], row[2*v+1] = -1, -rho
        rows.append(row)
    return params, sp.Matrix(rows)


def main():
    cases = families = rank_possible = 0
    rank_pairs = {}
    possible = []
    for positions in combinations_with_replacement(range(N), 3):
        topology = make_topology((0, 1, 2, 3), positions)
        if topology is None:
            continue
        params, A = generator_matrix(topology)
        coeff = topology[4]
        for fire in combinations(range(N), 4):
            rest = tuple(v for v in range(N) if v not in fire)
            # Kernel is all of Q^rest, so every cycle row must vanish there.
            M = sp.Matrix([[coeff[i][v][p] for p in range(8)]
                           for i in range(4) for v in rest])
            basis = M.nullspace()
            cases += 1
            if not basis:
                continue
            families += 1
            terminal_cols = tuple(2*v+j for v in rest for j in (0, 1))
            z = sp.symbols(f"z0:{len(basis)}")
            subs = {params[i]: sum(z[j]*basis[j][i] for j in range(len(basis)))
                    for i in range(8)}
            As = A.subs(subs)
            # DomainMatrix computes exact generic rank over Q(z), avoiding
            # heuristic expression-zero tests.
            pair = (As.to_DM().to_field().rank(),
                    As[:, terminal_cols].to_DM().to_field().rank())
            rank_pairs[pair] = rank_pairs.get(pair, 0)+1
            if pair == (10, 4):
                rank_possible += 1
                possible.append((positions, fire, len(basis)))
    print("LADDER_FULL_F4_SYMBOLIC_AUDIT")
    print(f"cases={cases};families={families};rank_possible={rank_possible}")
    print(f"rank_pairs={tuple(sorted(rank_pairs.items()))}")
    print(f"possible={tuple(possible)}")


if __name__ == "__main__":
    main()

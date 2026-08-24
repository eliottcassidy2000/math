#!/usr/bin/env python3
"""Rank-variety attack on the remaining saturated [2,3] score-5/3 traces.

This is an exploratory exact-Q eliminator.  It fixes the common rung-slope
gauge to zero and encodes every prescribed firing wave by exact ranks of the
cycle-slope matrices on successive live vertex sets.
"""

from itertools import combinations, combinations_with_replacement, product
import argparse
import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

from ladder_all_cost10_kernel_search import make_topology


PARAMS = sp.symbols("r0:8")
_CYCLE_CACHE = {}
_RANK_CACHE = {}
_SURVIVOR_UPPER_CACHE = {}


def cycle_matrix(topology, live):
    """Return C diag(rho) B on the current live vertex set."""
    live = tuple(sorted(live))
    cache_key = (topology[1], live)
    if cache_key in _CYCLE_CACHE:
        return _CYCLE_CACHE[cache_key]
    index = {v: i for i, v in enumerate(live)}
    active = []
    for u, v, p in topology[2]:
        row = [sp.Integer(0)] * len(live)
        if u in index:
            row[index[u]] += 1
        if v in index:
            row[index[v]] -= 1
        if any(row):
            active.append((row, p))
    B = sp.Matrix([x[0] for x in active])
    if B.rank() != len(live):
        _CYCLE_CACHE[cache_key] = None
        return None
    cycles = sp.Matrix.hstack(*B.T.nullspace()).T
    D = sp.diag(*(PARAMS[p] for _, p in active))
    result = sp.simplify(cycles * D * B)
    _CYCLE_CACHE[cache_key] = result
    return result


def minors(matrix, size):
    if size == 0:
        return (sp.Integer(1),)
    if size > min(matrix.shape):
        return ()
    result = []
    for rs in combinations(range(matrix.rows), size):
        for cs in combinations(range(matrix.cols), size):
            value = sp.factor(matrix.extract(rs, cs).det())
            if value != 0 and value not in result:
                result.append(value)
    return tuple(result)


def rank_conditions(matrix, exact_rank):
    if matrix is None or exact_rank > min(matrix.shape):
        return None
    lower = minors(matrix, exact_rank)
    if not lower:
        return None
    upper = minors(matrix, exact_rank + 1)
    gram = sp.expand(sum(x * x for x in lower))
    return upper, gram


def headroom_profiles(sizes):
    q = len(sizes)
    result = []
    for profile in product(range(5), repeat=q):
        if profile[0] > 4 or profile[-1] != sizes[-1]:
            continue
        if any(profile[j] < sizes[j] for j in range(q)):
            continue
        if any(not 0 <= profile[j] - profile[j + 1] <= sizes[j]
               for j in range(q - 1)):
            continue
        result.append(profile)
    return tuple(result)


def trace_families(kind):
    universe = set(range(6))
    if kind == "corner3":
        first = frozenset((0,))
        for mask in range(1, 1 << 5):
            second = frozenset(v for i, v in enumerate((1, 2, 3, 4, 5)) if mask & (1 << i))
            third = frozenset(universe - first - second)
            if third:
                yield (first, second, third)
    elif kind == "rung3":
        first = frozenset((0, 3))
        rest = (1, 2, 4, 5)
        for mask in range(1, 1 << 4):
            second = frozenset(v for i, v in enumerate(rest) if mask & (1 << i))
            third = frozenset(universe - first - second)
            if third:
                yield (first, second, third)
    elif kind == "corner4":
        first, second = frozenset((0,)), frozenset((3,))
        cell = (1, 2, 4, 5)
        for mask in range(1, 1 << 4):
            third = frozenset(v for i, v in enumerate(cell) if mask & (1 << i))
            fourth = frozenset(universe - first - second - third)
            if fourth:
                yield (first, second, third, fourth)
    else:
        raise ValueError(kind)


def audit_case(topology, waves, profile):
    equations = [PARAMS[0]]
    grams = []
    live = set(range(6))
    for wave, headroom in zip(waves, profile):
        live_tuple = tuple(sorted(live))
        W = cycle_matrix(topology, live_tuple).subs({PARAMS[0]: 0})
        rank_key = (topology[1], live_tuple, headroom)
        if rank_key not in _RANK_CACHE:
            _RANK_CACHE[rank_key] = rank_conditions(W, headroom)
        condition = _RANK_CACHE[rank_key]
        if condition is None:
            return "generic_rank_fail", None
        upper, gram = condition
        equations.extend(upper)
        grams.append(gram)
        survivors = live - set(wave)
        restricted = W[:, [i for i, v in enumerate(sorted(live)) if v in survivors]]
        survivor_bound = headroom - len(wave)
        if survivor_bound < 0:
            return "generic_rank_fail", None
        survivor_key = (topology[1], live_tuple, tuple(sorted(survivors)), survivor_bound)
        if survivor_key not in _SURVIVOR_UPPER_CACHE:
            _SURVIVOR_UPPER_CACHE[survivor_key] = minors(restricted, survivor_bound + 1)
        equations.extend(_SURVIVOR_UPPER_CACHE[survivor_key])
        live = survivors
    equations = tuple(dict.fromkeys(sp.factor(x) for x in equations if x != 0))
    # Gauge is already substituted in every matrix, so remove its standalone
    # equation and variable from the actual elimination ideal.
    equations = tuple(x for x in equations if x != PARAMS[0])
    nonzero = tuple(sp.factor(x) for x in grams)
    variables = sorted(set().union(
        *(x.free_symbols for x in equations + nonzero)), key=str)
    sat = sp.symbols(f"sat0:{len(nonzero)}")
    gb = sp.groebner(list(equations) +
                     [1 - z * witness for z, witness in zip(sat, nonzero)],
                     *sat, *variables, order="grevlex")
    if gb == [1]:
        return "unit", None
    signature = (len(equations), len(gb.polys),
                 tuple(str(x.as_expr()) for x in gb.polys[-3:]))
    return "potential", signature


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=("corner3", "rung3", "corner4"))
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--potential-limit", type=int, default=20)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()

    jobs = []
    for positions in combinations_with_replacement(range(6), 3):
        topology = make_topology((0, 1, 2, 3), positions)
        for waves in trace_families(args.kind):
            sizes = tuple(len(x) for x in waves)
            for profile in headroom_profiles(sizes):
                jobs.append((positions, topology, waves, profile))
    total_jobs = len(jobs)
    stop = total_jobs if args.stop is None else min(args.stop, total_jobs)
    jobs = jobs[args.start:stop]
    counts = {"generic_rank_fail": 0, "unit": 0, "potential": 0}
    potentials = []
    for local, (positions, topology, waves, profile) in enumerate(jobs):
        status, signature = audit_case(topology, waves, profile)
        counts[status] += 1
        if status == "potential" and len(potentials) < args.potential_limit:
            potentials.append((args.start + local, positions,
                               tuple(tuple(sorted(x)) for x in waves), profile, signature))
        if args.progress and (local + 1) % 100 == 0:
            print(f"progress={local+1}/{len(jobs)};counts={counts}", flush=True)
    print("LADDER_MULTIROUND_RANK_VARIETY")
    print(f"kind={args.kind};global_jobs={total_jobs};slice=({args.start},{stop});slice_jobs={len(jobs)}")
    print(f"counts={counts}")
    print(f"potentials={tuple(potentials)}")


if __name__ == "__main__":
    main()

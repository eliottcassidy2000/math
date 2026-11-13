"""Independent constructions for the root lane's elementary theorems."""
from itertools import combinations_with_replacement, product
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def grounded_theta(arms):
    # Build the actual multigraph first, including direct parallel hub edges.
    n = sum(arms) - 1
    lap = [[0] * n for _ in range(n)]
    edges = []
    nxt = 2
    for length in arms:
        path = [0] + list(range(nxt, nxt + length - 1)) + [1]
        nxt += length - 1
        edges.extend(zip(path, path[1:]))
    need(nxt == n, 'theta vertex count')
    for u, v in edges:
        lap[u][u] += 1
        lap[v][v] += 1
        lap[u][v] -= 1
        lap[v][u] -= 1
    keep = [0] + list(range(2, n))
    result = [[lap[i][j] for j in keep] for i in keep]
    result[0][0] -= 1
    return result, edges, n


def signed_maps(m):
    return {tuple((s * x + d) % m for x in range(m))
            for s in (-1, 1) for d in range(m)}


def independent_states(m, b):
    perms = signed_maps(m)
    # Two loop holonomies suffice for the hostile intersections and every
    # claimed normal form. This is a finite independent control, not its proof.
    fixed = {tuple(x for x in range(m) if p[x] == x and q[x] == x)
             for p in perms for q in perms}
    states = set()
    for fs in fixed:
        for transports in product(perms, repeat=b - 1):
            states.add(frozenset((x,) + tuple(p[x] for p in transports) for x in fs))
    return states


def main():
    spec = importlib.util.spec_from_file_location('root_triples', HERE.parent / 'triples' / 'check.py')
    root = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(root)
    count = 0
    for arms in combinations_with_replacement(range(1, 13), 3):
        literal, _, _ = grounded_theta(arms)
        star, _ = root.theta(arms)
        need(literal == star, ('literal grounded theta minus root unit', arms))
        count += 1

    completions = []
    for r in (4, 6, 8):
        _, edges, n = grounded_theta((2, 3, r))
        total = 0
        for mask in range(1 << len(edges)):
            deg = [0] * n
            for i, (u, v) in enumerate(edges):
                if mask >> i & 1:
                    deg[u] += 1
                    deg[v] += 1
            total += deg == [2, 2] + [1] * (n - 2)
        need(total == 2, ('theta demand hostile', r))
        completions.append(dict(arms=[2, 3, r], hub_demands=[2, 2],
                                internal_demands=1, degree_completions=total))

    state_counts = []
    for m in range(2, 9):
        # R0 and R1 are square-labelled reflections. Close their generated group.
        generators = [tuple(-x % m for x in range(m)),
                      tuple((1 - x) % m for x in range(m))]
        group = {tuple(range(m))}
        frontier = list(group)
        while frontier:
            p = frontier.pop()
            for g in generators:
                composite = tuple(g[p[x]] for x in range(m))
                if composite not in group:
                    group.add(composite)
                    frontier.append(composite)
        need(group == signed_maps(m), ('square reflection generation', m))
        for b in (1, 2, 3):
            states = independent_states(m, b)
            expected = (1 + 2 ** (b - 1) if m == 2 else
                        1 + (m ** b if m % 2 else m ** b // 2) + (2 * m) ** (b - 1))
            need(len(states) == expected, ('independent boundary-state count', m, b))
            state_counts.append(dict(modulus=m, boundaries=b, states=len(states)))

    report = dict(status='PASS', literal_grounded_theta_matrices=count,
                  theta_degree_hostiles=completions, independent_boundary_counts=state_counts,
                  square_reflection_generator_moduli=list(range(2, 9)),
                  scope='Finite independent controls supplement the all-parameter proofs; '
                        'no tournament or literal square-sum realization inference.')
    (HERE / 'audit_root.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', theta_matrices=count, boundary_counts=len(state_counts))))


if __name__ == '__main__':
    main()

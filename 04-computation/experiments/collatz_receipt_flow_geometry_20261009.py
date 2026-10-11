"""Exact finite policy geometry and ROOT-receipt compilation.

Production only authenticates supplied receipts; it never discovers an orbit.
The fixed 27 word is declared finite data. Abstract escape models are not
claimed to be Collatz trajectories. Standard library; -O changes no checks.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product

import paper_local_global_grounding_20261007 as inherited

CHECKS = 0
ROOT27 = (1, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3,
          1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4)


def require(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def odd(n):
    require(type(n) is int and n > 0 and n % 2 == 1, 'exact positive odd label')


@dataclass(frozen=True)
class Kernel:
    nodes: tuple
    rows: tuple


def audit_kernel(k):
    require(type(k) is Kernel and type(k.nodes) is tuple and type(k.rows) is tuple,
            'exact Kernel with tuple data')
    for n in k.nodes:
        odd(n)
    require(1 in k.nodes and len(set(k.nodes)) == len(k.nodes), 'distinct labels and ROOT1')
    require(len(k.rows) == len(k.nodes), 'complete rows; missing mass is not ROOT')
    for n, row in zip(k.nodes, k.rows):
        require(type(row) is tuple and len(row) == len(k.nodes), 'complete tuple row')
        require(all(type(p) is F and p >= 0 for p in row), 'nonnegative exact Fraction entries')
        require(sum(row) == (0 if n == 1 else 1), 'killed ROOT, stochastic nonroot row')


def kernel(nodes, weighted_arcs):
    """Convenience construction; duplicate arcs are rejected, never overwritten."""
    require(type(nodes) is tuple and type(weighted_arcs) is tuple, 'tuple inputs')
    for n in nodes:
        odd(n)
    require(len(set(nodes)) == len(nodes), 'distinct nodes')
    rows = [[F(0) for _ in nodes] for _ in nodes]
    seen = set()
    for item in weighted_arcs:
        require(type(item) is tuple and len(item) == 3, 'arc triple')
        a, b, p = item
        odd(a); odd(b)
        require(a in nodes and b in nodes and (a, b) not in seen, 'listed unique arc')
        require(type(p) is F and p > 0, 'positive exact Fraction arc')
        seen.add((a, b))
        rows[nodes.index(a)][nodes.index(b)] = p
    out = Kernel(nodes, tuple(tuple(row) for row in rows))
    audit_kernel(out)
    return out


def reachable(k, source):
    audit_kernel(k); odd(source)
    require(source in k.nodes, 'listed source')
    seen, todo = {source}, [source]
    while todo:
        n = todo.pop()
        for target, p in zip(k.nodes, k.rows[k.nodes.index(n)]):
            if p and target not in seen:
                seen.add(target); todo.append(target)
    return frozenset(seen)


def source_kernel(k, source):
    keep = reachable(k, source) | {1}
    ids = tuple(i for i, n in enumerate(k.nodes) if n in keep)
    out = Kernel(tuple(k.nodes[i] for i in ids),
                 tuple(tuple(k.rows[i][j] for j in ids) for i in ids))
    audit_kernel(out)
    return out


def audit_rank(k, h):
    audit_kernel(k)
    require(type(h) is tuple and len(h) == len(k.nodes)
            and all(type(x) is F and x >= 0 for x in h), 'complete nonnegative Fraction rank')
    require(h[k.nodes.index(1)] == 0, 'ROOT rank zero')
    for n, row, value in zip(k.nodes, k.rows, h):
        if n != 1:
            require(value - sum((p*x for p, x in zip(row, h)), F(0)) >= 1,
                    'unit expected rank decrease')
    return max((sum((p*x for p, x in zip(row, h)), F(0))/value
                for n, row, value in zip(k.nodes, k.rows, h) if n != 1), default=F(0))


def potential(k):
    """None means an ungrounded support component, not Collatz nonconvergence."""
    audit_kernel(k)
    ids = tuple(i for i, n in enumerate(k.nodes) if n != 1)
    matrix = [[F(i == j)-k.rows[i][j] for j in ids] for i in ids]
    try:
        values = inherited.solve(matrix, [F(1)]*len(ids))
    except ValueError:
        return None
    h = [F(0)]*len(k.nodes)
    for i, value in zip(ids, values):
        h[i] = value
    h = tuple(h)
    audit_rank(k, h)
    return h


def compile_root(receipts, k, source, supplied_rank=None):
    """Greedy rank path, then bounded inherited first-hit receipt export.

    Every positive arc is authenticated, including arcs unused by the path.
    A source can compile even if a disjoint component is ungrounded.
    """
    audit_kernel(k); odd(source)
    require(source in k.nodes and type(receipts) is tuple, 'listed source and receipt tuple')
    arcs = {}
    for rec in receipts:
        inherited.verify_receipt(rec)
        for a, b, direction in ((rec.left, rec.right, True), (rec.right, rec.left, False)):
            arcs.setdefault((a, b), (rec, direction))
    for a, row in zip(k.nodes, k.rows):
        for b, p in zip(k.nodes, row):
            require(not p or a == b or (a, b) in arcs, 'positive arc lacks authenticated receipt')
    local = source_kernel(k, source)
    h = potential(local) if supplied_rank is None else supplied_rank
    if h is None:
        return None
    theta = audit_rank(local, h)
    current, path, steps = source, [source], []
    while current != 1:
        i = local.nodes.index(current)
        choices = [(h[j], n) for j, n in enumerate(local.nodes)
                   if local.rows[i][j] and h[j] < h[i]]
        require(choices, 'rank must expose a strictly descending successor')
        target = min(choices)[1]
        steps.append(arcs[(current, target)])
        path.append(target); current = target
        require(len(path) <= len(local.nodes), 'strict rank has no repeated node')
    clocks, word = inherited.authenticate_path(tuple(steps), source)
    return {'path': tuple(path), 'word': word, 'clock': clocks,
            'nodes': local.nodes, 'rank': h, 'contraction': theta}


def audit_common_rank(nodes, vertices, h):
    """Each row is an independently selectable convex hull of listed vertices."""
    require(type(vertices) is tuple and len(vertices) == len(nodes), 'complete row polytopes')
    require(all(type(rowset) is tuple and rowset for rowset in vertices), 'nonempty vertex tuples')
    base = Kernel(nodes, tuple(rowset[0] for rowset in vertices))
    audit_rank(base, h)
    for i, rowset in enumerate(vertices):
        for row in rowset:
            rows = list(base.rows); rows[i] = row
            audit_rank(Kernel(nodes, tuple(rows)), h)
    return True


def audit_refinement(nodes, vertices, h, fine, projection, mixtures):
    """Verify a finite instance of the countable hidden-state refinement theorem.

    Labels in this API are abstract proof states, not authenticated integers.
    Supplied convex coefficients are exact witnesses for every projected row.
    Only the genuine terminal may project to ROOT. No seam is reconstructed
    from this quotient: compile_root separately requires actual receipts.
    """
    audit_common_rank(nodes, vertices, h)
    audit_kernel(fine)
    require(type(projection) is tuple and len(projection) == len(fine.nodes), 'complete projection')
    require(type(mixtures) is tuple and len(mixtures) == len(fine.nodes), 'complete convex witnesses')
    for n, coarse in zip(fine.nodes, projection):
        odd(coarse)
        require(coarse in nodes and ((n == 1) == (coarse == 1)), 'terminal fiber must be exact')
    for row, coarse, weights in zip(fine.rows, projection, mixtures):
        i = nodes.index(coarse)
        require(type(weights) is tuple and len(weights) == len(vertices[i])
                and all(type(w) is F and w >= 0 for w in weights) and sum(weights) == 1,
                'exact convex coefficients')
        pushed = tuple(sum((p for p, target in zip(row, projection) if target == n), F(0))
                       for n in nodes)
        expected = tuple(sum((w*vertex[j] for w, vertex in zip(weights, vertices[i])), F(0))
                         for j in range(len(nodes)))
        require(pushed == expected, 'refinement has changed the certified coarse policy')
    lifted = tuple(h[nodes.index(n)] for n in projection)
    return lifted, audit_rank(fine, lifted)


def worst_survival(nodes, vertices, horizon):
    require(type(horizon) is int and horizon >= 0, 'nonnegative exact horizon')
    # Authentication/geometry is independent; validate all possible rows first.
    base = Kernel(nodes, tuple(rowset[0] for rowset in vertices))
    audit_kernel(base)
    for i, rowset in enumerate(vertices):
        require(type(rowset) is tuple and rowset, 'nonempty tuple vertices')
        for row in rowset:
            rows = list(base.rows); rows[i] = row
            audit_kernel(Kernel(nodes, tuple(rows)))
    values = tuple(F(n != 1) for n in nodes)
    for _ in range(horizon):
        values = tuple(F(0) if n == 1 else
                       max(sum((p*x for p, x in zip(row, values)), F(0)) for row in rowset)
                       for n, rowset in zip(nodes, vertices))
    return values


def rejects(call):
    try:
        call()
    except ValueError:
        return True
    return False


def main():
    R = inherited.Receipt
    square_receipts = (R(7, 3, (1, 1, 2, 3), (1,)), R(3, 5, (1,), ()), R(5, 1, (4,), ()))
    nodes = (1, 3, 5, 7)
    square = kernel(nodes, ((7, 3, F(1)), (3, 7, F(1, 2)), (3, 5, F(1, 2)),
                            (5, 3, F(1, 2)), (5, 1, F(1, 2))))
    check(potential(square) == (F(0), F(8), F(5), F(9)), 'uniform square rank')
    compiled = compile_root(square_receipts, square, 7)
    check(compiled['path'] == (7, 3, 5, 1) and compiled['word'] == (1, 1, 2, 3, 4),
          'rank extracts actual strict ROOT word')
    check(inherited.replay(7, compiled['word']) == 1, 'independent exported replay')
    vertices = (((F(0),)*4,),
                ((F(0), F(0), F(1, 4), F(3, 4)), (F(0), F(0), F(3, 4), F(1, 4))),
                ((F(1, 4), F(3, 4), F(0), F(0)), (F(3, 4), F(1, 4), F(0), F(0))),
                ((F(0), F(1), F(0), F(0)),))
    robust = (F(0), F(32), F(25), F(33))
    check(audit_common_rank(nodes, vertices, robust), 'row-polytope common rank')
    for choices in product(*vertices):
        fixed = Kernel(nodes, choices)
        check(audit_rank(fixed, robust) < 1, 'all four policy cube vertices')
        check(compile_root(square_receipts, fixed, 7, robust)['word'] == compiled['word'],
              'common rank compiles each policy')
    for horizon in range(49):
        check(max(worst_survival(nodes, vertices, horizon)) <= F(63, 64)**(horizon//3),
              'uniform floor controls all history-dependent policies')
    # Hidden-mode rows need not be Markov after forgetting their modes.
    fine_nodes = (1, 3, 5, 7, 9, 11, 13)
    projection = (1, 3, 3, 5, 5, 7, 7)
    refinement_count = 0
    last_fine = last_mixtures = None
    for hidden_choices, mixing in product(range(64), (F(1, 4), F(1, 2), F(3, 4))):
        rows, mixtures = [], []
        for i, coarse in enumerate(projection):
            j = nodes.index(coarse)
            weights = (F(1),) if len(vertices[j]) == 1 else (mixing, 1-mixing)
            coarse_row = tuple(sum((w*v[k] for w, v in zip(weights, vertices[j])), F(0))
                               for k in range(len(nodes)))
            row = [F(0)]*len(fine_nodes)
            for target, p in zip(nodes, coarse_row):
                choices = [k for k, value in enumerate(projection) if value == target]
                selected = choices[((hidden_choices >> max(i-1, 0)) & 1) % len(choices)]
                row[selected] += p
            rows.append(tuple(row)); mixtures.append(weights)
        fine = Kernel(fine_nodes, tuple(rows))
        lifted, theta = audit_refinement(nodes, vertices, robust, fine, projection, tuple(mixtures))
        check(lifted == (F(0), F(32), F(32), F(25), F(25), F(33), F(33)) and theta <= F(32, 33),
              'exact hidden-mode refinement rank')
        refinement_count += 1
        last_fine, last_mixtures = fine, tuple(mixtures)
    check(rejects(lambda: audit_refinement(nodes, vertices, robust, last_fine,
                  (1, 1, 3, 5, 5, 7, 7), last_mixtures)), 'false ROOT fiber rejected')
    wrong_rows = list(last_fine.rows)
    wrong_rows[1] = (F(0), F(1), F(0), F(0), F(0), F(0), F(0))
    check(rejects(lambda: audit_refinement(nodes, vertices, robust,
                  Kernel(fine_nodes, tuple(wrong_rows)), projection, last_mixtures)),
          'new trapping edge is not certified refinement')
    for a, b in product(range(1, 10), repeat=2):
        fixed = kernel(nodes, ((7, 3, F(1)), (3, 7, F(a, 10)), (3, 5, F(10-a, 10)),
                               (5, 3, F(b, 10)), (5, 1, F(10-b, 10))))
        check(potential(fixed) is not None, '81 strict-interior cube controls')
    closed = kernel(nodes, ((7, 3, F(1)), (3, 7, F(1)), (5, 1, F(1))))
    check(potential(closed) is None and compile_root(square_receipts, closed, 7) is None,
          'boundary policy erases exit')
    check(compile_root(square_receipts, closed, 5)['word'] == (4,), 'source-local restriction')

    # Every support on three nonroot vertices, with empty rows retained as loops.
    graph_count = 0
    for mask in range(1 << 9):
        arcs = []
        for i, n in enumerate(nodes[1:]):
            targets = [m for j, m in enumerate(nodes) if m != n and mask & (1 << (3*i+(j if j < i+1 else j-1)))]
            if not targets:
                targets = [n]
            arcs.extend((n, m, F(1, len(targets))) for m in targets)
        k = kernel(nodes, tuple(arcs))
        h = potential(k)
        proper = all(1 in reachable(k, n) for n in nodes)
        check((h is not None) == proper, '512 directed support census')
        if h is not None:
            check(audit_rank(k, h) < 1, 'exact spectral-radius certificate')
        graph_count += 1

    rec27 = (R(31, 27, (1,), (1, 2, 1)), R(27, 1, ROOT27, ()))
    loop27 = kernel((1, 27, 31), ((27, 31, F(1)), (31, 27, F(1))))
    check(compile_root(rec27[:1], loop27, 31) is None, 'closed 27/31 component')
    with27 = kernel((1, 27, 31), ((27, 31, F(1, 2)), (27, 1, F(1, 2)), (31, 27, F(1))))
    check(potential(with27) == (F(0), F(3), F(4)), 'supplied terminal opens component')
    out31 = compile_root(rec27, with27, 31)
    check(inherited.replay(31, out31['word']) == 1 and len(ROOT27) == 41 and sum(ROOT27) == 70,
          'fixed terminal provenance and 31 export')
    check(rejects(lambda: compile_root(rec27[:1], with27, 31)), 'missing terminal receipt rejected')
    check(rejects(lambda: compile_root((replace(rec27[1], left_word=ROOT27+(2,)),), with27, 27)),
          'ROOT-padding receipt rejected')
    for malformed in (True, 27.0, 0, 2):
        check(rejects(lambda n=malformed: compile_root(rec27, with27, n)), 'exact source typing')
    check(rejects(lambda: audit_kernel(Kernel((1, 3), ((F(0), F(0)), (F(1, 2), F(0)))))),
          'lost mass is not authentication')
    check(rejects(lambda: audit_kernel(Kernel((1, 3), ((F(0), F(0)), (1.0, F(0)))))),
          'float weight rejected')
    check(rejects(lambda: audit_rank(square, (F(0), F(1), F(1), F(1)))), 'false common rank')

    survival = F(1)
    for N in range(1, 257):
        p = F(1, (N+1)**2)
        survival *= 1-p
        check(survival == F(N+2, 2*(N+1)) and survival > F(1, 2), 'temporal positivity hostile')
        check((1-p)*(1/p) == 1/p-1, 'every frozen policy has a finite rank')
    # Exact finite gambler's-ruin harmonic controls, with limit 2^-n.
    r = F(1, 2)
    for ceiling in range(2, 33):
        hit = tuple((r**n-r**ceiling)/(1-r**ceiling) for n in range(ceiling+1))
        check(hit[0] == 1 and hit[-1] == 0, 'finite escape boundaries')
        for n in range(1, ceiling):
            check(hit[n] == F(1, 3)*hit[n-1]+F(2, 3)*hit[n+1] and 0 < hit[n] < r**n,
                  'infinite-path escape limit')
    check(F(3, 5)**2+F(4, 5)**2 == 1, 'positive sphere coordinates')
    print('PROVED/FINITE-EXACT: finite authenticated policy geometry; no universal Collatz claim.')
    print('Directed support census:', graph_count, '; interior square policies:81; robust vertices:4.')
    print('Square h(ROOT,3,5,7)=(0,8,5,9); robust h=(0,32,25,33); ROOT7 word', compiled['word'])
    print('Closed27/31 cannot compile. Supplied27 rank41/cost70 gives h=(0,3,4); ROOT31 rank/cost',
          len(out31['word']), sum(out31['word']))
    print('Adaptive p_t=1/(t+2)^2: survival_N=(N+2)/(2N+2) ->1/2 despite each frozen policy proper.')
    print('Uniform floor1/4: all-policy tail <=(63/64)^floor(t/3); 49 exact horizon controls.')
    print('Hidden-mode refinement:', refinement_count,
          'exact pushforward certificates; common rank gives tail <=h(source)*(32/33)^t.')
    print('Infinite biased walk: every state has ROOT path, but hit probability at n is2^-n.')
    print('Checks:', CHECKS)


if __name__ == '__main__':
    main()

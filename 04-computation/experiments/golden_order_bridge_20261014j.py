"""Exact marked-cycle/fence bridge and linear-extension boundary flow.

No external dependencies, discovery oracle, or presumed Collatz convergence.
Run: python -B golden_order_bridge_20261014j.py
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import permutations, product
from math import prod
import json


def integer(x, lower=0):
    if type(x) is not int or x < lower:
        raise ValueError("exact integer outside domain")
    return x


def bits(w, minimum=3):
    if type(w) is not tuple or len(w) < minimum:
        raise ValueError("bit tuple of sufficient length required")
    if any(type(b) is not int or b not in (0, 1) for b in w):
        raise ValueError("exact bits required")
    return w


def positive(x):
    if type(x) not in (int, Fraction) or x <= 0:
        raise ValueError("positive exact rational required")
    return Fraction(x)


@dataclass(frozen=True)
class Poset:
    n: int
    predecessors: tuple


def poset(n, edges):
    integer(n)
    if type(edges) is not tuple:
        raise ValueError("edge tuple required")
    pred = [0] * n
    for edge in edges:
        if type(edge) is not tuple or len(edge) != 2:
            raise ValueError("exact edge pair required")
        a, b = edge
        integer(a); integer(b)
        if a >= n or b >= n or a == b:
            raise ValueError("invalid strict edge")
        pred[b] |= 1 << a
    for k in range(n):
        for v in range(n):
            if pred[v] >> k & 1:
                pred[v] |= pred[k]
    if any(pred[v] >> v & 1 for v in range(n)):
        raise ValueError("directed cycle")
    return Poset(n, tuple(pred))


def audit_poset(p):
    if type(p) is not Poset:
        raise ValueError("Poset required")
    integer(p.n)
    if type(p.predecessors) is not tuple or len(p.predecessors) != p.n:
        raise ValueError("predecessor tuple required")
    for v, row in enumerate(p.predecessors):
        integer(row)
        if row >= 1 << p.n or row >> v & 1:
            raise ValueError("invalid strict predecessor mask")
    edges = tuple((u, v) for v in range(p.n) for u in range(p.n)
                  if p.predecessors[v] >> u & 1)
    if poset(p.n, edges) != p:
        raise ValueError("predecessors not transitively closed")
    return p


def _ideal(p, mask):
    return all(not (mask >> v & 1) or p.predecessors[v] & mask == p.predecessors[v]
               for v in range(p.n))


def is_ideal(p, mask):
    audit_poset(p); integer(mask)
    if mask >= 1 << p.n:
        raise ValueError("mask outside source")
    return _ideal(p, mask)


def fence(n):
    integer(n, 3)
    return poset(n, tuple((i, j) if i % 2 else (j, i)
                         for i, j in ((i, i + 1) for i in range(n - 1))))


@dataclass(frozen=True)
class CutState:
    word: tuple
    cut: int
    ideal: int
    seam: int


def compile_cut(word, cut=0):
    """Cut edge (cut-1,cut); retain a path state and its actual closing defect."""
    bits(word); integer(cut)
    n = len(word)
    if cut >= n:
        raise ValueError("cut outside source cycle")
    local = tuple(word[(cut + i) % n] for i in range(n))
    if any(local[i] * local[i+1] for i in range(n-1)):
        raise ValueError("internal no-11 constraint failed")
    mask = sum((b ^ (i % 2)) << i for i, b in enumerate(local))
    return CutState(word, cut, mask, local[0] * local[-1])


def audit_cut(s):
    if type(s) is not CutState:
        raise ValueError("CutState required")
    bits(s.word); integer(s.cut); integer(s.ideal); integer(s.seam)
    if s.seam not in (0, 1) or compile_cut(s.word, s.cut) != s:
        raise ValueError("source/ideal/seam mismatch")
    return s


def from_ideal(n, cut, mask):
    integer(n, 3); integer(cut); integer(mask)
    if cut >= n or not is_ideal(fence(n), mask):
        raise ValueError("invalid chart or ideal")
    local = tuple(((mask >> i) & 1) ^ (i % 2) for i in range(n))
    w = [0] * n
    for i, b in enumerate(local):
        w[(cut + i) % n] = b
    return compile_cut(tuple(w), cut)


def rotate(s, steps=1):
    audit_cut(s); integer(steps)
    n = len(s.word)
    k = steps % n
    return compile_cut(s.word[k:] + s.word[:k], (s.cut-k) % n)


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def activities(z):
    if type(z) is not tuple:
        raise ValueError("activity tuple required")
    return tuple(positive(x) for x in z)


def path_partition(z):
    z = activities(z)
    end0, end1 = Fraction(1), Fraction(0)
    for q in z:
        end0, end1 = end0 + end1, q * end0
    return end0 + end1


def cycle_partition(z):
    z = activities(z)
    if len(z) < 3:
        raise ValueError("simple cycle length at least three")
    m = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))
    for q in z:
        m = matmul(m, ((1, q), (1, 0)))
    return m[0][0] + m[1][1]


def closure_ledger(z):
    z = activities(z)
    if len(z) < 3:
        raise ValueError("simple cycle length at least three")
    defect = z[0] * z[-1] * path_partition(z[2:-2])
    return path_partition(z), defect, cycle_partition(z)


@dataclass(frozen=True)
class ExtensionFlow:
    source: Poset
    forward: tuple
    backward: tuple


def _flow(p):
    size = 1 << p.n
    f, g = [0] * size, [0] * size
    f[0] = g[-1] = 1
    for mask in range(size):
        if not _ideal(p, mask):
            continue
        for x in range(p.n):
            if not (mask >> x & 1) and p.predecessors[x] & mask == p.predecessors[x]:
                f[mask | (1 << x)] += f[mask]
    for mask in reversed(range(size-1)):
        if _ideal(p, mask):
            g[mask] = sum(g[mask | (1 << x)] for x in range(p.n)
                          if not (mask >> x & 1) and
                          p.predecessors[x] & mask == p.predecessors[x])
    return ExtensionFlow(p, tuple(f), tuple(g))


def extension_flow(p):
    audit_poset(p)
    return _flow(p)


def audit_flow(f):
    if type(f) is not ExtensionFlow:
        raise ValueError("ExtensionFlow required")
    audit_poset(f.source)
    for row in (f.forward, f.backward):
        if type(row) is not tuple or len(row) != 1 << f.source.n:
            raise ValueError("exact complete count tuple required")
        for x in row:
            integer(x)
    if f != _flow(f.source):
        raise ValueError("source-bound counts failed")
    return f


def rank_masses(f, rank):
    audit_flow(f); integer(rank)
    if rank > f.source.n:
        raise ValueError("rank outside source")
    z = f.forward[-1]
    return tuple((mask, Fraction(f.forward[mask]*f.backward[mask], z))
                 for mask in range(1 << f.source.n)
                 if mask.bit_count() == rank and f.forward[mask])


def rank_transport(f, rank, allowed):
    """Return retained mass and its normalized law; never silently lose the defect."""
    audit_flow(f)
    if type(allowed) is not tuple:
        raise ValueError("exact allowed-ideal tuple required")
    for mask in allowed:
        if not is_ideal(f.source, mask):
            raise ValueError("allowed state is not an ideal of the source")
    if len(set(allowed)) != len(allowed):
        raise ValueError("duplicated allowed state")
    masses = tuple((mask, mass) for mask, mass in rank_masses(f, rank) if mask in allowed)
    retained = sum((mass for _, mass in masses), Fraction(0))
    return retained, tuple((mask, mass/retained) for mask, mass in masses) if retained else ()


def comparison(f, x, y):
    """Uniform-linear-extension Pr(x occurs before y), via marked adding-x edges."""
    audit_flow(f); integer(x); integer(y)
    p = f.source
    if x >= p.n or y >= p.n or x == y:
        raise ValueError("two distinct source vertices required")
    numerator = sum(f.forward[mask]*f.backward[mask | (1 << x)]
                    for mask in range(1 << p.n)
                    if not (mask >> x & 1) and not (mask >> y & 1)
                    and p.predecessors[x] & mask == p.predecessors[x])
    return Fraction(numerator, f.forward[-1])


def fibonacci(n):
    integer(n)
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a


def swap_poset(n):
    """The finite thin-order segment: only consecutive labels are incomparable."""
    integer(n, 2)
    return poset(n, tuple((i, j) for i in range(n) for j in range(i+2, n)))


def swap_extension(word, cyclic=False):
    bits(word, 1)
    if type(cyclic) is not bool:
        raise ValueError("exact Boolean cycle flag required")
    if cyclic and len(word) < 3:
        raise ValueError("simple cycle length at least three")
    if any(word[i]*word[i+1] for i in range(len(word)-1)) or (cyclic and word[0]*word[-1]):
        raise ValueError("swaps must be vertex-disjoint")
    n = len(word) if cyclic else len(word)+1
    perm = list(range(n))
    for i, b in enumerate(word):
        if b:
            j = (i+1) % n
            perm[i], perm[j] = perm[j], perm[i]
    return tuple(perm)


def unswap_extension(perm):
    if type(perm) is not tuple or len(perm) < 2:
        raise ValueError("exact permutation tuple required")
    for x in perm:
        integer(x)
    if set(perm) != set(range(len(perm))):
        raise ValueError("not a permutation of the source labels")
    positions = {x: i for i, x in enumerate(perm)}
    word = tuple(int(positions[i+1] < positions[i]) for i in range(len(perm)-1))
    if swap_extension(word) != perm:
        raise ValueError("not an extension of the marked thin-order segment")
    return word


def swap_probability(n, edge):
    integer(n, 2); integer(edge)
    if edge >= n-1:
        raise ValueError("edge outside source path")
    return Fraction(fibonacci(edge+1)*fibonacci(n-edge-1), fibonacci(n+1))


def thin_gold_table(n):
    """All-N exact Gold table: first 0<1, then 1<2, true rows/columns first."""
    integer(n, 3)
    a, b = fibonacci(n-1), fibonacci(n-2)
    return ((a, b), (a, 0))


def gold_budget(table):
    if type(table) is not tuple or len(table) != 2:
        raise ValueError("exact two-row table required")
    for row in table:
        if type(row) is not tuple or len(row) != 2:
            raise ValueError("exact two-column row required")
        for value in row:
            integer(value)
    total = sum(map(sum, table))
    if not total:
        raise ValueError("empty source law")
    heavy, largest = max(map(sum, table)), max(v for row in table for v in row)
    return total, heavy, largest, heavy + largest <= total


def comparison_table(p, first, second):
    audit_poset(p)
    for pair in (first, second):
        if type(pair) is not tuple or len(pair) != 2:
            raise ValueError("exact comparison pair required")
        for x in pair:
            integer(x)
        if pair[0] == pair[1] or max(pair) >= p.n:
            raise ValueError("invalid comparison source")
    base = tuple((u, v) for v in range(p.n) for u in range(p.n) if p.predecessors[v] >> u & 1)
    rows = []
    for flip1 in (False, True):
        row = []
        for flip2 in (False, True):
            a = first[::-1] if flip1 else first
            b = second[::-1] if flip2 else second
            try:
                augmented = poset(p.n, base + (a, b))
            except ValueError:
                row.append(0)
            else:
                row.append(extension_flow(augmented).forward[-1])
        rows.append(tuple(row))
    return tuple(rows)


def rational_cycle(word):
    bits(word)
    if any(word[i] * word[(i+1) % len(word)] for i in range(len(word))):
        raise ValueError("cyclic no-11 word required")
    a, b = Fraction(1), Fraction(0)
    for bit in word:
        a, b = (3*a, 3*b+1) if bit else (a/2, b/2)
    source = b / (1-a)
    x = source
    for bit in word:
        if x.denominator % 2 == 0 or x.numerator % 2 != bit:
            raise ValueError("rational parity authentication failed")
        x = 3*x+1 if bit else x/2
    if x != source:
        raise ValueError("cycle did not close")
    return source


def golden_label(word):
    bits(word)
    if len(word) != 5:
        raise ValueError("only the proved eleven-state chart is exported")
    rational_cycle(word)
    return sum(word[j]*pow(4, 4-j, 11) for j in range(5)) % 11


def encode(x):
    if isinstance(x, Fraction):
        return str(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def main():
    checks = 0
    def check(claim):
        nonlocal checks
        checks += 1
        if not claim:
            raise RuntimeError("check failed at " + str(checks))
    def rejects(fn):
        try:
            fn()
        except (ValueError, TypeError):
            check(True)
        else:
            check(False)

    cycle_counts = {}
    for n in range(3, 14):
        p = fence(n)
        source_words = tuple(product((0, 1), repeat=n))
        path_words = tuple(w for w in source_words if all(not w[i]*w[i+1] for i in range(n-1)))
        closed = tuple(w for w in path_words if not w[0]*w[-1])
        ideals = tuple(mask for mask in range(1 << n) if _ideal(p, mask))
        check(len(ideals) == len(path_words))
        for w in path_words:
            s = compile_cut(w)
            check(is_ideal(p, s.ideal))
            check(from_ideal(n, 0, s.ideal) == s)
            for k in range(n):
                rotated = rotate(s, k)
                check(rotated.ideal == s.ideal and rotated.seam == s.seam)
                check(rotate(rotated, n-k) == s)
        gauges = sum(all((g[i] ^ g[(i+1) % n]) == 1 for i in range(n))
                     for g in source_words)
        check(gauges == (2 if n % 2 == 0 else 0))
        fib = [0, 1]
        for _ in range(n):
            fib.append(fib[-1] + fib[-2])
        check(len(closed) == fib[n-1] + fib[n+1])
        check(all(sum(w[i] for w in closed) == fib[n-1] for i in range(n)))
        check(3*fib[n-1] < len(closed))
        for z in (tuple([1]*n), tuple(Fraction(i % 3 + 1, i % 2 + 1) for i in range(n))):
            literal = sum(prod((z[i] for i in range(n) if w[i]), start=Fraction(1)) for w in closed)
            path, bad, cycle = closure_ledger(z)
            check(cycle == literal == path-bad)
        cycle_counts[n] = len(closed)

    # Independent full permutations versus DAG counts for every naturally labelled poset, n<=5.
    total_posets = 0
    for n in range(6):
        forward_edges = tuple((a, b) for a in range(n) for b in range(a+1, n))
        unique = set()
        for selection in range(1 << len(forward_edges)):
            edges = tuple(e for j, e in enumerate(forward_edges) if selection >> j & 1)
            unique.add(poset(n, edges))
        total_posets += len(unique)
        for p in sorted(unique, key=lambda q: q.predecessors):
            f = extension_flow(p)
            exts = []
            for perm in permutations(range(n)):
                position = {x: i for i, x in enumerate(perm)}
                if all(position[a] < position[b] for b in range(n) for a in range(n)
                       if p.predecessors[b] >> a & 1):
                    exts.append(perm)
            check(f.forward[-1] == len(exts) == f.backward[0])
            for rank in range(n+1):
                masses = dict(rank_masses(f, rank))
                check(sum(masses.values()) == 1)
                actual = {}
                for perm in exts:
                    mask = sum(1 << x for x in perm[:rank])
                    actual[mask] = actual.get(mask, 0) + 1
                check(masses == {k: Fraction(v, len(exts)) for k, v in actual.items()})
            # Source-preserving conservation at every ideal, not just total extension count.
            for mask in range(1 << n):
                if not f.forward[mask]:
                    continue
                through = f.forward[mask]*f.backward[mask]
                if mask:
                    incoming = sum(f.forward[mask ^ (1 << x)]*f.backward[mask]
                                   for x in range(n) if mask >> x & 1 and
                                   _ideal(p, mask ^ (1 << x)))
                    check(incoming == through)
                if mask != (1 << n)-1:
                    outgoing = sum(f.forward[mask]*f.backward[mask | (1 << x)]
                                   for x in range(n) if not (mask >> x & 1) and
                                   p.predecessors[x] & mask == p.predecessors[x])
                    check(outgoing == through)
            for x in range(n):
                for y in range(n):
                    if x != y:
                        check(comparison(f, x, y) == Fraction(sum(e.index(x) < e.index(y) for e in exts), len(exts)))

    p, f = fence(5), extension_flow(fence(5))
    golden = tuple(compile_cut(w) for w in product((0, 1), repeat=5)
                   if all(not w[i]*w[(i+1) % 5] for i in range(5)))
    labels = {golden_label(s.word): rational_cycle(s.word) for s in golden}
    check(set(labels) == set(range(11)))
    for s in golden:
        check(golden_label(rotate(s).word) == 4*golden_label(s.word) % 11)
    membership = tuple(Fraction(sum(s.ideal >> i & 1 for s in golden), 11) for i in range(5))
    check(membership == (Fraction(3, 11), Fraction(8, 11), Fraction(3, 11), Fraction(8, 11), Fraction(3, 11)))
    check(not any(Fraction(1, 3) <= q <= Fraction(2, 3) for q in membership))
    path_masks = tuple(i for i in range(32) if _ideal(p, i))
    ordinary_membership = tuple(Fraction(sum(mask >> i & 1 for mask in path_masks), 13) for i in range(5))
    check(ordinary_membership[0] == ordinary_membership[-1] == Fraction(5, 13))
    good = {s.ideal for s in golden}
    bad = sorted(set(path_masks)-good)
    check(bad == [27, 31])
    rank_retention = tuple(sum(mass for mask, mass in rank_masses(f, k) if mask in good) for k in range(6))
    check(f.forward[-1] == 16)
    check(rank_retention == (1, 1, 1, 1, Fraction(5, 8), 0))
    for k in range(6):
        retained, law = rank_transport(f, k, tuple(sorted(good)))
        check(retained == rank_retention[k])
        check(sum(mass for _, mass in law) == (1 if retained else 0))
        if k <= 3:
            check(law == rank_masses(f, k))
    check(comparison(f, 0, 4) == Fraction(1, 2))
    # The omitted rank-four state's incoming flow is 6/16; it is not 1/13 or 1/11.
    check(Fraction(f.forward[27]*f.backward[27], 16) == Fraction(3, 8))
    check(cycle_partition((2,)*5) == 31)
    check(sum(2**s.ideal.bit_count() for s in golden) == 73)

    swap_counts = {}
    for n in range(2, 12):
        q = swap_poset(n)
        qflow = extension_flow(q)
        words = tuple(w for w in product((0, 1), repeat=n-1)
                      if not any(w[i]*w[i+1] for i in range(n-2)))
        family = {swap_extension(w) for w in words}
        check(len(family) == fibonacci(n+1) == qflow.forward[-1])
        for word in words:
            perm = swap_extension(word)
            check(unswap_extension(perm) == word)
            positions = {x: i for i, x in enumerate(perm)}
            check(all(positions[a] < positions[b] for a in range(n) for b in range(a+2, n)))
        if n <= 7:
            literal = {perm for perm in permutations(range(n))
                       if all(perm.index(a) < perm.index(b) for a in range(n) for b in range(a+2, n))}
            check(literal == family)
        for edge in range(n-1):
            check(swap_probability(n, edge) == comparison(qflow, edge+1, edge))
            check(swap_probability(n, edge) == Fraction(sum(w[edge] for w in words), len(words)))
        swap_counts[n] = len(family)
        if n >= 3:
            table = thin_gold_table(n)
            check(comparison_table(q, (0, 1), (1, 2)) == table)
            total, heavy, largest, valid = gold_budget(table)
            check(valid and total == fibonacci(n+1) and heavy == fibonacci(n)
                  and largest == fibonacci(n-1) and heavy + largest == total)
    selected_extensions = {swap_extension(s.word) for s in golden}
    forced = tuple((a, b) for a in range(6) for b in range(6) if a != b
                   and all(perm.index(a) < perm.index(b) for perm in selected_extensions))
    check(poset(6, forced) == swap_poset(6))
    check(extension_flow(poset(6, forced)).forward[-1] == 13 > len(selected_extensions) == 11)
    missing_extension = swap_extension((1, 0, 0, 0, 1))
    check(missing_extension not in selected_extensions)
    selected_table = tuple(tuple(sum((perm.index(0) < perm.index(1)) == first and
                                    (perm.index(1) < perm.index(2)) == second
                                    for perm in selected_extensions)
                                for second in (True, False)) for first in (True, False))
    check(selected_table == ((5, 3), (3, 0)))
    check(gold_budget(selected_table) == (11, 8, 5, False))
    for a in range(6):
        for b in range(6):
            if a != b:
                q = Fraction(sum(perm.index(a) < perm.index(b) for perm in selected_extensions), 11)
                check(not Fraction(1, 3) <= q <= Fraction(2, 3))
    for n in range(3, 10):
        selected_words = tuple(w for w in product((0, 1), repeat=n)
                               if all(not w[i]*w[(i+1) % n] for i in range(n)))
        selected = tuple(swap_extension(w) for w in selected_words)
        common = tuple((a, b) for a in range(n+1) for b in range(n+1) if a != b
                       and all(perm.index(a) < perm.index(b) for perm in selected))
        check(poset(n+1, common) == swap_poset(n+1))
        check(len(selected) < fibonacci(n+2))

    h = {pow(4, j, 11) for j in range(5)}
    group = tuple((a, b) for a in sorted(h) for b in range(11))
    check(h == {1, 3, 4, 5, 9} == {pow(3, j, 11) for j in range(5)})
    check(len(group) == 55)
    images = {tuple(sorted((b, (a+b) % 11))) for a, b in group}
    check(len(images) == 55)
    check(all((v-u) % 11 in h for u, v in ((0, 1), (1, 2), (2, 0))))
    check(labels[0] == 0)
    check({labels[k] for k in h} == {Fraction(2**j, 13) for j in range(5)})
    check({labels[k] for k in set(range(1, 11))-h} == {-5, -14, -7, -20, -10})

    s = golden[0]
    for fn in (lambda: compile_cut((False, 0, 0)), lambda: compile_cut((0, 0, 0), True),
               lambda: compile_cut((1, 1, 0)), lambda: audit_cut(replace(s, ideal=float(s.ideal))),
               lambda: audit_cut(replace(s, seam=True)), lambda: from_ideal(5, 0, 1),
               lambda: poset(2, ((0, 1), (1, 0))), lambda: poset(2, ((True, 0),)),
               lambda: audit_poset(Poset(2, (0, True))), lambda: audit_flow(replace(f, backward=list(f.backward))),
               lambda: audit_flow(replace(f, forward=(True,)+f.forward[1:])),
               lambda: comparison(f, True, 4), lambda: rank_masses(f, 6),
               lambda: rank_transport(f, 1, (False,)), lambda: rank_transport(f, 1, (0, 0)),
               lambda: rank_transport(f, 1, (1,)),
               lambda: cycle_partition((1.0, 1, 1)), lambda: cycle_partition((1, 0, 1)),
               lambda: swap_extension((1, 1)), lambda: swap_extension((0, 1), 1),
               lambda: unswap_extension((False, 1)), lambda: unswap_extension((2, 1, 0)),
               lambda: swap_probability(5, True),
               lambda: thin_gold_table(True), lambda: gold_budget(((True, 1), (1, 0))),
               lambda: comparison_table(p, (False, 1), (1, 2)),
               lambda: golden_label((0, 0, 0)), lambda: golden_label((1, 1, 0, 0, 0))):
        rejects(fn)

    return {"status": "PROVED scoped + FINITE-EXACT; no new universal poset/Collatz claim",
            "checks": checks, "cycle_counts": cycle_counts,
            "golden_labels_to_actual_rational_sources": labels,
            "fence5": {"ideal_count": len(path_masks), "golden_count": len(golden),
                       "omitted_ideal_masks": bad, "golden_uniform_membership": membership,
                       "unconditioned_ideal_membership": ordinary_membership,
                       "linear_extensions": f.forward[-1], "rank_retention": rank_retention,
                       "activity2_true_mass": 31, "naive_ideal_cardinality_mass": 73,
                       "comparison_0_before_4": comparison(f, 0, 4)},
            "group": {"size": len(group), "multipliers": sorted(h), "paley_cycle": [0, 1, 2, 0]},
            "thin_order_swap_bridge": {"counts": swap_counts, "selected_six_label_extensions": len(selected_extensions),
                                       "their_forced_relations": forced, "extension_hull_size": 13,
                                       "missing_hull_extension": missing_extension,
                                       "Q6_gold_table": thin_gold_table(6),
                                       "selected_gold_table": selected_table,
                                       "all_N_gold_budget": "F_N+F_(N-1)=F_(N+1), equality",
                                       "bulk_probability_limit": "(5-sqrt(5))/10",
                                       "boundary_probability_limit": "(3-sqrt(5))/2"},
            "universes": {"cycles_all_words": [3, 13], "all_naturally_labelled_posets_orders": [0, 5],
                          "distinct_posets_tested": total_posets, "activity_laws_per_cycle": 2},
            "scope": "Marked-cut conditioning is lossless. Uniform ideal, weighted rank-prefix, and linear-extension pair laws differ. Period five refers to standard C, not shortcut T. Golay numerical coincidences import no additional conclusion."}


if __name__ == "__main__":
    print(json.dumps(encode(main()), indent=2, sort_keys=True))

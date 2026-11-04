"""Exact all-tail inversion of a positive linear objective on a fixed funnel core."""
from pathlib import Path
import importlib.util
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / '05-knowledge/results/polarized_path_calculus_20261103'))
import recognize as R
import explore as X
spec = importlib.util.spec_from_file_location(
    'component_inherited_counts',
    ROOT / '05-knowledge/results/orbit_spectrum_session_20261103/coverage/check.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def tail_bits(z):
    need(z >= 1, 'nonempty strong tail encoding')
    return tuple((z >> i) & 1 for i in range(z.bit_length()))


def before_after(z):
    """B_tail(1), independently specialized from the inherited six-state automaton."""
    states = [1, 0, 0, 0, 0, 0]
    for bit in tail_bits(z):
        nxt = [0] * 6
        for before in range(3):
            for after in range(2):
                weight = states[2 * before + after]
                nxt[2 * (1 if bit else 2) + after] += weight
                if after or not bit:
                    nxt[2 * before + 1] += weight
        states = nxt
    return sum(states[:4])


def all_data(a):
    """(h, e_w, P(2), E_w(2), S_w(2)) for every physical mark."""
    h, ends = c.tables(a)
    full = (1 << len(a)) - 1
    reverse = tuple(sum(1 << j for j in range(len(a)) if a[j] >> i & 1)
                    for i in range(len(a)))
    hr, starts = c.tables(reverse)
    need(h == hr, 'path reversal preserves every subset count')
    p2 = sum(h[m] * h[full ^ m] for m in range(full + 1))
    return [(h[full], ends[full][r], p2,
             sum(ends[m][r] * h[full ^ m] for m in range(1, full + 1)),
             sum(starts[m][r] * h[full ^ m] for m in range(1, full + 1)))
            for r in range(len(a))]


def forward(values, word):
    h, e, p, f, s = values
    if not word:
        return s + h - e, p - f
    need(word[-1] == 1, 'strong nonempty funnel tail must end in one')
    z = sum(bit << i for i, bit in enumerate(word))
    return (z + 1) * s + before_after(z) * h - int(z % 2 == 0) * e, (z + 1) * p - f


def solve_weighted(values, alpha, beta, target, allow_empty=True):
    """Complete all-tail solver for alpha*h_new + beta*F_new == target.

    Returned tested count covers candidate nonempty integers; the empty-tail
    boundary is checked separately and reported explicitly.
    """
    need(alpha >= 1 and beta >= 1, 'positive objective coefficients')
    h, e, p, f, s = values
    hits = []
    if allow_empty:
        nh, nf = forward(values, ())
        if alpha * nh + beta * nf == target:
            hits.append({'tail': '', 'h': nh, 'F': nf})
    k = alpha * s + beta * p
    t = target + beta * f
    denominator = k + 2 * alpha * h
    lower = max(2, (t + 2 * alpha * h + denominator - 1) // denominator)
    upper = (t + alpha * e) // k
    for length_code in range(lower, upper + 1):
        word = tail_bits(length_code - 1)
        nh, nf = forward(values, word)
        if alpha * nh + beta * nf == target:
            hits.append({'tail': ''.join(map(str, word)), 'h': nh, 'F': nf})
    return {'candidate_integers': max(0, upper - lower + 1),
            'empty_tail_checked': allow_empty, 'hits': hits}


def solve_outer(values, target, allow_empty=True):
    """All transitive prefix clocks s>=1, with the same all-tail completeness."""
    h, _, p, f, _ = values
    need(p - f >= h >= 1, 'positive clock bound')
    tested = clocks = 0
    hits = []
    clock, multiplier = 1, 2
    while h + multiplier * (p - f) <= target:
        result = solve_weighted(values, 1, multiplier, target, allow_empty)
        tested += result['candidate_integers']
        clocks += 1
        hits.extend(dict(s=clock, **hit) for hit in result['hits'])
        clock += 1
        multiplier *= 2
    return {'candidate_integers': tested, 'clocks': clocks, 'hits': hits}

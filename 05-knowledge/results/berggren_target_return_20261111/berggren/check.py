"""Independent literal controls for apex charts and ready exits.

No imported counting module.  No assert statements, so -O preserves checks.
Run: python -B check.py
"""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import json, math

HERE = Path(__file__).resolve().parent


def need(x, message):
    if not x:
        raise RuntimeError(message)


def graph(w):
    m = len(w)
    rows = [0] * (m + 1)
    for i in range(m):
        for j in range(i + 1, m):
            rows[i] |= 1 << j
        if w[i] == 'O':
            rows[m] |= 1 << i
        else:
            rows[i] |= 1 << m
    return tuple(rows)


def induced(rows, vertices):
    return tuple(sum(((rows[u] >> v) & 1) << j for j, v in enumerate(vertices))
                 for u in vertices)


def path_tables(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, amount in enumerate(dp[mask]):
            if not amount:
                continue
            rest = rows[v] & ~mask
            while rest:
                bit = rest & -rest
                rest -= bit
                dp[mask | bit][bit.bit_length() - 1] += amount
    totals = [sum(x) for x in dp]
    totals[0] = 1
    return dp, totals


def H(rows):
    return path_tables(rows)[1][-1]


def root_profile(rows, root):
    ends, totals = path_tables(rows)
    full = len(totals) - 1
    converse = tuple(sum(((rows[u] >> v) & 1) << u for u in range(len(rows)))
                     for v in range(len(rows)))
    starts, _ = path_tables(converse)
    f_end = sum(ends[mask][root] * totals[full ^ mask]
                for mask in range(1, full + 1) if mask >> root & 1)
    return {'H': totals[-1], 'S': starts[-1][root], 'E': ends[-1][root],
            'd': totals[full ^ (1 << root)], 'F': f_end,
            'indegree': sum(row >> root & 1 for row in rows)}


def acyclic_order(rows):
    """Kahn elimination, independently of cut-word or degree formulas."""
    rem = set(range(len(rows)))
    order = []
    while rem:
        sources = [v for v in rem if all(v == u or not (rows[u] >> v & 1) for u in rem)]
        if len(sources) != 1:
            return None
        v = sources[0]
        order.append(v)
        rem.remove(v)
    return tuple(order)


def apices(rows):
    return [r for r in range(len(rows))
            if acyclic_order(induced(rows, [v for v in range(len(rows)) if v != r])) is not None]


def strong(rows):
    if not rows:
        return False
    for r in range(len(rows)):
        seen = 1 << r
        while True:
            more = seen
            for v, out in enumerate(rows):
                if seen >> v & 1:
                    more |= out
            if more == seen:
                break
            seen = more
        if seen != (1 << len(rows)) - 1:
            return False
    return True


def canonical(rows):
    """Independent degree-partition canonical adjacency string and automorphisms."""
    degrees = sorted(set(map(int.bit_count, rows)))
    groups = [[v for v, row in enumerate(rows) if row.bit_count() == d] for d in degrees]
    best = None
    multiplicity = 0
    for parts in product(*(permutations(g) for g in groups)):
        order = sum((tuple(part) for part in parts), ())
        key = tuple(sum(((rows[v] >> u) & 1) << i for i, u in enumerate(order)) for v in order)
        if best is None or key < best:
            best, multiplicity = key, 1
        elif key == best:
            multiplicity += 1
    return best, multiplicity


def birth(rows, parent):
    """The new vertex loses only to parent."""
    n = len(rows)
    return tuple(row | ((1 << n) if v == parent else 0) for v, row in enumerate(rows)) + (((1 << n) - 1) ^ (1 << parent),)


def state(w):
    h, p, e = 1, 0, 0
    for letter in w:
        if letter == 'O':
            p, e = 2 * p + 1, 1
        else:
            h, p, e = h + p + e, 2 * p, 0
    return h, p, e, 1


SUBS = {'A': {'O': 'OI', 'I': 'I'},
        'B': {'O': 'OI', 'I': 'OIO'},
        'C': {'O': 'O', 'I': 'OOI'}}


def branch(w, letter):
    return ''.join(SUBS[letter][c] for c in w)


def census():
    data = []
    for n in range(1, 7):
        classes, connected = {}, set()
        labeled = 0
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        for code in range(1 << len(pairs)):
            rows = [0] * n
            for k, (i, j) in enumerate(pairs):
                if code >> k & 1:
                    rows[i] |= 1 << j
                else:
                    rows[j] |= 1 << i
            if not apices(rows):
                continue
            labeled += 1
            key, aut = canonical(rows)
            classes[key] = aut
            if strong(rows):
                connected.add(key)
        expected = 1 if n == 1 else (1 << (n - 1)) - 1 - (n - 2) * (n - 1) // 2
        expected_strong = 1 if n in (1, 3) else (0 if n == 2 else (1 << (n - 3)) - 1)
        need(len(classes) == expected, 'ordinary apex class census')
        need(len(connected) == expected_strong, 'strong apex class census')
        need(sum(math.factorial(n) // a for a in classes.values()) == labeled, 'automorphism mass')
        if n >= 3:
            need(labeled * 3 == math.factorial(n) * (3 * expected - 2 * (n - 2)), 'labeled closed form')
        data.append({'n': n, 'universe': 1 << len(pairs), 'labeled_one_apex': labeled,
                     'ordinary': len(classes), 'strong_ordinary': len(connected),
                     'automorphism_distribution': dict(Counter(classes.values()))})
    return data


def controls():
    counts = Counter()
    for m in range(9):
        for letters in product('OI', repeat=m):
            w = ''.join(letters)
            rows = graph(w)
            h, p, e, _ = state(w)
            need(H(rows) == h, 'independent path DP versus four-state count')
            if m:
                profile = root_profile(rows, m)
                q = int(''.join('1' if c == 'I' else '0' for c in reversed(w)), 2)
                need(profile == {'H': h, 'S': int(w[0] == 'O'), 'E': int(w[-1] == 'I'),
                                 'd': 1, 'F': q + 1, 'indegree': w.count('I')}, 'apex rooted interface')
            if strong(rows) and m >= 2:
                need(w[0] == 'O' and w[-1] == 'I', 'strong cut boundary')
                expected_apices = 3 if m == 2 else (2 if min(w.count('O'), w.count('I')) == 1 else 1)
                need(len(apices(rows)) == expected_apices, 'complete admissible apex classification')
                profile = root_profile(rows, 0)
                need(profile['indegree'] == 1, 'first base vertex indegree one')
                need(profile['F'] - h == int(w[1] == 'O'), 'first-base ready defect')
                if w.startswith('OO') and m <= 6:
                    child = birth(birth(rows, 0), len(rows))
                    need(H(child) == 3 * h + 2, 'literal ready exit target')
                    need(not apices(child), 'reset leaves apex class')
                    old_apex = m
                    deleted = induced(child, [v for v in range(len(child)) if v != old_apex])
                    need(H(deleted) == 3, 'C3 plus transitive deletion core')
                    # Two explicitly vertex-disjoint cycles: (old_apex,1,last), (0,x,y).
                    for cycle in ((m, 1, m - 1), (0, m + 1, m + 2)):
                        need(all(child[a] >> b & 1 for a, b in zip(cycle, cycle[1:] + cycle[:1])), 'disjoint cycle witness')
                    counts['ready_exits'] += 1
                counts['strong_cuts'] += 1
            counts['all_cuts'] += 1
    witnesses = []
    for addr, w in [('B', 'OIOIO'), ('C', 'OOOI'), ('AA', 'OIII')]:
        children = {a: state(branch(w, a))[0] for a in 'ABC'}
        need(H(graph(w)) == 9, 'H9 source actual count')
        witnesses.append({'address': addr, 'word': w, 'rows': graph(w), 'state': state(w),
                          'apex_profile': root_profile(graph(w), len(w)),
                          'strong': strong(graph(w)), 'children': children})
    need([w['children']['A'] for w in witnesses] == [107, 73, 17], 'requested collision')
    need(canonical(graph('OOOI'))[0] == canonical(graph('OIII'))[0], 'actual C-AA isomorphism')
    need(canonical(graph('OOOI'))[0] != canonical(graph('OIOI'))[0], 'two distinct H9 cores')
    return dict(counts), witnesses


def reverse_complement_value(w):
    return sum((c == 'I') << i for i, c in enumerate(w))


def full_state(rows, root):
    p = root_profile(rows, root)
    converse = tuple(sum(((rows[u] >> v) & 1) << u for u in range(len(rows)))
                     for v in range(len(rows)))
    cp = root_profile(converse, root)
    base = induced(rows, [v for v in range(len(rows)) if v != root])
    totals = path_tables(base)[1]
    k = sum(totals[mask] * totals[(len(totals) - 1) ^ mask] for mask in range(len(totals)))
    return [p['H'], p['E'], p['S'], cp['F'], p['F'], p['d'], k]


def future_state(frame):
    h, e, s, u, v, d, k = frame
    return [h, u, u - d, h + u - e]


def rewrite_controls():
    witnesses = []
    for t in range(1, 5):
        for z, coefficient, next_defect in [('OOI', 75, 65), ('OII', 105, 91)]:
            suffix = 'I' * t + z
            w = 'OOOO' + suffix
            target_word = 'OIOOIO' + suffix
            source = graph(w)
            apex = len(w)
            x, y = len(source), len(source) + 1
            first = birth(source, 0)
            standard = birth(first, x)
            rewritten = list(standard)
            for v in (0, 1, 2):
                need(rewritten[y] >> v & 1, 'rewrite old orientation')
                rewritten[y] ^= 1 << v
                rewritten[v] |= 1 << y
            rewritten = tuple(rewritten)
            order = [0, x, 1, 2, y, 3, *range(4, apex), apex]
            need(induced(rewritten, order) == graph(target_word), 'literal returned cut chart')
            h = coefficient * (1 << t) - 9
            target = 3 * h + 2
            need(H(source) == h and H(standard) == target and H(rewritten) == target, 'exact three-arc neutral rewrite')
            need(strong(standard) and strong(rewritten), 'both outputs strong')
            need(not apices(standard) and apices(rewritten) == [apex], 'feedback grade changes from two to one')
            profile = root_profile(rewritten, x)
            need(profile['indegree'] == 1 and profile['F'] - target == next_defect * (1 << t) - 4,
                 'returned indegree-one root growing defect')
            before, after = full_state(standard, apex), full_state(rewritten, apex)
            p, q, power = state(w)[1], reverse_complement_value(w), 1 << len(w)
            need(before == [target, 3, 1, 3 * (p + 1), 6 * q + 9, 3, 6 * power], 'H3-core full old-apex interface')
            delta = [a - b for a, b in zip(future_state(before), future_state(after))]
            pv = state(suffix)[1]
            need(delta == [0, 2 * pv + 2, 2 * pv, 2 * pv], 'uniform future-state drop')
            need(root_profile(rewritten, 0)['F'] - target == 0, 'first base root becomes zero-ready')
            witnesses.append({'t': t, 'z': z, 'source_word': w, 'returned_word': target_word,
                              'source_rows': source, 'standard_rows': standard, 'returned_rows': rewritten,
                              'old_apex': apex, 'x': x, 'y': y, 'returned_base_order_then_apex': order,
                              'trace': [h, H(first), target, target], 'returned_x_profile': profile,
                              'old_apex_before': before, 'old_apex_after': after, 'future_drop': delta})
    # Complete exact eligibility through suffix length nine, including rejected near misses.
    accepted = []
    for length in range(1, 10):
        for bits in product('OI', repeat=length - 1):
            suffix = ''.join(bits) + 'I'
            old = state('OOOO' + suffix)[0]
            new = state('OIOOIO' + suffix)[0]
            rule = suffix.startswith('I') and state(suffix)[0] == 5
            need((new == 3 * old + 2) == rule, 'all-suffix exact eligibility')
            if rule:
                core = suffix.lstrip('I')
                need(core in ('OOI', 'OII'), 'complete H5 suffix core')
                accepted.append(suffix)
    # The first source C(OI)=OOOI exits to H29, outside the entire apex class.
    # The inherited strong-word bound H >= 2^(m-2)+1 reduces H29 to m<=6.
    h29 = []
    for m in range(2, 7):
        for body in product('OI', repeat=m - 2):
            w = 'O' + ''.join(body) + 'I'
            if state(w)[0] == 29:
                h29.append(w)
    need(not h29, 'all-order H29 one-apex obstruction with inherited length bound')
    # A repeat using the same suffix and P -> 3P would require a negative suffix count.
    need(state('OOOO') == (1, 15, 1, 1), 'first prefix')
    need(state('OIOOIO') == (15, 45, 1, 1), 'returned prefix')
    need(state('OIIIIOOO') == (17, 135, 1, 1), 'second multiplied prefix')
    return {'witnesses': witnesses, 'accepted_suffixes_through_length9': accepted,
            'H29_strong_cut_bound': 6, 'H29_witnesses': h29,
            'second_prefix_H_difference': 17 - 3 * 15}


def main():
    small, witnesses = controls()
    result = {'status': 'PASS', 'universe': 'All labeled tournaments through order6; all binary cuts through length8; all strong OO-start cuts through length6 for literal two-birth controls.',
              'census': census(), 'controls': small, 'H9_sources': witnesses,
              'neutral_three_arc_returns': rewrite_controls()}
    (HERE / 'results.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'census': result['census'], 'controls': small}, indent=2))


if __name__ == '__main__':
    main()

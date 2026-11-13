"""Literal same-target wrapper preparation and accepted return clocks.

No repository implementation is imported.  The universal proof is in PROOF.md;
the finite controls here use subset endpoint DP and, at orders five and six,
independent permutation/cut enumeration.
"""
from itertools import permutations
from math import gcd
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def require(test, message):
    if not test:
        raise RuntimeError(message)


def wrapper(n, position=0):
    """Physical root 0 sits at the declared position in the TT_n module."""
    order = list(range(1, position + 1)) + [0] + list(range(position + 1, n))
    rows = [0] * (n + 2)
    for i, u in enumerate(order):
        for v in order[i + 1:]:
            rows[u] |= 1 << v
        rows[u] |= 1 << n
    rows[n] = 1 << (n + 1)
    rows[n + 1] = (1 << n) - 1
    return tuple(rows), tuple(order)


def endpoint_dp(rows):
    n = len(rows)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, full + 1):
        for last, count in enumerate(dp[mask]):
            if count:
                choices = rows[last] & (full ^ mask)
                while choices:
                    bit = choices & -choices
                    choices ^= bit
                    dp[mask | bit][bit.bit_length() - 1] += count
    return dp


def state(rows, root=0):
    full = (1 << len(rows)) - 1
    end = endpoint_dp(rows)
    start = endpoint_dp(tuple(full ^ (1 << i) ^ r for i, r in enumerate(rows)))
    counts = [sum(a) for a in end]
    counts[0] = 1
    return (counts[full], start[full][root], counts[full ^ (1 << root)],
            sum(end[m][root] * counts[full ^ m] for m in range(1, full + 1)))


def legal(rows, word):
    return all(rows[u] >> v & 1 for u, v in zip(word, word[1:]))


def permutation_state(rows, root=0):
    n = len(rows)
    hp = [p for p in permutations(range(n)) if legal(rows, p)]
    other = [v for v in range(n) if v != root]
    q = sum(legal(rows, p) for p in permutations(other))
    covers = 0
    for p in permutations(range(n)):
        cut = p.index(root) + 1
        covers += legal(rows, p[:cut]) and legal(rows, p[cut:])
    return len(hp), sum(p[0] == root for p in hp), q, covers


def strong(rows):
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(r == (1 << len(rows)) - 1 for r in reach)


def birth(rows, root, new_root=False):
    n = len(rows)
    out = list(rows)
    out[root] |= 1 << n
    out.append(((1 << n) - 1) ^ (1 << root))
    return tuple(out), n if new_root else root


def native(rows, path):
    require(legal(rows, path), 'illegal native word')
    retained = {frozenset((u, v)) for u, v in zip(path, path[1:])}
    out = list(rows)
    for u in range(len(rows)):
        for v in range(u + 1, len(rows)):
            if frozenset((u, v)) not in retained:
                out[u] ^= 1 << v
                out[v] ^= 1 << u
    return tuple(out)


def native_prepare(n):
    rows, _ = wrapper(n)
    temporary = len(rows)
    # The new vertex beats only root 0.  These are the inherited compiler's
    # actual temporary cut and two legal Hamilton words, on successive inputs.
    ext = [r | (1 << temporary) if i != 0 else r for i, r in enumerate(rows)]
    ext.append(1)
    first = tuple(range(n + 2)) + (temporary,)
    second = (0, temporary) + tuple(range(1, n + 2))
    middle = native(tuple(ext), first)
    final = native(middle, second)
    mask = (1 << temporary) - 1
    endpoint = tuple(r & mask for r in final[:temporary])
    require(endpoint == wrapper(n, 1)[0], 'native compiler changed extra old arcs')
    return {'temporary_vertex': temporary, 'temporary_outmask': 1,
            'first_word': first, 'second_word': second, 'endpoint_rows': endpoint}


def formula(n, j):
    return (2**n + 1, 2**(n - 1 - j), 2**(n - 1) + 1,
            2**(n - j) * 3**j + 2**(j + 1) - (j == n - 1))


def returned(n, k):
    L = 2**(n - 2)
    t = (4 * L + 3) * 2**k + L - 1
    b = (12 * L + 9) * 2**k + 5 * L + 1
    return b, t


def literal_return(n, k):
    rows, _ = wrapper(n, 1)
    root = 0
    for _ in range(k):
        rows, root = birth(rows, root)
    rows, root = birth(rows, root, True)
    before = state(rows, root)
    rows, root = birth(rows, root)
    z = state(rows, root)
    b, t = returned(n, k)
    require(z == (b, t, t, b), 'literal return disagrees with formula')
    require(gcd(*before) == 1 and strong(rows), 'pair-content or strongness failure')
    L = 2**(n - 2)
    other = state(rows, len(rows) - 1)
    require(other == (b, t + L + 2, 2 * t - L + 1, 4 * t + 3)
            and gcd(*other) == 1, 'other retained port is not primitive')
    return {'n': n, 'k': k, 'rows': rows, 'root': root,
            'state_h_s_q_F': z, 'pre_last_R_state': before,
            'other_retained_root': len(rows) - 1, 'other_port_state': other}


def main():
    profiles = []
    for n in range(3, 8):
        for j in range(n):
            rows, order = wrapper(n, j)
            z = state(rows)
            require(z == formula(n, j), 'wrapper profile formula failed')
            require(strong(rows), 'wrapper is not strong')
            if n <= 4:
                require(z == permutation_state(rows), 'independent enumeration failed')
            profiles.append({'n': n, 'position': j, 'rows': rows, 'module_order': order,
                             'root': 0, 'indegree': j + 1, 'state_h_s_q_F': z})
    compilers = [{'n': n, **native_prepare(n)} for n in range(3, 8)]
    cases = [(3, 1), (3, 2), (4, 1), (4, 2), (4, 3), (5, 2), (6, 3), (7, 4)]
    returns = [literal_return(n, k) for n, k in cases]
    for r in returns:
        b, t, _, _ = r['state_h_s_q_F']
        require(gcd(b, t) == 1, 'accepted return is not primitive')
        c, u = 61 * b - 2 * t, 9 * b
        require(gcd(c, u) == 1 and c % 3 and u % 3 == 0, 'C3 entry failed')
        r['C3_profile_parameters'] = [c, u]
    # This bounded scalar check exercises the all-parameter certificate only;
    # it is not a new census or a substitute for its modular proof.
    for n in range(3, 41):
        for multiplier in (1, 2, 3):
            k = multiplier * max(1, n - 3)
            require(gcd(*returned(n, k)) == 1, 'factor-free clock certificate failed')
    hostile = returned(7, 6)
    require(hostile == (25313, 8415) and gcd(*hostile) == 17, 'wrong unrestricted-clock hostile')
    out = {'status': 'FINITE-EXACT controls; all-parameter theorem in PROOF.md',
           'universe': 'W_n root positions n3..7; eight literal returns through order15',
           'profiles': profiles, 'native_preparations': compilers, 'literal_returns': returns,
           'off_clock_hostile': {'n': 7, 'k': 6, 'balanced_b_t': hostile,
                                'returned_root_gcd': 17,
                                'scope': 'proved state formula, not a 19-vertex DP recount'}}
    (HERE / 'wrapper_repair_receipt.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'passed': True, 'position_profiles': len(profiles),
                      'literal_return_counts': [r['state_h_s_q_F'][0] for r in returns],
                      'off_clock_gcd': 17}))


if __name__ == '__main__':
    main()

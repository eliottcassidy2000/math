"""Independent audit of the sibling all-k target return identities."""
from functools import lru_cache
from pathlib import Path
import hashlib, json
from check import need, H, strong, induced, birth, root_profile, full_state

HERE = Path(__file__).resolve().parent
RESET = HERE.parent / 'reset'


def inflate(rows, module, k):
    names = [i for i in range(len(rows)) for _ in range(k if i == module else 1)]
    out = []
    for i, u in enumerate(names):
        out.append(sum((int(i < j) if u == v else (rows[u] >> v & 1)) << j
                       for j, v in enumerate(names) if i != j))
    return tuple(out)


def coefficient_polynomial(rows, module):
    """Backward suffix DP on external subsets and an indistinguishable run token.

    Each visit to X is one nonempty run, adjacent X's are prohibited, and the
    other vertices occur once.  Coefficient j is subsequently multiplied by
    j! S(k,j). This does not import the forward inherited module-run engine.
    """
    n = len(rows)
    ext = [v for v in range(n) if v != module]
    X, sentinel = n, n + 1
    full = (1 << len(ext)) - 1

    @lru_cache(None)
    def suffix(mask, last):
        poly = [int(mask == full)]

        def add(other, shift):
            while len(poly) < len(other) + shift:
                poly.append(0)
            for i, value in enumerate(other):
                poly[i + shift] += value

        for i, vertex in enumerate(ext):
            if mask >> i & 1:
                continue
            if last == sentinel or (rows[module if last == X else last] >> vertex & 1):
                add(suffix(mask | (1 << i), vertex), 0)
        if last != X and (last == sentinel or rows[last] >> module & 1):
            add(suffix(mask, X), 1)
        return tuple(poly)

    result = list(suffix(0, sentinel))
    while len(result) < n + 1:
        result.append(0)
    return result[1:]


def marked_polynomials(P, module):
    n = len(P)
    sep = tuple(row | (1 << n) for row in P) + ((1 << n) - 1,)
    sources = {
        'parent': P,
        'child': birth(P, 0),
        'two_covers': sep,
        'start_covers': tuple(row & ~1 for row in sep),
        'end_paths': tuple(0 if i == 0 else row for i, row in enumerate(P)),
        'start_paths': tuple(row & ~1 for row in P),
    }
    return {key: coefficient_polynomial(T, module) for key, T in sources.items()}


def main():
    data = json.loads((RESET / 'target_ready_return.json').read_text())
    polys = []
    direct = []
    for stage, (P, module) in enumerate(zip(data['parent_templates'], data['module_vertices'])):
        p = marked_polynomials(tuple(P), module)
        for key, expected in data['marked_runs'][stage].items():
            observed = p[key]
            # Independent engines may retain a different number of trailing zeros.
            length = max(len(observed), len(expected))
            need(observed + [0] * (length - len(observed)) == expected + [0] * (length - len(expected)),
                 'complete independent placeholder polynomial: ' + key)
        polys.append(p)
    for rec in data['controls'][:3]:
        k = rec['k']
        for stage, key in enumerate(('ready221', 'ready665')):
            P = inflate(data['parent_templates'][stage], data['module_vertices'][stage], k)
            R = birth(P, 0)
            need(list(R) == rec[key], 'independent whole-module inflation and birth')
            need(strong(P) and strong(R), 'independent strongness')
            need(root_profile(P, 0)['S'] == 1, 'independent unique-start parent')
            need(full_state(R, len(R) - 1) == rec['states'][stage], 'independent seven-state literal DP')
            need(H(R) == ((99 * (1 << k) + 23) if stage == 0 else (297 * (1 << k) + 71)), 'absolute target')
            direct.append({'k': k, 'stage': stage, 'order': len(R), 'H': H(R)})
        S = tuple(rec['source'])
        R, Q = tuple(rec['ready221']), tuple(rec['ready665'])
        need(induced(S, rec['first_kept']) == induced(R, rec['first_map']), 'same-gauge retained map')
        converse = tuple(sum((R[u] >> v & 1) << u for u in range(len(R))) for v in range(len(R)))
        need(induced(S, rec['converse_kept']) == induced(converse, rec['converse_map']), 'converse retained map')
        F = birth(birth(R, len(R) - 1), len(R))
        need(H(F) == H(Q), 'second same-count target')
        need(induced(F, rec['second_kept']) == induced(Q, rec['second_map']), 'second retained map')
    inverse = json.loads((RESET / 'ready_two_inverse.json').read_text())
    need(inverse['bank_rows'] == 32605, 'declared two-entry bank size')
    for target in ('1997', '2005'):
        need(not inverse['hits'][target], 'declared two-ready exclusion')
    # Independent arithmetic for the analytic order-tail gate.
    minimum16 = 5 ** 5
    need(2 * minimum16 > 5993 + 2, 'strong parent order >=16 excluded for N<=5993,b<=2')
    need(min(9545, 18905) > 5993, 'inherited sharp ready-order15,16 minima exclude remaining tail')
    pins = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (RESET / 'target_ready_return.py', RESET / 'target_ready_return.json',
             RESET / 'ready_two_inverse.py', RESET / 'ready_two_inverse.json')}
    result = {'status': 'PASS', 'scope': 'Independent complete placeholder coefficient polynomials; literal k1..3 graph/profile/map audit; analytic all-order tail arithmetic. The full 32605-row inverse census is inherited, not independently replayed here.',
              'complete_polynomials': polys, 'literal_controls': direct, 'source_sha256': pins}
    (HERE / 'reset_audit.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print('PASS: independent complete run polynomials, six literal ready graphs, both target-preserving map gauges, and order-tail arithmetic')


if __name__ == '__main__':
    main()

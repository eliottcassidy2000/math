"""Independent literal small-cover controls and pinned large-response audit."""
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = HERE.parents[3]
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def decode(n, code):
    rows = [0] * n
    bit = 0
    for u in range(n):
        for v in range(u + 1, n):
            if code >> bit & 1:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
            bit += 1
    return tuple(rows)


def opposite(rows):
    return tuple(((1 << len(rows)) - 1) ^ (1 << u) ^ row for u, row in enumerate(rows))


def inflate(q, root, m):
    blocks = []
    offset = 0
    for z in range(len(q)):
        size = len(m) if z == root else 1
        blocks.append(tuple(range(offset, offset + size)))
        offset += size
    rows = [0] * offset
    for z, block in enumerate(blocks):
        if z == root:
            for i, u in enumerate(block):
                for j, v in enumerate(block):
                    if m[i] >> j & 1:
                        rows[u] |= 1 << v
        for w, other in enumerate(blocks):
            if q[z] >> w & 1:
                for u in block:
                    for v in other:
                        rows[u] |= 1 << v
    return tuple(rows), blocks


@lru_cache(maxsize=None)
def literal(rows):
    n = len(rows)
    s, e, b2, e2 = ([0] * n for _ in range(4))
    c = [[0] * n for _ in range(n)]
    for p in permutations(range(n)):
        bad = [k for k in range(n - 1) if not rows[p[k]] >> p[k + 1] & 1]
        if not bad:
            s[p[0]] += 1
            e[p[-1]] += 1
        if len(bad) > 1:
            continue
        cuts = bad if bad else range(n - 1)
        for k in cuts:
            c[p[k]][p[k + 1]] += 1
            b2[p[0]] += 1
            e2[p[k]] += 1
    return {'S': s, 'E': e, 'C': c, 'B2': b2, 'E2': e2, 'H': sum(s)}


def hamilton(rows):
    n = len(rows)
    pred = opposite(rows)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for mask in range(1, full + 1):
        vs = mask
        while vs:
            bit = vs & -vs
            vs -= bit
            v = bit.bit_length() - 1
            rest = mask ^ bit
            before = rest & pred[v]
            value = int(not rest)
            while before:
                ubit = before & -before
                before -= ubit
                value += dp[rest][ubit.bit_length() - 1]
            dp[mask][v] = value
    return sum(dp[full])


def actual_query(rows, a, c):
    n = len(rows)
    new = list(rows) + [0, 0]
    for i in range(n):
        for port, mask in ((n, a), (n + 1, c)):
            if mask >> i & 1:
                new[i] |= 1 << port
            else:
                new[port] |= 1 << i
    new[n] |= 1 << (n + 1)
    return hamilton(tuple(new))


def literal_controls():
    twins, crosses, inherit, triangles = 0, 0, 0, 0
    for n in range(1, 5):
        for code in range(1 << (n * (n - 1) // 2)):
            q = decode(n, code)
            d = literal(q)
            for root in range(n):
                t, blocks = inflate(q, root, (2, 0))
                f, r = literal(t), literal(opposite(t))
                a, b = blocks[root]
                need(f['S'][a] - f['S'][b] == d['S'][root], 'literal twin start difference')
                need(f['E'][b] - f['E'][a] == d['E'][root], 'literal twin end difference')
                twins += 1
                for z in range(n):
                    if z == root:
                        continue
                    w = blocks[z][0]
                    first = f['C'][b][w] - f['C'][a][w]
                    second = f['C'][w][a] - f['C'][w][b]
                    need(first == d['C'][root][z] and second == d['C'][z][root], 'literal twin two-cover contractions')
                    reverse = r['C'][b][w] - r['C'][a][w] + r['C'][w][a] - r['C'][w][b]
                    need(first + second > 0 and reverse < 0, 'opposite twin contrast signs')
                    crosses += 1
    for qn in range(1, 4):
        for qc in range(1 << (qn * (qn - 1) // 2)):
            q = decode(qn, qc)
            qd = literal(q)
            for root in range(qn):
                for mn in range(1, 4):
                    for mc in range(1 << (mn * (mn - 1) // 2)):
                        m = decode(mn, mc)
                        md = literal(m)
                        t, blocks = inflate(q, root, m)
                        td = literal(t)
                        for u in range(mn):
                            for v in range(mn):
                                a, b = blocks[root][u], blocks[root][v]
                                need((td['C'][a][b] - td['C'][b][a] - qd['H'] * (md['C'][u][v] - md['C'][v][u])) % 4 == 0, 'module C-skew inheritance identity')
                        for z in range(qn):
                            if z != root:
                                w = blocks[z][0]
                                need((td['S'][w] - qd['S'][z]) % 2 == 0 and (td['E'][w] - qd['E'][z]) % 2 == 0, 'outside endpoint parity inheritance')
                        if all((s - e) % 2 == 0 for s, e in zip(td['S'], td['E'])):
                            need(all((s - e) % 2 == 0 for s, e in zip(qd['S'], qd['E'])), 'balanced quotient implication')
                        inherit += 1
    small_modules = [decode(n, code) for n in range(1, 5)
                     for code in range(1 << (n * (n - 1) // 2))]
    small_modules += [decode(5, 0), decode(5, 741), (6, 12, 24, 17, 3)]
    for m in small_modules:
        md = literal(m)
        t, blocks = inflate((2, 4, 1), 0, m)
        td = literal(t)
        for u, w in enumerate(blocks[0]):
            need(td['S'][w] == md['S'][u] + md['B2'][u], 'exact C3 quotient start count')
            need(td['E'][w] == md['E'][u] + md['E2'][u], 'exact C3 quotient end count')
            need((td['S'][w] - md['E'][u]) % 2 == 0 and (td['E'][w] - md['S'][u]) % 2 == 0, 'C3 quotient swapped parity')
        triangles += 1
    hostile = (10, 20, 33, 22, 37, 11)
    hd = literal(hostile)
    need(hd['H'] == 45 and hd['S'] == [7, 7, 7, 8, 8, 8] and hd['E'] == [8, 8, 8, 7, 7, 7], 'unique6 literal endpoints')
    need(all((hd['C'][u][v] - hd['C'][v][u]) % 4 == 0 for u in range(6) for v in range(6)), 'unique6 literal C residues')
    for a, b in ((5, 0), (3, 1), (4, 2)):
        need(hostile[a] >> b & 1 and all((hostile[a] >> w & 1) == (hostile[b] >> w & 1) for w in range(6) if w not in (a, b)), 'unique6 actual twin pair')
    return {'rooted_twin_quotients': twins, 'twin_cross_signs': crosses,
            'module_inheritance_cases': inherit, 'triangle_exact_cases': triangles}


def main():
    names = ['recovery/algebra_review.py', 'recovery/ALGEBRA_RESPONSE_REVIEW.md',
             'algebra/PROOF.md', 'algebra/scout.cpp', 'algebra/run.py',
             'algebra/module_checks.cpp', 'algebra/module_check.py',
             'algebra/normal.json', 'algebra/optimized.json',
             'algebra/modules.json', 'algebra/modules.optimized.json',
             'response/check.cpp', 'response/run.py', 'response/normal.json', 'response/optimized.json']
    paths = [PACKET / name for name in names]
    residual = ROOT / '05-knowledge/results/two_port_seven_core_20261103/residual'
    paths += [residual / 'full10.json', residual / 'INDEPENDENT_REVIEW.md',
              ROOT / '05-knowledge/results/session176_tournament_information_20261104/circuits/PROOF.md']
    pins = {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}
    for lane, ordinary, optimized in [('algebra', 'normal.json', 'optimized.json'),
                                       ('algebra', 'modules.json', 'modules.optimized.json'),
                                       ('response', 'normal.json', 'optimized.json')]:
        folder = PACKET / lane
        receipt = json.loads((folder / ordinary).read_text())
        need(receipt == json.loads((folder / optimized).read_text()), 'normal/optimized producer equality')
        for name, digest in receipt['source_sha256'].items():
            need(sha256((folder / name).read_bytes()).hexdigest() == digest, 'producer source pin')
    scout = json.loads((PACKET / 'algebra/normal.json').read_text())
    bank = ROOT / '05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    need(sha256(bank.read_bytes()).hexdigest() == scout['bank_sha256'], 'complete small class bank pin')
    need(scout['small_core_C_symmetric_mod4'] == {'1': 1, '2': 1, '3': 0, '4': 0, '5': 0, '6': 1, '7': 0, '8': 0}, 'small complete residue table')
    full10 = json.loads((residual / 'full10.json').read_text())
    need(full10['complete_bank'] and full10['classes'] == 9733056 and full10['endpoint_parity_balanced'] == 12670 and full10['combined_residual_addresses'] == [], 'complete order10 joint exclusion')
    for name, digest in full10['source_sha256'].items():
        need(sha256((ROOT / name).read_bytes()).hexdigest() == digest, 'order10 source pin')
    flags = ROOT / '05-knowledge/results/consecutive_flags_order10_20261028'
    for name, digest in full10['bank_sha256'].items():
        matches = list(flags.rglob(name))
        need(len(matches) == 1 and sha256(matches[0].read_bytes()).hexdigest() == digest, 'complete order10 compressed bank pin')
    controls = literal_controls()
    modules = json.loads((PACKET / 'algebra/modules.json').read_text())
    prime_subsets = 0
    for control in modules['prime_controls']:
        rows = control['rows']
        tested = 0
        for mask in range(1, 4095):
            if mask.bit_count() < 2:
                continue
            # Independent direction: every exterior vertex must beat all
            # or none of the subset, rather than comparing interior rows.
            uniform = all((rows[w] & mask) in (0, mask)
                          for w in range(12) if not mask >> w & 1)
            need(not uniform, 'saved12 prime control has no proper module')
            tested += 1
        need(tested == control['tested_subsets'] == 4082 and control['proper_nontrivial_modules'] == 0,
             'complete independent prime-control subset universe')
        prime_subsets += tested
    response = json.loads((PACKET / 'response/normal.json').read_text())
    a, b = response['controls']
    n = 12
    need(a['n'] == b['n'] == n and tuple(b['rows']) == opposite(a['rows']), 'two actual row-encoded converse controls')
    need(a['S'] == b['E'] and a['E'] == b['S'], 'full endpoint reversal')
    for i in range(n):
        for j in range(n):
            need(a['M'][i*n+j] == b['M'][j*n+i] and a['C'][i*n+j] == b['C'][j*n+i], 'full ordered pair reversal')
            for k in range(n):
                need(a['U'][(i*n+j)*n+k] == b['V'][(k*n+j)*n+i] and a['V'][(i*n+j)*n+k] == b['U'][(k*n+j)*n+i], 'full two-path reversal')
                for l in range(n):
                    need(a['D'][((i*n+j)*n+k)*n+l] == b['D'][((l*n+k)*n+j)*n+i], 'full three-path reversal')
    need(len(a['battery']) == len(b['battery']) == 6241 and a['queries'] == b['queries'] == 6250, 'complete low-weight and high-weight query counts')
    differences = sum(x != y for x, y in zip(a['battery'], b['battery']))
    need(differences == response['distinguishing_battery_entries'] == 6240, 'full battery differences')
    actual = [actual_query(tuple(row['rows']), 0, 1) for row in (a, b)]
    need(actual == [627597, 1139039], 'independent actual fourteen-vertex query')
    need(pins == {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}, 'review sources unchanged')
    print(json.dumps({'status': 'PASS INDEPENDENT ALGEBRA AND RESPONSE AUDIT',
                      'checks': CHECKS, 'literal_controls': controls,
                      'actual_query_counts': actual, 'response_differences': differences,
                      'prime_controls_subsets': prime_subsets,
                      'source_sha256': pins}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

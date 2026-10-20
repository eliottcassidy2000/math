"""Exact centre/deficit audits; run from any directory, also under python -O.

The complete n3..9 class universe is inherited, with source hashes recorded.
Low-deficit witnesses and the infinite-family samples are independently
checked using literal Hamiltonian words rather than the census canonicalizer.
No imported computation is promoted from finite to universal.
"""
from collections import Counter
import csv
import gzip
import hashlib
from itertools import permutations, product
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def decode(code, n):
    adj = [0] * n
    bit = n * (n - 1) // 2 - 1
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if code >> bit & 1 else (j, i)
            adj[u] |= 1 << v
            bit -= 1
    return tuple(adj)


def ordered_sum(parts):
    answer = []
    offset, total = 0, sum(map(len, parts))
    for part in parts:
        outside = ((1 << total) - 1) ^ ((1 << (offset + len(part))) - 1)
        answer.extend((row << offset) | outside for row in part)
        offset += len(part)
    return tuple(answer)


ONE = (0,)
C3 = (2, 4, 1)
C5 = tuple((1 << ((i + 1) % 5)) | (1 << ((i + 2) % 5)) for i in range(5))
R4 = decode(8, 4)


def quiet(word):
    return ordered_sum([ONE if x == 1 else C3 for x in word])


def cap(u):
    """x -> U -> y -> x."""
    m = len(u)
    return (((1 << (m + 1)) - 2),
            *[(row << 1) | (1 << (m + 1)) for row in u], 1)


def paths(adj):
    n = len(adj)
    full = (1 << n) - 1

    def dfs(path, used):
        if used == full:
            yield tuple(path)
            return
        available = full if not path else adj[path[-1]] & ~used
        while available:
            bit = available & -available
            available -= bit
            yield from dfs(path + [bit.bit_length() - 1], used | bit)
    yield from dfs([], 0)


def word_code(adj, path):
    code = 0
    for i, u in enumerate(path):
        for v in path[i + 1:]:
            code = (code << 1) | (adj[u] >> v & 1)
    return code


def blue(adj, path):
    n = len(path)
    return all((adj[path[i]] >> path[j] & 1) ==
               (adj[path[n - 1 - j]] >> path[n - 1 - i] & 1)
               for i in range(n) for j in range(i + 1, n))


def literal_masses(adj):
    all_words, blue_words = set(), set()
    hp = bp = 0
    for p in paths(adj):
        hp += 1
        w = word_code(adj, p)
        all_words.add(w)
        if blue(adj, p):
            bp += 1
            blue_words.add(w)
    return dict(H=hp, Hblue=bp, h=len(all_words), b=len(blue_words),
                beta=len(all_words) - len(blue_words))


def automorphisms(adj):
    groups = {}
    for v, row in enumerate(adj):
        groups.setdefault(row.bit_count(), []).append(v)
    groups = list(groups.values())
    vertices = sum(groups, [])
    answer = []
    for images in product(*(permutations(g) for g in groups)):
        p = [0] * len(adj)
        for u, v in zip(vertices, sum((list(x) for x in images), [])):
            p[u] = v
        if all((adj[i] >> j & 1) == (adj[p[i]] >> p[j] & 1)
               for i in range(len(adj)) for j in range(i + 1, len(adj))):
            answer.append(tuple(p))
    return answer


def feedback(adj):
    n = len(adj)
    dp = [0] * (1 << n)
    for mask in range(1, 1 << n):
        avail = mask
        best = n * n
        while avail:
            bit = avail & -avail
            avail -= bit
            v = bit.bit_length() - 1
            rest = mask ^ bit
            best = min(best, dp[rest] + (adj[v] & rest).bit_count())
        dp[mask] = best
    return dp[-1]


def strong(adj):
    n = len(adj)
    full = (1 << n) - 1
    for start in range(n):
        seen = 1 << start
        while True:
            nxt = seen
            for v in range(n):
                if seen >> v & 1:
                    nxt |= adj[v]
            if nxt == seen:
                break
            seen = nxt
        if seen != full:
            return False
    return True


def quiet_words(max_size, prefix=()):
    yield prefix
    for x in (1, 3):
        if sum(prefix) + x <= max_size:
            yield from quiet_words(max_size, prefix + (x,))


def finite_census():
    output, sources = {}, {}
    for n in range(3, 10):
        relative = ('05-knowledge/results/blue_quotient_20261020b/n9_classes.tsv.gz'
                    if n == 9 else
                    f'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n{n}.tsv')
        path = ROOT / relative
        sources[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
        opener = gzip.open if n == 9 else open
        with opener(path, 'rt', encoding='utf-8') as stream:
            table = list(csv.DictReader(stream, delimiter='\t'))
        hist, shist = Counter(), Counter()
        witnesses = []
        count_sc = 0
        for row in table:
            if n == 9:
                code = int(row['code'])
                sc = code == int(row['converse_code'])
                h = int(row['SC_H']) // int(row['Aut'])
                b = int(row['blue_mass'])
            else:
                code, sc = int(row['canon']), int(row['SC'])
                h, b = int(row['tc']), int(row['bm'])
            if not sc:
                continue
            count_sc += 1
            adj = decode(code, n)
            delta = h - b
            need(h % 2 == b % 2 == 1 and delta >= 0, 'mass parity')
            hist[delta] += 1
            isstrong = strong(adj)
            if isstrong:
                shist[delta] += 1
            if delta <= 20:
                check = literal_masses(adj)
                need((check['h'], check['b']) == (h, b), f'word check n{n}:{code}')
                aut = automorphisms(adj)
                need(len(aut) == int(row['Aut']), 'independent automorphism check')
                need(check['H'] == h * len(aut), 'free orbit identity')
                moved = sum(any(p[i] != i for p in aut) for i in range(n))
                f = feedback(adj)
                witnesses.append(dict(code=code, strong=isstrong, h=h, b=b,
                                      beta=delta, f=f, moved=moved,
                                      symmetry_slack=3 * f - moved))
        output[str(n)] = dict(ordinary_classes=len(table), self_converse=count_sc,
                             all_beta_le20=dict(sorted((d, c) for d, c in hist.items() if d <= 20)),
                             strong_beta_le20=dict(sorted((d, c) for d, c in shist.items() if d <= 20)),
                             witnesses_beta_le20=witnesses)
    need([output[str(n)]['strong_beta_le20'].get(2, 0) for n in range(3, 10)] ==
         [0, 1, 0, 0, 0, 0, 0], 'finite deficit-two centre inventory')
    return output, sources


def family_audit():
    output = []
    for word in quiet_words(9):
        if sum(word) < 2:
            continue
        adj = cap(quiet(word))
        r, t = word.count(1), word.count(3)
        check = literal_masses(adj)
        expected_h = 1 + 2 ** (r + 2 * t)
        pal = word == word[::-1]
        expected_b = (0 if not pal else 1 if sum(word) % 2 else
                      1 + 2 ** ((r + 2 * t) // 2))
        need(check['h'] == expected_h and check['b'] == expected_b,
             f'cap formula {word}: {check}')
        aut = automorphisms(adj)
        need(len(aut) == 3 ** t, 'cap automorphism identity')
        moved = sum(any(p[i] != i for p in aut) for i in range(len(adj)))
        f = feedback(adj)
        need(f == t + 1 and moved == 3 * t, 'slack-three independent check')
        need(strong(adj), 'cap must be strong')
        output.append(dict(word=word, n=len(adj), **check,
                           Aut=len(aut), f=f, moved=moved, symmetry_slack=3))
    return output


def polynomial_audit():
    # p_l(z): count W by |W| at fixed effective weight #1+2#C3=l.
    p = [Counter({0: 1}), Counter({1: 1})]
    for ell in range(2, 13):
        c = Counter({k + 1: v for k, v in p[-1].items()})
        c.update({k + 3: v for k, v in p[-2].items()})
        p.append(c)
    output = []
    for ell, c in enumerate(p):
        direct = Counter({ell + t: comb(ell - t, t) for t in range(ell // 2 + 1)})
        need(c == direct, 'Fibonacci polynomial coefficient identity')
        output.append(dict(ell=ell, coefficients=dict(sorted(c.items())),
                           even_center_beta=2 ** ell * (2 ** ell - 1),
                           singleton_center_beta=2 ** (2 * ell + 1),
                           triangle_center_beta=2 ** (2 * ell + 2)))
    return output


def stack_controls():
    examples = []
    for name, parts, target in [
        ('C5,C5', [C5, C5], (9, 3, 6)),
        ('C5,1,C5', [C5, ONE, C5], (9, 3, 6)),
        ('C5,C3,C5', [C5, C3, C5], (9, 3, 6)),
        ('C5,R4,C5', [C5, R4, C5], (45, 9, 36)),
        ('R4,R4', [R4, R4], (25, 5, 20)),
        ('R4,R4,R4', [R4] * 3, (125, 15, 110)),
        ('R4,R4,R4,R4', [R4] * 4, (625, 25, 600)),
    ]:
        data = literal_masses(ordered_sum(parts))
        need(tuple(data[k] for k in ('h', 'b', 'beta')) == target, 'stack control')
        examples.append(dict(name=name, n=sum(map(len, parts)), **data))
    # Arithmetic classifications of the first outer-mass corrections.
    choices = {}
    for beta in range(2, 43, 2):
        candidates = []
        for q in range(3, beta + 1, 2):
            if beta % q or q * (q - 1) > beta:
                continue
            for c in range(1, beta + 1, 2):
                d = q * c - beta // q
                if 1 <= d <= c:
                    candidates.append((q, c, d))
        choices[beta] = candidates
    need(choices[6] == [(3, 1, 1)] and choices[12] == [] and
         choices[18] == [(3, 3, 3)] and choices[20] == [(5, 1, 1)] and
         choices[24] == [(3, 3, 1)] and choices[36] == [(3, 5, 3)],
         'first deficit correction table')
    return examples, choices


def main():
    census, sources = finite_census()
    family = family_audit()
    polynomials = polynomial_audit()
    stacks, choices = stack_controls()
    receipt = dict(status='PASS', universe='complete inherited classes n3..9; independent low-deficit checks; all quiet words of size at most9',
                   sources_sha256=sources, census=census, cap_family=family,
                   fibonacci_polynomials=polynomials, stack_controls=stacks,
                   possible_nonrigid_outer_corrections=choices)
    (HERE / 'centers.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print('PASS complete inherited self-converse census n3..9; all beta<=20 witnesses independently checked')
    print('Strong beta<=10:', {n: {d: c for d, c in v['strong_beta_le20'].items() if d <= 10} for n, v in census.items()})
    print('PASS', len(family), 'quiet-word caps of inner order<=9: literal h,b, Aut, feedback, moved support')
    print('PASS Fibonacci polynomial identities ell0..12 and seven literal SCC-stack controls')
    print('FINITE-EXACT D2(z)=z^4 through order9; PROVED within the entire symmetry-slack3 strong-centre family')
    print('All-order D2(z)=z^4 beyond that family remains OPEN.')


if __name__ == '__main__':
    main()

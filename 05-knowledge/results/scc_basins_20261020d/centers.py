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
from itertools import combinations, permutations, product
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


def reverse_path(adj, path):
    answer = list(adj)
    for u, v in zip(path, path[1:]):
        need(answer[u] >> v & 1, 'path reversal input arc')
        answer[u] ^= 1 << v
        answer[v] ^= 1 << u
    return tuple(answer)


def triangles(adj):
    total = 0
    for vertices in combinations(range(len(adj)), 3):
        if all(sum(adj[u] >> v & 1 for v in vertices) == 1 for u in vertices):
            total += 1
    return total


def neighbour_controls():
    output = []
    for word in quiet_words(4):
        left = quiet(word)
        tournament = ordered_sum([left, R4, quiet(word[::-1])])
        p = next(p for p in paths(tournament) if not blue(tournament, p))
        neighbour = reverse_path(tournament, p)
        masses = literal_masses(neighbour)
        need(masses['b'] == 0, 'mixed-leaf neighbour must be nonselfconverse')
        aut = automorphisms(neighbour)
        centre = tuple(range(len(left), len(left) + 4))
        scores = [row.bit_count() for row in neighbour]
        if word:
            shell = tuple(i for i, score in enumerate(scores)
                          if score in (len(left) + 1, len(left) + 2))
            need(shell == centre, 'intrinsic central score shell')
            need(strong(neighbour) and len(aut) == 1, 'strong asymmetric neighbour')
            expected_triangles = len(tournament) + 2 * word.count(3) - 2 * (word[0] == 3)
            need(triangles(neighbour) == expected_triangles, 'neighbour triangle current')
            for score in set(scores):
                vertices = [i for i, s in enumerate(scores) if s == score]
                sub = tuple(sum((neighbour[u] >> v & 1) << j for j, v in enumerate(vertices))
                            for u in vertices)
                need(triangles(sub) == 0, 'score class is transitive')
        else:
            need(not strong(neighbour) and len(aut) == 3 and masses['h'] == 1,
                 'empty-padding exception')
        output.append(dict(word=word, n=len(tournament), original_rows=tournament,
                           black_path=p, neighbour_rows=neighbour, scores=scores,
                           strong=strong(neighbour), Aut=len(aut),
                           directed_triangles=triangles(neighbour), **masses))
    return output


def recover_left_word(neighbour):
    """Decode the canonical orientation whose central shell is C3 -> 1."""
    n = len(neighbour)
    if n == 4:
        return ()
    r = (n - 4) // 2
    scores = [row.bit_count() for row in neighbour]
    left = [i for i, s in enumerate(scores) if s >= r + 3]
    need(len(left) == r, 'left shell size')
    top_score = max(scores[i] for i in left)
    top = [i for i in left if scores[i] == top_score]
    if len(top) == 1:
        prefix, removed = (1,), set(top)
    else:
        need(len(top) == 2, 'top shell size')
        following = [i for i in left if scores[i] == top_score - 1]
        need(len(following) <= 1, 'next shell size')
        vertices = top + following
        is_triangle = (len(vertices) == 3 and
                       all(sum(neighbour[u] >> v & 1 for v in vertices) == 1
                           for u in vertices))
        if following and not is_triangle:
            prefix, removed = (3,), set(vertices)
        else:
            prefix, removed = (1, 1), set(top)
    rest = Counter(scores[i] for i in left if i not in removed)
    need(all(size in (1, 3) for size in rest.values()), 'remaining atom shells')
    return prefix + tuple(rest[s] for s in sorted(rest, reverse=True))


def large_neighbour_inverse():
    addresses, counts = {}, Counter()
    digest = hashlib.sha256()
    for word in quiet_words(18):
        r = sum(word)
        tournament = ordered_sum([quiet(word), R4, quiet(word[::-1])])
        p = tuple(range(r)) + (r, r + 3, r + 2, r + 1) + tuple(range(r + 4, 2 * r + 4))
        neighbour = reverse_path(tournament, p)
        need(recover_left_word(neighbour) == word, 'intrinsic inverse word decoder')
        if word:
            scores = [row.bit_count() for row in neighbour]
            groups = [[i for i, s in enumerate(scores) if s == score] for score in sorted(set(scores))]
            ordering = []
            for group in groups:
                degrees = {i: sum(neighbour[i] >> j & 1 for j in group) for i in group}
                need(sorted(degrees.values()) == list(range(len(group))), 'score group transitive')
                ordering.extend(sorted(group, key=degrees.get))
            address = (len(neighbour), word_code(neighbour, ordering))
            need(address not in addresses, 'distinct mixed leaves share black neighbour')
            addresses[address] = word
            # Literal closed reversed path certifies strongness.
            cycle = p[::-1]
            need(all(neighbour[u] >> v & 1 for u, v in zip(cycle, cycle[1:] + cycle[:1])),
                 'neighbour Hamilton cycle')
        else:
            address = (4, 0)
        digest.update(repr((word, address)).encode('ascii'))
        digest.update(b'\n')
        counts[r] += 1
    return dict(max_half_order=18, max_total_order=40,
                checked_words=sum(counts.values()), words_by_half_order=dict(sorted(counts.items())),
                actual_neighbour_collisions=0, exact_inverse_failures=0,
                word_address_sha256=digest.hexdigest())


def main():
    census, sources = finite_census()
    family = family_audit()
    polynomials = polynomial_audit()
    stacks, choices = stack_controls()
    neighbours = neighbour_controls()
    inverse = large_neighbour_inverse()
    receipt = dict(status='PASS', universe='complete inherited classes n3..9; independent low-deficit checks; all quiet words of size at most9',
                   sources_sha256=sources, census=census, cap_family=family,
                   fibonacci_polynomials=polynomials, stack_controls=stacks,
                   exact_black_neighbours=neighbours,
                   black_neighbour_inverse=inverse,
                   possible_nonrigid_outer_corrections=choices)
    (HERE / 'centers.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print('PASS complete inherited self-converse census n3..9; all beta<=20 witnesses independently checked')
    print('Strong beta<=10:', {n: {d: c for d, c in v['strong_beta_le20'].items() if d <= 10} for n, v in census.items()})
    print('PASS', len(family), 'quiet-word caps of inner order<=9: literal h,b, Aut, feedback, moved support')
    print('PASS Fibonacci polynomial identities ell0..12 and seven literal SCC-stack controls')
    print('PASS', len(neighbours), 'explicit mixed-leaf neighbours: score shell, nonselfconverseness, Aut, triangle current')
    print('PASS', inverse['checked_words'], 'mixed-leaf words of half-order<=18: intrinsic inverse and unique neighbour addresses')
    print('FINITE-EXACT D2(z)=z^4 through order9; PROVED within the entire symmetry-slack3 strong-centre family')
    print('All-order D2(z)=z^4 is imported from deficit_two_centers_20261020; its residual census is not replayed by this script.')


if __name__ == '__main__':
    main()

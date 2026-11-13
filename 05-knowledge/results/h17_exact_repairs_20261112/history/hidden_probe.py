"""Exact finite physical-source probe for the one-birth response fibre."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent


def require(test, message):
    if not test:
        raise ValueError(message)


def rows_from_mask(n, mask):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if mask >> k & 1 else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def literal_profile(rows, root=0):
    n = len(rows)

    def path(word):
        return all(rows[u] >> v & 1 for u, v in zip(word, word[1:]))

    h = s = f = 0
    for word in permutations(range(n)):
        good = path(word)
        h += good
        s += good and word[0] == root
        cut = word.index(root) + 1
        f += path(word[:cut]) and path(word[cut:])
    q = sum(path(word) for word in permutations(v for v in range(n) if v != root))
    return h, s, q, f


def subset_profile(rows, root=0):
    n = len(rows)
    full = (1 << n) - 1

    @lru_cache(None)
    def end(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(end(rest, u) for u in range(n)
                   if rest >> u & 1 and rows[u] >> v & 1)

    @lru_cache(None)
    def start(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(start(rest, u) for u in range(n)
                   if rest >> u & 1 and rows[v] >> u & 1)

    @lru_cache(None)
    def h(mask):
        return sum(end(mask, v) for v in range(n) if mask >> v & 1) if mask else 1

    return (h(full), start(full, root), h(full ^ (1 << root)),
            sum(end(mask, root) * h(full ^ mask)
                for mask in range(1 << n) if mask >> root & 1))


def strong(rows):
    n = len(rows)
    for v in range(n):
        reached = 1 << v
        while True:
            expanded = reached
            for u in range(n):
                if reached >> u & 1:
                    expanded |= rows[u]
            if expanded == reached:
                break
            reached = expanded
        if reached != (1 << n) - 1:
            return False
    return True


def birth(rows, root, letter):
    n = len(rows)
    out = list(rows)
    out[root] |= 1 << n
    out.append(((1 << n) - 1) ^ (1 << root))
    return tuple(out), n if letter == 'N' else root


def step(z, letter):
    h, s, q, f = z
    if letter == 'R':
        return h + f - s, q, q, 2 * f - q
    return h + f - s, h - s, h, h + f


def cyclic_blowup(a, b, c):
    groups = []
    pos = 0
    for size in (a, b, c):
        groups.append(list(range(pos, pos + size)))
        pos += size
    rows = [0] * (pos + 2)
    for i, group in enumerate(groups):
        for j, v in enumerate(group):
            rows[v] = sum(1 << w for w in group[j + 1:] + groups[(i + 1) % 3])
            rows[v] |= 1 << pos
    rows[pos] = 1 << (pos + 1)
    rows[pos + 1] = (1 << pos) - 1
    return tuple(rows), pos


def main():
    universes = []
    transitions = {}
    for n in range(1, 6):
        pairs = tuple(combinations(range(n), 2))
        changes = [sum(1 << i for i in choice)
                   for size in (1, 2) for choice in combinations(range(len(pairs)), size)]
        profiles = [literal_profile(rows_from_mask(n, mask))
                    for mask in range(1 << len(pairs))]
        tested = Counter()
        retained = Counter()
        for mask, z in enumerate(profiles):
            h, s, q, f = z
            hr = step(z, 'R')[0]
            hrr = step(step(z, 'R'), 'R')[0]
            hnr = step(step(z, 'N'), 'R')[0]
            require((hnr - h) == 2 * f and h + f - hr == s
                    and 2 * q == 2 * f + hr - hrr, 'four-observation inversion')
            for change in changes:
                other = mask ^ change
                if other <= mask:
                    continue
                radius = change.bit_count()
                tested[radius] += 1
                w = profiles[other]
                if h != w[0] or f - s != w[3] - w[1] or z == w:
                    continue
                retained[radius] += 1
                key = (n, radius, tuple(sorted((z, w))))
                if key not in transitions:
                    left = rows_from_mask(n, mask)
                    right = rows_from_mask(n, other)
                    transitions[key] = {'order': n, 'changed_arcs': radius,
                                       'left_rows': left, 'right_rows': right,
                                       'left_state': z, 'right_state': w,
                                       'both_strong': strong(left) and strong(right),
                                       'labelled_pair_count': 0}
                transitions[key]['labelled_pair_count'] += 1
        universes.append({'order': n, 'root': 0, 'labelled_sources': len(profiles),
                          'unordered_pairs_tested': dict(tested),
                          'hidden_pairs': dict(retained)})

    old = (16, 9, 3, 5, 14)
    new = (18, 24, 3, 5, 12)
    z, w = subset_profile(old), subset_profile(new)
    require(z == (15, 3, 3, 21) and w == (15, 3, 5, 21), 'selected source')
    require(strong(old) and strong(new), 'selected source strongness')
    require(z == literal_profile(old) and w == literal_profile(new), 'independent source recount')
    first_children = []
    for rows in (old, new):
        child, root = birth(rows, 0, 'R')
        first_children.append({'rows': child, 'root': root,
                               'profile': subset_profile(child, root)})
    require([x['profile'] for x in first_children] == [(33, 3, 3, 39), (33, 5, 5, 37)],
            'first actual target')
    word_controls = 0
    for length in range(9):
        for letters in product('RN', repeat=length):
            word = ''.join(letters)
            zz, ww = z, w
            for letter in word:
                zz, ww = step(zz, letter), step(ww, letter)
            equal = not word or word == 'R' or word.startswith('N')
            require((zz[0] == ww[0]) == equal, ('equality language', word))
            if word.startswith('R') and len(word) >= 2:
                require(ww[0] < zz[0], ('strict sign', word))
            word_controls += 1
    returns = []
    for k in (1, 2, 3):
        results = []
        for rows in (old, new):
            root = 0
            for letter in 'R' * k + 'NR':
                rows, root = birth(rows, root, letter)
            results.append(subset_profile(rows, root))
        require(gcd(*results[0]) == 3 and gcd(*results[1]) == 1, 'return primitivity')
        require(results[1] == (48 * 2 ** k + 11, 16 * 2 ** k + 1,
                               16 * 2 ** k + 1, 48 * 2 ** k + 11), 'inherited return formula')
        returns.append({'k': k, 'old_profile': results[0], 'repaired_profile': results[1]})

    h9_rows = ((22, 8, 2, 5, 14), (22, 8, 2, 21, 6))
    h9_controls = []
    for index, rows in enumerate(h9_rows):
        def path(word):
            return all(rows[u] >> v & 1 for u, v in zip(word, word[1:]))

        state = subset_profile(rows)
        expected = ((9, 3, 3, 9), (9, 4, 5, 10))[index]
        require(state == expected == literal_profile(rows), 'H9 two-source profile')
        require(strong(rows), 'H9 strong source')
        triangles = [word for word in combinations(range(5), 3)
                     if all(sum(rows[u] >> v & 1 for v in word) == 1 for u in word)]
        overlap_edges = [(i, j) for i, j in combinations(range(len(triangles)), 2)
                         if len(set(triangles[i]) & set(triangles[j])) == 2]
        pentagons = [(0,) + word for word in permutations(range(1, 5))
                     if path((0,) + word + (0,))]
        require(len(triangles) == 3 and len(pentagons) == 1, 'H9 odd-cycle inventory')
        require(len(overlap_edges) == 2 + index, 'H9 P3/K3 overlap distinction')
        hpaths = [word for word in permutations(range(5)) if path(word)]
        deleted_paths = [word for word in permutations(range(1, 5)) if path(word)]
        predecessor_paths = [word for word in deleted_paths if word[0] == 3]
        require(len(predecessor_paths) == index, 'H9 unit readiness')
        require(not index or predecessor_paths == [(3, 4, 2, 1)], 'readiness carrier')
        child, root = birth(rows, 0, 'N')
        n_state = subset_profile(child, root)
        require(n_state == ((15, 6, 9, 18), (15, 5, 9, 19))[index], 'H9 actual N target')
        exit_rows, root = birth(child, root, 'R')
        nr_state = subset_profile(exit_rows, root)
        require(nr_state == ((27, 9, 9, 27), (29, 9, 9, 29))[index], 'H9 actual NR exit')
        require(gcd(*n_state) == (3, 1)[index], 'same-target primitive repair')
        require(state[3] == state[0] + state[2] - state[1], 'indegree-one identity')
        h9_controls.append({'rows': rows, 'profile': state, 'Hamilton_paths': hpaths,
                            'deleted_paths': deleted_paths,
                            'predecessor_start_paths': predecessor_paths,
                            'triangles': triangles, 'overlap_edges': overlap_edges,
                            'directed_pentagons_up_to_rotation': pentagons,
                            'N_rows': child, 'N_root': 5, 'N_profile': n_state,
                            'NR_rows': exit_rows, 'NR_root': root, 'NR_profile': nr_state})
    require(h9_rows[0][0] == h9_rows[1][0] == 22, 'whole root star retained')
    changed = [(i, j) for i, j in combinations(range(5), 2)
               if (h9_rows[0][i] ^ h9_rows[1][i]) >> j & 1]
    require(changed == [(3, 4)] and h9_rows[0][4] >> 3 & 1,
            'only old arc 4 to 3 reverses')

    blowups = []
    for abc in product((1, 2), repeat=3):
        rows, root = cyclic_blowup(*abc)
        z0 = subset_profile(rows, root)
        variations = []
        for v in range(root):
            altered = list(rows)
            for i, j in ((v, root), (root + 1, v)):
                altered[i] ^= 1 << j
                altered[j] ^= 1 << i
            zz = subset_profile(tuple(altered), root)
            good = z0[0] == zz[0] and z0[3] - z0[1] == zz[3] - zz[1]
            require(good == (abc == (1, 1, 1)), 'bounded blowup boundary')
            variations.append({'selected_vertex': v, 'delta_h': zz[0] - z0[0],
                               'same_h_and_first_target': good})
        blowups.append({'block_sizes': abc, 'source_profile': z0, 'variations': variations})

    receipt = {'status': 'FINITE-EXACT; all-word conclusions proved separately',
               'universe': universes, 'profile_transition_classes': list(transitions.values()),
               'selected_witness': {'root': 0, 'old_rows': old, 'new_rows': new,
                                    'old_profile': z, 'new_profile': w,
                                    'first_R_children': first_children},
               'return_controls': returns, 'word_controls': word_controls,
               'H9_root_star_ready_repair': h9_controls,
               'bounded_cyclic_blowup_probe': blowups,
               'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output = json.dumps(receipt, indent=2) + '\n'
    destination = HERE / 'hidden_receipt.json'
    if '--write' in sys.argv:
        destination.write_text(output)
    else:
        require(destination.read_text() == output, 'receipt mismatch')
    print(json.dumps({'passed': True, 'labelled_sources': sum(u['labelled_sources'] for u in universes),
                      'transition_classes': len(transitions), 'word_controls': word_controls,
                      'receipt_sha256': hashlib.sha256(output.encode()).hexdigest()}))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact finite controls for the proved saturated near-twin lemma; stdlib only."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import gzip
import hashlib
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def closure(pred):
    pred = list(pred)
    for k in range(len(pred)):
        for j in range(len(pred)):
            if pred[j] >> k & 1:
                pred[j] |= pred[k]
    return tuple(pred)


def natural_posets(n):
    edges = tuple(combinations(range(n), 2))
    seen = set()
    for bits in range(1 << len(edges)):
        pred = [0] * n
        for i, (a, b) in enumerate(edges):
            if bits >> i & 1:
                pred[b] |= 1 << a
        pred = closure(pred)
        if pred not in seen:
            seen.add(pred)
            yield pred


def successors(pred):
    return tuple(sum(1 << j for j in range(len(pred)) if pred[j] >> i & 1)
                 for i in range(len(pred)))


def extensions(pred):
    n = len(pred)
    def visit(word, used):
        if len(word) == n:
            yield tuple(word)
            return
        for v in range(n):
            if not (used >> v & 1) and not (pred[v] & ~used):
                yield from visit(word + [v], used | (1 << v))
    return list(visit([], 0))


def relabel(word, mapping):
    return tuple(mapping.get(v, v) for v in word)


def check_extension_maps():
    counts = Counter()
    size_counts = []
    max_probability = Fraction(0)
    sharp_example = None
    for n in range(1, 7):
        count = 0
        for pred in natural_posets(n):
            count += 1
            succ = successors(pred)
            near = []
            twins = []
            for x, z in product(range(n), repeat=2):
                if x == z or pred[x] != pred[z]:
                    continue
                if succ[x] == succ[z] and x < z:
                    twins.append((x, z))
                extra = succ[x] & ~succ[z]
                if not (succ[z] & ~succ[x]) and extra.bit_count() == 1:
                    near.append((x, z, extra.bit_length() - 1))
            if not near and not twins:
                continue
            words = extensions(pred)
            wordset = set(words)
            ranks = {w: {v: i for i, v in enumerate(w)} for w in words}
            for x, z in twins:
                left = {w for w in words if ranks[w][x] < ranks[w][z]}
                right = wordset - left
                images = {relabel(w, {x:z, z:x}) for w in left}
                check(images == right, 'Twin swap is not a bijection')
                counts['twin_pairs'] += 1
                counts['twin_extensions_checked'] += len(words)
            for x, z, a in near:
                classes = [set(), set(), set()]
                for w in words:
                    r = ranks[w]
                    if r[z] < r[x] < r[a]:
                        classes[0].add(w)
                    elif r[x] < r[z] < r[a]:
                        classes[1].add(w)
                    elif r[x] < r[a] < r[z]:
                        classes[2].add(w)
                    else:
                        raise RuntimeError('Unexpected relative order')
                E0, E1, E2 = classes
                swapped = {relabel(w, {x:z, z:x}) for w in E0}
                rotated = {relabel(w, {x:z, a:x, z:a}) for w in E2}
                check(swapped == E1, 'Near-twin swap is not a bijection')
                check(len(rotated) == len(E2) and rotated <= E0,
                      'Near-twin rotation is not an injection into E0')
                p = Fraction(len(E1) + len(E2), len(words))
                check(Fraction(1,2) <= p <= Fraction(2,3), 'Probability bound failed')
                if p > max_probability:
                    max_probability = p
                    sharp_example = {'pred': pred, 'x':x, 'z':z, 'a':a, 'p':str(p)}
                counts['near_twin_pairs'] += 1
                counts['near_twin_extensions_checked'] += len(words)
                counts['rotation_images_checked'] += len(E2)
        size_counts.append(count)
    check(size_counts == [1,2,7,40,357,4824], 'Natural-poset census mismatch')
    check(max_probability == Fraction(2,3), 'Sharp upper endpoint not observed')
    return {'natural_poset_counts_n_1_to_6':size_counts, **dict(counts),
            'maximum_observed_probability':str(max_probability), 'sharp_example':sharp_example}


def c4_cores():
    c = 4
    full = (1 << c) - 1
    for zp in natural_posets(c):
        ideals = [s for s in range(1 << c)
                  if all(not (s >> i & 1) or not (zp[i] & ~s) for i in range(c))]
        upsets = [full ^ s for s in ideals]
        for lower, d1, d2 in product(upsets, ideals, ideals):
            for chain, ba1, ba2 in product((0,1), repeat=3):
                if chain and (d1 & ~d2 or ba1 > ba2):
                    continue
                if lower & d1 and not ba1 or lower & d2 and not ba2:
                    continue
                p = [0,4,0] + [(zp[i] << 3) | (4 if lower >> i & 1 else 0)
                               for i in range(c)]
                p += [1 | (d1 << 3) | (4 if ba1 else 0),
                      1 | (d2 << 3) | (4 if ba2 else 0) | (128 if chain else 0)]
                yield tuple(p)


def check_census():
    seen = set()
    counts = Counter()
    coverage = []
    catalogue = hashlib.sha256()
    for p in c4_cores():
        check(p not in seen, 'Duplicate morphology')
        check(p == closure(p), 'Nontransitive morphology')
        seen.add(p)
        succ = successors(p)
        degrees = [8 - (p[i] | succ[i]).bit_count() for i in range(9)]
        counts['all_canonical_cores'] += 1
        if max(degrees) > 7:
            continue
        core_index = counts['degree_valid_cores']
        counts['degree_valid_cores'] += 1
        catalogue.update((json.dumps(p, separators=(',',':')) + '\n').encode())
        low = any((p[z] | succ[z]) == 4 for z in range(3,7))
        upper = any((p[z] | succ[z]) in (128,256) for z in range(3,7))
        counts['twin_y_certificate'] += low
        counts['one_extra_upper_certificate'] += upper
        counts['both_structural_certificates'] += low and upper
        counts['either_structural_certificate'] += low or upper
        counts['remaining_after_both'] += not (low or upper)
        if low or upper:
            choices = [(z, 'twin_y' if (p[z] | succ[z]) == 4 else 'one_extra_upper')
                       for z in range(3,7) if (p[z] | succ[z]) in (4,128,256)]
            z, kind = choices[0]
            coverage.append({'core_index':core_index, 'predecessor_masks':p,
                             'kind':kind, 'common_vertex':z})
        for z in range(3,7):
            nbr = p[z] | succ[z]
            if nbr == 4:
                check(p[z] == p[1] and succ[z] == succ[1], 'Core twin sets differ')
            elif nbr in (128,256):
                check(p[z] == p[0] and succ[0] == (succ[z] | (384 ^ nbr)),
                      'Core near-twin sets differ')
    expected = {'all_canonical_cores':96428, 'degree_valid_cores':79121,
                'twin_y_certificate':13885, 'one_extra_upper_certificate':21893,
                'both_structural_certificates':4210,
                'either_structural_certificate':31568, 'remaining_after_both':47553}
    check(dict(counts) == expected, 'Structural pruning census mismatch')
    return dict(counts), coverage, catalogue.hexdigest()


def main():
    census, records, catalogue_hash = check_census()
    coverage = {'scope':'Exact ambient balanced-pair reductions, independent of majority order.',
                'degree_valid_catalogue_sha256':catalogue_hash,
                'covered_cores':len(records), 'certificates':records}
    coverage_data = gzip.compress((json.dumps(coverage, sort_keys=True, separators=(',',':')) + '\n').encode(), mtime=0)
    coverage_path = Path(__file__).with_name('c4_saturated_near_twin_coverage_20261019.json.gz')
    coverage_path.write_bytes(coverage_data)
    result = {'status':'PASS',
              'scope':'Finite controls for separately proved structural lemmas; no complete c=4 theorem.',
              'extension_map_controls':check_extension_maps(),
              'c4_structural_census':census,
              'coverage_file':coverage_path.name,
              'coverage_sha256':hashlib.sha256(coverage_data).hexdigest(),
              'degree_valid_catalogue_sha256':catalogue_hash,
              'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    target = Path(__file__).with_name('saturated_near_twin_lemma_controls_20261019.json')
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

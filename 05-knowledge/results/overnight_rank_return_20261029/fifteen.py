"""Complete fifteen-Lagrangian carrier and quadratic-incidence census."""
from functools import lru_cache
from itertools import combinations, product
from fractions import Fraction
from pathlib import Path
import argparse
import json


def need(ok, why):
    if not ok:
        raise ValueError(why)


def omega(u, v):
    return (((u & 3) & (v >> 2)).bit_count() +
            ((u >> 2) & (v & 3)).bit_count()) % 2


def Q(v):
    return ((v & 3) & (v >> 2)).bit_count() % 2


def plucker(line):
    u, v = line[:2]
    def p(i, j):
        return ((u >> i & 1) * (v >> j & 1)) ^ ((u >> j & 1) * (v >> i & 1))
    a, b, c, d, t = p(0, 1), p(0, 3), p(1, 2), p(2, 3), p(0, 2)
    need(p(1, 3) == t, 'isotropic Plucker constraint')
    need(t == (a*d ^ b*c), 'quadratic free bit')
    return a | b << 1 | c << 2 | d << 3, t


def rank2(rows):
    basis = {}
    for row in rows:
        while row:
            j = row.bit_length() - 1
            if j in basis:
                row ^= basis[j]
            else:
                basis[j] = row
                break
    return len(basis)


def rankq(rows):
    matrix = [[Fraction((row >> i) & 1) for i in range(15)] for row in rows]
    r = 0
    for j in range(15):
        pivot = next((i for i in range(r, len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            continue
        matrix[r], matrix[pivot] = matrix[pivot], matrix[r]
        v = matrix[r][j]
        matrix[r] = [x/v for x in matrix[r]]
        for i in range(len(matrix)):
            if i != r:
                v = matrix[i][j]
                matrix[i] = [a-v*b for a, b in zip(matrix[i], matrix[r])]
        r += 1
    return r


def apply(columns, vector):
    answer = 0
    for i, col in enumerate(columns):
        if vector >> i & 1:
            answer ^= col
    return answer


def census():
    lagrangians = sorted(set(tuple(sorted((u, v, u ^ v)))
                            for u in range(1, 16) for v in range(u+1, 16)
                            if omega(u, v) == 0))
    codes = [plucker(line)[0] for line in lagrangians]
    need(set(codes) == set(range(1, 16)), 'four-bit bijection')
    forms = [v for v in range(16) if Q(v) == 0]
    rows = [sum(int(all((Q(x) ^ omega(v, x)) == 0 for x in line)) << i
                for i, line in enumerate(lagrangians)) for v in forms]
    need(len(forms) == 10 and all(row.bit_count() == 6 for row in rows), 'even forms')
    need(all(sum(row >> i & 1 for row in rows) == 4 for i in range(15)), 'incidence degrees')
    weights = []
    for row in rows:
        coeff = [c for c in range(32) if all(
            (((c & 15) & code).bit_count() ^ ((c >> 4) * plucker(line)[1])) % 2
            == (row >> i & 1) for i, (code, line) in enumerate(zip(codes, lagrangians)))]
        need(len(coeff) == 1 and coeff[0] & 16, 'unique five-coordinate response')
        weights.append(coeff[0])
    separators = []
    for k in range(1, 6):
        found = [choice for choice in combinations(range(10), k)
                 if len({tuple(rows[j] >> i & 1 for j in choice)
                         for i in range(15)}) == 15]
        if found:
            separators = found
            need(k == 5, 'minimum fixed query count')
            break
    @lru_cache(None)
    def depth(state):
        if state.bit_count() <= 1:
            return 0
        return 1 + min(max(depth(state & row), depth(state & ~row))
                       for row in rows if state & row and state & ~row)
    need(depth((1 << 15)-1) == 5, 'adaptive optimum')
    need(rank2(rows) == 5 and rankq(rows) == 10, 'field-sensitive response ranks')
    need(all((a & b).bit_count() == (6 if i == j else 2)
             for i, a in enumerate(rows) for j, b in enumerate(rows)), 'rational Gram matrix')
    # Exhaust every symplectic ordered basis, and verify the transported
    # four-coordinate action without retaining the quadratic fifth bit.
    group = []
    standard = [1, 2, 4, 8]
    for cols in product(range(1, 16), repeat=4):
        if not all(omega(cols[i], cols[j]) == omega(standard[i], standard[j])
                   for i in range(4) for j in range(i)):
            continue
        image = {}
        for line, code in zip(lagrangians, codes):
            new_line = tuple(sorted(apply(cols, x) for x in line))
            image[code] = plucker(new_line)[0]
        induced = tuple(image[x] for x in standard)
        need(all(apply(induced, code) == target for code, target in image.items()),
             'linear four-bit transport')
        group.append((cols, induced))
    need(len(group) == 720 and len({b for a, b in group}) == 720, 'faithful transported group')
    @lru_cache(None)
    def fixed_cost(state):
        indices = [i for i in range(15) if state >> i & 1]
        for k in range(6):
            if any(len({tuple(rows[j] >> i & 1 for j in choice) for i in indices})
                   == len(indices) for choice in combinations(range(10), k)):
                return k
        raise ValueError('missing separator')
    phase_histogram = {}
    for _, induced in group:
        fixed = sum((apply(induced, code) == code) << i for i, code in enumerate(codes))
        key = (fixed.bit_count(), fixed_cost(fixed), depth(fixed))
        phase_histogram[key] = phase_histogram.get(key, 0) + 1
    need(phase_histogram == {(15, 5, 5): 1, (7, 3, 3): 15, (3, 2, 2): 100,
                             (1, 0, 0): 300, (0, 0, 0): 304}, 'Frobenius query census')
    def quadratic(x):
        return ((x & 1)*(x >> 3 & 1)) ^ ((x >> 1 & 1)*(x >> 2 & 1))
    def pairing(x, y):
        return quadratic(x ^ y) ^ quadratic(x) ^ quadratic(y)
    for x in range(16):
        for y in range(16):
            for z in range(16):
                need(pairing(x ^ y, z) == pairing(x, z) ^ pairing(y, z), 'LCA bilinearity')
                need(pairing(x, y) ^ pairing(x ^ y, z) ==
                     pairing(y, z) ^ pairing(x, y ^ z), 'LCA tree reassociation')
    for _, induced in group:
        eta = sum((quadratic(apply(induced, e)) ^ quadratic(e)) << i
                  for i, e in enumerate(standard))
        for x in range(16):
            need(quadratic(apply(induced, x)) ^ quadratic(x) == (eta & x).bit_count() % 2,
                 'linear transport of the fifth parity coordinate')
    need(quadratic(1) == quadratic(8) == 0 and quadratic(1 ^ 8) == 1,
         'same aggregate vector but different joint response')
    # Save all matrices: a Frobenius matrix can be looked up exactly.
    return dict(lagrangians=[dict(nonzero_vectors=l, four_bit_code=c,
                                  fifth_bit=plucker(l)[1]) for l, c in zip(lagrangians, codes)],
                form_labels=forms, incidence_rows=rows, linear_response_coefficients=weights,
                fixed_query_minimum=5, adaptive_query_minimum=5,
                number_of_five_query_separators=len(separators),
                example_separator=[forms[j] for j in separators[0]],
                binary_response_rank=5, rational_response_rank=10,
                lca_merge_triples=4096, realization_defect_hostile=[1, 8, 9],
                frobenius_query_census=[dict(fixed_kernels=k[0], fixed_queries=k[1],
                                             adaptive_queries=k[2], matrices=v)
                                        for k, v in sorted(phase_histogram.items())],
                symplectic_transport=[dict(input_columns=a, output_columns=b) for a, b in group])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    result = census()
    path = Path(__file__).with_suffix('.json')
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.write:
        path.write_text(encoded, encoding='utf-8')
    else:
        need(json.loads(path.read_text(encoding='utf-8')) == json.loads(encoded), 'stored census')
    print(json.dumps({k: v for k, v in result.items() if k.endswith('minimum') or
                     k.endswith('rank') or k.endswith('separators') or k == 'example_separator'}))


if __name__ == '__main__':
    main()

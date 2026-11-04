"""Small independent controls for the arithmetic note; no bank census."""
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
SOURCES = {
    'reader': RESULTS / 'zero_ready_component_returns_20261103/inverse/reader.py',
    'sections': HERE.parent / 'selector/inverse.py',
    's159': RESULTS / 'session159_constructions_20261103/automatic/NOTE.md',
}


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


reader = module('inherited_reader', SOURCES['reader'])
sections = module('section_inverse', SOURCES['sections'])


def excess(r):
    return 0 if r == 0 else reader.word_h(reader.decode_r(r)) - 1


def one_apex(word):
    n = len(word)
    rows = [sum(1 << j for j in range(i + 1, n)) for i in range(n)] + [0]
    for i, bit in enumerate(word):
        if bit == '1':
            rows[n] |= 1 << i
        else:
            rows[i] |= 1 << n
    return tuple(rows)


def ordinal(left, right):
    n, m = len(left), len(right)
    return tuple(x | (((1 << m) - 1) << n) for x in left) + tuple(x << n for x in right)


def transitive(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def wrap(seed, mark, word, clock):
    block = reader.graph(seed, mark, word)
    rows = ordinal(transitive(clock), block)
    return reader.birth(rows, ((1 << len(rows)) - 1) ^ (1 << (clock + mark)))


def determinant(rows):
    a = [list(map(Fraction, row)) for row in rows]
    answer = Fraction(1)
    for i in range(len(a)):
        pivot = next(j for j in range(i, len(a)) if a[j][i])
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            answer *= -1
        z = a[i][i]
        answer *= z
        for j in range(i + 1, len(a)):
            ratio = a[j][i] / z
            a[j] = [x - ratio * y for x, y in zip(a[j], a[i])]
    need(answer.denominator == 1, 'Integral determinant')
    return int(answer)


def main():
    literal = 0
    for r in range(1, 128):
        rows = one_apex(reader.decode_r(r))
        count = reader.tables(rows)[0][-1]
        need(count == excess(r) + 1, 'Literal one-apex Hamilton count')
        literal += 1
    for n in range(256):
        a, e = excess(n), n & 1
        need(a % 2 == 0 and 0 <= a <= n, 'Bounds and parity')
        need(excess(2*n) == a+n+e and excess(2*n+1) == a, 'Digital sections')
        need(excess(4*n)-excess(2*n) == 2*n, 'Third integral kernel basis')
        need(excess(4*n+2)-a-2*n == 2, 'Fourth integral kernel basis')
        need(sections.excess(n) == a, 'Independent word accumulator')
    det = determinant([[excess(n), n+(n&1), 2*n, 2] for n in range(4)])
    need(abs(det) == 16, 'Kernel basis independence')
    for m in range(2, 17):
        need(excess((1 << m)-2) == 1 << (m-1), 'Global monotonicity hostile left')
        need(excess((1 << m)-1) == 0, 'Global monotonicity hostile right')
    K,q,d = 1620,11,2
    descending = [K*r+q*excess(r)+d*(r&1) for r in (510,511)]
    need(descending == [829016,827822], 'Actual-seed global-monotonicity failure')

    hostile = []
    for m in range(2, 9):
        for k in range(2, 9):
            r = ((1 << m)-1) << k
            pairs = ['11']*(m-1)+['10']+['00']*(k-2)+['01','00']
            x = int(''.join(z[0] for z in pairs), 2)
            y = int(''.join(z[1] for z in pairs), 2)
            need(x == r and y == r-(1 << k)+2, 'Regular pattern values')
            need(excess(r) == r-(1 << m)+2, 'All-block formula control')
            need((y == excess(r)) == (m == k), 'Equal-block hostile intersection')
            hostile.append((m, k, r, excess(r), y))

    inversions = 0
    for K in range(1, 5):
        for q in range(1, 5):
            for d in range(q+1):
                for residual in range(-2, 81):
                    expected = [r for r in range(max(0, residual//K)+1)
                                if K*r+q*excess(r)+d*(r&1) == residual]
                    got, _ = sections.sections(K, q, d, residual)
                    need(got == expected, 'Exact section inverse versus bounded enumeration')
                    inversions += 1

    C3, R4 = (2,4,1), (2,12,9,1)
    prefixes = (C3, R4)
    seeds = ((0,), (2,0), C3)
    words = ('', '0', '10')
    grades = []
    for prefix in prefixes:
        pf = reader.all_fields(prefix)[0]
        H, P2 = pf[0], pf[6]
        for seed in seeds:
            for mark in range(len(seed)):
                for word in words:
                    block = reader.graph(seed, mark, word)
                    old = reader.all_fields(block)[mark]
                    joined = ordinal(prefix, block)
                    actual = reader.all_fields(joined)[len(prefix)+mark]
                    scales = (H,H,H,P2,P2,H,P2,P2)
                    need(actual == tuple(a*b for a,b in zip(old,scales)), 'Literal eight-field grades')
                    need(joined == reader.graph(ordinal(prefix,seed),len(prefix)+mark,word), 'Reader commutes with prefix')
                    grades.append((len(prefix),len(seed),mark,word,H,P2))
    need(reader.all_fields(C3)[0][6] == 12, 'C3 grade')
    need((reader.all_fields(R4)[0][0],reader.all_fields(R4)[0][6]) == (5,32), 'Primitive R4 grade')
    strong_controls = 0
    for seed in seeds:
        for mark in range(len(seed)):
            for word in ('0','10','00','110'):
                rows = wrap(seed,mark,word,1)
                need(reader.strong(rows), 'Weak-seed nonempty zero-ending wrapper')
                f = reader.run(reader.all_fields(seed)[mark],word)
                need(reader.tables(rows)[0][-1] == f[0]+2*f[3], 'Literal wrapper count')
                strong_controls += 1
    need(not reader.strong(wrap((2,0),0,'',1)), 'Empty-word strongness hostile')

    moduli = {}
    for modulus in (49,343):
        checked = 0
        # This checks scalar modular closure only, not all eight fields.
        for r in range(256):
            state = (0,0,0,1)
            for bit in bin(r)[2:]:
                a,n,e,one = state
                state = ((a+n+e,2*n,0,one) if bit == '0'
                         else (a,2*n+one,one,one))
                state = tuple(x%modulus for x in state)
            need(state == tuple(x%modulus for x in (excess(r),r,r&1,1)), 'Scalar modular weighted carrier')
            checked += 1
        moduli[str(modulus)] = checked
    return {
        'status':'FINITE-EXACT controls; all-order arguments are in NOTE.md',
        'source_sha256':{name:sha256(path.read_bytes()).hexdigest() for name,path in SOURCES.items()},
        'literal_one_apex_graphs':literal,
        'digital_recurrence_inputs':256,
        'kernel_basis_determinant':det,
        'actual_seed_adjacent_descent':descending,
        'regular_pattern_pairs':len(hostile),
        'regular_pattern_diagonal_hits':sum(a==b for a,b,*_ in hostile),
        'inverse_parameter_targets':inversions,
        'literal_prefix_grade_cases':len(grades),
        'literal_weak_seed_strong_wrappers':strong_controls,
        'empty_word_hostile':'seed rows[2,0],mark0,clock1 is not strong',
        'scalar_modular_inputs':moduli,
        'limitations':'No fixed-bank census, exact spectrum claim, or arbitrary joint-attachment reduction.',
    }


if __name__ == '__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))

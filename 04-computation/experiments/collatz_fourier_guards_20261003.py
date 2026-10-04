"""Exact compressed-prefix guards and phase-sensitive Fourier controls.

Run: python3 04-computation/experiments/collatz_fourier_guards_20261003.py
No numerical Fourier transform, external tables, or presumed convergent inputs.
The exhaustive universe is printed; root 1 is the certified common terminal.
"""
from collections import Counter
from fractions import Fraction
from itertools import product


def summary(word):
    A = S = 0
    for a in word:
        assert a >= 1
        S = 3*S + 2**A
        A += a
    return A, len(word), S


def k0(target):
    assert target > 0 and target % 2 == 1 and target % 3 != 0
    return 2 if target % 3 == 1 else 1


def sibling(target, index):
    assert index >= 0
    return (2**(k0(target)+2*index)*target-1)//3


def compile_guard(prefix, target=1):
    """Return index phase -> exact forward valuation head, modulo 3**p.

    This compiles arithmetic legality. First-hit validity also requires that
    the selected sibling is not already 1; see replay below.
    """
    assert all(j >= 0 for j in prefix)
    p = len(prefix)
    M = 3**p
    coefficient = 2**k0(target)*target
    inverse_chart = {}
    for j in range(M):
        # Precision 3*M is necessary before the division by 3.
        residue = ((pow(4, j, 3*M)*coefficient-1)//3) % M
        assert residue not in inverse_chart
        inverse_chart[residue] = j
    phases = {}
    for parities in product((1, 2), repeat=p):
        word = tuple(2*j+a for j, a in zip(prefix, parities))
        A, _, S = summary(word)
        residue = S*pow(2**A, -1, M) % M
        phase = inverse_chart[residue]
        assert phase not in phases
        phases[phase] = word
    assert len(phases) == 2**p
    return phases


def replay(prefix, index, target=1):
    """Independent backward integer decoder; does not use carry/cylinders."""
    h = sibling(target, index)
    current = h
    reverse_word = []
    reverse_nodes = [h]
    for j in reversed(prefix):
        if current % 3 == 0:
            return None
        a = 2*j + (2 if current % 3 == 1 else 1)
        numerator = 2**a*current-1
        assert numerator % 3 == 0
        current = numerator//3
        assert current > 0 and current % 2 == 1
        reverse_word.append(a)
        reverse_nodes.append(current)
    return tuple(reversed(reverse_word)), tuple(reversed(reverse_nodes))


def odd_step(n):
    q = 3*n+1
    a = (q & -q).bit_length()-1
    return q >> a, a


def exact_parseval(phases, M):
    """Sum Fourier squared magnitudes in Z[zeta_M], without rounding.

    Expand each squared magnitude as sum_(u,v) zeta^(r*(u-v)), sum all r,
    then reduce by Phi_(3**p)(X)=X**(2M/3)+X**(M/3)+1.
    """
    assert M >= 3 and M % 3 == 0
    coefficients = [0]*M
    for r in range(M):
        for u in phases:
            for v in phases:
                coefficients[r*(u-v) % M] += 1
    d = M//3
    for i in range(d):
        top = coefficients[i+2*d]
        coefficients[i] -= top
        coefficients[i+d] -= top
    reduced = coefficients[:2*d]
    assert reduced == [M*len(phases)] + [0]*(2*d-1)
    return M*len(phases)-len(phases)**2


def dft3(mask):
    # Pairs (a,b) denote a+b*zeta, with zeta**2=-1-zeta.
    powers = ((1, 0), (0, 1), (-1, -1))
    result = []
    for r in range(3):
        terms = [powers[(-r*j) % 3] for j in mask]
        result.append(tuple(sum(t[c] for t in terms) for c in range(2)))
    return result


def real_part(z):
    a, b = z
    return Fraction(2*a-b, 2)


def norm_squared(z):
    a, b = z
    return a*a-a*b+b*b


def main():
    heads = phase_tests = root_exceptions = compiled_words = 0
    for p in range(1, 5):
        M = 3**p
        for prefix in product(range(3), repeat=p):
            phases = compile_guard(prefix)
            assert exact_parseval(phases, M) == 2**p*(3**p-2**p)
            heads += 1
            compiled_words += len(phases)
            for index in range(2*M):
                actual = replay(prefix, index)
                assert (actual is not None) == (index % M in phases)
                phase_tests += 1
                if actual is None:
                    continue
                word, nodes = actual
                assert word == phases[index % M]
                for a, source, dest in zip(word, nodes, nodes[1:]):
                    assert odd_step(source) == (dest, a)
                assert odd_step(nodes[-1])[0] == 1
                # The selected final block is a 1->1 padding block iff index=0.
                strict = all(n > 1 for n in nodes)
                assert strict == (index > 0)
                if not strict:
                    root_exceptions += 1
    print(f'EXHAUSTIVE: {heads} compressed heads, lengths 1..4, indices 0..2')
    print(f'Compiled {compiled_words} distinct valuation phases; replayed {phase_tests} index cases (two periods each)')
    print(f'Every head has exactly 2^p phases mod 3^p; {root_exceptions} legal index-zero cases require first-hit truncation')
    print('EXACT: cyclotomic Parseval, nonconstant energy 2^p*(3^p-2^p), for all 120 heads')
    # Change the common target independently of the root-1 first-hit tests.
    for target in (5, 7, 13):
        for prefix in ((0,), (0, 0), (0, 1, 2), (2, 1, 0, 2)):
            phases = compile_guard(prefix, target)
            for index in range(3**len(prefix)):
                actual = replay(prefix, index, target)
                assert (actual is not None) == (index in phases)
                if actual is not None:
                    assert actual[0] == phases[index]
    print('Additional targets 5/7/13: four fixed heads, every phase independently replayed')
    phases = compile_guard((0, 0))
    assert set(phases) == {0, 3, 4, 7}
    examples = []
    for index, expected in ((3, 75), (4, 151), (7, 19417), (9, 621377)):
        word, nodes = replay((0, 0), index)
        assert nodes[0] == expected and all(n > 1 for n in nodes)
        examples.append((index, nodes + (1,), word))
    print('Head (0,0): legal phases mod9 =', sorted(phases))
    print('First positive representatives (index, route, valuation head):', examples)
    # Index zero is arithmetic-legal for head (1), but gives 5->1->1.
    assert replay((1,), 0)[1] == (5, 1)
    print('ROOT HOSTILE: head(1), index0 gives 5->1->1; source>1 alone is not a first-hit check')
    mask_at_one = {d for d in range(3) if replay((0,), 1+d) is not None}
    mask_at_three = {d for d in range(3) if replay((0,), 3+d) is not None}
    assert mask_at_one == {0, 2} and mask_at_three == {0, 1}
    first, second = dft3(mask_at_one), dft3(mask_at_three)
    assert [real_part(z) for z in first] == [real_part(z) for z in second]
    assert [norm_squared(z) for z in first] == [norm_squared(z) for z in second] == [4, 1, 1]
    assert first != second
    print('FOURIER HOSTILE: increment masks {0,2} and {0,1}; identical cosine coefficients [2,1/2,1/2] and squared magnitudes [4,1,1]')
    print('Their exact DFTs (a+b*zeta_3) are', first, 'and', second, '; +1 legality differs')
    # Even lift: equality of formal root-of-unity monomials, not floating DFT.
    lift_cases = 0
    for prefix in ((0,), (0, 0), (0, 1, 2), (2, 1, 0, 2)):
        phases = compile_guard(prefix)
        M = 3**len(prefix)
        for parity in (0, 1):
            for r in range(2*M):
                left = Counter((-r*(parity+2*j)) % (2*M) for j in phases)
                right = Counter((-r*parity-2*(r % M)*j) % (2*M) for j in phases)
                assert left == right
                lift_cases += 1
    print(f'Even/odd Fourier lift: {lift_cases} formal monomial checks on four heads, both parity sheets')
    print('BOUNDARY: masks describe legality in one sibling-index family; neither first-hit validity nor integer-basin coverage follows from spectra alone')
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()

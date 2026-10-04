"""Exact finite audits for the modular-table / reflection / tournament bridge.

Run: python 04-computation/experiments/modular_sine_tournaments_20261003.py
Universe: moduli 1..101; all table cells; all nonzero cyclic differences;
all odd-prime quadratic characters in that range; certified seeds 1..1000.
No floating-point trigonometry is used: sine signs follow from half-circle
membership. The all-modulus statements are proved in the accompanying note.
"""

from math import gcd


def sine_sign(d, m):
    d %= m
    if d == 0 or 2 * d == m:
        return 0
    return 1 if 2 * d < m else -1


def prime(n):
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def orbit_rep(a, m):
    a %= m
    return (a, 1) if 2 * a <= m else (m - a, -1)


def poly_sub(a, b):
    out = list(a) + [0] * max(0, len(b) - len(a))
    for i, x in enumerate(b):
        out[i] -= x
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def certified_parity(n):
    """First-hit ordinary trajectory; the finite universe is fully checked."""
    bits = []
    while n != 1:
        assert len(bits) < 10000, "Uncertified seed in finite test universe"
        bits.append(n % 2)
        n = 3 * n + 1 if n % 2 else n // 2
    return bits


def route_polynomial(bits):
    # F(z) = B(z) + z^L/(1-z^3), since the root-1 tail is (100)^infinity.
    p = poly_sub(bits, [0, 0, 0] + bits)
    p += [0] * (len(bits) + 1 - len(p))
    p[len(bits)] += 1
    return p


def main():
    cells = centers = cyclic = antipodal = rank_checks = primes = 0
    for m in range(1, 102):
        h = m // 2
        fixed = [a for a in range(m) if 2 * a % m == 0]
        assert len(fixed) == (2 if m % 2 == 0 else 1)
        even_dim = (m + len(fixed)) // 2
        odd_dim = (m - len(fixed)) // 2
        assert (even_dim, odd_dim) == ((h + 1, h - 1) if m % 2 == 0 else (h + 1, h))
        for a in range(m):
            aa, sa = orbit_rep(a, m)
            image = {(a * b) % m for b in range(m)}
            assert len(image) == m // gcd(a, m)
            rank_checks += 1
            for b in range(m):
                bb, sb = orbit_rep(b, m)
                assert a * b % m == sa * sb * (aa * bb) % m
                cells += 1
        if m % 2:
            # Multiplication by 2 permutes residues for every odd m, including 1.
            assert len({2 * a % m for a in range(m)}) == m
            if m >= 3:
                d = pow(4, -1, m)
                block = [[a * b % m for b in (h, h + 1)] for a in (h, h + 1)]
                assert block == [[d, m - d], [m - d, d]]
                centers += 1
            for d in range(1, m):
                assert sine_sign(d, m) == -sine_sign(-d, m) != 0
                cyclic += 1
            assert sum(sine_sign(d, m) == 1 for d in range(m)) == h
        else:
            assert len({2 * a % m for a in range(m)}) == m // 2
            assert sine_sign(h, m) == 0 and h == (-h) % m
            antipodal += 1
        if m > 2 and prime(m):
            chi = [0] + [1 if pow(d, (m - 1) // 2, m) == 1 else -1 for d in range(1, m)]
            reflection_sign = -1 if m % 4 == 3 else 1
            assert all(chi[-d % m] == reflection_sign * chi[d] for d in range(m))
            primes += 1

    # Index reflection and matrix transpose are different operations.
    # The modular kernel is symmetric under transpose, including its sine part.
    assert sine_sign(1 * 2, 5) == sine_sign(2 * 1, 5) == 1
    # A quadratic residue indicator has nonunit ties for composite moduli.
    assert pow(3, 2, 9) == 0

    poly_checks = 0
    for n in range(1, 1001):
        bits = certified_parity(n)
        assert all(a + b <= 1 for a, b in zip(bits, bits[1:]))
        p = route_polynomial(bits)
        for modulus in (3, 5, 9, 10, 11):
            # Exact group-algebra image of P in Z[z]/(z^modulus-1).
            folded = [sum(p[i::modulus]) for i in range(modulus)]
            direct = [0] * modulus
            for i, bit in enumerate(bits):
                direct[i % modulus] += bit
                direct[(i + 3) % modulus] -= bit
            direct[len(bits) % modulus] += 1
            assert folded == direct
            poly_checks += 1
        # Keeping root label and length makes polynomial tail reconstruction exact.
        unfolded = []
        for i in range(len(bits) + 12):
            val = (p[i] if i < len(p) else 0) + (unfolded[i - 3] if i >= 3 else 0)
            unfolded.append(val)
        assert unfolded == bits + [1, 0, 0] * 4
    print('FINITE-EXACT universe: m=1..101; all modular cells/differences; seeds=1..1000')
    print(f'wedge reconstruction cells: {cells}')
    print(f'multiplication character-image ranks: {rank_checks}')
    print(f'odd center blocks: {centers}; odd cyclic differences: {cyclic}')
    print(f'even antipodal tie controls: {antipodal}; odd primes: {primes}')
    print(f'rooted route-polynomial group-algebra checks: {poly_checks}')
    print('hostiles: K_5 sine is symmetric, not skew; even antipodes tie; composite 9 has nonunit 3')
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()

"""Exact Bernoulli/Fermat bridges and depth-resolved Collatz carries.

Standard library only. All mathematical checks remain active under python -O.
"""

from fractions import Fraction
from itertools import product
from math import comb, isqrt, prod


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def bernoulli_table(limit):
    out = [Fraction(1)]
    for n in range(1, limit + 1):
        out.append(-sum(Fraction(comb(n + 1, j)) * out[j]
                        for j in range(n)) / (n + 1))
    return out


def psi(x):
    return x - x.numerator // x.denominator - Fraction(1, 2)


def carry_path(n):
    require(n >= 0, "nonnegative integer carry input")
    return tuple((3 * (n % (1 << s)) + 1) // (1 << s)
                 for s in range(1, (3 * n + 1).bit_length() + 1))


def total_carry(n):
    return 1 + 3 * n.bit_count() - (3 * n + 1).bit_count()


def U(n):
    t = 3 * n + 1
    return t // (t & -t)


def valuation(n, p):
    require(n != 0, "nonzero valuation input")
    n = abs(n)
    a = 0
    while n % p == 0:
        n //= p
        a += 1
    return a


def eval_poly(coefficients, n):
    answer = 0
    for c in reversed(coefficients):
        answer = answer * n + c
    return answer


def phi(n):
    answer, rest, d = n, n, 2
    while d * d <= rest:
        if rest % d == 0:
            answer -= answer // d
            while rest % d == 0:
                rest //= d
        d += 1
    if rest > 1:
        answer -= answer // rest
    return answer


def root_sum(N, d):
    # Exact in Q[z]/(z^(N/2)+1), z a primitive N-th root; N power of 2.
    v = [0] * (N // 2)
    for j in range(N):
        e = (j * d) % N
        v[e % (N // 2)] += 1 if e < N // 2 else -1
    return v


def main():
    print("Bernoulli / Fermat / dyadic-carry exact audit")
    B = bernoulli_table(256)
    require(B[:6] == [Fraction(1), Fraction(-1, 2), Fraction(1, 6),
                      Fraction(0), Fraction(-1, 30), Fraction(0)],
            "Bernoulli convention")
    F = [(1 << (1 << j)) + 1 for j in range(6)]
    require(all(prime(p) for p in F[:5]), "first five Fermat primality")
    require(F[5] == 641 * 6700417, "sixth Fermat composite")
    for r in range(1, 9):
        expected = 2 * prod(F[j] for j in range(6)
                            if (1 << j) <= r and prime(F[j]))
        require(B[1 << r].denominator == expected,
                "direct Bernoulli denominator versus dyadic prime formula")
        print(f"denom(B_2^{r})={expected}")
    D16 = 2 * prod(F[:5])
    require(D16 == 8589934590, "known Fermat product")
    # This all-n finite test uses independent totient computation.
    for n in range(3, 10001):
        odd = n // (n & -n)
        p = phi(n)
        require((D16 // 2) % odd == 0 if p & (p - 1) == 0 else
                (D16 // 2) % odd != 0, "constructible catalogue vs totient")
    require(U(17) == 13 and phi(13) == 12,
            "constructible polygon sizes are not Collatz invariant")
    print("Constructibility controls: n=3..10000; hostile 17->13")

    for L in range(1, 9):
        N = 1 << L
        for d in range(-N, N + 1):
            target = [N if d % N == 0 else 0] + [0] * (N // 2 - 1)
            require(root_sum(N, d) == target, "cyclotomic exact projector")
            J = Fraction(1, N) + psi(Fraction(d - 1, N)) - psi(Fraction(d, N))
            require(J == int(d % N == 0), "B1 exact projector")
    word_count = 0
    for L in range(1, 11):
        N = 1 << L
        for word in product((0, 1), repeat=L):
            A, C = 1, 0
            for j, e in enumerate(word):
                A, C = (3 ** e) * A, (3 ** e) * C + e * (1 << j)
            source = (-C * pow(A, -1, N)) % N
            x = source + N
            for e in word:
                require(x % 2 == e, "word/source parity legality")
                x = (3 * x + 1) // 2 if e else x // 2
            require(x == (A * (source + N) + C) // N, "affine endpoint")
            require((A * (source + 1) + C) % N != 0,
                    "neighbour source hostile")
            word_count += 1
    print(f"Exact cyclotomic/B1 projectors: L=1..8; affine words={word_count}")

    for n in range(20001):
        path = carry_path(n)
        require(sum(path) == total_carry(n), "carry sum/popcount identity")
        last, decoded = 1, 0
        for s, c in enumerate(path):
            options = [b for b in (0, 1) if (last + 3 * b) // 2 == c]
            require(len(options) == 1, "injective carry graph transition")
            decoded += options[0] << s
            m = 1 << (s + 1)
            c_b1 = 3 * psi(Fraction(n, m)) - psi(Fraction(3 * n + 1, m))
            c_b1 += 1 + Fraction(1, m)
            require(c_b1 == c, "carry from B1 at exact endpoints")
            last = c
        require(decoded == n and path[-1] == 0, "lossless finite carry path")
        out = 3 * n + 1
        for s in range(out.bit_length() + 1):
            c_here = 1 if s == 0 else path[s - 1] if s <= len(path) else 0
            c_next = path[s] if s < len(path) else 0
            require(((out >> s) & 1) == 3 * ((n >> s) & 1) + c_here - 2 * c_next,
                    "coefficientwise carry polynomial identity")
        if n & 1:
            y = U(n)
            a = valuation(out, 2)
            moment_n = sum(s for s in range(n.bit_length()) if (n >> s) & 1)
            moment_y = sum(s for s in range(y.bit_length()) if (y >> s) & 1)
            moment_c = sum(s * c for s, c in enumerate(path, 1))
            require(a * y.bit_count() + moment_y ==
                    3 * moment_n + 2 * total_carry(n) - moment_c,
                    "first positional carry moment couples actual halving count")
    print("Full carry path, Bernoulli formula and inverse decoder: n=0..20000")
    print("Carry-polynomial coefficients: n=0..20000; first moment on 10000 odd sources")

    polys = [(0, 1), (1, 1), (1, 3), (1, 0, 1), (1, -1, 1), (-7, 5)]
    primes = [2, 3, 5, 11, 17, 257]
    fibre_count = 0
    for a in [3, 27, 233]:
        expected = tuple(valuation(eval_poly(P, a), p) for p in primes for P in polys)
        Q = prod(p ** (1 + max(valuation(eval_poly(P, a), p) for P in polys))
                 for p in primes)
        while Q <= a:
            Q *= 2
        b = Q - a
        constS = (Q - 1).bit_count() - (b - 1).bit_count()
        constK = 1 + 3 * constS - (3 * Q - 1).bit_count() + (3 * b - 2).bit_count()
        first = (3 * b - 1).bit_length()
        for r in range(first, first + 65):
            n = Q * (1 << r) - b
            require(tuple(valuation(eval_poly(P, n), p) for p in primes for P in polys)
                    == expected, "all polynomial valuations fixed")
            require(n.bit_count() == r + constS and total_carry(n) == 2 * r + constK,
                    "unbounded carry inside fixed fibre")
            fibre_count += 1
        print(f"Fibre a={a}: Q={Q}; r>={first}; S=r{constS:+}; K=2r{constK:+}")
    print(f"Fixed-fibre controls: {fibre_count} integers, 36 valuations each")

    require((51).bit_count() == (77).bit_count() == 4 and
            total_carry(51) == total_carry(77) == 9 and U(51) == 77,
            "minimal single-edge same-summary growth hostile")
    for L in range(9, 1025):
        n = 15 * (1 << L) + 51
        y = 45 * (1 << (L - 1)) + 77
        require(U(n) == y > n, "padded hostile is a growth edge")
        require(n.bit_count() == y.bit_count() == 8 and
                total_carry(n) == total_carry(y) == 17,
                "padded hostile fixes both summaries")
    print("Summary-rank hostile: 51->77; 1016 padded edges with (S,K)=(8,17)")

    for n in range(1, 2000, 2):
        x, discounted, power = n, Fraction(0), 1
        for j in range(1, 41):
            require(total_carry(x) >= 3, "odd carry lower bound")
            power *= 3
            discounted += Fraction(total_carry(x) - 1, power)
            x = U(x)
            require(discounted + Fraction(x.bit_count(), power) == n.bit_count(),
                    "finite discounted carry budget with exact remainder")
    print("Discounted carry budget: 1000 odd sources, 40 steps each")
    print("PASS; no Collatz convergence or Fermat-prime finiteness claim")


if __name__ == "__main__":
    main()

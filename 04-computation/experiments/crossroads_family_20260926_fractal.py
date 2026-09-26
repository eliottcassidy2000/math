"""Exact controls for parity-language covers and a two-completion IFS hostile.

Stdlib only. All checks survive python -O. No finite check is a Collatz proof.
"""
from fractions import Fraction
from itertools import product
from math import comb, log2, log


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def affine(word):
    a, c = 1, 0
    for j, bit in enumerate(word):
        if bit:
            a, c = 3*a, 3*c + 2**j
    return a, c


def positive(word):
    a = 1
    for j, bit in enumerate(word, 1):
        if bit:
            a *= 3
        if a <= 2**j:
            return False
    return True


def counts(depth):
    row, answer = {0: 1}, [1]
    for j in range(1, depth+1):
        nxt = {}
        for k, value in row.items():
            for b in (0, 1):
                if 3**(k+b) > 2**j:
                    nxt[k+b] = nxt.get(k+b, 0) + value
        row = nxt
        answer.append(sum(row.values()))
    return answer


def v2(n):
    n = abs(n)
    if not n:
        raise ValueError("zero valuation is infinite")
    return (n & -n).bit_length()-1


def main():
    w = counts(512)
    for k in range(1, 17):
        brute = sum(positive(bits) for bits in product((0, 1), repeat=k))
        check(brute == w[k], "brute/DP count")
    # Independent Spitzer recurrence, using exact binomial tails.
    tails = [0] + [sum(comb(k, j) for j in range(k+1) if 3**j > 2**k)
                   for k in range(1, 129)]
    for k in range(1, 129):
        check(k*w[k] == sum(tails[j]*w[k-j] for j in range(1, k+1)),
              "Spitzer identity")
    rho = log(2)/log(3)
    h = -rho*log2(rho)-(1-rho)*log2(1-rho)
    print("EXACT: ballot DP to 512; brute binary words to 16; Spitzer to 128")
    print("h (display approximation only):", format(h, ".12f"))
    for k in (8, 16, 32, 64, 128, 256, 512):
        print("block", k, "W", w[k], "dimension approx", format(log2(w[k])/k, ".12f"))
    fixed = 0
    for k in range(1, 9):
        good = [bits for bits in product((0, 1), repeat=k) if positive(bits)]
        for bits in good:
            a, c = affine(bits)
            x = Fraction(c, 2**k-a)
            check(x < 0 and x.denominator % 2 == 1, "negative periodic point")
            modulus = 2**k
            residue = (-c*pow(a, -1, modulus)) % modulus
            n = residue
            got = []
            for _ in range(k):
                b = n % 2
                got.append(b)
                n = (3*n+1)//2 if b else n//2
            check(tuple(got) == bits, "affine inverse branch parity")
            check(positive(bits+(1,)*k), "extendible positive cylinder")
            fixed += 1
        # Every pair of allowed blocks concatenates, for this small universe.
        for left in good:
            for right in good:
                check(positive(left+right), "positive block concatenation")
    print("EXACT: inverse residues, negative periodic points and block concatenation;", fixed, "words")
    # The same rational sequence converges to a negative real and to +1 in Z_2.
    x, a, c = 1, 1, 0
    saved = {}
    for j in range(1, 81):
        digit = 1 if x % 2 else 2
        x = (3*x+digit)//2
        # Prefix composition phi_d1 o ... o phi_dj at zero is -c/3^j.
        c = 3*c + 2**(j-1)*digit
        a *= 3
        approximant = Fraction(-c, a)
        error = approximant-1
        check(error == -Fraction(2**j*x, 3**j), "two-place exact identity")
        check(v2(error.numerator)-v2(error.denominator) >= j, "2-adic convergence")
        check(-2 <= approximant < 0, "real invariant interval with initial zero")
        if j in (10, 20, 40, 80):
            saved[j] = (float(approximant), v2(error.numerator), x)
    for j, (real, precision, orbit) in saved.items():
        print("two-place", j, "real approx", format(real, ".12f"),
              "2-adic agreement bits", precision, "positive inverse orbit", orbit)
    for k in range(1, 13):
        modulus = 2**k
        inv3 = pow(3, -1, modulus)
        images = [{((2*y-d)*inv3) % modulus for y in range(modulus)} for d in (1, 2)]
        check(not images[0] & images[1], "disjoint 2-adic branches")
        check(images[0] | images[1] == set(range(modulus)), "full 2-adic cover")
    print("EXACT: two-loop IFS partitions every dyadic quotient through 2^12")
    print("PASS. Scope: cover counts and exact controls, not positive-integer exclusion.")


if __name__ == "__main__":
    main()

"""Exact controls for cubic depth, missing phases, and golden ideal branches.

Universe: six characteristic-two fields, four factored binary cyclotomic
layers, and eight golden norm layers. Standard library only; checks survive -O.
Run from the repository root with python -X utf8 -B this_file.py.
"""

from math import gcd, isqrt, prod


CHECKS = 0

# Exact Pocklington witnesses (q, exponent in F, witness). Every supporting
# q is small enough for exhaustive trial division; no probable-prime test.
POCKLINGTON = {
    16753783618801: [(2, 4, 7), (3, 5, 2), (5, 2, 3),
                     (19, 1, 3), (9071791, 1, 3)],
    192971705688577: [(2, 9, 5), (3, 5, 2), (2207, 1, 3)],
    3712990163251158343: [(2, 1, 3), (3, 6, 5), (7, 1, 3), (382021, 1, 3)],
}


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(message)


def prime(n):
    if n in POCKLINGTON:
        certificate = POCKLINGTON[n]
        F = prod(q**e for q, e, _ in certificate)
        check((n-1) % F == 0 and F*F > n, "Pocklington factor bound")
        for q, e, a in certificate:
            check(prime(q) and e > 0, "Pocklington support prime")
            check(pow(a, n-1, n) == 1, "Pocklington Fermat return")
            check(gcd(pow(a, (n-1)//q, n)-1, n) == 1,
                  "Pocklington order exclusion")
        return True
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def certify_factorization(n, factors):
    check(prod(p**e for p, e in factors.items()) == n, "factor product")
    for p, e in factors.items():
        check(prime(p) and e > 0, "certified prime factor")


def polymod(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def polygcd(a, b):
    while b:
        a, b = b, polymod(a, b)
    return a


class BinaryField:
    def __init__(self, a):
        self.a = a
        self.m = 3 ** (a - 1)
        self.d = 2 * self.m
        self.f = (1 << self.d) | (1 << self.m) | 1

    def mul(self, x, y):
        z = 0
        while y:
            if y & 1:
                z ^= x
            y >>= 1
            x <<= 1
            if x >> self.d:
                x ^= self.f
        return z

    def power(self, x, e):
        z = 1
        while e:
            if e & 1:
                z = self.mul(z, x)
            x = self.mul(x, x)
            e >>= 1
        return z

    def irreducibility_check(self):
        # Rabin's criterion, independent of the number-theoretic proof.
        x = 2
        squares = [x]
        for _ in range(self.d):
            x = self.mul(x, x)
            squares.append(x)
        check(squares[self.d] == 2, "Frobenius closes")
        for p in (2, 3):
            if self.d % p == 0:
                check(polygcd(squares[self.d // p] ^ 2, self.f) == 1,
                      "Rabin proper subfield exclusion")

    def order(self, x, factors):
        value = 2**self.d - 1
        for p in factors:
            while value % p == 0 and self.power(x, value // p) == 1:
                value //= p
        check(self.power(x, value) == 1, "claimed order returns")
        for p in factors:
            if value % p == 0:
                check(self.power(x, value // p) != 1, "order minimality")
        return value


FIELD_FACTORS = [
    {3: 1},
    {3: 2, 7: 1},
    {3: 3, 7: 1, 19: 1, 73: 1},
    {3: 4, 7: 1, 19: 1, 73: 1, 262657: 1, 87211: 1},
    {3: 5, 7: 1, 19: 1, 73: 1, 163: 1, 2593: 1, 135433: 1,
     262657: 1, 87211: 1, 71119: 1, 97685839: 1, 272010961: 1},
]


def binary_layer(k):
    x = 2 ** (3 ** (k - 1))
    return x*x - x + 1


def finite_fields():
    print("FINITE FIELDS: a, dimension, order(1+zeta), index, order(normalized)")
    for a, factors in enumerate(FIELD_FACTORS, 1):
        field = BinaryField(a)
        field.irreducibility_check()
        size = 2**field.d - 1
        certify_factorization(size, factors)
        n, q, zeta, b = 3**a, 2**field.m, 2, 3
        check(field.order(zeta, factors) == n, "cyclotomic root order")
        check(field.power(b, q - 1) == field.power(zeta, n - 1),
              "half-field Frobenius relation")
        c = field.mul(b, field.power(zeta, (n - 1) // 2))
        check(field.power(c, q) == c, "normalized subfield element")
        ob, oc = field.order(b, factors), field.order(c, factors)
        check(ob == n * oc, "coprime phase/subfield order decomposition")
        compulsory = (q + 1) // n
        check((q + 1) % n == 0, "integral missing cofactor")
        check(compulsory == prod(binary_layer(j)//3 for j in range(2, a)),
              "telescoping missing cyclotomic layers")
        check(size // ob >= compulsory, "unavoidable missed phases")
        check(oc == q - 1, "finite-only normalized primitivity, a<=5")
        # Canonical multiplicative splitting; exhaustive at dimensions2,6.
        def split(x):
            r = field.power(field.power(x, q + 1), q // 2)
            v = field.mul(x, field.power(r, q - 2))
            if compulsory == 1:
                vn, vd = v, 1
            else:
                en = compulsory * pow(compulsory, -1, n)
                ed = n * pow(n, -1, compulsory)
                vn, vd = field.power(v, en), field.power(v, ed)
            check(field.power(r, q - 1) == 1, "subfield component")
            check(field.power(vn, n) == 1, "ternary phase component")
            check(field.power(vd, compulsory) == 1, "missing phase component")
            check(field.mul(r, field.mul(vn, vd)) == x, "three-factor decoder")
            return r, vn, vd

        samples = range(1, size + 1) if a <= 2 else range(1, 17)
        for x in samples:
            split(x)
        check(split(b) == (c, field.power(zeta, (n + 1)//2), 1),
              "b has trivial missing-phase coordinate")
        if a >= 2:
            previous_degree = 2 * 3 ** (a - 2)
            norm_exponent = size // (2**previous_degree - 1)
            check(field.power(b, norm_exponent) == (1 ^ field.power(zeta, 3)),
                  "relative norm compatibility")
            lower_zeta = field.power(zeta, 3)
            lower_c = field.mul(1 ^ lower_zeta,
                                field.power(lower_zeta, (3**(a-1)-1)//2))
            check(field.power(c, norm_exponent) == lower_c,
                  "normalized relative norm compatibility")
        print(a, field.d, ob, size // ob, oc)
    field = BinaryField(3)
    check(field.power(3, (2**18 - 1)//19) == 1, "hostile: generator misses19")
    check(field.power(3, (2**18 - 1)//73) != 1, "control: retains73")
    print("Hostile at dimension18: order13797, index19; 73 is retained.")


def modular_order(a, modulus, candidate, factors):
    check(pow(a, candidate, modulus) == 1, "modular order return")
    for p in factors:
        check(candidate % p == 0 and pow(a, candidate // p, modulus) != 1,
              "modular order minimality")


def binary_layers():
    factors = [{3: 1}, {3: 1, 19: 1}, {3: 1, 87211: 1},
               {3: 1, 163: 1, 135433: 1, 272010961: 1}]
    print("BINARY CYCLOTOMIC LAYERS: k, Phi_(2*3^k)(2), new factors")
    for k, fac in enumerate(factors, 1):
        value = binary_layer(k)
        certify_factorization(value, fac)
        check(value % 3 == 0 and value % 9 != 0, "one inherited3")
        for p in fac:
            if p != 3:
                modular_order(2, p, 2 * 3**k, (2, 3))
        print(k, value, {p: e for p, e in fac.items() if p != 3})


def sixth_field():
    """A complete order certificate without factoring the unused q+1."""
    field = BinaryField(6)
    field.irreducibility_check()
    q, n = 2**243, 729
    factors = {p: 1 for p in (7, 73, 487, 2593, 71119, 262657, 97685839,
                              16753783618801, 192971705688577,
                              3712990163251158343)}
    certify_factorization(q-1, factors)
    c = field.mul(3, field.power(2, (n-1)//2))
    check(field.power(c, q) == c, "sixth normalized subfield")
    check(field.power(c, q-1) == 1, "sixth normalized return")
    for p in factors:
        check(field.power(c, (q-1)//p) != 1, "sixth normalized minimality")
    order_b = n*(q-1)
    check(field.power(3, order_b) == 1, "sixth b return")
    for p in (3, *factors):
        check(field.power(3, order_b//p) != 1, "sixth b minimality")
    check(field.power(3, q-1) == field.power(2, n-1), "sixth Frobenius")
    missing = (q+1)//n
    check(missing == prod(binary_layer(j)//3 for j in range(2, 6)),
          "sixth missing-factor telescope")
    check((q*q-1)//order_b == missing, "sixth exact index")
    lower_zeta = field.power(2, 3)
    lower_c = field.mul(1 ^ lower_zeta, field.power(lower_zeta, (243-1)//2))
    relative_exponent = (q*q-1)//(2**162-1)
    check(field.power(c, relative_exponent) == lower_c, "sixth normalized norm")
    check(field.power(3, relative_exponent) == (1 ^ lower_zeta), "sixth b norm")
    print("SIXTH FIELD: dimension486; c primitive of exact order2^243-1.")
    print("Subfield factorization:", factors)
    print("b order:", order_b, "index:", missing)
    print("Large-prime certificates: three Pocklington proofs with trial-prime supports.")


def golden_mul(x, y):
    a, b = x
    c, d = y
    return a*c + b*d, a*d + b*c + b*d


def golden_power(x, exponent):
    value = (1, 0)
    while exponent:
        if exponent & 1:
            value = golden_mul(value, x)
        x = golden_mul(x, x)
        exponent >>= 1
    return value


def golden_norm(x):
    a, b = x
    return a*a + a*b - b*b


def golden_conjugate(x):
    a, b = x
    return a+b, -b


def golden_layers():
    values = []
    print("GOLDEN NORM LAYERS: k, N (or decimal digits), two root orders")
    for k in range(1, 9):
        d = 3**(k-1)
        u = golden_power((0, 1), d)
        u2 = golden_mul(u, u)
        A = (u2[0] + u[0] + 1, u2[1] + u[1])
        B = (u2[0] - u[0] + 1, u2[1] - u[1])
        N = golden_norm(A)
        lucas = 2*u[0] + u[1]
        check(N == lucas*lucas + 3 == golden_norm(B), "equal golden norms")
        check(golden_mul(A, B) == (N*u2[0], N*u2[1]), "golden product")
        check(golden_mul(u2, golden_conjugate(A)) == B, "conjugate branch")
        if values:
            old = values[-1]
            check(N == old**3 - 3*old**2 + 3, "golden cubic recurrence")
        for old in values:
            check(gcd(old, N) == 1 and N % old == 3 % old,
                  "all earlier norms coprime")
        check(N % 3 == 1 and N % 5 == 4, "unramified at3 and5")
        values.append(N)
        if k >= 2:
            check(N % 2 == 1, "odd golden branch norm")
            check(gcd(A[1], N) == gcd(B[1], N) == 1, "linear root decoders")
            r = (-A[0] * pow(A[1], -1, N)) % N
            s = (-B[0] * pow(B[1], -1, N)) % N
            check((r*r-r-1) % N == (s*s-s-1) % N == 0,
                  "golden quotient roots")
            check((r+s) % N == 1 and gcd(r-s, N) == 1, "CRT branch separation")
            modular_order(r, N, 3**k, (3,))
            modular_order(s, N, 2*3**k, (2, 3))
            print(k, N if k <= 4 else str(len(str(N))) + " digits", 3**k, 2*3**k)
        else:
            print(k, N, "coincident dyadic branches; excluded from odd branch statement")
    certify_factorization(values[3], {3079: 1, 62650261: 1})
    for p in (19, 5779, 3079, 62650261):
        check(prime(p), "golden factor primality")
    modular_order(2, 5779, 5778, (2, 3, 107))
    check(pow(2, 54, 5779) == 2944, "hostile to binary/golden clock identification")
    check(binary_layer(2)//3 == values[1] == 19, "shared19 control")
    check(binary_layer(3)//3 == 87211 != values[2] == 5779,
          "next-stage divergence")
    print("Hostile: N4=3079*62650261; norm layers need not be prime.")
    print("Hostile: ord_5779(2)=5778; golden branch orders27 and54.")


def main():
    finite_fields()
    sixth_field()
    binary_layers()
    golden_layers()
    print("PASS:", CHECKS, "exact checks; no randomized or floating-point tests.")


if __name__ == "__main__":
    main()

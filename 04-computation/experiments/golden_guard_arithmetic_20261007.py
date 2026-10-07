"""Exact joint golden/binary carrier and authenticated odd-word receipts.

See the matching proof note. All assertions remain active under python -O.
No floating point, third-party dependencies, or imported convergence claim.
"""

from dataclasses import dataclass
from itertools import product
from math import gcd, lcm


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(n, minimum=None):
    need(type(n) is int and (minimum is None or n >= minimum), "exact integer required")
    return n


@dataclass(frozen=True)
class State:
    a: int
    b: int
    debt: int

    def __post_init__(self):
        for n in (self.a, self.b, self.debt):
            integer(n)

    @property
    def source(self):
        return self.a + 2*self.b + self.debt


def checked(s):
    need(type(s) is State, "exact State required")
    return s


def append(s, digit):
    checked(s)
    integer(digit)
    return State(s.b + digit, s.a + s.b, 2*s.debt + s.b)


def add(s, t):
    checked(s); checked(t)
    return State(s.a+t.a, s.b+t.b, s.debt+t.debt)


def multiply(s, t):
    checked(s); checked(t)
    a, b, d = s.a, s.b, s.debt
    c, e, f = t.a, t.b, t.debt
    return State(a*c+b*e, a*e+b*c+b*e,
                 b*e+d*(c+2*e)+f*(a+2*b)+d*f)


def digits_state(digits):
    need(type(digits) is tuple and all(type(d) is int for d in digits), "integer digit tuple required")
    s = State(0, 0, 0)
    for d in digits:
        s = append(s, d)
    return s


def binary(n):
    integer(n, 0)
    return tuple(int(d) for d in bin(n)[2:])


def bits_checked(bits):
    need(type(bits) is tuple and bool(bits)
         and all(type(d) is int and d in (0, 1) for d in bits), "nonempty exact bit tuple required")
    return digits_state(bits)


def polynomial(p):
    need(type(p) is tuple and all(type(c) is int for c in p), "integer coefficient tuple required")
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def padd(p, q, scale=1):
    p, q = polynomial(p), polynomial(q)
    integer(scale)
    out = [0]*max(len(p), len(q))
    for i, c in enumerate(p): out[i] += c
    for i, c in enumerate(q): out[i] += scale*c
    return polynomial(tuple(out))


def pmul(p, q):
    p, q = polynomial(p), polynomial(q)
    if not p or not q: return ()
    out = [0]*(len(p)+len(q)-1)
    for i, c in enumerate(p):
        for j, d in enumerate(q): out[i+j] += c*d
    return polynomial(tuple(out))


def pscale(p, c):
    integer(c)
    return polynomial(tuple(c*x for x in polynomial(p)))


def shift(p, n):
    integer(n, 0)
    p = polynomial(p)
    return (0,)*n+p if p else ()


def peval(p, x):
    out = 0
    for c in reversed(polynomial(p)): out = out*x+c
    return out


def pread(p):
    return digits_state(tuple(reversed(polynomial(p))))


def remainder3(s):
    checked(s)
    # a+b*z+debt*(z*z-z-1)
    return polynomial((s.a-s.debt, s.b-s.debt, s.debt))


def divide_monic(p, q):
    p, q = list(polynomial(p)), polynomial(q)
    need(q and q[-1] == 1, "nonzero monic divisor required")
    out = [0]*max(0, len(p)-len(q)+1)
    while p and len(p) >= len(q):
        k, c = len(p)-len(q), p[-1]
        out[k] = c
        for j, d in enumerate(q): p[k+j] -= c*d
        while p and p[-1] == 0: p.pop()
    return polynomial(tuple(out)), tuple(p)


F = (-1, -1, 1)
JOINT = (2, 1, -3, 1)


def valuation(n):
    integer(n, 1)
    return (n & -n).bit_length()-1


def word_checked(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word), "positive exponent tuple required")
    return word


def carrier(word):
    """Ordered carry C(z), scalar P=3^r, cost A; C(2) is the usual B."""
    word_checked(word)
    p, a, c = 1, 0, ()
    for exponent in word:
        c = padd(pscale(c, 3), (0,)*a+(1,))
        p *= 3
        a += exponent
    return p, a, c


def replay(n, word, first_hit=False):
    integer(n, 1)
    need(n % 2 == 1, "positive odd source required")
    need(type(first_hit) is bool, "boolean first-hit flag required")
    word_checked(word)
    for exponent in word:
        need(not first_hit or n != 1, "premature ROOT")
        numerator = 3*n+1
        need(valuation(numerator) == exponent, "wrong actual valuation")
        n = numerator >> exponent
    if first_hit: need(n == 1, "ROOT endpoint required")
    return n


@dataclass(frozen=True)
class Receipt:
    source_bits: tuple
    target_bits: tuple
    word: tuple
    carry: tuple


def compile_receipt(source, word):
    target = replay(source, word)
    src, dst = binary(source), binary(target)
    p, a, c = carrier(word)
    numerator = padd(padd(pscale(tuple(reversed(src)), p), c),
                     shift(tuple(reversed(dst)), a), -1)
    quotient, remainder = divide_monic(numerator, (-2, 1))
    need(not remainder, "source/target identity failed")
    return Receipt(src, dst, word, quotient)


def verify_receipt(receipt, root=False):
    """Polynomial identity + positive odd endpoints implies exact word legality."""
    need(type(receipt) is Receipt and type(root) is bool, "typed receipt and flag required")
    src, dst = bits_checked(receipt.source_bits), bits_checked(receipt.target_bits)
    need(src.source > 0 and src.source % 2 == 1 and dst.source > 0 and dst.source % 2 == 1,
         "positive odd endpoints required")
    p, a, c = carrier(receipt.word)
    polynomial(receipt.carry)
    left = padd(padd(pscale(tuple(reversed(receipt.source_bits)), p), c),
                shift(tuple(reversed(receipt.target_bits)), a), -1)
    need(left == pmul((-2, 1), receipt.carry), "false carry identity")
    if root:
        # Extra chronological boundary, not inferred from an endpoint value.
        replay(src.source, receipt.word, first_hit=True)
    return src.source, dst.source


def concatenate(left, right):
    _, middle = verify_receipt(left)
    source, _ = verify_receipt(right)
    need(middle == source, "source seam mismatch")
    # Bit strings may have leading zeros; canonicalize the shared polynomial.
    need(polynomial(tuple(reversed(left.target_bits))) ==
         polynomial(tuple(reversed(right.source_bits))), "polynomial seam mismatch")
    p, _, _ = carrier(right.word)
    _, a, _ = carrier(left.word)
    carry = padd(pscale(left.carry, p), shift(right.carry, a))
    out = Receipt(left.source_bits, right.target_bits, left.word+right.word, carry)
    verify_receipt(out)
    return out


def golden_order(m):
    integer(m, 1)
    if m == 1: return 1
    a, b, k = 1, 0, 0
    while True:
        a, b = b, (a+b) % m
        k += 1
        if (a, b) == (1, 0): return k


def two_order(m):
    integer(m, 1)
    need(m % 2 == 1, "2 is not a unit modulo an even modulus")
    if m == 1: return 1
    k, x = 0, 1
    while True:
        x = 2*x % m
        k += 1
        if x == 1: return k


def joint_schedule(m):
    """Zero-padding: worst-case preperiod and exact eventual period."""
    integer(m, 1)
    k = valuation(m)
    odd = m >> k
    return k, lcm(golden_order(m), two_order(odd))


def orbit_to(n, endpoint=1):
    integer(n, 1); integer(endpoint, 1)
    need(n % 2 and endpoint % 2, "odd endpoints required")
    word = []
    for _ in range(10000):
        if n == endpoint: return tuple(word)
        a = valuation(3*n+1)
        word.append(a)
        n = (3*n+1) >> a
    raise ValueError("declared finite control budget exhausted")


def main():
    count = 0
    def check(ok, label):
        nonlocal count
        need(ok, label)
        count += 1

    # Complete small signed polynomial universe, independent monic division.
    polynomials = [p for length in range(6) for p in product((-1, 0, 1), repeat=length)]
    for p in polynomials:
        s = pread(p)
        _, rem = divide_monic(p, JOINT)
        check(rem == remainder3(s), "joint canonical remainder")
        check(s.source == peval(p, 2), "binary evaluation")
        check(pread(padd(p, pmul(JOINT, (2, -3, 1)))) == s, "kernel relation")
    for coords in product(range(-2, 3), repeat=3):
        s = State(*coords)
        check(pread(remainder3(s)) == s, "integral CRT surjectivity")
    short = [p for length in range(4) for p in product((0, 1), repeat=length)]
    for p in short:
        for q in short:
            check(multiply(pread(p), pread(q)) == pread(pmul(p, q)), "bit-polynomial product")
    z = State(0, 1, 0)
    check(multiply(z, z) == State(1, 1, 1), "square retains carry")
    check(multiply(State(0, 0, 1), State(0, 0, 1)) == State(0, 0, 1), "binary idempotent")
    check(pread((0, 0, 0, 1)) == pread((-2, -1, 3)), "cubic relation")
    square3 = multiply(bits_checked(binary(3)), bits_checked(binary(3)))
    normal9 = bits_checked(binary(9))
    check(square3 == State(2,3,1) and normal9 == State(2,2,3)
          and square3.source == normal9.source == 9,
          "binary normalization preserves source but changes golden value")

    # Every source 0..4095; all K=1..10 guards, plus canonical bit reading.
    for n in range(4096):
        bits = binary(n)
        s = bits_checked(bits)
        fib0, fib1, ga, gb = 1, 0, 0, 0
        for digit in reversed(bits):
            ga += digit*fib0; gb += digit*fib1
            fib0, fib1 = fib1, fib0+fib1
        check((s.a, s.b) == (ga, gb) and s.source == n, "independent power sum")
        for k in range(1, 11):
            modulus = 1 << k
            check((s.a+2*s.b+s.debt) % modulus == n % modulus, "source guard recovery")
    for n, m in ((3, 4), (7, 9), (13, 17)):
        s, t = bits_checked(binary(n)), bits_checked(binary(m))
        check((s.a,s.b) == (t.a,t.b) and s.source != t.source, "golden source collision")
    check(valuation(3*7+1) == 1 and valuation(3*9+1) == 2, "collision changes actual step")
    difference = padd(tuple(reversed(binary(233))), tuple(reversed(binary(223))), -1)
    check(difference == pmul(F, (0,1,0,1)), "named 223/233 golden collision")
    for t in range(128):
        n, m = 256*t+223, 256*t+233
        s, u = bits_checked(binary(n)), bits_checked(binary(m))
        check((s.a,s.b) == (u.a,u.b) and u.debt-s.debt == 10,
              "all-height named collision control")
        check(valuation(3*n+1) == 1 and valuation(3*m+1) == 2,
              "equal golden source values have opposite first drift")
    for k in range(2, 13):
        n, m = 2**(k+1)+1, 3*2**(k-1)+1
        s, t = bits_checked(binary(n)), bits_checked(binary(m))
        check((s.a,s.b) == (t.a,t.b) and (s.debt-t.debt) % 2**(k-1) == 0
              and (s.debt-t.debt) % 2**k != 0, "one fewer debt bit loses a source guard")

    clocks = []
    for m in (3, 5, 7, 11, 19, 105, 223, 233, 425):
        g, b = golden_order(m), two_order(m)
        joint = lcm(g, b)
        clocks.append((m, g, b, joint))
        state = State(1, 0, 0)
        first = None
        for k in range(1, joint+1):
            state = append(state, 0)
            state = State(state.a % m, state.b % m, state.debt % m)
            if state == State(1, 0, 0) and first is None: first = k
        check(first == joint, "joint clock independent state replay")
    s = digits_state((1,)+(0,)*80)
    check((s.a % 105, s.b % 105) == (1, 0) and s.source % 105 == 46,
          "golden period 80 is not a binary guard period")
    for k in range(1, 9):
        m = 2**k
        state = State(1, 0, 0)
        for j in range(k+1):
            check((state.source % m == 0) == (j >= k), "exact binary preperiod")
            state = append(state, 0)
        check(golden_order(m) == 3*2**(k-1), "golden clock continues after binary loss")
    mixed = []
    for m in (18, 332):
        k, period = joint_schedule(m)
        mixed.append((m, golden_order(m), k, period))
        state = State(1,0,0)
        seen = {state:0}
        for j in range(1, k+period+1):
            state = append(state,0)
            state = State(state.a % m, state.b % m, state.debt % m)
            if state in seen:
                check(seen[state] == k and j-seen[state] == period, "exact mixed-modulus schedule")
                break
            seen[state] = j
        else: raise ValueError("mixed schedule failed to return")

    # All actual source prefixes through six odd steps; split composition independently.
    for n in range(1, 1024, 2):
        word, value = [], n
        for _ in range(6):
            a = valuation(3*value+1)
            word.append(a); value = (3*value+1) >> a
            receipt = compile_receipt(n, tuple(word))
            check(verify_receipt(receipt) == (n, value), "verified polynomial receipt")
            if len(word) >= 2:
                split = len(word)//2
                left = compile_receipt(n, tuple(word[:split]))
                right = compile_receipt(replay(n, tuple(word[:split])), tuple(word[split:]))
                check(concatenate(left, right) == receipt, "ordered carry composition")

    rows = []
    for original in (105, 223, 233, 322, 332, 425):
        e = valuation(original)
        n = original >> e
        word = orbit_to(n)
        check(verify_receipt(compile_receipt(n, word), root=True) == (n, 1), "numeric first-hit ROOT control")
        state = bits_checked(binary(original))
        rows.append((original, e, (state.a,state.b,state.debt), len(word), sum(word)))
    word83 = orbit_to(83, 425)
    check(replay(83, word83) == 425, "332 odd core reaches common hub")
    check(2+sum(word83) == 33, "332 to 425 Terras clock")

    good = compile_receipt(7, (1,))
    bads = [lambda: State(True,0,0), lambda: binary(1.0),
            lambda: bits_checked((True,)), lambda: carrier((False,)),
            lambda: replay(True,()), lambda: two_order(18),
            lambda: verify_receipt(Receipt(binary(9),good.target_bits,good.word,good.carry)),
            lambda: verify_receipt(Receipt(good.source_bits,good.target_bits,good.word,(1,))),
            lambda: verify_receipt(compile_receipt(1,(2,)),root=True),
            lambda: concatenate(good,compile_receipt(13,(3,)))]
    for action in bads:
        try: action()
        except ValueError: count += 1
        else: raise ValueError("malformed/false receipt accepted")
    print("JOINT CLOCKS (modulus, golden, binary, joint):", clocks)
    print("MIXED CLOCKS (modulus, golden, preperiod, joint period):", mixed)
    print("NUMERIC ROWS (n, initial halvings, (A,B,D), odd rank, valuation cost):", rows)
    print("332->425: initial halvings2, odd word", word83, "; Terras clock33")
    print("UNIVERSES: signed polynomials length0..5; binary products length0..3; sources0..4095; guards1..10 bits")
    print("RECEIPTS: all512 odd sources<1024, prefixes1..6; all nontrivial midpoint splits;10 hostile APIs")
    print("CHECKS", count, "; PROVED identities plus FINITE-EXACT controls; no new universal coverage")


if __name__ == "__main__":
    main()

"""Finite-word signed Collatz machine and independently expandable run passports.

Pure standard library. No claim that repeated outer rewrites terminate.
Run: python 04-computation/experiments/creative_transducer_20260925.py
"""
from fractions import Fraction


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def bits(n):
    require(n >= 0, "nonnegative binary input")
    return bin(n)[2:][::-1] if n else ""


def value(word):
    return sum((b == "1") << j for j, b in enumerate(word))


def multiply_three(word, carry, audit=False):
    """Read least-significant bit first; terminal marker flushes final carry."""
    require(carry in (0, 1, 2), "invalid finite control state")
    require(all(b in "01" for b in word), "nonbinary input")
    initial = 3 * value(word) + carry
    output = []
    for j, bit in enumerate(word):
        total = 3 * int(bit) + carry
        output.append(str(total % 2))
        carry = total // 2
        require(carry in (0, 1, 2), "carry escapes three-state control")
        if audit:
            require(initial == value("".join(output))
                    + 2 ** (j + 1) * (3 * value(word[j + 1:]) + carry),
                    "local arithmetic invariant")
    # Distinct terminal rules: [0]# -> #, [1]# -> 1#, [2]# -> 01#.
    return ("".join(output) + bits(carry)).rstrip("0")


def step_word(word, sigma):
    require(sigma in (-1, 1), "unknown sign")
    require(word and word[-1] == "1" and all(b in "01" for b in word),
            "finite positive canonical input required")
    if word[0] == "0":
        return word[1:]
    return multiply_three(word[1:], 2 if sigma == 1 else 1)


def step_integer(n, sigma):
    return (3 * n + sigma) // 2 if n % 2 else n // 2


def v2(n):
    require(n > 0, "v2 needs nonzero positive integer")
    return (n & -n).bit_length() - 1


def passport(n, sigma):
    require(n > 0 and n % 2 == 1, "passport source must be positive odd")
    require(sigma in (-1, 1), "unknown sign")
    if sigma == -1 and n == 1:
        return None  # A genuine fixed point; no finite odd-run length.
    a = v2(n + sigma)
    u = (n + sigma) >> a
    z = 3 ** a * u - sigma
    b = v2(z)
    return a, u, b, z >> b


def validate(n, sigma, record):
    if sigma not in (-1, 1) or record is None:
        return False
    a, u, b, m = record
    return (n > 0 and n % 2 == 1 and a >= 1 and b >= 1
            and u > 0 and u % 2 == 1 and m > 0 and m % 2 == 1
            and n + sigma == 2 ** a * u
            and 3 ** a * u - sigma == 2 ** b * m)


def expand(n, sigma, record, seek=None):
    require(validate(n, sigma, record), "invalid run passport")
    a, _, b, m = record
    x, first_hit = n, 0 if n == seek else None
    for j in range(a + b):
        require(x % 2 == int(j < a), "claimed maximal run has wrong parity")
        x = step_integer(x, sigma)
        if x == seek and first_hit is None:
            first_hit = j + 1
    require(x == m, "expanded endpoint mismatch")
    return x, first_hit


def main():
    print("STATUS: exact transducer/passport rules; outer termination remains OPEN")
    print("SCOPE: positive finite binary words; shortcut odd step(3n+sigma)/2")
    print("TRANSITIONS (carry,input)->(output,nextcarry)")
    for carry in range(3):
        print(carry, [(b, (3 * b + carry) % 2, (3 * b + carry) // 2) for b in range(2)])
    print("END marker: carry0 emits empty;carry1 emits1;carry2 emits01 (LSB first)")
    for n in range(4097):
        for carry in range(3):
            out = multiply_three(bits(n), carry, audit=n <= 128)
            require(value(out) == 3 * n + carry, "finite word multiplication")
    for n in range(1, 16385):
        for sigma in (-1, 1):
            out = step_word(bits(n), sigma)
            require(out == bits(step_integer(n, sigma)), "signed machine arithmetic")
    print("THREE-STATE controls: all n0..4096,all3carries;local invariant n0..128")
    print("SIGNED wrapper controls: all n1..16384,both signs PASS")
    for n in range(1, 16385, 2):
        for sigma in (-1, 1):
            rec = passport(n, sigma)
            if rec is None:
                require((n, sigma) == (1, -1), "unhandled passport source")
                continue
            require(validate(n, sigma, rec), "generated passport rejected")
            m, _ = expand(n, sigma, rec)
            a, _, b, _ = rec
            if b >= a:
                require(m < n or (n, sigma, m) == (1, 1, 1), "safe-bank descent")
    print("MAXIMAL RUN controls: odd n1..16383,both signs;independent expansion PASS")
    # Compile exact maximal-run cylinders before following any trajectories.
    compiled = 0
    for sigma in (-1, 1):
        for a in range(1, 13):
            for b in range(1, 13):
                modulus = 2 ** (b + 1)
                u0 = pow(3 ** a, -1, modulus) * (sigma + 2 ** b) % modulus
                for t in (0, 1, 3):
                    u = u0 + modulus * t
                    n = 2 ** a * u - sigma
                    rec = passport(n, sigma)
                    require(rec is not None and rec[:3] == (a, u, b), "cylinder compiler")
                    expand(n, sigma, rec)
                    compiled += 1
    print("EXACT CYLINDER compiler: a,b1..12,both signs,t=0,1,3:", compiled, "PASS")
    for sigma in (1, -1):
        rows = []
        for a in range(1, 7):
            u0 = sigma * pow(3 ** a, -1, 2 ** a) % (2 ** a)
            residue, modulus = 2 ** a * u0 - sigma, 4 ** a
            rows.append((a, residue, modulus))
            for t in (0, 1, 3):
                n = residue + modulus * t
                x = n
                for j in range(2 * a):
                    require(x % 2 == int(j < a), "generated 1^a0^a word")
                    x = step_integer(x, sigma)
                require(x < n or (sigma, n) == (1, 1), "simple compiler bank")
        print("INHERITED COROLLARY safe-bank rows (a,residue,modulus),sign", sigma, rows)
    # Finite union density has a closed remainder bound; no orbit-frequency claim.
    require(sum(Fraction(1, 4 ** a) for a in range(1, 13))
            == Fraction(1, 3) * (1 - Fraction(1, 4 ** 12)), "density partial sum")
    print("SIMPLE BANK density:1/3 of all integers in odd cylinders;plus evens gives5/6")
    print("This is a simpler inherited sub-bank,not an improvement on guards-valves density")
    plus_root = passport(1, 1)
    minus5, minus17, minus41 = passport(5, -1), passport(17, -1), passport(41, -1)
    require(plus_root == (1, 1, 1, 1), "plus root cycle")
    require(minus5 == (2, 1, 1, 5), "minus5 cycle")
    require(minus17 == (4, 1, 1, 41) and minus41 == (3, 5, 3, 17), "minus17 cycle")
    print("CYCLE passports:plus1", plus_root, ";minus1 fixed;minus5", minus5)
    print("CYCLE passports:minus17", minus17, ";minus41", minus41)
    for k in range(1, 129):
        for sigma in (-1, 1):
            n = 2 ** (k + 1) - sigma
            require(step_integer(n, sigma) % (2 ** k) == n % (2 ** k), "finite-prefix rank hostile")
    print("PREFIX-RANK hostile:k1..128,both signs,n=2^(k+1)-sigma shares firstkbits with nextstate")
    # Known one-odd-step family: a generated positive root certificate control.
    for k in list(range(2, 33)) + [4096]:
        n = (4 ** k - 1) // 3
        rec = passport(n, 1)
        require(rec[0] == 1 and rec[2] == 2 * k - 1 and rec[3] == 1, "alternating-bit family")
        _, hit = expand(n, 1, rec, seek=4)
        require(hit == 2 * k - 2, "generated root4 control")
    print("GENERATED ROOT family n=(4^k-1)/3:k2..32and4096;root4 reached after2k-2steps")
    print("LARGE control:8191 input bits;passport counts(1,8191);8190-step root4 certificate")
    original = passport(3, 1)
    a, u, b, m = original
    require(not validate(3, 1, (a, u, b + 1, m)), "bad length not rejected")
    require(not validate(3, -1, original), "sign transplant not rejected")
    print("HOSTILES:wrong run length and unmodified sign-transplant rejected")
    print("PASS; no all-start generator or termination theorem inferred")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
verify_eta_dictionary_agentr.py   (agentr, 2026-06-11)

Exact verification of the Construction-A theta dictionary for Type II
(doubly-even self-dual) binary codes, with every q-convention pinned.

================================ CONVENTIONS ================================
Formal nome notation.  q is a formal variable (if you want tau: q = e^{pi i tau},
so q^2 = e^{2 pi i tau} is the "classical" modular-forms nome).

  theta_2(q) = sum_{m in Z} q^{(m+1/2)^2}            (exponents 1/4, 9/4, ...)
  theta_3(q) = sum_{m in Z} q^{m^2}
  theta_4(q) = sum_{m in Z} (-1)^m q^{m^2}
  eta(q)     = q^{1/24} prod_{n>=1} (1 - q^n)        (FORMAL eta; the classical
               Dedekind eta(tau) equals eta(q^2) in this notation, because
               eta(tau) = e^{pi i tau/12} prod (1 - e^{2 pi i n tau}))
  Delta(q)   = eta(q)^24 = q prod_{n>=1} (1-q^n)^24
  E4(q)      = 1 + 240 sum_{n>=1} sigma_3(n) q^n

Construction A.  For C a doubly-even self-dual binary code of length n,
  L_C = { v / sqrt(2) : v in Z^n, (v mod 2) in C }      (even unimodular)
  Theta_{L_C}(q) := sum_{x in L_C} q^{<x,x>}
                  = W_C( theta_3(q^2), theta_2(q^2) ),
because a coordinate with codeword bit 0 contributes sum_m q^{(2m)^2/2}
 = theta_3(q^2), and bit 1 contributes sum_m q^{(2m+1)^2/2} = theta_2(q^2).
NOTE the half-integer exponents of theta_2(q^2) = 2 sum_{m>=0} q^{(2m+1)^2/2};
they recombine into integer (indeed even) exponents because Type II weights
are divisible by 4 (sum of four odd squares is 4 mod 8).

================================ CLAIMS VERIFIED ============================
(1) W_e8hat(x,y) = x^8 + 14 x^4 y^4 + y^8  maps to  E4(q^2)
    -- the theta series of E8 (checked exactly through q^100, far past q^10).
(2) The "code discriminant"  P24(x,y) = x^4 y^4 (x^4 - y^4)^4  maps to
        16 * eta(q^2)^24 = 16 * Delta(q^2)        [formal notation]
    i.e. with q = e^{pi i tau}: image = 16 * Delta(tau), the classical
    normalized weight-12 cusp form in e^{2 pi i tau}.  NOT eta(q)^24 and NOT
    eta(q^4)^24 (both shown to mismatch).  Derivation chain, each step
    verified in exact series here:
      (a) Jacobi:        theta_3^4 = theta_2^4 + theta_4^4
                         => x^4 - y^4 |--> theta_4(q^2)^4
      (b) triple product: theta_2(q) theta_3(q) theta_4(q) = 2 eta(q^2)^3
                         (this is "theta2 theta3 theta4 = 2 eta^3" in the
                          tau-convention where all functions share one tau)
      (c) eta quotient:   theta_4(q) = eta(q)^2 / eta(q^2)
      Then with Q = q^2:
        theta_2(Q)^4 theta_3(Q)^4 theta_4(Q)^16
          = (theta_2 theta_3 theta_4)^4(Q) * theta_4(Q)^12
          = 16 eta(Q^2)^12 * eta(Q)^24 / eta(Q^2)^12 = 16 eta(Q)^24
          = 16 eta(q^2)^24.
    Contrast: x^8 y^8 (x^4-y^4)^2 |--> 256 eta(q^4)^24 = 256 Delta(2 tau),
    a level-2 form -- the exponent pattern (4,4,16) is what lands on SL2(Z).
(3) Extremal length-24 (Golay) enumerator in the basis {W_e8hat^3, P24}:
        W_g24 = W_e8hat^3 - 42 * P24          (exact integer coefficients)
    hence Theta_{A(g24)} = E4(q^2)^3 - 672 Delta(q^2)
    (48 roots: the Niemeier lattice N(A1^24); Leech would be -720 Delta).

Independent cross-checks (no theta identities used):
  * brute-force lattice-point counts for Construction A on the [8,4,4]
    extended Hamming code (all |v|^2 <= 20) and on the [24,12,8] extended
    Golay code (all |v|^2 <= 8), compared coefficient-by-coefficient;
  * code weight enumerators computed by enumerating all codewords.

Output is teed to 05-knowledge/results/verify_eta_dictionary_agentr.out
"""

import itertools
from fractions import Fraction
from math import comb

# =========================================================================
# Exact truncated power series in t = q^(1/2)  (so theta_2(q^2) has integer
# exponents).  Lists of ints, index = exponent of t, exact through t^N.
# =========================================================================
N = 200          # exact through t^200 = q^100

def zero():
    return [0] * (N + 1)

def smul(a, b):
    c = [0] * (N + 1)
    for i in range(N + 1):
        ai = a[i]
        if ai:
            lim = N - i
            for j in range(lim + 1):
                bj = b[j]
                if bj:
                    c[i + j] += ai * bj
    return c

def spow(a, k):
    r = zero(); r[0] = 1
    base = a
    while k:
        if k & 1:
            r = smul(r, base)
        k >>= 1
        if k:
            base = smul(base, base)
    return r

def ssub(a, b):
    return [x - y for x, y in zip(a, b)]

def sscale(a, c):
    return [c * x for x in a]

def lincomb(*pairs):
    out = zero()
    for c, s in pairs:
        for i in range(N + 1):
            out[i] += c * s[i]
    return out

# ---- theta blocks at argument q^2, written in t = q^(1/2) ---------------
def theta3_q2():                     # sum q^{2 m^2} = sum t^{4 m^2}
    s = zero(); m = 0
    while 4 * m * m <= N:
        s[4 * m * m] += 1 if m == 0 else 2
        m += 1
    return s

def theta2_q2():                     # 2 sum_{m>=0} q^{(2m+1)^2/2} = 2 sum t^{(2m+1)^2}
    s = zero(); m = 0
    while (2 * m + 1) ** 2 <= N:
        s[(2 * m + 1) ** 2] += 2
        m += 1
    return s

def theta4_q2():                     # sum (-1)^m q^{2 m^2} = sum (-1)^m t^{4 m^2}
    s = zero(); m = 0
    while 4 * m * m <= N:
        s[4 * m * m] += 1 if m == 0 else (2 if m % 2 == 0 else -2)
        m += 1
    return s

def eta_like(step, power, shift):
    """t^shift * prod_{n>=1} (1 - t^(step*n))^power, exact through t^N."""
    s = zero(); s[0] = 1
    n = 1
    while step * n <= N:
        f = zero()
        k = 0
        while step * n * k <= N and k <= power:
            f[step * n * k] = (-1) ** k * comb(power, k)
            k += 1
        s = smul(s, f)
        n += 1
    out = zero()
    for i in range(N + 1 - shift):
        out[i + shift] = s[i]
    return out

def E4_q2():                         # E4(q^2) = 1 + 240 sum sigma_3(n) q^{2n} -> t^{4n}
    s = zero(); s[0] = 1
    for n in range(1, N // 4 + 1):
        sig3 = sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
        s[4 * n] = 240 * sig3
    return s

def q_expansion(s, kmax, label):
    """Pretty-print a t-series whose support is even t-powers, as a q-series."""
    assert all(s[i] == 0 for i in range(1, N + 1, 2)), label + ": odd t-power present!"
    terms = []
    for k in range(0, kmax + 1):
        c = s[2 * k]
        if c:
            terms.append(f"{c}*q^{k}" if k else f"{c}")
    return label + " = " + (" + ".join(terms) if terms else "0") + " + O(q^%d)" % (kmax + 1)

# =========================================================================
# Bivariate polynomial helpers (dicts {(i,j): int} for x^i y^j)
# =========================================================================
def pmul(A, B):
    C = {}
    for (i, j), a in A.items():
        for (k, l), b in B.items():
            key = (i + k, j + l)
            C[key] = C.get(key, 0) + a * b
    return {k: v for k, v in C.items() if v}

def plin(*pairs):
    C = {}
    for c, A in pairs:
        for k, v in A.items():
            C[k] = C.get(k, 0) + c * v
    return {k: v for k, v in C.items() if v}

def poly_to_series(P, x_pows, y_pows):
    out = zero()
    for (i, j), c in P.items():
        s = smul(x_pows[i], y_pows[j])
        for idx in range(N + 1):
            out[idx] += c * s[idx]
    return out

# =========================================================================
# Codes
# =========================================================================
def hamming_e8_codewords():
    """Extended Hamming [8,4,4]: G = [I4 | A], A = J - I (rows 0111,1011,1101,1110)."""
    rows_bits = [
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 1, 1, 1, 0],
    ]
    rows = [sum(b << i for i, b in enumerate(r)) for r in rows_bits]
    words = []
    for s in range(16):
        w = 0
        for i in range(4):
            if (s >> i) & 1:
                w ^= rows[i]
        words.append(w)
    return words

def golay_g24_codewords():
    """Extended binary Golay [24,12,8] from the cyclic [23,12] QR code with
    generator g(x) = x^11+x^10+x^6+x^5+x^4+x^2+1, plus overall parity bit."""
    g = 0
    for e in (0, 2, 4, 5, 6, 10, 11):
        g |= 1 << e
    # sanity: g divides x^23 - 1 over GF(2)
    rem = 1 | (1 << 23)               # x^23 + 1
    while rem.bit_length() >= g.bit_length():
        rem ^= g << (rem.bit_length() - g.bit_length())
    assert rem == 0, "g(x) does not divide x^23+1 over GF(2)"
    rows = []
    for i in range(12):
        r23 = g << i                  # degree <= 22, fits in 23 bits
        parity = bin(r23).count("1") & 1
        rows.append(r23 | (parity << 23))
    words = []
    for s in range(4096):
        w = 0
        for i in range(12):
            if (s >> i) & 1:
                w ^= rows[i]
        words.append(w)
    return words

def weight_distribution(words, n):
    dist = {}
    for w in words:
        wt = bin(w).count("1")
        dist[wt] = dist.get(wt, 0) + 1
    return dist

# =========================================================================
# Independent lattice-point enumeration for Construction A
#   L_C = {v/sqrt(2) : v in Z^n, v mod 2 in C};  Theta coeff of q^k counts
#   v with |v|^2 = 2k.  Returns counts[s] = #{v : |v|^2 = s}, s <= cap.
# =========================================================================
def lattice_counts_blocks(words, n, block, rng, cap):
    """Meet-in-the-middle: enumerate Z^block vectors with entries in rng and
    |v|^2 <= cap, bucketed by parity mask; then convolve along each codeword."""
    blocks = n // block
    dist = {}
    for v in itertools.product(rng, repeat=block):
        s2 = sum(c * c for c in v)
        if s2 <= cap:
            mask = 0
            for i, c in enumerate(v):
                mask |= (c & 1) << i
            dist.setdefault(mask, [0] * (cap + 1))[s2] += 1
    empty = [0] * (cap + 1)
    total = [0] * (cap + 1)
    mask_block = (1 << block) - 1
    for w in words:
        acc = [0] * (cap + 1); acc[0] = 1
        for b in range(blocks):
            part = dist.get((w >> (b * block)) & mask_block, empty)
            nxt = [0] * (cap + 1)
            for s1 in range(cap + 1):
                a1 = acc[s1]
                if a1:
                    for s2 in range(cap + 1 - s1):
                        p = part[s2]
                        if p:
                            nxt[s1 + s2] += a1 * p
            acc = nxt
        for s in range(cap + 1):
            total[s] += acc[s]
    return total

# =========================================================================
# MAIN
# =========================================================================
def main():
    print(__doc__)
    print("=" * 74)
    print("STEP 0: build exact series (t = q^(1/2), exact through q^%d)" % (N // 2))
    x = theta3_q2()       # x = theta_3(q^2)
    y = theta2_q2()       # y = theta_2(q^2)
    t4 = theta4_q2()      # theta_4(q^2)
    x_pows = {0: None}; y_pows = {0: None}
    one = zero(); one[0] = 1
    for d in (4, 8, 12, 16, 20, 24):
        x_pows[d] = spow(x, d)
        y_pows[d] = spow(y, d)
    x_pows[0] = one; y_pows[0] = one

    # ---------------------------------------------------------------------
    print()
    print("=" * 74)
    print("CHECK A (preliminaries): the three classical identities, exact series")
    # (a) Jacobi at argument q^2:  theta3^4 = theta2^4 + theta4^4
    jac = ssub(x_pows[4], y_pows[4])
    th4_4 = spow(t4, 4)
    assert jac == th4_4, "Jacobi theta3^4 - theta2^4 = theta4^4 FAILED"
    print("  (a) theta_3(q^2)^4 - theta_2(q^2)^4 == theta_4(q^2)^4      OK (to q^100)")

    # (b) triple product, at nome q, in the variable u = q^(1/4):
    #     theta_2(q) theta_3(q) theta_4(q) = 2 eta(q^2)^3
    #     (u-series shapes coincide with the q^2-at-t shapes used above)
    tp = smul(smul(theta2_q2(), theta3_q2()), theta4_q2())   # u-series of LHS
    eta_q2_cubed = eta_like(8, 3, 1)                         # u * prod(1-u^{8n})^3
    assert tp == sscale(eta_q2_cubed, 2), "triple product FAILED"
    print("  (b) theta_2(q)theta_3(q)theta_4(q) == 2*eta(q^2)^3         OK (to q^50)")
    print("      [tau-convention: theta_2(tau)theta_3(tau)theta_4(tau) = 2 eta(tau)^3,")
    print("       since classical eta(tau) = formal eta(q^2) with q = e^{pi i tau}]")

    # (c) eta quotient: theta_4(q) = eta(q)^2/eta(q^2), i.e.
    #     theta_4(q) * prod(1-q^{2n}) = prod(1-q^n)^2   (prefactors q^{1/12} cancel)
    th4_plain = zero(); m = 0
    while m * m <= N:
        th4_plain[m * m] += 1 if m == 0 else (2 if m % 2 == 0 else -2)
        m += 1
    lhs = smul(th4_plain, eta_like(2, 1, 0))
    rhs = eta_like(1, 2, 0)
    assert lhs == rhs, "theta_4 eta-quotient FAILED"
    print("  (c) theta_4(q) == eta(q)^2/eta(q^2)  (as q-series)         OK (to q^200)")

    # ---------------------------------------------------------------------
    print()
    print("=" * 74)
    print("CHECK 1: W_e8hat(x,y) = x^8 + 14 x^4 y^4 + y^8  |-->  E4(q^2)")
    We8_poly = {(8, 0): 1, (4, 4): 14, (0, 8): 1}
    img_e8 = poly_to_series(We8_poly, x_pows, y_pows)
    e4 = E4_q2()
    assert img_e8 == e4, "W_e8hat image != E4(q^2)"
    print("  EXACT equality of all coefficients through q^100.")
    print("  " + q_expansion(img_e8, 20, "Theta_E8(q) = W_e8hat(theta3(q^2),theta2(q^2))"))
    print("  (= E4(q^2) with E4(q) = 1 + 240 sum sigma_3(n) q^n;  240*sigma_3(10) = %d)"
          % (240 * sum(d ** 3 for d in (1, 2, 5, 10))))

    # independent: code weight enumerator + brute-force lattice count
    words8 = hamming_e8_codewords()
    wd8 = weight_distribution(words8, 8)
    assert wd8 == {0: 1, 4: 14, 8: 1}, "e8hat weight distribution wrong: %s" % wd8
    print("  Code check: [8,4] ext. Hamming weight distribution = {0:1, 4:14, 8:1}  OK")
    counts8 = lattice_counts_blocks(words8, 8, 4, range(-4, 5), 20)
    for s in range(21):
        assert counts8[s] == img_e8[s], "E8 lattice count mismatch at |v|^2=%d" % s
    assert all(counts8[s] == 0 for s in range(21) if s % 4 != 0)
    print("  Brute-force Construction-A point count over Z^8 (all |v|^2 <= 20):")
    print("    norms q^k, k=0..10: %s   == series coefficients  OK"
          % [counts8[2 * k] for k in range(11)])

    # ---------------------------------------------------------------------
    print()
    print("=" * 74)
    print("CHECK 2: P24 = x^4 y^4 (x^4-y^4)^4  |-->  16 * eta(q^2)^24 = 16 * Delta(q^2)")
    diff4 = ssub(x_pows[4], y_pows[4])
    img_P24 = smul(smul(x_pows[4], y_pows[4]), spow(diff4, 4))
    delta_q2 = eta_like(4, 24, 4)     # eta(q^2)^24 = q^2 prod(1-q^{2n})^24 = t^4 prod(1-t^{4n})^24
    assert img_P24 == sscale(delta_q2, 16), "P24 image != 16 eta(q^2)^24"
    print("  EXACT equality of all coefficients through q^100.  Constant c = 16.")
    print("  " + q_expansion(img_P24, 14, "P24 image"))
    print("  " + q_expansion(delta_q2, 14, "eta(q^2)^24"))
    # disambiguate the argument: eta(q)^24 and eta(q^4)^24 both FAIL
    delta_q1 = eta_like(2, 24, 2)     # eta(q)^24 = q prod(1-q^n)^24 = t^2 prod(1-t^{2n})^24
    delta_q4 = eta_like(8, 24, 8)     # eta(q^4)^24 = q^4 prod(1-q^{4n})^24
    assert img_P24 != sscale(delta_q1, 16) and all(
        img_P24 != sscale(delta_q1, c) for c in (1, 2, 4, 8, 32, 256))
    assert all(img_P24 != sscale(delta_q4, c) for c in (1, 16, 256))
    print("  Disambiguation: eta(q)^24 starts at q^1, eta(q^4)^24 starts at q^4;")
    print("    the image starts at q^2  ==>  argument is q^2, period.")
    print("  [tau-language, q = e^{pi i tau}: image = 16 * Delta(tau), the classical")
    print("   normalized weight-12 cusp form q - 24q^2 + 252q^3 - ... in e^{2 pi i tau}.]")
    # contrast candidate with the wrong exponent pattern
    img_alt = smul(smul(x_pows[8], y_pows[8]), spow(diff4, 2))
    assert img_alt == sscale(delta_q4, 256), "x^8y^8(x^4-y^4)^2 != 256 eta(q^4)^24"
    print("  Contrast: x^8 y^8 (x^4-y^4)^2 |--> 256 * eta(q^4)^24 = 256 Delta(2 tau)")
    print("    (level-2 form, NOT the SL2(Z) discriminant -- exponents (4,4,16) matter).")

    # ---------------------------------------------------------------------
    print()
    print("=" * 74)
    print("CHECK 3: W_g24 in the basis {W_e8hat^3, P24}")
    Wg24_poly = {(24, 0): 1, (16, 8): 759, (12, 12): 2576, (8, 16): 759, (0, 24): 1}
    We8_cubed = pmul(pmul(We8_poly, We8_poly), We8_poly)
    xy4 = {(4, 4): 1}
    d4 = {(4, 0): 1, (0, 4): -1}
    P24_poly = pmul(xy4, pmul(pmul(d4, d4), pmul(d4, d4)))
    # solve W_g24 = a * We8^3 + b * P24 from two monomials, then verify all
    a = Fraction(Wg24_poly[(24, 0)], We8_cubed[(24, 0)])
    b = Fraction(Wg24_poly.get((20, 4), 0) - a * We8_cubed.get((20, 4), 0),
                 P24_poly[(20, 4)])
    assert a == 1 and b == -42, "unexpected coefficients (a,b)=(%s,%s)" % (a, b)
    combo = plin((int(a), We8_cubed), (int(b), P24_poly))
    assert combo == Wg24_poly, "W_g24 != W_e8^3 - 42 P24 (full monomial check)"
    print("  Exact identity of polynomials in Z[x,y]:")
    print("      W_g24 = W_e8hat^3 - 42 * P24")
    print("      i.e.  x^24+759x^16y^8+2576x^12y^12+759x^8y^16+y^24")
    print("            = (x^8+14x^4y^4+y^8)^3 - 42 * x^4y^4(x^4-y^4)^4")
    print("  W_e8hat^3 expanded: %s" % {k: v for k, v in sorted(We8_cubed.items(), reverse=True)})
    print("  P24 expanded:       %s" % {k: v for k, v in sorted(P24_poly.items(), reverse=True)})

    # theta side: Theta_{A(g24)} = E4(q^2)^3 - 672 Delta(q^2)
    img_g24 = poly_to_series(Wg24_poly, x_pows, y_pows)
    e4_cubed = smul(smul(e4, e4), e4)
    rhs = lincomb((1, e4_cubed), (-672, delta_q2))
    assert img_g24 == rhs, "Theta_{A(g24)} != E4(q^2)^3 - 672 Delta(q^2)"
    print("  Theta side (exact through q^100):")
    print("      Theta_{A(g24)} = E4(q^2)^3 - 42*16*Delta(q^2) = E4(q^2)^3 - 672 Delta(q^2)")
    print("  " + q_expansion(img_g24, 6, "Theta_{A(g24)}"))
    print("  48 roots ==> A(g24) is the Niemeier lattice N(A_1^24), NOT Leech")
    print("  (Theta_Leech = E4^3 - 720 Delta has 0 roots; extremal ENUMERATOR does not")
    print("   give the extremal LATTICE under plain Construction A).")

    # independent: Golay code + brute-force lattice count
    words24 = golay_g24_codewords()
    wd24 = weight_distribution(words24, 24)
    assert wd24 == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, "Golay weights wrong: %s" % wd24
    print("  Code check: [24,12] ext. Golay weight distribution = "
          "{0:1, 8:759, 12:2576, 16:759, 24:1}  OK")
    counts24 = lattice_counts_blocks(words24, 24, 8, range(-2, 3), 8)
    for s in range(9):
        assert counts24[s] == img_g24[s], "Golay lattice count mismatch at |v|^2=%d" % s
    print("  Brute-force Construction-A point count over Z^24 (all |v|^2 <= 8):")
    print("    q^0,q^1,q^2,q^3,q^4 coefficients: %s == series  OK"
          % [counts24[2 * k] for k in range(5)])
    print("    (48 = 2*24 doubled unit vectors; 195408 = 4*C(24,2) + 759*2^8 = 1104+194304)")

    # ---------------------------------------------------------------------
    print()
    print("=" * 74)
    print("SUMMARY -- ALL ASSERTIONS PASSED.  Exact identities (formal nome q,")
    print("x = theta_3(q^2), y = theta_2(q^2), eta(q) = q^{1/24} prod(1-q^n)):")
    print()
    print("  (1)  x^8 + 14 x^4 y^4 + y^8                  |-->  E4(q^2)")
    print("  (2)  x^4 y^4 (x^4 - y^4)^4                   |-->  16 eta(q^2)^24 = 16 Delta(q^2)")
    print("       [x^8 y^8 (x^4 - y^4)^2                  |-->  256 eta(q^4)^24, level 2]")
    print("  (3)  W_g24 = W_e8hat^3 - 42 * x^4y^4(x^4-y^4)^4")
    print("       Theta_{A(g24)} = E4(q^2)^3 - 672 Delta(q^2)   (Niemeier A_1^24)")
    print()
    print("  In tau-language (q = e^{pi i tau}, classical nome e^{2 pi i tau}):")
    print("  Theta_E8 = E4(tau),  P24-image = 16 Delta(tau),  with Delta = eta(tau)^24.")

if __name__ == "__main__":
    main()

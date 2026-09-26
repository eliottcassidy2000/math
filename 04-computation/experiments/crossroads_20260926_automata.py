"""Exact probes: affine parity sections and Thue--Morse block-swap certificates.

Run from repository root with python3 (standard library only). No floating
point computation enters a verdict. Proofs and full scope are in the note.
"""
from fractions import Fraction


def require(ok, label):
    if not ok:
        raise AssertionError(label)


def step(n, q=3, sigma=1):
    return (q*n+sigma)//2 if n & 1 else n//2


def parity(n, depth, q=3, sigma=1):
    word = 0
    ones = 0
    for j in range(depth):
        bit = n & 1
        word |= bit << j
        ones += bit
        n = step(n, q, sigma)
    return word, ones, n


def v2(n):
    require(n != 0, 'nonzero valuation argument')
    return (abs(n) & -abs(n)).bit_length()-1


def residue(f, precision):
    f = Fraction(f)
    require(f.denominator & 1, '2-adic integral denominator')
    mod = 1 << precision
    return f.numerator * pow(f.denominator, -1, mod) % mod


def tm(n):
    return n.bit_count() & 1


def rs(n):
    return (n & (n >> 1)).bit_count() & 1


def poly_eval(coeff, x):
    out = 0
    for c in reversed(coeff):
        out = out*x+c
    return out


RS_P = [1, 1, 0, -2, -1, 1, -4]
RS_Q = [1, 0, -1, 0, -1, 0, -1]


def rs_series_mod(p, b, precision):
    mod = 1 << precision
    x = p*pow(b, -1, mod) % mod
    power = 1
    out = 0
    n = 0
    while power:
        out = (out+(1-2*rs(n))*power) % mod
        power = power*x % mod
        n += 1
    return out


def rs_pade(p, b, k):
    rho = Fraction(p, b)
    n = 1 << k
    x = rho**n
    aa = sum(((1-2*rs(j))*rho**j for j in range(n//2)), Fraction())
    bb = sum(((1-2*rs(j))*rho**j for j in range(n//2, n)), Fraction())
    return (aa*poly_eval(RS_P, x)+bb*poly_eval(RS_P, -x))/poly_eval(RS_Q, x)


def even_pade_contact(d):
    """Full-rank even-denominator ansatz only; return None if singular."""
    coeff = [Fraction(1-2*rs(n)) for n in range(160)]
    rows = [[coeff[n-2*j] for j in range(1, d+1)]+[-coeff[n]]
            for n in range(2*d+1, 3*d+1)]
    for col in range(d):
        pivot = next((r for r in range(col, d) if rows[r][col]), None)
        if pivot is None:
            return None
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = rows[col][col]
        rows[col] = [x/scale for x in rows[col]]
        for r in range(d):
            if r != col:
                scale = rows[r][col]
                if scale:
                    rows[r] = [x-scale*y for x, y in zip(rows[r], rows[col])]
    q = [Fraction(1)]+[r[-1] for r in rows]
    return next(n for n in range(2*d+1, 160)
                if sum(q[j]*coeff[n-2*j] for j in range(d+1)))


def tm_product_mod(p, b, precision):
    mod = 1 << precision
    x = p * pow(b, -1, mod) % mod
    out = 1
    while x:
        out = out * (1-x) % mod
        x = x*x % mod
    return out


def pade(p, b, k):
    rho = Fraction(p, b)
    prod = Fraction(1)
    for j in range(k):
        prod *= 1-rho**(1 << j)
    x = rho**(1 << k)
    return prod * (1-2*x*x)/(1+x)


def bernstein_mod(bits, q, precision):
    mod = 1 << precision
    total = 0
    power_q = 1
    for j, bit in enumerate(bits):
        if bit:
            power_q *= q
            total -= (1 << j) * pow(power_q, -1, mod)
    return total % mod


def block_word(q, length, blocks, selector=tm):
    # q is intentionally not used: the same word can live under different maps.
    a = [1]*(length-1)+[0]
    b = [1]*(length-2)+[0, 1]
    return [bit for n in range(blocks) for bit in (b if selector(n) else a)]


def main():
    section_checks = 0
    for q, sigma in ((3, 1), (3, -1), (5, 1)):
        for depth in range(1, 11):
            outputs = set()
            for r in range(1 << depth):
                w, k, u = parity(r, depth, q, sigma)
                outputs.add(w)
                for z in (0, 1, 2, 7):
                    w2, k2, u2 = parity(r+(1 << depth)*z, depth, q, sigma)
                    require((w2, k2, u2) == (w, k, u+q**k*z), 'section law')
                    section_checks += 1
            require(len(outputs) == 1 << depth, 'parity tree bijection')
    print('section_law_checks', section_checks)
    states = {(0, 0)}
    census = []
    for depth in range(1, 19):
        states = {(k+e, (3**e*(u+3**k*b)+e)//2)
                  for k, u in states for b in (0, 1) for e in ((u+b)%2,)}
        census.append(len(states))
    print('exact_distinct_sections_depth_1_to_18', census)
    for depth in range(1, 201):
        w, k, u = parity((1 << depth)-1, depth)
        require((w, k, u) == ((1 << depth)-1, depth, 3**depth-1), 'Mersenne section')
    # Prefixes r=1 and r=4 at depth 3 share u=2 but have different slopes.
    one = parity(1, 3)
    two = parity(4, 3)
    require(one[2] == two[2] == 2 and one[1] != two[1], 'slope-loss hostile')
    require(parity(one[2]+3**one[1], 3)[0] != parity(two[2]+3**two[1], 3)[0],
            'slope affects future outputs')
    print('state_carry_only_hostile', {'depth': 3, 'sources': [1, 4],
          'states': [(one[1], one[2]), (two[1], two[2])]})

    coeff = [1-2*tm(n) for n in range(128)]
    defect = [coeff[n]+(coeff[n-1] if n else 0)
              -(1 if n == 0 else -2 if n == 2 else 0) for n in range(128)]
    require(defect[:6] == [0]*6 and defect[6] == 2, 'Pade contact order')
    require(all(c % 2 == 0 for c in defect), 'defect even coefficients')
    print('Pade_defect_first_nonzero', [(n, c) for n, c in enumerate(defect[:17]) if c])
    valuation_checks = 0
    for p, b in ((2, 1), (4, 3), (8, 9), (1024, 19683), (32, 625), (16, 81)):
        require(max(p, b) < (1 << (2*v2(p))), 'irrationality inequality')
        for k in range(9):
            expected = 6*v2(p)*(1 << k)+1
            precision = expected+7
            approximation = pade(p, b, k)
            difference = (tm_product_mod(p, b, precision)-residue(approximation, precision)) % (1 << precision)
            require(v2(difference) == expected, 'exact Pade error valuation')
            require(max(abs(approximation.numerator), approximation.denominator)
                    <= 2*max(p, b)**(3*(1 << k)-1), 'height upper bound')
            valuation_checks += 1
    print('Pade_value_and_height_checks', valuation_checks)

    bridge_checks = 0
    for q in (3, 5, 7, 9, 11):
        for length in range(3, 13):
            rho = Fraction(1 << length, q**(length-1))
            ga = -sum((Fraction(2**j, q**(j+1)) for j in range(length-1)), Fraction())
            delta = -Fraction(2**(length-2), q**(length-1))
            for blocks in (8, 32, 128):
                precision = length*blocks
                mod = 1 << precision
                bits = block_word(q, length, blocks)
                direct = bernstein_mod(bits, q, precision)
                formula = (residue((ga+delta/2)/(1-rho), precision)
                           -residue(delta/2, precision)*tm_product_mod(rho.numerator, rho.denominator, precision)) % mod
                require(direct == formula, 'block Bernstein Mahler bridge')
                # Completely independent dynamics on the least residue realizes the word.
                require(parity(direct, precision, q)[0] == sum(bit << j for j, bit in enumerate(bits)),
                        'actual parity reconstruction')
                bridge_checks += 1
    print('Bernstein_product_bridges_and_independent_dynamics', bridge_checks)
    for q, length in ((3, 10), (5, 5), (5, 10), (7, 3)):
        p, b = 1 << length, q**(length-1)
        print('block_certificate', {'q': q, 'length': length, 'ones': length-1,
             'rho': str(Fraction(p, b)), 'sufficient_inequality': max(p, b) < p*p})

    rs_coeff = [1-2*rs(n) for n in range(128)]
    rs_defect = [sum(RS_Q[j]*rs_coeff[n-j] for j in range(min(6, n)+1))
                 -(RS_P[n] if n < len(RS_P) else 0) for n in range(128)]
    require(rs_defect[:11] == [0]*11 and rs_defect[11] == -4, 'RS Pade contact')
    require(all(c % 2 == 0 for c in rs_defect), 'RS even defect')
    print('Rudin_Shapiro_Pade_defect_first_nonzero',
          [(n, c) for n, c in enumerate(rs_defect[:23]) if c])
    rs_checks = 0
    for p, b in ((2, 1), (4, 3), (8, 25), (1024, 19683)):
        require(max(p, b)**7 < (1 << (11*v2(p))), 'RS irrationality inequality')
        for k in range(1, 8):
            n = 1 << k
            expected = 11*v2(p)*n+2
            precision = expected+5
            approximation = rs_pade(p, b, k)
            difference = (rs_series_mod(p, b, precision)-residue(approximation, precision)) % (1 << precision)
            require(v2(difference) == expected, 'RS exact Pade valuation')
            require(max(abs(approximation.numerator), approximation.denominator)
                    <= 10*n*max(p, b)**(7*n-1), 'RS Pade height')
            rs_checks += 1
    print('Rudin_Shapiro_value_and_height_checks', rs_checks)
    for k in range(1, 8):
        n = 1 << k
        for j in range(5*n):
            low = j % n
            high = j//n
            sign = (1-2*rs(low))*(1-2*rs(high))
            if low >= n//2:
                sign *= (-1)**high
            require(sign == 1-2*rs(j), 'RS two-channel factorization')
    rs_bridge_checks = 0
    for q, length in ((3, 3), (3, 10), (3, 117), (3, 118), (5, 3)):
        blocks = 16
        precision = length*blocks
        mod = 1 << precision
        rho = Fraction(1 << length, q**(length-1))
        ga = -sum((Fraction(2**j, q**(j+1)) for j in range(length-1)), Fraction())
        delta = -Fraction(2**(length-2), q**(length-1))
        bits = block_word(q, length, blocks, rs)
        direct = bernstein_mod(bits, q, precision)
        formula = (residue((ga+delta/2)/(1-rho), precision)
                   -residue(delta/2, precision)*rs_series_mod(rho.numerator, rho.denominator, precision)) % mod
        require(direct == formula, 'RS Bernstein bridge')
        require(parity(direct, precision, q)[0] == sum(bit << j for j, bit in enumerate(bits)),
                'RS independent parity reconstruction')
        rs_bridge_checks += 1
        print('Rudin_Shapiro_block_certificate', {'q': q, 'length': length,
              'sufficient_inequality': max(rho.numerator, rho.denominator)**7 < 2**(11*length)})
    print('Rudin_Shapiro_bridges_and_independent_dynamics', rs_bridge_checks)
    contacts = [(d, even_pade_contact(d)) for d in range(1, 33)]
    ranked = [(Fraction(h, 2*d+1), d, h) for d, h in contacts if h is not None]
    best = max(ranked)
    require(best == (Fraction(11, 7), 3, 11), 'bounded Pade ansatz census')
    print('Rudin_Shapiro_full_rank_even_Pade_d_le_32_best', str(best))
    require(all(3**(7*(length-1)) < 2**(11*length) for length in range(3, 118)),
            'RS q3 certified block lengths')
    require(3**(7*117) >= 2**(11*118), 'RS first method miss')
    print('Rudin_Shapiro_first_q3_block_length_not_certified_by_this_form', 118)
    print('PASS: finite exact checks; infinite irrationality and nonfinite-state claims use proofs in note.')


if __name__ == '__main__':
    main()

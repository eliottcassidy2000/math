# Fixed source: no calls, arrays, recursion, or sequence tables.
def count(n):
    order = n
    phase = 0
    accumulated = 0
    while phase < 2:
        n = order
        parity = 0
        if phase == 1:
            parity = n % 2
            n = n // 2
        factorial = 1
        d = 1
        while d <= n:
            factorial = factorial * d
            d += 1
        limit = 1
        d = 1
        while d <= n:
            radix = n // d + 1
            limit = limit * radix
            d += 2
        total = 0
        code = 0
        while code < limit:
            scan = code
            weight = 0
            d = 1
            while d <= n:
                radix = n // d + 1
                multiplicity = scan % radix
                scan = scan // radix
                weight += d * multiplicity
                d += 2
            if weight == n:
                exponent = 0
                length = 0
                denominator = 1
                scan = code
                d = 1
                while d <= n:
                    radix = n // d + 1
                    multiplicity = scan % radix
                    scan = scan // radix
                    length += multiplicity
                    exponent += multiplicity * ((d - 1) // 2)
                    if multiplicity > 0:
                        exponent += d * ((multiplicity * (multiplicity - 1)) // 2)
                    t = 1
                    while t <= multiplicity:
                        denominator = denominator * d * t
                        t += 1
                    previous = code
                    e = 1
                    while e < d:
                        old_radix = n // e + 1
                        old_multiplicity = previous % old_radix
                        previous = previous // old_radix
                        x = d
                        y = e
                        while y != 0:
                            r = x % y
                            x = y
                            y = r
                        exponent += x * multiplicity * old_multiplicity
                        e += 2
                    d += 2
                summand = (factorial // denominator) * (2 ** ((1 + phase) * exponent + parity * length))
                total += summand
            code += 1
        accumulated += total // factorial
        phase += 1
    return accumulated // 2

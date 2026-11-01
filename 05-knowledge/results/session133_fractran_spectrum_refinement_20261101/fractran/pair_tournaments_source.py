# Fixed source: no calls, arrays, recursion, or sequence tables.
def count(n):
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
            denominator = 1
            scan = code
            d = 1
            while d <= n:
                radix = n // d + 1
                multiplicity = scan % radix
                scan = scan // radix
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
            summand = (factorial // denominator) * (4 ** exponent)
            total += summand
        code += 1
    return total // factorial

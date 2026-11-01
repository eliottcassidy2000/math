# Fixed source: no calls, arrays, recursion, or sequence tables.
def count(n):
    if n == 0:
        return 0
    order = n
    positive = 0
    negative = 0
    cut_code = 0
    cut_limit = 2 ** (n - 1)
    while cut_code < cut_limit:
        cut_scan = cut_code
        position = 1
        component_size = 0
        components = 0
        component_product = 1
        while position <= order:
            component_size += 1
            cut = 1
            if position < order:
                cut = cut_scan % 2
                cut_scan = cut_scan // 2
            if cut == 1:
                n = component_size
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
                        summand = (factorial // denominator) * (2 ** exponent)
                        total += summand
                    code += 1
                component_product = component_product * (total // factorial)
                component_size = 0
                components += 1
            position += 1
        if components % 2 == 1:
            positive += component_product
        else:
            negative += component_product
        cut_code += 1
    return positive - negative

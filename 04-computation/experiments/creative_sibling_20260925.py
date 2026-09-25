"""Guarded sibling normal forms and a strictly decreasing certificate grammar.

Pure Python, exact arithmetic. Never assumes universal Collatz convergence.
"""
from collections import Counter


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def valuation(n):
    return (n & -n).bit_length() - 1


def odd_step(n, sign=1):
    value = 3 * n + sign
    return value >> valuation(value)


def sibling(n, sign=1):
    return 4 * n + sign


def normal(n, sign=1):
    target = 5 if sign == 1 else 3
    depth = 0
    while n % 8 == target:
        n = (n - sign) // 4
        depth += 1
    return n, depth


def quotient(n, sign=1):
    return normal(odd_step(n, sign), sign)[0]


def ancestors(n, sign=1):
    result = [(n, 0)]
    while n % 8 == (5 if sign == 1 else 3):
        n = (n - sign) // 4
        result.append((n, len(result)))
    return result


def reverse_descent(n, sign=1):
    if n % 3 != (-sign) % 3:
        return None
    m = (2 * n - sign) // 3
    return m if 0 < m < n else None


def validate_odd_path(path, sign=1):
    check(all(n > 0 and n % 2 for n in path), "odd path domain")
    check(all(odd_step(x, sign) == y for x, y in zip(path, path[1:])),
          ("invalid odd path", path, sign))


def transport_join(n, m, depth, target_path, sign=1):
    value = m
    for _ in range(depth):
        value = sibling(value, sign)
    check(value == odd_step(n, sign), "join witness")
    check(target_path[0] == m and target_path[-1] == 1, "certificate endpoints")
    validate_odd_path(target_path, sign)
    if m == 1:
        result = [n, odd_step(n, sign), 1]
    else:
        result = [n, odd_step(n, sign)] + target_path[1:]
    validate_odd_path(result, sign)
    return result


def transport_quotient(n, target_path, sign=1):
    m, depth = normal(odd_step(n, sign), sign)
    return transport_join(n, m, depth, target_path, sign)


def transport_reverse(n, target_path, sign=1):
    m = reverse_descent(n, sign)
    check(m is not None and target_path[0] == m and target_path[-1] == 1,
          "reverse certificate endpoints")
    validate_odd_path(target_path, sign)
    check(len(target_path) > 1 and target_path[1] == n, "reverse forced prefix")
    result = target_path[1:]
    validate_odd_path(result, sign)
    return result


def ordinary_root4(odd_path):
    """Expand only an already checked finite certificate; append 1->4.

    No orbit search occurs here, and the path is cut at its first visit to4.
    """
    validate_odd_path(odd_path)
    check(odd_path[-1] == 1, "odd root")
    raw = [odd_path[0]]
    for start, end in zip(odd_path, odd_path[1:]):
        value = 3 * start + 1
        raw.append(value)
        while value % 2 == 0:
            value //= 2
            raw.append(value)
        check(value == end, "ordinary expansion")
    raw.append(4)
    raw = raw[:raw.index(4) + 1]
    check(all((3*x+1 if x % 2 else x//2) == y for x, y in zip(raw, raw[1:])),
          "root4 raw edge")
    return raw


def main():
    print("creative_sibling_20260925: guarded common-future descent")
    for sign in (1, -1):
        for n in range(1, 20001, 2):
            m, depth = normal(n, sign)
            a = valuation(3 * n + sign)
            check(depth == (a - 1) // 2, ("normal depth", sign, n))
            check(valuation(3 * m + sign) in (1, 2), "canonical valuation")
            check(odd_step(m, sign) == odd_step(n, sign), "same fibre")
            check(odd_step(sibling(n, sign), sign) == odd_step(n, sign), "sibling")
            reconstructed = m
            for _ in range(depth):
                reconstructed = sibling(reconstructed, sign)
            check(reconstructed == n, "normal reconstruction")
            q = quotient(n, sign)
            expected = n > 1 and (n % 4 == (1 if sign == 1 else 3)
                                  or n % 16 == (3 if sign == 1 else 13))
            check((q < n) == expected, ("quotient descent residues", sign, n, q))
            p = reverse_descent(n, sign)
            if p is not None:
                check(odd_step(p, sign) == n, "reverse step")
            left, right = odd_step(n, sign), n
            for _ in range(8):
                check(left == odd_step(right, sign), "clock intertwining")
                left, right = odd_step(left, sign), quotient(right, sign)
        print("sign", sign, "odd inputs<=20000:normal form,valuation,descent,8 clocks PASS")
    residual = [n for n in range(1, 49, 2)
                if n % 16 in (7, 11, 15) and n % 3 != 2]
    check(residual == [7, 15, 27, 31, 39, 43], "combined residual classes")
    print("Plus residual classes mod48:", residual, "; covered18/24 oddclasses=3/4")

    # Full strong-induction closure of ONLY the listed strict-descent rules.
    # Target labels are smaller, so no certificate depends on a later source.
    cap = 100000
    baseline, normalized, both, combined = ({1: [1]} for _ in range(4))
    ladder, ladder_ports = {1: [1]}, {1: [1]}
    reverse_uses = 0
    for n in range(3, cap + 1, 2):
        b, q, p = odd_step(n), quotient(n), reverse_descent(n)
        if b < n and b in baseline:
            baseline[n] = [n] + baseline[b]
        if q < n and q in normalized:
            normalized[n] = transport_quotient(n, normalized[q])
        both_candidates = []
        if b < n and b in both:
            both_candidates.append([n] + both[b])
        if q < n and q in both:
            both_candidates.append(transport_quotient(n, both[q]))
        if both_candidates:
            both[n] = min(both_candidates, key=len)
        candidates = []
        if b < n and b in combined:
            candidates.append([n] + combined[b])
        if q < n and q in combined:
            candidates.append(transport_quotient(n, combined[q]))
        if p is not None and p in combined:
            candidates.append(transport_reverse(n, combined[p]))
            reverse_uses += 1
        if candidates:
            combined[n] = min(candidates, key=len)
        for family, ports in ((ladder, False), (ladder_ports, True)):
            joins = [transport_join(n, m, j, family[m])
                     for m, j in ancestors(b) if m < n and m in family]
            if ports and p is not None and p in family:
                joins.append(transport_reverse(n, family[p]))
            if ports and n % 18 == 13:
                y, middle = (8*n-5)//9, (4*n-1)//3
                check(odd_step(y) == middle and odd_step(middle) == n,
                      "two-step inverse port")
                if y in family:
                    check(family[y][:3] == [y, middle, n], "two-step forced prefix")
                    joins.append(family[y][2:])
            if joins:
                family[n] = min(joins, key=len)
    check(set(baseline) <= set(both) == set(combined), "closure inclusion/equality")
    check(tuple(map(len, (baseline, normalized, both, combined))) == (252, 274, 541, 541),
          "closure census")
    check(set(both) <= set(ladder) == set(ladder_ports), "full ladder/ports inclusion")
    check(len(ladder) == 640, "full ladder count")
    for certificate in ladder_ports.values():
        ordinary_root4(certificate)
    for n in ladder:
        if odd_step(n) <= cap:
            check(odd_step(n) in ladder, "finite forward closure")
    extra = sorted(set(both) - set(baseline))
    lost = sorted(set(baseline) - set(normalized))
    check(lost[0] == 241, "canonicalization loses certificate")
    print("Exact odd source universe<=100000; baseline/Q/(F+Q)/(F+Q+reverse) counts:",
          len(baseline), len(normalized), len(both), len(combined))
    print("New sources beyond every-step forward descent:", len(extra), "first30", extra[:30])
    print("Q-only loses forward certificate:241 path", baseline[241], "Q241", quotient(241))
    print("Full sibling ladder / ladder plus inverse1,2 counts:", len(ladder), len(ladder_ports))
    print("First ladder-interior gain:483 ancestors", ancestors(odd_step(483)),
          "certificate", ladder[483])
    print("All640 generated certificates compile to ordinary root4:PASS")
    print("Reverse rule usable on already certified sources:", reverse_uses)
    print("First source outside both generated families:", next(n for n in range(1, cap, 2) if n not in combined))

    # Explicit infinite schema, tested on finite parameters without needing
    # arbitrary input orbits. Each input comes with an existing certificate.
    family_rows = []
    for a in sorted(ladder)[:20]:
        j0 = (1 - a) % 3
        for ell in range(8):
            j = j0 + 3 * ell
            m = (4 ** j * (3 * a + 1) - 1) // 3
            n = (8 * m + 1) // 3
            check(n % 16 == 3 and 3*n == 8*m+1, "family guard")
            check(odd_step(n) == sibling(m), "family first edge")
            check(odd_step(odd_step(n)) == odd_step(a), "family merge")
            check(valuation(3*n+1) == 1, "family first clock")
            check(valuation(3*odd_step(n)+1) == valuation(3*a+1)+2*j+2,
                  "family second clock")
            cert = [n, odd_step(n), 1] if a == 1 else [n, odd_step(n)] + ladder[a][1:]
            ordinary_root4(cert)
            if a == 1:
                check(n == (32 * 64 ** ell - 5) // 9, "base family")
                family_rows.append((ell, n, odd_step(n)))
    print("Family20 certifiedbases x8 parameters:160 root4 certificates PASS")
    print("Base1 two-odd-step family (ell,n,first image):", family_rows)

    # Negative controls retain known distinct basins rather than root-certify.
    check(quotient(5, -1) == 7 and quotient(7, -1) == 5, "minus2cycle")
    check(reverse_descent(7, -1) == 5 and reverse_descent(5, -1) is None,
          "minus local minimum")
    check(odd_step(13, -1) == 19 and odd_step(19, -1) == 7
          and odd_step(7, -1) == 5 and odd_step(5, -1) == 7, "minus family basin")
    print("Minus hostiles:Q has5<->7; reverse7->5 stops at5; family13->19->7->5 stays offroot1")
    # The classical long exponent-one prefix survives sibling normalization.
    for length in range(1, 21):
        n = 2 ** (length + 2) - 1
        q = n
        for j in range(1, length + 1):
            q = quotient(q)
            expected = 3 ** j * 2 ** (length + 2 - j) - 1
            check(q == expected > n, ("long rising Q prefix", length, j))
    print("Long growth hostile L=1..20:2^(L+2)-1 has L uncompacted increasingQsteps PASS")
    print("ALL CHECKS PASSED; universal termination remains OPEN")


if __name__ == "__main__":
    main()

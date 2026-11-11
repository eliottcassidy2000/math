"""Complete small-order census and boundary controls for extension gcds."""
from collections import Counter
from itertools import combinations, permutations, product
from math import factorial, gcd
from functools import reduce
import tournament_checks as tc


def canonical(rows):
    n = len(rows)
    degrees = [r.bit_count() for r in rows]
    cells = [[i for i in range(n) if degrees[i] == d] for d in sorted(set(degrees))]
    best = None
    for blocks in product(*(permutations(cell) for cell in cells)):
        order = sum((tuple(block) for block in blocks), ())
        code, bit = 0, 0
        for i in range(n):
            for j in range(i + 1, n):
                code |= ((rows[order[i]] >> order[j]) & 1) << bit
                bit += 1
        if best is None or code < best:
            best = code
    return best if best is not None else 0


def decode(n, code):
    rows, bit = [0] * n, 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> bit & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            bit += 1
    return rows


def hp(rows):
    n = len(rows)
    if not n:
        return 1
    end = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        end[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            amount = end[mask][v]
            if not amount:
                continue
            available = rows[v] & ~mask
            while available:
                bit = available & -available
                available -= bit
                end[mask | bit][bit.bit_length() - 1] += amount
    return sum(end[-1])


def dfs_hp(rows):
    n, full = len(rows), (1 << len(rows)) - 1
    def visit(mask, v):
        if mask == full:
            return 1
        choices, count = rows[v] & ~mask, 0
        while choices:
            bit = choices & -choices
            choices -= bit
            count += visit(mask | bit, bit.bit_length() - 1)
        return count
    return sum(visit(1 << v, v) for v in range(n))


def child(rows, mask):
    n = len(rows)
    return [rows[i] | (0 if mask >> i & 1 else 1 << n) for i in range(n)] + [mask]


def one_content(rows):
    h, _, _, q, starts, ends = tc.path_tables(rows)
    n = len(rows)
    def value(mask):
        result = sum(starts[i] if mask >> i & 1 else ends[i] for i in range(n))
        result += sum(q[i][j] for i in range(n) for j in range(n)
                      if not mask >> i & 1 and mask >> j & 1)
        return result
    small = [mask for mask in range(1 << n) if mask.bit_count() <= 2]
    small_gcd = reduce(gcd, (value(mask) for mask in small))
    full_gcd = reduce(gcd, (value(mask) for mask in range(1 << n)))
    tc.require(small_gcd == full_gcd and h % small_gcd == 0 and small_gcd % 2 == 1, "extension content")
    return h, small_gcd


def alternating_coefficient(rows):
    n, k = len(rows), len(rows) - 1
    vertices = []
    for i in range(n):
        vertices.append(i)
        if i < n - 1:
            vertices.append(n + i)
    edges = [tuple(sorted((a, b))) for a, b in zip(vertices, vertices[1:])]
    coefficient = 0
    for mask in range(1 << len(edges)):
        # Baseline: all free orientation bits zero, so higher label beats lower.
        extension = list(rows) + [sum(1 << i for i in range(j)) for j in range(n, n + k)]
        for index, (a, b) in enumerate(edges):
            if mask >> index & 1:
                extension[b] &= ~(1 << a)
                extension[a] |= 1 << b
        coefficient += (-1) ** (len(edges) - mask.bit_count()) * hp(extension)
    tc.require(abs(coefficient) == 2, "alternating top coefficient")
    return coefficient


def run():
    expected = [1, 1, 1, 2, 4, 12, 56, 456]
    bank, census, coefficients = [[]], [], []
    for n in range(1, 8):
        codes = {canonical(child(rows, mask)) for rows in bank for mask in range(1 << (n - 1))}
        bank = [decode(n, code) for code in sorted(codes)]
        tc.require(len(bank) == expected[n], "A000568 class control")
        mass = sum(factorial(n) // tc.automorphism_count(rows) for rows in bank)
        tc.require(mass == 2 ** (n * (n - 1) // 2), "independent labelled mass")
        distribution, nonunit = Counter(), []
        lifetimes = Counter()
        for code, rows in zip(sorted(codes), bank):
            h, content = one_content(rows)
            distribution[content] += 1
            if content > 1:
                _, _, _, q, starts, ends = tc.path_tables(rows)
                cylinder_content = reduce(gcd, starts + ends + [a for row in q for a in row])
                next_gcd, certificate = content, []
                for mask in range(1 << n):
                    child_h, child_content = one_content(child(rows, mask))
                    next_gcd = gcd(next_gcd, child_content)
                    certificate.append({"mask": mask, "child_H": child_h, "child_C1": child_content})
                    if next_gcd == 1:
                        break
                tc.require(next_gcd == 1, "C2 terminal in small census")
                nonunit.append({"code": code, "H": h, "content": content, "rows": rows,
                                "cylinder_content": cylinder_content, "C2": next_gcd,
                                "C2_certificate": certificate})
                lifetimes[2] += 1
            else:
                lifetimes[0 if h == 1 else 1] += 1
            if 2 <= n <= 4:
                coefficients.append({"n": n, "code": code, "coefficient": alternating_coefficient(rows)})
        census.append({"n": n, "ordinary_classes": len(bank), "labelled_mass": mass,
                       "C1_distribution": dict(sorted(distribution.items())),
                       "lifetime_distribution": dict(sorted(lifetimes.items())), "nonunit_classes": nonunit})
    p7 = tc.recursive_rows(7, 1)
    deletion = [row & 63 for row in p7[:6]]
    tc.require(canonical(deletion) == 8850, "first content exception is Paley7 deletion")
    base = tc.recursive_rows(3, 2)
    first = child(base, 15)
    second = child(first, 3)
    h = hp(second)
    tc.require(h == 51085 and dfs_hp(second) == h, "two-vertex escape witness")
    _, c1_base = one_content(base)
    _, c1_first = one_content(first)
    tc.require(c1_base == 3 and c1_first == 1, "content boundary")
    # A sink after the first child has count15033; gcd(15033,51085)=1.
    tc.require(gcd(hp(child(first, 0)), h) == 1, "two-vertex content one")
    return {"complete_ordinary_census": census, "alternating_coefficient_checks": coefficients,
            "escape_witness": {"rows": second, "H": h, "first_mask": 15, "second_mask": 3,
                               "C0_base": 3159, "C1_base": 3, "C2_base": 1,
                               "C1_asymmetric_child": c1_first}}

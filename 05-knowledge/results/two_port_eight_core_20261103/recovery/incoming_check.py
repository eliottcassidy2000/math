"""Independent complete short component-word acceptance controls."""
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
CHECKS = 0


def need(condition, description):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(description)


def join(a, b):
    h, p, f, g, j, e, s = a
    H, P, F, G, J, E, S = b
    return h * H, p * P, f * H + p * F, g * P + h * G, \
        j * H + h * J + e * S + g * F, h * E, s * H


def atom(radix, digit):
    if radix == 2:
        return 1, 2, digit, digit, 0, digit, digit
    return 3, 12, 3 * digit, 3 * digit, digit * (digit - 1), digit, digit


def words(f):
    cmax, remaining = 0, f
    while remaining % 3 == 0:
        cmax += 1
        remaining //= 3
    accepted = []
    tested = 0
    for c in range(cmax + 1):
        q = f // (3 ** c)
        # Nonzero final digit has weight at least 2^(length-1),
        # so every valid word has length <= bit_length(q).
        for length in range(max(1, c), q.bit_length() + 1):
            for triangle_positions in combinations(range(length), c):
                radices = [4 if i in triangle_positions else 2 for i in range(length)]
                for digits in product(*(range(r) for r in radices)):
                    if not digits[-1]:
                        continue
                    tested += 1
                    value, scale = 0, 1
                    for radix, digit in zip(radices, digits):
                        value += scale * digit
                        scale *= radix
                    if value != q:
                        continue
                    signature = (1, 1, 0, 0, 0, 0, 0)
                    for radix, digit in zip(radices, digits):
                        signature = join(signature, atom(radix, digit))
                    h, p, actual_f, g, j, e, s = signature
                    need(h == 3 ** c and actual_f == f and p == h * scale,
                         "atom law matches literal mixed-radix representation")
                    accepted.append({"word": [[r, d] for r, d in zip(radices, digits)],
                                     "signature": list(signature), "intercept": h - j})
    return accepted, tested


def j_transitive(mask):
    selected = [i for i in range(mask.bit_length()) if mask >> i & 1]
    return sum(2 ** (k - i - 1) + int(k == i + 1)
               for i, k in combinations(selected, 2))


def main():
    source_paths = [HERE / "incoming_check.py", HERE / "INCOMING.md",
                    RESULTS / "multiport_defect_grammar_20261103/COMPONENT_LANGUAGE.md",
                    RESULTS / "multiport_defect_grammar_20261103/TRANSITIVE_CHART.md",
                    RESULTS / "multiport_defect_grammar_20261103/transport/TWO_DIGIT_ENTRY.md",
                    RESULTS / "recursive_seed_switching_20261103/transport/FOLLOWUP.md"]
    pins = {str(p.relative_to(RESULTS)): sha256(p.read_bytes()).hexdigest() for p in source_paths}
    need(594839 + 16 * 827538 == 13835447, "named seed repair count")
    need(37 * 179 * 2089 == 13835447, "bank obstruction is a composite count")
    need(j_transitive(7) == 6 and j_transitive(21) == 12, "first transitive chart joint")
    need(3 * j_transitive(7) - j_transitive(21) - 2 == 4, "first chart digit")
    need(3 * j_transitive(21) - j_transitive(63) - 2 == -28, "output chart hostile")
    for clock in range(10):
        source = -5 + 7 * 2 ** clock
        need(-13 + 21 * 2 ** clock == 3 * source + 2, "common-source +2 law")
        need(-11 + 21 * 2 ** clock == 3 * source + 4, "common-source +4 law")
    expected_intercepts = {
        21: [-15, -13, -11],
        63: [-61, -43, -39, -37, -33, -31, -19],
        189: [-183, -181, -139, -137, -129, -91, -69, -57, -49],
    }
    expected_words = {21: 4, 63: 7, 189: 13}
    expected_digits = {21: ["2", "4"], 63: ["22", "42"], 189: ["422"]}
    rows = []
    for depth, f in enumerate((21, 63, 189), 1):
        accepted, tested = words(f)
        intercepts = sorted({w["intercept"] for w in accepted})
        need(len(accepted) == expected_words[f], "complete atom word count")
        need(intercepts == expected_intercepts[f], "complete intercept spectrum")
        admitted_digits = []
        for digits in product((2, 4), repeat=depth):
            intercept = -5
            for d in digits:
                intercept = 3 * intercept + d
            if intercept in intercepts:
                admitted_digits.append("".join(map(str, digits)))
        need(admitted_digits == expected_digits[f], "exact fixed-clock digit acceptance")
        rows.append({"f": f, "bounded_words_tested": tested, "accepted": accepted,
                     "intercepts": intercepts, "admitted_digits": admitted_digits})
    need(pins == {str(p.relative_to(RESULTS)): sha256(p.read_bytes()).hexdigest()
                  for p in source_paths}, "source pins remain stable")
    print(json.dumps({"status": "FINITE-EXACT FIXED COMPONENT GRAMMAR AND COMMON CLOCK",
                      "checks": CHECKS, "rows": rows, "source_sha256": pins},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

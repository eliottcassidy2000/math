"""Explicit checks remain enabled under python -O; no assert statements."""

import itertools
import json
import random
import sys
from pathlib import Path

import port_numeration as p


checks = 0


def check(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def direct(coefficients):
    return sum(d * w for d, w in zip(coefficients, p.weights(len(coefficients))))


def raises(exception, fn, *args):
    try:
        fn(*args)
    except exception:
        check(True, "expected exception")
        return
    raise RuntimeError(f"{fn.__name__}{args} did not reject malformed input")


def run():
    ws = p.weights(1002)
    check(ws[:8] == [1, 3, 8, 21, 55, 144, 377, 987], "seed weights")
    for i in range(1000):
        a, b, c = ws[i:i + 3]
        check(a == p.fibonacci(2 * i + 2), "even Fibonacci identity")
        check(a * b + 1 == (b - a) ** 2, "adjacent D(1) pair")
        check(a * c + 1 == b * b, "outer D(1) pair")
        check(b * c + 1 == (c - b) ** 2, "second adjacent D(1) pair")
        check(c == a + b + 2 * (b - a), "regular extension")
        check(a % 2 == [1, 1, 0][i % 3], "weight parity cycle")
        check(p.encode(a) == "1" + "0" * i, "large power boundary")
        if i:
            check(p.encode(a - 1) == "2" + "1" * (i - 1), "large maximal word")

    enumerated = 0
    legal_total = 0
    for length in range(11):
        legal, locked_legal = [], []
        for chars in itertools.product("012", repeat=length):
            word = "".join(chars)
            enumerated += 1
            val = direct([int(x) for x in reversed(word)])
            check(p.value(word) == val, "Horner vs independent place values")
            if p.valid(word):
                legal.append(val)
                check(p.parity(word) == val % 2, "automaton parity")
            if p.valid(word, locked=True):
                locked_legal.append(val)
        previous = ws[length - 1] if length else 0
        check(legal == list(range(ws[length])), "free-state interval and lex rank")
        check(locked_legal == list(range(ws[length] - previous)),
              "locked-state interval and lex rank")
        legal_total += len(legal)

    for n in range(100001):
        word = p.encode(n)
        check(p.valid(word), "greedy legality")
        check(p.value(word) == n, "greedy round trip")
        check(p.successor(word) == p.encode(n + 1), "lex successor")
        check(p.parity(word) == n % 2, "successor parity")
        z = p.zeckendorf(n)
        check("11" not in z, "Zeckendorf nonadjacency")
        check(p.zeckendorf_value(z) == n, "Zeckendorf round trip")
        zl = [int(d) for d in reversed(z)] + [0, 0]
        cs = [zl[2 * j] + (zl[2 * j - 1] if j else 0) - zl[2 * j + 1]
              for j in range((len(z) + 2) // 2)]
        check(all(c in (-1, 0, 1) for c in cs), "signed even coefficients")
        check(direct(cs) == n, "signed Zeckendorf to even-weight identity")
        x, y = p.split_into_zeckendorf(word)
        check("11" not in x and "11" not in y, "split words are Zeckendorf")
        check(p.zeckendorf_value(x) + p.zeckendorf_value(y) == n,
              "split bridge preserves value")

    normalized = 0
    rewrite_steps = 0
    for length in range(9):
        for ds in itertools.product(range(3), repeat=length):
            word, trace = p.normalize(ds, True)
            check(word == p.encode(direct(ds)), "local normal form vs greedy")
            check(all(t["after"] < t["before"] for t in trace), "strict potential")
            check(len(trace) <= sum(ds), "termination bound")
            normalized += 1
            rewrite_steps += len(trace)

    rng = random.Random(125_138)
    for _ in range(1000):
        ds = [rng.randrange(10) for _ in range(rng.randrange(1, 65))]
        word, trace = p.normalize(ds, True)
        check(word == p.encode(direct(ds)), "large-coefficient normalization")
        check(all(t["after"] < t["before"] for t in trace), "random strict potential")
        normalized += 1
        rewrite_steps += len(trace)
        n = rng.getrandbits(512)
        cg = p.encode(n)
        check(p.valid(cg) and p.value(cg) == n, "512-bit round trip")
        check(p.successor(cg) == p.encode(n + 1), "512-bit successor")

    for low in range(30):
        for high in range(low + 1, 40):
            left = 2 * ws[low] + 2 * ws[high] + sum(ws[low + 1:high])
            right = ws[high + 1] + (ws[low - 1] if low else 0)
            check(left == right, "forbidden-block carry identity")

    # Independent unbounded-coin dynamic program: optimize summand count,
    # and count unordered coefficient arrays attaining that minimum.
    limit = 100000
    best, ways = [limit + 1] * (limit + 1), [0] * (limit + 1)
    best[0], ways[0] = 0, 1
    for coin in ws:
        if coin > limit:
            break
        for n in range(coin, limit + 1):
            score = best[n - coin] + 1
            if score < best[n]:
                best[n], ways[n] = score, ways[n - coin]
            elif score == best[n]:
                ways[n] += ways[n - coin]
    for n in range(limit + 1):
        check(ways[n] == 1, "unique minimum-summand coefficient array")
        check(best[n] == sum(map(int, p.encode(n))), "canonical summand minimum")

    states = p.parity_states()
    check(len(states) == 8, "eight live legality/parity states reachable")
    for (locked, n, shifted), word in states.items():
        check(p.valid(word) and p.parity(word) == n, "state witness")
        check(p.value(word + "0") % 2 == shifted, "shifted register witness")

    for arg in (-1, -10):
        raises(ValueError, p.encode, arg)
    for arg in (True, 1.5, "1"):
        raises(TypeError, p.encode, arg)
    for word in ("22", "212", "21112", "3", "x", "10212"):
        check(not p.valid(word), "invalid word rejected")
        raises(ValueError, p.successor, word)
        raises(ValueError, p.parity, word)
    for ds in ([-1], [1, -2]):
        raises(ValueError, p.normalize, ds)
    raises(TypeError, p.normalize, [True])
    check(p.normalize([]) == "0", "empty normalization")
    check(p.successor("") == "1", "empty successor")
    check(p.value("") == 0, "empty value")
    check(p.successor("00021") == "100", "leading-zero successor")
    check(1 + 3 == 4 and 1 + 8 == 9 and 3 + 8 == 11,
          "D(1) does not assert a square-sum triangle")
    report = {
        "status": "PASS", "optimized_python": not __debug__, "explicit_checks": checks,
        "all_ternary_words_through_length": 10, "ternary_words_enumerated": enumerated,
        "legal_words_including_leading_zeros_enumerated": legal_total,
        "consecutive_integer_roundtrips": 100001,
        "D1_triples_checked": 1000,
        "local_normalizations": normalized, "local_rewrite_steps": rewrite_steps,
        "random_512_bit_cases": 1000,
        "live_legality_parity_states": len(states),
        "independent_minimum_summand_DP_cases": limit + 1,
        "proof_boundary": "finite tests supplement universal proofs; no Hamiltonicity claim",
    }
    print(json.dumps(report, sort_keys=True, indent=2))
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    run()

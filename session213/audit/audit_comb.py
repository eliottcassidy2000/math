"""Independent comb formula and density-certificate checks."""
import itertools
import json
from independent_counts import apex, bits, brute_paths, check, count_paths


def comb(m, q, tau, forward):
    a = apex([0, tau] + bits(q, m))
    b = len(a) - 1
    for row in a:
        row.append(True)
    a.append([False] * len(a[0]))
    a[1][-1] = False
    a[-1][1] = True
    if forward:
        a[b][-1] = False
        a[-1][b] = True
    return a


def quadratic(m, q, tau, forward):
    e = bits(q, m)
    value = (2 ** m + 1 + forward * (4 * 2 ** m - 2)
             + tau * (2 ** (m + 1) + 2) - tau * forward * 2 ** m)
    for i in range(m):
        u = 2 ** (m - 1 - i) + 2 ** i
        d = 2 ** (m + 1 - i) * 3 ** i + 2 ** (i + 2) - 2 * (i == m - 1)
        linear = d - forward * (u - (i == 0) + (i == m - 1)) - tau * (u + (i == 0) - (i == m - 1))
        value += e[i] * linear
        for j in range(i + 1, m):
            gap = j - i
            weight = 2 ** (m - gap) * (3 ** (gap - 1) + 1) + 2 ** (gap - 1) + (gap == 1)
            value -= e[i] * e[j] * weight
    return value


def path_cover_two(a):
    n = len(a)
    value = 0
    for mask in range(1 << n):
        left = [i for i in range(n) if mask >> i & 1]
        right = [i for i in range(n) if not mask >> i & 1]
        value += count_paths([[a[i][j] for j in left] for i in left]) * count_paths([[a[i][j] for j in right] for i in right])
    return value


def context_readout(m, q, tau):
    e = bits(q, m)
    cover = path_cover_two(apex(e))
    reverse = 1 + sum((1 - e[i]) * 2 ** (m - 1 - i) for i in range(m))
    mixed = sum(e[i] * 3 ** i * 2 ** (m - 1 - i) for i in range(m))
    return (count_paths(apex([0, tau] + e)) - (2 - tau) * e[-1]
            + (1 + tau) * 2 ** m + 2 * mixed + cover - (2 - tau) * reverse)


def valid_cycle(a, word):
    return len(set(word)) == len(word) and all(a[word[i]][word[(i + 1) % len(word)]] for i in range(len(word)))


def main():
    checks = {"quadratic_DP": 0, "cycle_DP": 0, "permutation": 0, "context_DP": 0,
              "odd_cycle_pair_certificates": 0, "density_lower_bounds": 0}
    spectra = {}
    for m in range(1, 7):
        values = set()
        for q in range(1 << m):
            for tau in range(2):
                for forward in range(2):
                    a = comb(m, q, tau, forward)
                    actual = count_paths(a)
                    expected = quadratic(m, q, tau, forward)
                    check(actual == expected, ("quadratic formula", m, q, tau, forward, actual, expected))
                    checks["quadratic_DP"] += 1
                    check(actual >= 1 + 2 ** m, ("one-cycle lower bound", m, q, tau, forward, actual))
                    c = count_paths(a, start=0, close=True)
                    check(c == q + tau + forward, ("cycle count", m, q, tau, forward, c))
                    check(actual >= len(a) * c, ("cycle-cut lower bound", m, q, tau, forward))
                    checks["cycle_DP"] += 1
                    if m <= 3:
                        check(brute_paths(a) == actual, ("permutation", m, q, tau, forward))
                        checks["permutation"] += 1
                    if not forward and m <= 4:
                        context = context_readout(m, q, tau)
                        check(context == actual, ("context readout", m, q, tau, context, actual))
                        checks["context_DP"] += 1
                    if q:
                        j = q.bit_length() - 1
                        t = m - j - 1
                        four_p = 2 ** t * (3 ** j + 1) - (1 + (-1) ** j) * (t == 0)
                        check(actual >= 1 + 2 ** m + four_p, ("pair lower bound", m, q, tau, forward, actual))
                        check(actual >= 2 ** t * 3 ** j, ("density lower bound", m, q, tau, forward, actual))
                        checks["density_lower_bounds"] += 1
                    values.add(actual)
        spectra[m] = {"distinct": len(values), "minimum": min(values), "maximum": max(values)}
        for j in range(m):
            a = comb(m, 1 << j, 0, False)
            t = m - j - 1
            certificate_count = 0
            for mask_b in range(1 << j):
                subset_b = [i for i in range(j) if mask_b >> i & 1]
                if len(subset_b) % 2:
                    continue
                remaining = [i for i in range(m) if i != j and i not in subset_b]
                for mask_a in range(1 << len(remaining)):
                    subset_a = [i for k, i in enumerate(remaining) if mask_a >> k & 1]
                    if len(subset_a) % 2 == 0:
                        continue
                    cycle_a = [m + 3, 1] + [i + 2 for i in subset_a]
                    cycle_b = [m + 2, 0] + [i + 2 for i in subset_b] + [j + 2]
                    check(valid_cycle(a, cycle_a) and valid_cycle(a, cycle_b), ("invalid cycle certificate", m, j, cycle_a, cycle_b))
                    check(not set(cycle_a).intersection(cycle_b), "cycles not disjoint")
                    certificate_count += 1
            four_p = 2 ** t * (3 ** j + 1) - (1 + (-1) ** j) * (t == 0)
            check(4 * certificate_count == four_p, ("pair count", m, j, certificate_count, four_p))
            checks["odd_cycle_pair_certificates"] += certificate_count
    print(json.dumps({"status": "PASS", "checks": checks, "spectra": spectra,
                      "scope": "Independent finite controls, not a substitute for the all-order proofs."}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

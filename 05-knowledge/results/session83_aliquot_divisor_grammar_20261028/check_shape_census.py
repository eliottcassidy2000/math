#!/usr/bin/env python3
"""Complete target-window aliquot shape census, independently source-sieved.

Standard-library only.  The inverse DFS is adapted from Session82's audited
check_aliquot_fibers.py; the direct divisor-addition sieve shares no pruning.
The report is stdout, with no raw fiber dump or generated JSON.
"""
from array import array
from collections import Counter, defaultdict
import argparse
import hashlib
import math


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime_table(limit):
    prime = bytearray(b"\x01") * (limit + 1)
    prime[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if prime[p]:
            prime[p*p:limit+1:p] = b"\x00" * ((limit-p*p)//p + 1)
    return prime, [p for p in range(2, limit + 1) if prime[p]]


def inverse_fiber(t, prime, primes):
    """Increasing-prime-prefix DFS; all terminal sources have s(n)=t."""
    need(t >= 2, "finite fibers require t >= 2")
    out = set()

    def visit(m, sig, start):
        sm = sig - m
        if sm:
            num = t - sig
            if num > 0 and num % sm == 0:
                q = num // sm
                if q < len(prime) and prime[q] and (start == 0 or q > primes[start-1]):
                    out.add(m*q)
        for i in range(start, len(primes)):
            p = primes[i]
            terminal2 = sm*p*p + sig*(p+1)
            continuation1 = sm*p*p + 2*sig*p + sig
            if min(terminal2, continuation1) > t:
                break
            pp, ps, exponent = p, 1+p, 1
            while True:
                n, sn = m*pp, sig*ps
                value = sn-n
                if value > t:
                    break
                if exponent >= 2 and value == t:
                    out.add(n)
                if i+1 < len(primes) and value*primes[i+1]+sn <= t:
                    visit(n, sn, i+1)
                pp *= p
                ps += pp
                exponent += 1

    visit(1, 1, 0)
    return sorted(out)


def source_sieve_fibers(limit):
    """Independent oracle over the proved full source universe n<=(V-1)^2."""
    bound = (limit-1)**2
    # Each stored s(n) is below n^2; use unsigned 64-bit cells explicitly.
    need(bound*bound < 2**64, "requested source sieve exceeds 64-bit safety guard")
    aliquot = array("Q", [0]) * (bound+1)
    for d in range(1, bound//2+1):
        for n in range(2*d, bound+1, d):
            aliquot[n] += d
    fibers = [[] for _ in range(limit+1)]
    for n, t in enumerate(aliquot):
        if 2 <= t <= limit:
            fibers[t].append(n)
    return fibers


def factor(n, primes):
    out = []
    for p in primes:
        if p*p > n:
            break
        if n % p == 0:
            a = 0
            while n % p == 0:
                a += 1
                n //= p
            out.append((p, a))
    if n > 1:
        out.append((n, 1))
    return out


def kind(fs):
    if len(fs) == 1:
        return "I"
    if len(fs) == 2 and all(a == 1 for _, a in fs):
        return "II"
    return "III"


def square(n):
    return math.isqrt(n)**2 == n


def canonical_digest(fibers):
    h = hashlib.sha256()
    for t in range(2, len(fibers)):
        h.update((str(t)+":"+",".join(map(str, fibers[t]))+"\n").encode("ascii"))
    return h.hexdigest()


def census(limit):
    prime, primes = prime_table(limit)
    fibers = [[] for _ in range(limit+1)]
    for t in range(2, limit+1):
        fibers[t] = inverse_fiber(t, prime, primes)
    oracle = source_sieve_fibers(limit)
    for t in range(2, limit+1):
        need(fibers[t] == oracle[t], f"independent full-source mismatch at target {t}")

    classes = {}
    incidences, masks, singletons, shapes, parity = (Counter() for _ in range(5))
    square_edges = defaultdict(list)
    for t in range(2, limit+1):
        cs = defaultdict(list)
        for n in fibers[t]:
            fs = factor(n, primes)
            exponents = tuple(sorted((a for _, a in fs), reverse=True))
            c = kind(fs)
            tau = math.prod(a+1 for _, a in fs)
            sig = math.prod((p**(a+1)-1)//(p-1) for p, a in fs)
            need(sig-n == t, f"factorized re-encoding failed: {n}")
            need(n <= (t-1)**2, "universal source bound")
            if not (len(fs) == 1 and fs[0][1] == 2):
                need(4*n < (t-1)**2, "non-prime-square strict source bound")
            need((tau-1 == 3) == (c == "II" or exponents == (3,)), "exact three-term rule")
            need(c != "III" or tau-1 >= 5, "class III minimum term count")
            need((sig % 2 == 1) == (square(n) or (n % 2 == 0 and square(n//2))), "sigma parity")
            if c == "II":
                p, q = fs[0][0], fs[1][0]
                need(t == 1+p+q, "distinct-prime pair formula")
                need((t-1)**2-4*n == (q-p)**2, "prime-gap square deficit")
                if t % 2 == 0:
                    need(p == 2 and q == t-3, "even-target semiprime rigidity")
                else:
                    need(p % 2 == q % 2 == 1, "odd-target semiprime parity")
            if square(sig):
                need(n != t, "square-sigma edge must not be a loop")
                square_edges[c].append((t, n))
                if c == "II" and t % 2 == 0:
                    k = math.isqrt((t-2)//12)
                    need(t == 12*k*k+2 and n == 24*k*k-2, "even square-sigma semiprime formula")
                if c == "II" and t % 2 == 1:
                    h, r = (t-1)//2, (fs[1][0]-fs[0][0])//2
                    need(r*r+sig == (h+1)**2, "Pythagorean prime-pair guard")
                    need(any(p % 4 == 1 for p, _ in factor(h+1, primes)), "hypotenuse prime obstruction")
            cs[c].append(n)
            incidences[c] += 1
            shapes[c, exponents, tau-1, len(fs)] += 1
            parity[c, t % 2, n % 2] += 1
        classes[t] = dict(cs)
        masks[tuple(sorted(cs))] += 1
        if len(fibers[t]) == 1:
            singletons[next(iter(cs))] += 1

    controls = {2: [], 3: [4], 5: [], 6: [6, 25], 7: [8], 12: [121], 13: [27, 35],
                14: [22, 169],
                16: [12, 26], 28: [28], 31: [32, 125, 161, 209, 221],
                40: [44, 74, 81], 64: [56, 76, 122]}
    for t, expected in controls.items():
        if t <= limit:
            need(fibers[t] == expected, f"hostile/positive control F_{t}")

    predicates = [
        ("all three classes", lambda cs: len(cs) == 3),
        ("only III, nonempty", lambda cs: set(cs) == {"III"}),
        ("I and III present", lambda cs: "I" in cs and "III" in cs),
        ("II and III present", lambda cs: "II" in cs and "III" in cs),
        ("multiple I sources", lambda cs: len(cs.get("I", [])) > 1),
        ("multiple II sources", lambda cs: len(cs.get("II", [])) > 1),
        ("multiple III sources", lambda cs: len(cs.get("III", [])) > 1),
    ]
    print(f"# Reproduced aliquot shape census: targets 2..{limit}\n")
    print(f"Independent divisor-addition sieve through n={(limit-1)**2}: PASS.\n")
    print(f"Full sorted fiber SHA-256: `{canonical_digest(fibers)}`\n")
    print("## Source incidences and singleton fibers\n")
    print("| Class | Source incidences | Singleton fibers | Square-sigma incidences |")
    print("|---|---:|---:|---:|")
    for c in ("I", "II", "III"):
        print(f"| {c} | {incidences[c]} | {singletons[c]} | {len(square_edges[c])} |")
    print(f"\nTotal source incidences: {sum(incidences.values())}.\n")
    print("## Per-class multiplicity\n\n| Class | Absent | Exactly one | At least two | Maximum | Maximizing targets |\n|---|---:|---:|---:|---:|---|")
    for c in ("I", "II", "III"):
        histogram = Counter(len(cs.get(c, [])) for cs in classes.values())
        maximum = max(histogram)
        maximizing = [t for t, cs in classes.items() if len(cs.get(c, [])) == maximum]
        print(f"| {c} | {histogram[0]} | {histogram[1]} | {sum(v for k, v in histogram.items() if k >= 2)} | {maximum} | {maximizing} |")
    print("## Class presence masks\n\n| Present classes | Target count |\n|---|---:|")
    for mask in [(), ("I",), ("II",), ("III",), ("I", "II"), ("I", "III"), ("II", "III"), ("I", "II", "III")]:
        print(f"| {','.join(mask) or 'none'} | {masks[mask]} |")
    print("\n## First witnesses\n")
    for label, predicate in predicates:
        t = next((t for t in range(2, limit+1) if predicate(classes[t])), None)
        print(f"- {label}: " + (f"t={t}, classes={classes[t]}." if t is not None else "not present in this window."))
    max_size = max(map(len, fibers[2:]))
    print(f"- Maximum full fiber size: {max_size}; targets={[t for t in range(2, limit+1) if len(fibers[t]) == max_size]}.")
    print(f"- Odd targets with only III: {[t for t in range(3, limit+1, 2) if set(classes[t]) == {'III'}]}.")
    print("\n## Exact source shape stratification\n")
    print("Exponents are an unordered decreasing tuple; the primes are distinct.\n")
    print("| Class | Exponent shape | Proper-divisor terms | omega | Incidences |")
    print("|---|---|---:|---:|---:|")
    for (c, shape, terms, omega), count in sorted(shapes.items(), key=lambda x: (x[0][0], x[0][2], x[0][3], x[0][1])):
        print(f"| {c} | {shape} | {terms} | {omega} | {count} |")
    print("\n## Parity stratification\n\n| Class | Target parity | Source parity | Incidences |\n|---|---|---|---:|")
    for (c, tp, np), count in sorted(parity.items()):
        print(f"| {c} | {'odd' if tp else 'even'} | {'odd' if np else 'even'} | {count} |")
    print("\n## Earliest square-sigma source-target pairs\n")
    for c in ("I", "II", "III"):
        print(f"- {c}: {square_edges[c][:15]} (pairs are (target, source)).")
    if limit == 2000:
        need(incidences == Counter(I=332, II=26647, III=1791), "frozen incidence totals")
        need(sum(masks.values()) == 1999 and masks[()] == 196, "frozen target totals")
        need(singletons == Counter(I=148, II=87, III=170), "frozen singleton totals")
        need({c: len(square_edges[c]) for c in ("I", "II", "III")} == {"I": 2, "II": 283, "III": 74}, "frozen square-edge totals")
        need(max_size == 92 and [t for t in range(2, limit+1) if len(fibers[t]) == max_size] == [1891], "frozen maximum")
        need(square_edges["I"] == [(40, 81), (57, 343)], "frozen prime-power square edges")
        need(square_edges["II"][0] == (14, 22) and square_edges["III"][0] == (74, 70), "minimal square-edge targets")
        need(canonical_digest(fibers) == "01c598ee1831b6c1b255950803fd0271e66934fc408a7cec6c752663d83adfe9", "frozen full-fiber digest")
    print("\nPASS: complete fibers, independent oracle, factor re-encoding, shape, parity, gap, square-edge, and hostile controls.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=2000)
    args = parser.parse_args()
    need(args.limit >= 2, "limit must be at least 2")
    census(args.limit)

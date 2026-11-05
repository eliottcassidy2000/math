#!/usr/bin/env python3
"""Self-contained physical controls for unit-determinant synthesis."""
from functools import lru_cache
from itertools import permutations
import json


def require(value, message):
    if not value:
        raise RuntimeError(message)


def validate(rows):
    for i, row in enumerate(rows):
        require(not (row >> i) & 1, "loop")
        for j in range(i):
            require(((row >> j) & 1) + ((rows[j] >> i) & 1) == 1,
                    "missing or doubled arc")


def table(rows):
    n = len(rows)
    size = 1 << n
    total = [0] * size
    total[0] = 1
    end = [[0] * n for _ in range(size)]
    for mask in range(1, size):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            old = mask ^ (1 << v)
            end[mask][v] = (1 if not old else
                           sum(end[old][u] for u in range(n)
                               if ((old >> u) & 1) and ((rows[u] >> v) & 1)))
        total[mask] = sum(end[mask])
    return total, end


def profile(rows, root):
    validate(rows)
    n = len(rows)
    full = (1 << n) - 1
    total, end = table(rows)
    converse = [full ^ (1 << i) ^ rows[i] for i in range(n)]
    rev, start = table(converse)
    require(rev == total, "converse counts")
    R = sum(start[m][root] * total[full ^ m]
            for m in range(1, full + 1) if (m >> root) & 1)
    L = sum(end[m][root] * total[full ^ m]
            for m in range(1, full + 1) if (m >> root) & 1)
    return total[full], end[full][root], R, L, total[full ^ (1 << root)]


def literal_profile(rows, root):
    h = e = R = L = 0
    def legal(p):
        return all((rows[a] >> b) & 1 for a, b in zip(p, p[1:]))
    for p in permutations(range(len(rows))):
        j = p.index(root)
        if legal(p):
            h += 1
            e += p[-1] == root
        R += legal(p[:j]) and legal(p[j:])
        L += legal(p[:j+1]) and legal(p[j+1:])
    rest = [x for x in range(len(rows)) if x != root]
    d = sum(legal(p) for p in permutations(rest))
    return h, e, R, L, d


def tournament(n, code):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> k) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            k += 1
    return rows


def inflate(rows, root, m):
    require(m >= 1, "positive clone multiplicity")
    old = [x for x in range(len(rows)) if x != root]
    clones = list(range(len(old), len(old) + m))
    out = [0] * (len(old) + m)
    for i, a in enumerate(old):
        for j, b in enumerate(old):
            if (rows[a] >> b) & 1:
                out[i] |= 1 << j
        for c in clones:
            if (rows[a] >> root) & 1:
                out[i] |= 1 << c
            else:
                out[c] |= 1 << i
    for i, c in enumerate(clones):
        for d in clones[i+1:]:
            out[c] |= 1 << d
    validate(out)
    return out, clones


def synthesis(x, rx, y, ry, m):
    require(m >= 1, "positive synthesis multiplicity")
    A = [i for i in range(len(x)) if i != rx]
    B = [i for i in range(len(y)) if i != ry]
    mx = {rx: 0, **{v: i+1 for i, v in enumerate(A)}}
    my = {ry: 0, **{v: len(A)+1+i for i, v in enumerate(B)}}
    first = len(A) + len(B) + 1
    Z = list(range(first, first + m))
    out = [0] * (first + m)
    for rows, mapping in ((x, mx), (y, my)):
        for a, aa in mapping.items():
            for b, bb in mapping.items():
                if (rows[a] >> b) & 1:
                    out[aa] |= 1 << bb
    for a in A:
        for b in B:
            out[mx[a]] |= 1 << my[b]
        for z in Z:
            out[mx[a]] |= 1 << z
    for i, z in enumerate(Z):
        out[z] |= 1
        for b in B:
            out[z] |= 1 << my[b]
        for zz in Z[i+1:]:
            out[z] |= 1 << zz
    validate(out)
    return out, Z


def strong(rows):
    n = len(rows)
    full = (1 << n) - 1
    for start in range(n):
        seen = 1 << start
        while True:
            new = seen
            for v in range(n):
                if (seen >> v) & 1:
                    new |= rows[v]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def apex(n, q):
    rows = [0] * (n + 1)
    root = n
    for i in range(n):
        for j in range(i + 1, n):
            rows[i] |= 1 << j
        if (q >> (n - 1 - i)) & 1:
            rows[root] |= 1 << i
        else:
            rows[i] |= 1 << root
    return rows, root


@lru_cache(None)
def f(q):
    if not q:
        return 1
    n = q // 2
    return f(n) if q % 2 else f(n) + n + n % 2


@lru_cache(None)
def slack(q):
    if not q:
        return 0
    n = q // 2
    return slack(n) if q % 2 == 0 else n + slack(n) + n % 2


def isprime(n):
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def switch(rows, root, b, which):
    H = profile(rows, root)[0]
    y, clones = inflate(rows, root, 2)
    c = clones[which]
    g, ey, Ry, L, t = profile(y, c)
    require(t == H and g == 2*H-b, "clone determinant input")
    out, Z = synthesis([6, 4, 0], 0, y, c, 1)
    require(len(out) == len(rows)+4, "switch vertex count")
    N = profile(out, Z[0])[0]
    require(N == 6*L-b and N >= 12*H-7*b, "switch count and lower bound")
    require(not strong(rows) or strong(out), "strongness preservation")
    for m in (2, 3):
        inflated, _ = inflate(out, Z[0], m)
        got = table(inflated)[0][-1]
        require(got == 3*L*2**m-b, "new marked clone ray")
    return out, Z[0], {"input_H": H, "old_offset": b, "marked_Y_clone": which,
                        "Y_profile": [g, ey, Ry, L, t], "new_amplitude": 3*L,
                        "new_offset": -b, "new_H": N,
                        "rows": out, "root": Z[0]}


def main():
    profiles = set()
    rooted_count = 0
    U1_checks = 0
    for n in range(1, 6):
        for code in range(1 << (n*(n-1)//2)):
            rows = tournament(n, code)
            for r in range(n):
                p = profile(rows, r)
                rooted_count += 1
                if n <= 4:
                    require(p == literal_profile(rows, r), "literal profile disagreement")
                h, e, R, L, d = p
                U = R-d
                sigma = R-h-d+e
                require(sigma >= 0 and sigma % 2 == 0, "physical slack")
                if U == 1:
                    require((h, e, R, d, sigma) == (1, 0, 2, 1, 0), "U=1 lemma")
                    U1_checks += 1
                profiles.add(p)

    prime_cases = set()
    checked_pairs = 0
    for h, e, R, LX, d in profiles:
        U = R-d
        sigma = R-h-d+e
        for g, ey, Ry, L, t in profiles:
            alpha, beta = U*L, d*g-sigma*t
            if alpha % 2 == 0 or not isprime(alpha) or beta not in (-1, 1):
                continue
            checked_pairs += 1
            require((U, L, g, t, sigma) == (alpha, 1, 1, 1, d-beta),
                    "general prime-amplitude reduction")
            if d == 1:
                require(alpha in (3, 5) and beta == -1, "d=1 prime-amplitude barrier")
            prime_cases.add((alpha, beta, d, sigma))

    physical_words = 0
    for n in range(8):
        for q in range(1 << n):
            rows, r = apex(n, q)
            h, e, R, L, d = profile(rows, r)
            require((h, e, R, d) == (f(q), 1-q%2, q+1, 1), "apex profile recurrence")
            require(R-h-d+e == slack(q), "apex slack recurrence")
            physical_words += 1
    arithmetic_controls = 65536
    for q in range(arithmetic_controls):
        s = slack(q)
        require(s == q+1-q%2-f(q), "arithmetic sigma identity")
        require(s % 2 == 0 and s >= 0, "sigma parity")
        if q:
            v = q // (q & -q)
            require(v <= 2*s+1 and slack(v) == s, "fixed-slack odd-core theorem")
            require((s == 0) == ((q & (q-1)) == 0), "zero-slack exact classification")
            require((s == 2) == (v in (3, 5)), "slack-two exact classification")

    missing_slack_odd_cores = [(v, slack(v)) for v in range(1, 22, 2)]
    require(all(s != 10 for v, s in missing_slack_odd_cores), "slack10 exhaustive odd cores")
    require([s for v, s in missing_slack_odd_cores] == [0,2,2,6,4,8,8,14,8,14,12],
            "displayed missing-slack certificate")
    require(profile([14,12,24,16,3], 4) == (11,1,13,13,1), "blocked physical partner")

    # Positive prime-amplitude examples with singleton Y.
    prime_witnesses = []
    for q in (3, 5):
        x, rx = apex(q.bit_length(), q)
        ray = []
        for m in range(1, 5):
            out, Z = synthesis(x, rx, [0], 0, m)
            H = table(out)[0][-1]
            require(H == q*2**m-1, "prime-amplitude positive synthesis")
            ray.append(H)
        prime_witnesses.append({"amplitude": q, "offset": -1, "counts": ray})

    x, rx = [16, 1, 3, 7, 14], 2
    require(profile(x, rx) == (9, 2, 16, 16, 5), "larger-deletion prime-amplitude input")
    positive_contrast = []
    for m in range(1, 6):
        out, Z = synthesis(x, rx, [0], 0, m)
        H = table(out)[0][-1]
        require(H == 11*2**m+1 and strong(out), "larger-deletion prime-amplitude ray")
        if m == 1:
            require(out == [6, 48, 34, 39, 45, 1] and Z == [5], "displayed alpha11 rows")
        positive_contrast.append(H)

    finite_partner_checks = 0
    for g, ey, Ry, L, t in profiles:
        for beta in (-1, 1):
            if (g-beta) % t:
                continue
            s = (g-beta)//t
            if s < 0 or s % 2:
                continue
            odd_cores = {v for v in range(1, 2*s+2, 2) if slack(v) == s}
            for q in range(1, 512):
                v = q // (q & -q)
                require((g-slack(q)*t == beta) == (v in odd_cores),
                        "finite-partner dyadic classification")
                finite_partner_checks += 1

    seeds = [([2, 4, 1], 0, 1), ([8, 1, 3, 6], 0, -1),
             ([8, 1, 3, 6], 1, 1), ([2, 12, 9, 17, 7], 0, 1)]
    switches = []
    for rows, r, b in seeds:
        for which in (0, 1):
            _, _, receipt = switch(rows, r, b, which)
            switches.append(receipt)
    rows, r, b = [2, 4, 1], 0, 1
    iterated = []
    for _ in range(2):
        rows, r, receipt = switch(rows, r, b, 0)
        b = -b
        iterated.append(receipt)

    print(json.dumps({"status": "PASS", "rooted_graphs_through_five": rooted_count,
                      "distinct_profiles": len(profiles), "U1_controls": U1_checks,
                      "prime_amplitude_profile_pairs": checked_pairs,
                      "observed_prime_cases_not_a_census_theorem": sorted(prime_cases),
                      "physical_apex_words": physical_words,
                      "arithmetic_slack_controls": arithmetic_controls,
                      "slack10_complete_odd_core_certificate": missing_slack_odd_cores,
                      "prime_amplitude_witnesses": prime_witnesses,
                      "larger_deletion_alpha11_counts": positive_contrast,
                      "finite_partner_checks": finite_partner_checks,
                      "sign_switch_witnesses": switches,
                      "iterated_sign_switches": iterated,
                      "scope": "All-order proofs plus finite controls; no full-spectrum claim"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

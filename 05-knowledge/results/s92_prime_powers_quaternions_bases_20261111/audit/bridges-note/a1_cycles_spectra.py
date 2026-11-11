#!/usr/bin/env python3
"""Auditor bridges-note, script A1 (independent code; written before reading the session scripts).

Sections 7 and 8 of the S92 note:
  - the five known cycles of 3n+1 on Z, their arcs and labels |4n+1|
  - H-spectra of ALL labelled tournaments of orders 3..7 (bitmask DP)
  - classes (isomorphism types) at orders 1..6: counts, strong classes, self-converse classes, |Aut|, H
  - order 4: arc-usage profiles, self-converse iff H = 1 mod 4 iff even number of 3-cycles
  - order 5: the strong classes with H, |Aut|, self-converse flag
  - even order 6: self-converse => H = 1 (mod 4)?  converse direction?
  - explicit witnesses (tournaments) for every label of the five cycles other than 7
"""
import itertools
import sys
import numpy as np
from numba import njit, prange


# ---------------------------------------------------------------- 3n+1 on Z
def step(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def cycle_of(seed):
    cyc = [seed]
    n = step(seed)
    while n != seed:
        cyc.append(n)
        n = step(n)
        if len(cyc) > 1000:
            raise RuntimeError("no cycle")
    return cyc


def label(n):
    return abs(4 * n + 1)


def n_of_count(m):
    chi = 1 if m % 4 == 1 else -1
    assert (chi * m - 1) % 4 == 0
    return (chi * m - 1) // 4


# ---------------------------------------------------------------- tournaments
@njit(cache=True)
def hcount(adj, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            out = adj[v] & ~mask
            w = 0
            while out:
                if out & 1:
                    dp[mask | (1 << w), w] += c
                out >>= 1
                w += 1
    s = 0
    for v in range(n):
        s += dp[full, v]
    return s


@njit(cache=True)
def decode(code, n, adj):
    # pairs (i<j) in lexicographic order; bit set means i -> j
    for v in range(n):
        adj[v] = 0
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> b) & 1:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
            b += 1


@njit(parallel=True, cache=True)
def all_H(n):
    m = n * (n - 1) // 2
    total = 1 << m
    out = np.zeros(total, np.int64)
    for code in prange(total):
        adj = np.zeros(n, np.int64)
        decode(code, n, adj)
        out[code] = hcount(adj, n)
    return out


@njit(cache=True)
def permuted_code(code, n, perm, pairidx):
    # new tournament: perm[i] -> perm[j] iff i -> j
    res = 0
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            bit = (code >> b) & 1
            a = perm[i]
            c = perm[j]
            if a < c:
                if bit:
                    res |= (1 << pairidx[a, c])
            else:
                if not bit:
                    res |= (1 << pairidx[c, a])
            b += 1
    return res


@njit(cache=True)
def canon_one(code, n, perms, pairidx):
    best = -1
    na = 0
    for k in range(perms.shape[0]):
        c = permuted_code(code, n, perms[k], pairidx)
        if best < 0 or c < best:
            best = c
        if c == code:
            na += 1
    return best, na


@njit(parallel=True, cache=True)
def canon_all(n, perms, pairidx):
    m = n * (n - 1) // 2
    total = 1 << m
    can = np.zeros(total, np.int64)
    aut = np.zeros(total, np.int64)
    for code in prange(total):
        b, na = canon_one(code, n, perms, pairidx)
        can[code] = b
        aut[code] = na
    return can, aut


def pair_index(n):
    idx = np.zeros((n, n), np.int64)
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            idx[i, j] = b
            b += 1
    return idx


def adj_of(code, n):
    adj = np.zeros(n, np.int64)
    decode(code, n, adj)
    return adj


def is_strong(code, n):
    adj = adj_of(code, n)
    full = (1 << n) - 1
    for mask in range(1, full):
        ok = True
        for v in range(n):
            if (mask >> v) & 1:
                # all arcs between v and the complement must leave v
                if (adj[v] & (full ^ mask)) != (full ^ mask):
                    ok = False
                    break
        if ok:
            return False
    return True


def c3_count(code, n):
    adj = adj_of(code, n)
    t = 0
    for a, b, c in itertools.combinations(range(n), 3):
        ab = (adj[a] >> b) & 1
        bc = (adj[b] >> c) & 1
        ca = (adj[c] >> a) & 1
        if ab == bc == ca:
            t += 1
    return t


def arc_usage(code, n):
    adj = adj_of(code, n)
    use = {}
    H = 0
    for p in itertools.permutations(range(n)):
        if all((adj[p[i]] >> p[i + 1]) & 1 for i in range(n - 1)):
            H += 1
            for i in range(n - 1):
                use[(p[i], p[i + 1])] = use.get((p[i], p[i + 1]), 0) + 1
    return H, sorted(use.values(), reverse=True), n * (n - 1) // 2 - len(use)


def scores(code, n):
    adj = adj_of(code, n)
    return tuple(sorted(bin(int(a)).count("1") for a in adj))


def main():
    print("=== A1.1 the five known cycles of 3n+1 on Z ===")
    seeds = [0, 1, -1, -5, -17]
    all_labels = []
    for sd in seeds:
        cyc = cycle_of(sd)
        labs = [label(n) for n in cyc]
        all_labels += labs
        print("cycle of %4d: length %2d: %s" % (sd, len(cyc), cyc))
        print("      labels |4n+1|: %s" % labs)
        print("      arcs: %s" % ", ".join("%d->%d" % (n, step(n)) for n in cyc))
    print("total arcs:", len(all_labels), " distinct labels:", len(set(all_labels)))
    print("sorted labels:", sorted(all_labels))
    print("7 in labels:", 7 in all_labels, " 21 in labels:", 21 in all_labels)
    for m in all_labels:
        assert label(n_of_count(m)) == m
    print("bijection check n = (chi_4(m) m - 1)/4 for odd m < 20001:",
          all(label(n_of_count(m)) == m for m in range(1, 20001, 2)),
          "; distinct n:", len({n_of_count(m) for m in range(1, 20001, 2)}) == 10000)

    print()
    print("=== A1.2 H-spectra of all labelled tournaments, orders 3..7 ===")
    spectra = {}
    Hs = {}
    for n in range(3, 8):
        h = all_H(n)
        Hs[n] = h
        vals = sorted(set(int(x) for x in h))
        spectra[n] = vals
        print("order %d: %d labelled tournaments; spectrum (%d values): %s" % (n, len(h), len(vals), vals))
        ns = sorted(n_of_count(v) for v in vals)
        lo, hi = ns[0], ns[-1]
        missing = [k for k in range(lo, hi + 1) if k not in ns]
        print("         as integers n: range [%d, %d], missing %s" % (lo, hi, missing))
        print("         missing odd counts below max: %s" % [m for m in range(1, vals[-1] + 1, 2) if m not in vals])

    print()
    print("=== A1.3 classes at orders 1..6 ===")
    info = {}
    for n in range(1, 7):
        if n == 1:
            print("order 1: 1 class, strong 1, self-converse 1")
            info[1] = (1, 1, 1)
            continue
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
        pidx = pair_index(n)
        can, aut = canon_all(n, perms, pidx)
        m = n * (n - 1) // 2
        full = (1 << m) - 1
        reps = sorted(set(int(c) for c in can))
        nstrong = 0
        nsc = 0
        rows = []
        for r in reps:
            st = is_strong(r, n) if n >= 2 else True
            conv = int(can[full ^ r])
            sc = (conv == r)
            hh = int(hcount(adj_of(r, n), n))
            rows.append((r, hh, int(aut[r]), st, sc, conv, scores(r, n), c3_count(r, n)))
            nstrong += st
            nsc += sc
        info[n] = (len(reps), nstrong, nsc)
        print("order %d: %d classes, strong %d, self-converse %d" % (n, len(reps), nstrong, nsc))
        if n in (4, 5):
            for row in sorted(rows, key=lambda t: (t[1], t[2])):
                print("    rep %5d  H=%3d |Aut|=%d strong=%-5s self-converse=%-5s converse-rep=%5d scores=%s c3=%d"
                      % row)
        if n == 4:
            print("  order 4 arc-usage profiles (nonzero usages; number of unused arcs):")
            for row in sorted(rows, key=lambda t: (t[1], t[0])):
                H, prof, unused = arc_usage(row[0], 4)
                print("    rep %3d H=%d profile %s unused arcs %d" % (row[0], H, prof, unused))
            # labelled census: self-converse labelled tournaments and their H
            sc_lab = {}
            for code in range(64):
                r = int(can[code])
                sc = int(can[63 ^ r]) == r
                key = (sc, int(Hs[4][code]))
                sc_lab[key] = sc_lab.get(key, 0) + 1
            print("  labelled 4-tournaments by (self-converse, H):", sorted(sc_lab.items()))
            ok = True
            for row in rows:
                H, c3, sc = row[1], row[7], row[4]
                nn = n_of_count(H)
                if not (sc == (H % 4 == 1) == (nn >= 0) == (c3 % 2 == 0)):
                    ok = False
            print("  order 4: self-converse iff H=1 mod 4 iff n>=0 iff c3 even:", ok)
        if n == 5:
            print("  order 5 strong classes: H list", sorted(r[1] for r in rows if r[3]),
                  " |Aut| list (same order)", [r[2] for r in sorted(rows, key=lambda t: (t[1], t[2])) if r[3]])
            print("  order 5 strong classes all self-converse:", all(r[4] for r in rows if r[3]))
            print("  order 5 non-self-converse classes: ",
                  [(r[0], r[1], r[3], r[5]) for r in rows if not r[4]], "(rep, H, strong, converse rep)")
        if n in (4, 6):
            viol = [r for r in rows if r[4] and r[1] % 4 != 1]
            conv_fail = [r for r in rows if (not r[4]) and r[1] % 4 == 1]
            print("  order %d: self-converse classes with H != 1 mod 4: %d ; non-self-converse classes with H = 1 mod 4: %d"
                  % (n, len(viol), len(conv_fail)))
            if conv_fail:
                print("    examples (rep, H, |Aut|, strong):", [(r[0], r[1], r[2], r[3]) for r in conv_fail[:6]])
        if n in (3, 5):
            viol = [r for r in rows if r[4] and r[1] % 4 != 1]
            print("  order %d (odd): self-converse classes with H = 3 mod 4: %d  e.g. %s"
                  % (n, len(viol), [(r[0], r[1]) for r in viol[:5]]))

    print()
    print("=== A1.4 generating function 1/(1 - S(x)) ===")
    strong = [info[n][1] for n in range(1, 7)]
    classes = [info[n][0] for n in range(1, 7)]
    print("strong classes by order 1..6:", strong)
    print("all classes by order 1..6:   ", classes)
    for trunc in (5, 6):
        S = [0] + strong[:trunc]
        coef = [1]
        for k in range(1, 8):
            coef.append(sum(S[j] * coef[k - j] for j in range(1, min(k, trunc) + 1)))
        print("1/(1 - S), S truncated at x^%d: coefficients orders 1..7: %s" % (trunc, coef[1:]))

    print()
    print("=== A1.5 witnesses: every label of the five cycles except 7 is a path count ===")
    labs = sorted(set(all_labels))
    small = set()
    for n in range(1, 8):
        if n in spectra:
            small |= set(spectra[n])
    small.add(1)
    wit = {}
    for m in labs:
        if m in small:
            o = min(n for n in spectra if m in spectra[n]) if m > 1 else 1
            wit[m] = "in the exhaustive spectrum of order %d" % o
    # products of two exhaustive-spectrum values (ordered join: H multiplies)
    sm = sorted(small)
    for m in labs:
        if m in wit:
            continue
        for a in sm:
            if a > 1 and m % a == 0 and (m // a) in small:
                wit[m] = "product %d * %d of exhaustive-spectrum values (ordered join)" % (a, m // a)
                break
    rest = [m for m in labs if m not in wit]
    print("labels needing a direct witness of order >= 8:", rest)
    # random search at orders 8, 9, 10
    found = search_witnesses(np.array(rest, dtype=np.int64))
    for m in rest:
        if m in found:
            n, code = found[m]
            adj = adj_of(code, n)
            assert int(hcount(adj, n)) == m
            wit[m] = "explicit tournament of order %d, code %d (bits over pairs i<j lexicographic, bit=1: i->j), H recomputed = %d" % (
                n, code, m)
    for m in labs:
        print("  label %5d (n = %5d): %s" % (m, n_of_count(m), wit.get(m, "NO WITNESS FOUND" if m != 7 else "not in any spectrum of order <= 7 (THM-1370: never)")))
    # sanity: ordered join multiplies H (check on a few pairs)
    chk = True
    rng = np.random.default_rng(1)
    for _ in range(200):
        n1, n2 = int(rng.integers(2, 6)), int(rng.integers(2, 5))
        c1 = int(rng.integers(0, 1 << (n1 * (n1 - 1) // 2)))
        c2 = int(rng.integers(0, 1 << (n2 * (n2 - 1) // 2)))
        a1, a2 = adj_of(c1, n1), adj_of(c2, n2)
        n = n1 + n2
        adj = np.zeros(n, np.int64)
        for v in range(n1):
            adj[v] = int(a1[v]) | (((1 << n2) - 1) << n1)
        for v in range(n2):
            adj[n1 + v] = int(a2[v]) << n1
        if int(hcount(adj, n)) != int(hcount(a1, n1)) * int(hcount(a2, n2)):
            chk = False
    print("ordered join multiplies H (200 random pairs):", chk)


@njit(cache=True)
def _search(targets, n, tries, seed):
    np.random.seed(seed)
    m = n * (n - 1) // 2
    found = np.full(targets.shape[0], -1, np.int64)
    adj = np.zeros(n, np.int64)
    left = targets.shape[0]
    for t in range(tries):
        code = 0
        for b in range(m):
            if np.random.random() < 0.5:
                code |= (1 << b)
        decode(code, n, adj)
        h = hcount(adj, n)
        for k in range(targets.shape[0]):
            if found[k] < 0 and targets[k] == h:
                found[k] = code
                left -= 1
        if left == 0:
            break
    return found


def search_witnesses(targets):
    found = {}
    for n, tries in ((8, 400000), (9, 400000), (10, 200000)):
        rem = np.array([t for t in targets if int(t) not in found], dtype=np.int64)
        if len(rem) == 0:
            break
        res = _search(rem, n, tries, 12345 + n)
        for t, c in zip(rem, res):
            if c >= 0:
                found[int(t)] = (n, int(c))
    return found


if __name__ == "__main__":
    main()

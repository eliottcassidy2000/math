"""Auditor's check of THM-4807 (d): the explicit separating word, built from the stated
recipe (cases alpha / beta / gamma) and tested for legality in both tournaments.

  * every ordered pair T != T' of labelled tournaments, n = 1..5 (n = 5: 1,047,552 pairs);
  * biased random pairs n = 6..12 (planted untraced arcs, many components, flips of
    untraced arcs only, flips of non-consecutive cross arcs only, random flips).
Also measures what goes wrong with the two tempting alternative choices in case gamma
(largest head index replaced by the smallest; least tail index replaced by the largest).
No code from the lane is imported.

usage: python3 audit_sepword.py [max_exhaustive_n=5] [random_pairs_per_n=1500] [seed=4807]
"""

import random
import sys
from collections import Counter

from audit_lib import (
    bits_iter, comp_index, from_mask, hp_tables, in_masks, is_ham_word, is_legal,
    ordered_join, path_from_tables, popcount, random_size, random_strong_on,
    random_tournament, relabel, reverse_arcs, scc_ordered, traced_from_tables,
)


class Info:
    __slots__ = ("n", "out", "inn", "F", "B", "tr", "comps", "idx")

    def __init__(self, n, out):
        self.n = n
        self.out = out
        self.inn = in_masks(n, out)
        self.F, self.B = hp_tables(n, out, self.inn)
        self.tr = traced_from_tables(n, out, self.F, self.B)
        self.comps, chain = scc_ordered(n, out)
        assert chain
        self.idx = comp_index(n, self.comps)

    def path(self, S, start=None, end=None):
        return path_from_tables(self.out, self.inn, self.F, self.B, S, start=start, end=end)


def word_alpha(T, a, b):
    """A Hamilton word of T through the traced arc a -> b."""
    n = T.n
    full = (1 << n) - 1
    for S in range(1, full):
        if (S >> a) & 1 and not (S >> b) & 1 and (T.F[S] >> a) & 1 and (T.B[full ^ S] >> b) & 1:
            return T.path(S, end=a) + T.path(full ^ S, start=b)
    raise AssertionError("traced arc without a word")


def chord_path(T, Bmask, u, v):
    """(b) applied to T[B] and the untraced arc u -> v: the path P'_B from c to v."""
    cu, chain = scc_ordered(T.n, T.out, Bmask & ~(1 << u))
    assert chain and len(cu) >= 2 and (cu[-1] >> v) & 1
    c = next(bits_iter(T.out[u] & cu[0]))
    P = []
    for t, C in enumerate(cu):
        if t == 0:
            P += T.path(C, start=c)
        elif t == len(cu) - 1:
            P += T.path(C, end=v)
        else:
            P += T.path(C)
    return P


def separating_word(T, Tp, variant="thm"):
    """Returns (case, word, legal_in_T, legal_in_Tp)."""
    n = T.n
    diff = [(a, b) for a in range(n) for b in bits_iter(T.out[a]) if (Tp.out[b] >> a) & 1]
    assert diff
    tracedd = [(a, b) for (a, b) in diff if (T.tr[a] >> b) & 1]
    if tracedd:
        a, b = tracedd[0]
        w = word_alpha(T, a, b)
        return "alpha", w, is_ham_word(n, T.out, tuple(w)), is_ham_word(n, Tp.out, tuple(w))
    intra = [(a, b) for (a, b) in diff if T.idx[a] == T.idx[b]]
    if intra:
        u, v = intra[0]
        Bi = T.idx[u]
        Bmask = T.comps[Bi]
        P = chord_path(T, Bmask, u, v)
        w = []
        for t, C in enumerate(T.comps):
            if t == Bi:
                w += P + [u]
            else:
                w += Tp.path(C)  # Redei path of T'[C]
        return "beta", w, is_ham_word(n, T.out, tuple(w)), is_ham_word(n, Tp.out, tuple(w))
    # gamma: every differing arc joins non-consecutive components of T
    for (a, b) in diff:
        if T.idx[b] - T.idx[a] < 2:
            raise AssertionError("case gamma with a consecutive or intra differing arc")
    if variant == "thm":
        i = min(T.idx[a] for a, b in diff)
        j = max(T.idx[b] for a, b in diff if T.idx[a] == i)
    elif variant == "smallest_head":
        i = min(T.idx[a] for a, b in diff)
        j = min(T.idx[b] for a, b in diff if T.idx[a] == i)
    elif variant == "largest_tail":
        i = max(T.idx[a] for a, b in diff)
        j = max(T.idx[b] for a, b in diff if T.idx[a] == i)
    else:
        raise ValueError(variant)
    u, v = next((a, b) for a, b in diff if T.idx[a] == i and T.idx[b] == j)
    C = T.comps
    k = len(C)
    w = []
    for t in range(0, i):
        w += T.path(C[t])
    for t in range(i + 1, j):
        w += T.path(C[t])
    w += T.path(C[j], end=v)
    w += T.path(C[i], start=u)
    for t in range(j + 1, k):
        w += T.path(C[t])
    return "gamma", w, is_ham_word(n, T.out, tuple(w)), is_ham_word(n, Tp.out, tuple(w))


def check_pair(T, Tp, cases, fails, alt):
    case, w, lt, ltp = separating_word(T, Tp)
    cases[case] += 1
    ok = (lt and not ltp) if case == "alpha" else (ltp and not lt)
    if not ok:
        fails[case] += 1
        if fails[case] <= 3:
            print("  FAILURE", case, T.out, Tp.out, w, lt, ltp)
    if case == "gamma":
        for variant in ("smallest_head", "largest_tail"):
            _, w2, lt2, ltp2 = separating_word(T, Tp, variant)
            alt[variant] += 1
            if not (ltp2 and not lt2):
                alt[variant + "_fails"] += 1
                key = variant + "_example"
                if key not in alt:
                    alt[key] = (T.n, T.out, Tp.out, w2)


def exhaustive(n):
    m = n * (n - 1) // 2
    infos = [Info(n, from_mask(n, mask)) for mask in range(1 << m)]
    cases, fails, alt = Counter(), Counter(), {}
    alt = Counter()
    alt_ex = {}
    for a in range(len(infos)):
        for b in range(len(infos)):
            if a != b:
                check_pair(infos[a], infos[b], cases, fails, alt)
    pairs = len(infos) * (len(infos) - 1)
    print(f"labelled n = {n}: {pairs} ordered pairs, cases {dict(sorted(cases.items()))}, "
          f"failures {dict(fails) if fails else 0}")
    print(f"    gamma alternatives: smallest head index fails {alt['smallest_head_fails']}/"
          f"{alt['smallest_head']}, largest tail index fails {alt['largest_tail_fails']}/{alt['largest_tail']}")
    for key in ("smallest_head_example", "largest_tail_example"):
        if key in alt:
            print(f"    first {key}: {alt[key]}")
    sys.stdout.flush()


def planted_untraced(rng, nmax):
    """Busch Cor. 4.14 structure: T - z = C_1 => ... => C_k (k >= 4), z -> x only in C_1,
    y -> z only in C_k, C_2 => z, z => C_(k-1); then x -> y is untraced."""
    while True:
        k = rng.randint(4, 5)
        sizes = [random_size(rng, 1, 4) for _ in range(k)]
        if sum(sizes) + 1 <= nmax:
            break
    n0, out0, comps = ordered_join(sizes, rng)
    n = n0 + 1
    z = n0
    out = out0 + [0]
    x = rng.choice(comps[0])
    y = rng.choice(comps[-1])
    for t, C in enumerate(comps):
        for a in C:
            if t == 0:
                to_z = (a != x)
            elif t == k - 1:
                to_z = (a == y)
            elif t == 1:
                to_z = True
            elif t == k - 2:
                to_z = False
            else:
                to_z = rng.random() < 0.5
            if to_z:
                out[a] |= 1 << z
            else:
                out[z] |= 1 << a
    return n, out


def random_T(rng, n):
    r = rng.random()
    if r < 0.25:
        return random_tournament(n, rng)
    if r < 0.55:
        # ordered join of random strong components
        sizes = []
        while sum(sizes) < n:
            s = random_size(rng, 1, 5)
            if sum(sizes) + s > n:
                s = 1
            sizes.append(s)
        _, out, _ = ordered_join(sizes, rng)
        return out
    # a planted strong block with an untraced arc, inside an ordered join
    while True:
        m, blk = planted_untraced(rng, n)
        if m <= n:
            break
    rest = n - m
    before = rng.randint(0, rest)
    after = rest - before
    out = [0] * n
    order = list(range(before)) + ["B"] + list(range(before, before + after))
    # vertices 0..before-1 singletons/transitive before, block at before..before+m-1, rest after
    blockv = list(range(before, before + m))
    tailv = list(range(before + m, n))
    headv = list(range(before))
    for a in headv:
        for b in blockv + tailv:
            out[a] |= 1 << b
    for a in blockv:
        for b in tailv:
            out[a] |= 1 << b
    for i, a in enumerate(headv):
        for b in headv[i + 1:]:
            out[a] |= 1 << b
    for i, a in enumerate(tailv):
        for b in tailv[i + 1:]:
            out[a] |= 1 << b
    for i in range(m):
        for j in bits_iter(blk[i]):
            out[before + i] |= 1 << (before + j)
    return out


def random_pairs(n, count, rng):
    cases, fails, alt = Counter(), Counter(), Counter()
    made = 0
    attempts = 0
    while made < count and attempts < 50 * count:
        attempts += 1
        out = random_T(rng, n)
        perm = list(range(n))
        rng.shuffle(perm)
        out = relabel(n, out, perm)
        T = Info(n, out)
        untr = [(a, b) for a in range(n) for b in bits_iter(out[a]) if not (T.tr[a] >> b) & 1]
        mode = rng.random()
        if mode < 0.6 and untr:
            k = rng.randint(1, len(untr))
            flips = rng.sample(untr, k)
        elif mode < 0.8:
            allarcs = [(a, b) for a in range(n) for b in bits_iter(out[a])]
            flips = rng.sample(allarcs, rng.randint(1, 3))
        else:
            if not untr:
                continue
            flips = [rng.choice(untr)]
        Tp = Info(n, reverse_arcs(n, out, flips))
        check_pair(T, Tp, cases, fails, alt)
        made += 1
    print(f"random n = {n}: {made} pairs, cases {dict(sorted(cases.items()))}, failures "
          f"{dict(fails) if fails else 0}; gamma alternatives: smallest head fails "
          f"{alt['smallest_head_fails']}/{alt['smallest_head']}, largest tail fails "
          f"{alt['largest_tail_fails']}/{alt['largest_tail']}")
    sys.stdout.flush()


def main():
    maxe = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cnt = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 4807
    for n in range(1, maxe + 1):
        exhaustive(n)
    rng = random.Random(seed)
    for n in range(6, 13):
        random_pairs(n, cnt if n <= 10 else cnt // 3, rng)


if __name__ == "__main__":
    main()

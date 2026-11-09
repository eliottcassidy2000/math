"""The carry-free layer: strong intercyclic tournaments (no two disjoint cycles), their odd-cycle counts A, and their
multiplicities.  In base b > a the b-coloured packing count 1 + a b is realised exactly by the tournaments whose unique
nontrivial strong component is strong intercyclic with a odd cycles (THM-4940 block lemma), so these multiplicities are
the 'Carmichael' multiplicities of the first block of every large base.

Completeness below 1024: by Session 16's classification (feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md)
a strong intercyclic tournament of order n is one-apex (reduced word of length n - 1 starting with 1 and ending with 0,
so alpha_1 >= 2^(n-4) from the first-last pairs) or lies on the two dyadic rays (alpha_1 = 3 * 2^(n-4) or
9 * 2^(n-5) - 2).  Hence alpha_1 < 1024 forces n <= 13, and intercyclic_strong_le13.txt (intercyclic_gen.c 13,
exhaustive by Moon's deletion lemma) is complete there.
Checks: (1) A below 1024 from the generator equals Session 16's formula (one-apex word values plus the two rays);
(2) the generator's counts at orders <= 9 equal the full census (1, 1, 6, 11, 20, 37, 70).
usage: python carry_free_layer.py intercyclic_strong_le13.txt
"""
import sys
import itertools
import collections


def one_apex_alpha1(w):
    m = len(w)
    s = 0
    for i in range(m):
        if w[i] != 1:
            continue
        for j in range(i + 1, m):
            if w[j] == 0:
                r = j - i - 1
                s += 1 if r == 0 else 1 << (r - 1)
    return s


def session16_A(X):
    vals = {0}
    for m in range(2, 15):                           # reduced words 1 u 0; length <= 14 suffices below 2^11
        for mid in itertools.product((0, 1), repeat=m - 2):
            v = one_apex_alpha1((1,) + mid + (0,))
            if v < X:
                vals.add(v)
    t = 1
    while 3 * 2 ** t < X:
        vals.add(3 * 2 ** t)
        t += 1
    t = 0
    while 9 * 2 ** t - 2 < X:
        vals.add(9 * 2 ** t - 2)
        t += 1
    return vals


def main():
    path = sys.argv[1]
    X = 1024
    mult = collections.Counter()
    sc = collections.Counter()
    per_order = collections.Counter()
    for line in open(path):
        f = line.split()
        n, a1, H, tau, isc = int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4])
        per_order[n] += 1
        if H:
            assert H == 1 + 2 * a1
        if a1 < X:
            mult[a1] += 1
            sc[a1] += isc
    assert [per_order[n] for n in range(3, 10)] == [1, 1, 6, 11, 20, 37, 70]
    A = {0} | set(mult)
    S16 = session16_A(X)
    assert A == S16, (sorted(A ^ S16)[:20])
    print(f"A below {X} from the exhaustive generator equals Session 16's formula: {len(A)} values")
    holes = [a for a in range(X) if a not in A]
    print(f"holes of A below {X}: {len(holes)} (density {len(holes) / X:.3f}); first: {holes[:40]}")
    ones = sorted(a for a in mult if mult[a] == 1)
    print(f"carry-free singletons (exactly one strong intercyclic class): {len(ones)}")
    print(f"  a = {ones}")
    print(f"  H = 1 + 2a = {[1 + 2 * a for a in ones]}")
    print(f"  all self-converse: {all(sc[a] == 1 for a in ones)}")
    hist = collections.Counter(mult.values())
    print(f"multiplicity histogram (multiplicity: number of values below {X}): {sorted(hist.items())}")
    print("binary comparison: a strong tournament with two disjoint cycles has order >= 6, so alpha_1 >= 6 (Moon) and")
    print("  alpha_2 >= 1, and its binary alpha-sum alpha_1 + 2 alpha_2 + ... is >= 8; hence for a <= 7 every strong")
    print("  tournament with H = 1 + 2a is intercyclic (carry-free), and the binary strong singletons with H <= 15 are")
    print("  exactly the carry-free singletons with a <= 7.  Below 625 the binary strong singletons are H = 3, 5, 11, 13")
    print("  (THM-4884 / THM-4810), i.e. a = 1, 2, 5, 6; every carry-free singleton with 8 <= a < 312 acquires a carried twin.")
    print(f"  the carry-free singletons with a <= 7 are {[a for a in ones if a <= 7]}")


if __name__ == "__main__":
    main()

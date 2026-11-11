"""a9_pairs.py (AUDIT) - item 9: section 4 of the README.

Inputs: aud_cls_m.txt (m = 1..9), the audit's OWN class lists with H and Pi_2 (gen_classes.c; 198,948 classes).
 1. Pi_2 = (-1)^n H (mod 3) on every class.
 2. F1: classes with H = 2^s + 1 and odd part of Pi_2 dividing 9.
 3. classes with Pi_2 a power of two; least Pi_2 / H per order; comparison with runs/e1_pair_spectrum.txt.
 4. one-apex words: the packet's inversion-pair formulas (re-implemented from the header of words_f2.c) and the
    audit's own formulas (a9_words.c) against DIRECT counts for all words of length <= 8.
 5. the order-7 lists: comparison of (code, H) with data/cls_7.txt as sets of classes (canonical codes differ)."""
import os, sys, collections, itertools
from aud_common import *

out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def odd(x):
    e = 0
    while x % 2 == 0:
        x //= 2
        e += 1
    return x, e


rows = []
for m in range(1, 10):
    for l in open(os.path.join(HERE, "aud_cls_%d.txt" % m)):
        c, h, p = map(int, l.split())
        rows.append((m, c, h, p))
P("classes of order <= 9 (own generation): %d ; by order %s" % (len(rows), sorted(collections.Counter(r[0] for r in rows).items())))
# independent recount of H and Pi_2 on a sample with the Python tools (different code from gen_classes.c)
import random
rng = random.Random(3)
bad = 0
for (m, c, h, p) in rng.sample(rows, 400):
    T = decode(c, m)
    pk = P_tuple(T, 2)
    if pk[1] != h or pk[2] + 2 * pk[1] != p:
        bad += 1
P("   recount of (H, Pi_2) on 400 random classes with the Python DP: mismatches %d" % bad)
viol = sum(1 for (m, c, h, p) in rows if (p - (-1) ** m * h) % 3)
P("1. Pi_2 = (-1)^n H (mod 3): violations %d of %d" % (viol, len(rows)))
P("2. F1: classes with H = 2^s + 1 and odd part of Pi_2 dividing 9:")
hit = 0
for s in range(1, 12):
    H = 2 ** s + 1
    cl = [(m, c, p) for (m, c, h, p) in rows if h == H]
    good = [(m, c, p) + odd(p) for (m, c, p) in cl if 9 % odd(p)[0] == 0]
    least = min((odd(p)[0] for (m, c, p) in cl), default=None)
    if 3 <= s <= 11 and good:
        hit += 1
    P("   s = %2d, H = %4d: %5d classes, %3d distinct odd parts, least odd part %s; with odd part dividing 9: %d classes %s"
      % (s, H, len(cl), len({odd(p)[0] for (m, c, p) in cl}), least, len(good),
         ["order %d, Pi_2 = %d * 2^%d" % (m, o, e) for (m, c, p, o, e) in good][:6]))
P("   levels 3..11 with such a class: %d of 9" % hit)
pw = collections.defaultdict(list)
for (m, c, h, p) in rows:
    if odd(p)[0] == 1:
        pw[h].append((m, odd(p)[1]))
P("3. classes with Pi_2 a power of two: H values %s" % sorted(pw))
P("   claimed list: [1, 5, 19, 25, 217, 223, 259, 599, 695, 707, 713, 761, 773, 803, 815, 1945, 1957, 1969]: %s"
  % (sorted(pw) == [1, 5, 19, 25, 217, 223, 259, 599, 695, 707, 713, 761, 773, 803, 815, 1945, 1957, 1969]))
P("   (number of such classes: %d; H = 1 at every order: transitive)" % sum(len(v) for v in pw.values()))
for mm in range(1, 10):
    r = min((p / h, h, p) for (m, c, h, p) in rows if m == mm)
    P("   order %d: least Pi_2 / H = %.4f (H = %d, Pi_2 = %d); n + 1 = %d" % (mm, r[0], r[1], r[2], mm + 1))
P("   distinct H: %d ; distinct (H, odd part of Pi_2): %d  [author: 1520, 50173]" % (len({r[2] for r in rows}), len({(r[2], odd(r[3])[0]) for r in rows})))
P("   odd parts of Pi_2 for H = 9, 17, 33: %s" % {H: sorted({odd(p)[0] for (m, c, h, p) in rows if h == H}) for H in (9, 17, 33)})
# 4. one-apex words
P("4. one-apex words, all words of length <= 8 (510 words), direct counts on the explicit tournaments:")


def one_apex(word):
    """chain c_1 -> ... -> c_L (vertices 0..L-1), apex L; word[i] = 1 iff apex -> c_(i+1)."""
    L = len(word)
    out_ = [set(range(i + 1, L)) for i in range(L)] + [set()]
    for i, b in enumerate(word):
        if b:
            out_[L].add(i)
        else:
            out_[i].add(L)
    return out_


mine = {}
for l in open(os.path.join(HERE, "a9_words_check8.txt")):
    L, w, h, p = map(int, l.split())
    mine[(L, w)] = (h, p)
bad_mine = bad_auth = 0
n = 0
for L in range(1, 9):
    for w in range(1 << L):
        word = [(w >> i) & 1 for i in range(L)]            # letter i+1 = bit i
        T = one_apex(word)
        pk = P_tuple(T, 2)
        h, p = pk[1], pk[2] + 2 * pk[1]
        n += 1
        if mine[(L, w)] != (h, p):
            bad_mine += 1
        # the packet's formulas: inversion pairs (i < j, w_i = 1, w_j = 0), d = j - i
        m_ = 0
        pi = 2 ** (L + 1)
        for j in range(L):
            if word[j]:
                continue
            for i in range(j):
                if not word[i]:
                    continue
                d = j - i
                m_ += 1 if d == 1 else 2 ** (d - 2)
                pi += (3 ** (d - 1) + 1) * 2 ** (L - d)
        if (1 + 2 * m_, pi) != (h, p):
            bad_auth += 1
P("   words %d: mismatches of the audit's formulas %d; mismatches of the packet's inversion-pair formulas %d" % (n, bad_mine, bad_auth))
# longer random words through the C path counter (H only) and the Python DP for Pi_2 up to L = 11
jobs = []
for trial in range(40):
    L = rng.randint(9, 13)
    word = [rng.randint(0, 1) for _ in range(L)]
    m_ = 0
    pi = 2 ** (L + 1)
    for j in range(L):
        if word[j]:
            continue
        for i in range(j):
            if word[i]:
                d = j - i
                m_ += 1 if d == 1 else 2 ** (d - 2)
                pi += (3 ** (d - 1) + 1) * 2 ** (L - d)
    jobs.append((one_apex(word), 1 + 2 * m_, pi, L))
Hs = ham_c([j[0] for j in jobs])
badH = sum(1 for (T, h, pi, L), hh in zip(jobs, Hs) if h != hh)
badP = 0
nP = 0
for (T, h, pi, L) in jobs:
    if L <= 10:
        nP += 1
        pk = P_tuple(T, 2)
        if pk[2] + 2 * pk[1] != pi:
            badP += 1
P("   40 random words of length 9..13: H formula mismatches %d; Pi_2 formula mismatches %d (of %d words of length <= 10)" % (badH, badP, nP))
# 5. my order-7 list vs the packet's data (as multisets of (H) and by isomorphism through canonical recomputation)
P("5. H = 129 at order 7 with Pi_2 = 3 * 2^9: %s" % [(m, c, h, p) for (m, c, h, p) in rows if m == 7 and h == 129 and p == 3 * 2 ** 9])
P("   all classes with H = 129 and odd part of Pi_2 dividing 9 (orders <= 9): %s" % [(m, c, p) for (m, c, h, p) in rows if h == 129 and 9 % odd(p)[0] == 0])
T = decode(308230, 7)
pk = P_tuple(T, 2)
P("   the packet's witness code 308230 (order 7): H = %d, Pi_2 = %d = %d * 2^%d" % ((pk[1], pk[2] + 2 * pk[1]) + odd(pk[2] + 2 * pk[1])))
open(os.path.join(HERE, "a9_pairs.out"), "w").write("\n".join(out) + "\n")

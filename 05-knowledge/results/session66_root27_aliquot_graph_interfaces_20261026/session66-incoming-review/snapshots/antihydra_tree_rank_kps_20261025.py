#!/usr/bin/env python3
"""
antihydra_tree_rank_kps_20261025.py  (kind-pasteur-2026-10-25-S4)

The source-preserving rank for arbitrarily shaped ternary-tree words on the root-27 (Antihydra) ray.

Setting (inherited from 05-knowledge/results/root27_machine_bridges_20261026/anchor/README.md):
  T(z) = ceil(3z/2) = (3z + c)/2, c = z mod 2 (the carry bit); credit q -> q + 2 if c = 1, q - 1 if c = 0,
  starting at q = 1 for the blank source z = 7.  A halt certificate is a finite full ternary-tree word
  w in {0,1}^n (Tree = 0 | 1 Tree Tree Tree) with  v_2(Q_z(w)) >= |w|,  Q_z(w) = 3^n z + C_w,
  C_w = sum_{i<n} c_i 2^i 3^{n-1-i}.  The open root-27 target is v_2(Q_7(w)) < |w| for every such w.

Theorem checked here (proof in the note):
  (K) there is one infinite carry word kappa(z) in {0,1}^N, independent of n, with
      3^n z + sum_{i<j} kappa_i 2^i 3^{n-1-i} = 0 (mod 2^j) for all j <= n and all n;
      kappa(z) is exactly the sequence of carry bits c_0, c_1, ... of the orbit z, T(z), T^2(z), ...
  (R) for every binary word w of length n (tree word or not), v_2(Q_z(w)) = LCP(w, kappa(z)) whenever
      that common-prefix length is < n, and v_2(Q_z(w)) >= n iff w is the length-n prefix of kappa(z).
  Hence the rank rho_z(w) := |w| - v_2(Q_z(w)) depends on the shape of w only through the position of
  its first disagreement with kappa(z): the root-27 inequality is equivalent to "no prefix of kappa(7)
  is a complete ternary-tree word", i.e. the credit walk of kappa(7) never returns to 0.
Controls: all full ternary-tree words with <= KMAX internal nodes, several sources including hostile
halting sources; the inherited leaf-expansion identity; frontier rank sequence; halting sources <= ZMAX.
"""
import sys
from math import gcd

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("" if detail == "" else "  [" + str(detail) + "]"))


def v2(x):
    if x == 0:
        return 10 ** 9
    return (x & -x).bit_length() - 1


def Q(z, w):
    n = len(w)
    return 3 ** n * z + sum((1 if w[i] == "1" else 0) * 2 ** i * 3 ** (n - 1 - i) for i in range(n))


def kappa(z, length):
    """Carry bits of the orbit of z under T(z) = ceil(3z/2)."""
    bits = []
    for _ in range(length):
        c = z & 1
        bits.append("1" if c else "0")
        z = (3 * z + c) // 2
    return "".join(bits)


def lcp(a, b):
    k = 0
    for x, y in zip(a, b):
        if x != y:
            break
        k += 1
    return k


def tree_words(k):
    """All full ternary-tree words with exactly k internal nodes (preorder: 1 T T T | 0)."""
    memo = {0: ["0"]}

    def gen(k):
        if k in memo:
            return memo[k]
        out = []
        for a in range(k):
            for b in range(k - a):
                c = k - 1 - a - b
                for wa in gen(a):
                    for wb in gen(b):
                        for wc in gen(c):
                            out.append("1" + wa + wb + wc)
        memo[k] = out
        return out

    return gen(k)


def is_tree_word(w):
    q = 1
    for i, ch in enumerate(w):
        q += 2 if ch == "1" else -1
        if q == 0:
            return i == len(w) - 1
        if q < 0:
            return False
    return False


def credit_walk_hits_zero(z, cap):
    """Return the length of the first prefix of kappa(z) that is a complete tree word, or None within cap."""
    q = 1
    for i in range(cap):
        c = z & 1
        q += 2 if c else -1
        if q == 0:
            return i + 1
        z = (3 * z + c) // 2
    return None


# ---------------------------------------------------------------- (K): the carry word is n-independent
for z in (7, 27, 5, 1, 3, 11, 39):
    kap = kappa(z, 60)
    ok = True
    for n in range(1, 41):
        for j in range(0, n + 1):
            val = 3 ** n * z + sum((1 if kap[i] == "1" else 0) * 2 ** i * 3 ** (n - 1 - i) for i in range(j))
            ok &= (val % (2 ** j) == 0)
    check("K1 z=%d: the orbit carry bits satisfy the fixed-source congruences for every n <= 40 and j <= n" % z, ok)
print("kappa(7)[:40] =", kappa(7, 40))
print("kappa(27)[:40] =", kappa(27, 40))

# ---------------------------------------------------------------- (R): v_2(Q_z(w)) = LCP(w, kappa(z))
KMAX = 7
words = []
for k in range(0, KMAX + 1):
    words += tree_words(k)
counts = [len(tree_words(k)) for k in range(0, KMAX + 1)]
check("R0 tree-word counts (1/(2k+1))C(3k,k) for k = 0..%d" % KMAX,
      counts == [1, 1, 3, 12, 55, 273, 1428, 7752][: KMAX + 1], counts)
for z in (7, 27, 5, 1, 3, 11, 39):
    kap = kappa(z, 200)
    ok_eq = True
    ok_full = True
    viol = 0
    for w in words:
        h = v2(Q(z, w))
        p = lcp(w, kap)
        if p < len(w):
            ok_eq &= (h == p)
        else:
            ok_full &= (h >= len(w))
            viol += 1
    check("R1 z=%d: v_2(Q_z(w)) = LCP(w, kappa(z)) for all %d tree words with LCP < |w|" % (z, len(words)), ok_eq)
    check("R2 z=%d: words that are prefixes of kappa(z) have v_2 >= |w| (count %d)" % (z, viol), ok_full, viol)
# random non-tree words too
import random
random.seed(5)
ok = True
for _ in range(20000):
    z = random.choice([7, 27, 5, 1, 3, 11, 39, 101, 2 ** 20 + 1])
    n = random.randint(1, 40)
    w = "".join(random.choice("01") for _ in range(n))
    h = v2(Q(z, w))
    p = lcp(w, kappa(z, n))
    ok &= (h == p) if p < n else (h >= n)
check("R3 the LCP law holds on 20000 random binary words (not necessarily tree words), several sources", ok)

# ---------------------------------------------------------------- leaf expansion identity and the rank law
ok_id = True
ok_rank = True
for z in (7, 27, 5):
    kap = kappa(z, 300)
    for w in words:
        if len(w) > 16:
            continue
        for p in range(len(w)):
            if w[p] != "0":
                continue
            u, v = w[:p], w[p + 1:]
            w2 = u + "1000" + v
            s = len(v)
            Cv = sum((1 if v[i] == "1" else 0) * 2 ** i * 3 ** (s - 1 - i) for i in range(s))
            ok_id &= (Q(z, w2) == 27 * Q(z, w) + 2 ** p * (3 ** (s + 3) - 38 * Cv))
            h = v2(Q(z, w))
            h2 = v2(Q(z, w2))
            if h < len(w):
                if p != h:
                    ok_rank &= (h2 == min(p, h))
                else:
                    ok_rank &= (h2 >= h + 1)
check("E1 inherited leaf-expansion identity Q(u1000v) = 27 Q(u0v) + 2^p (3^(s+3) - 38 C_v) on all expansions, |w| <= 16", ok_id)
check("E2 rank law: expanding at p != h gives precision min(p,h); at p = h precision rises (so rank changes by +3 - delta only there)", ok_rank)

# ---------------------------------------------------------------- frontier of z = 7: credit, precision, rank
z = 7
kap = kappa(z, 400)
q = 1
rows = []
word = ""
for i in range(60):
    word += kap[i]
    q += 2 if kap[i] == "1" else -1
    w = word + "0" * q              # frontier tree word: prefix plus pending leaves
    assert is_tree_word(w)
    h = v2(Q(z, w))
    rows.append((i + 1, kap[i], q, len(w), h, len(w) - h))
print("frontier of source 7: (symbols read, bit, credit q, |w|, precision, rank) for the first 20 symbols:")
for r in rows[:20]:
    print("   ", r)
check("F1 along the frontier the rank |w| - v_2(Q_7(w)) stays >= 1 for 60 symbols (q never hits 0)", all(r[5] >= 1 for r in rows) and all(r[2] >= 1 for r in rows))
check("F2 frontier rank equals q - LCP(0^q, kappa after the prefix) at every step",
      all(r[5] == r[2] - lcp("0" * r[2], kap[r[0]:]) for r in rows))

# ---------------------------------------------------------------- halting sources (hostile positive controls)
ZMAX = 4001
CAP = 20000
halting = {}
for z in range(1, ZMAX + 1, 2):
    t = credit_walk_hits_zero(z, CAP)
    if t is not None:
        halting[z] = t
print("odd sources z <= %d whose credit walk returns to 0 within %d symbols: %d of %d" % (ZMAX, CAP, len(halting), (ZMAX + 1) // 2))
print("first halting sources and halting lengths:", sorted(halting.items())[:15])
# positive control: for each halting source the frontier prefix is a complete tree word with v_2(Q) >= |w|
ok = True
for z, t in list(halting.items())[:200]:
    w = kappa(z, t)
    ok &= is_tree_word(w) and v2(Q(z, w)) >= len(w)
check("H1 for halting sources the kappa prefix is a complete tree word violating the inequality (positive control, 200 sources)", ok)
check("H2 the blank source 7 and 27 do not halt within the cap", 7 not in halting and 27 not in halting)
# larger census against the i.i.d. model: an odd source starts with bit 1, so its credit is 3 after one
# symbol and the model's halting probability is phi^-3 (klein A5 gives phi^-6 for the blank source after 1110).
PHI = (1 + 5 ** 0.5) / 2
ZMAX2, CAP2 = 40001, 500
halt2 = 0
tot2 = 0
lengths = {}
for z in range(1, ZMAX2 + 1, 2):
    tot2 += 1
    t = credit_walk_hits_zero(z, CAP2)
    if t is not None:
        halt2 += 1
        lengths[t] = lengths.get(t, 0) + 1
share = halt2 / tot2
print("odd sources z <= %d halting within %d symbols: %d of %d = %.4f; model phi^-3 = %.4f; phi^-6 = %.4f"
      % (ZMAX2, CAP2, halt2, tot2, share, PHI ** -3, PHI ** -6))
print("halting-length distribution (length: count), first 12 lengths:", sorted(lengths.items())[:12])
check("H3 the halting share of odd sources <= %d lies within 0.01 of phi^-3 (model comparison, FINITE-EXACT share)" % ZMAX2,
      abs(share - PHI ** -3) < 0.01, round(share, 4))
# how long do the walks from 7 and 27 run: report the minimum credit seen
for z0 in (7, 27):
    z = z0
    q = 1
    qmin = 10 ** 9
    for i in range(200000):
        c = z & 1
        q += 2 if c else -1
        qmin = min(qmin, q)
        z = (3 * z + c) // 2
    print("source %d: credit after 200000 symbols = %d, minimum credit seen = %d" % (z0, q, qmin))

bad = [n for n, ok in CHECKS if not ok]
print("%d checks, %d failed" % (len(CHECKS), len(bad)))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL CHECKS PASSED")

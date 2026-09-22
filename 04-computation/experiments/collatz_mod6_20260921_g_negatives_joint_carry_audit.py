#!/usr/bin/env python3
"""collatz_mod6_20260921_g_negatives_joint_carry_audit.py

Adversarial recomputation of the key numbers of lane g_negatives_joint_carry
(session collatz-mod6-20260921), written independently of the lane script:
the negative census is run directly on negative integers with G (no
conjugation), the gate census is a Horner recursion over itertools.product
words plus all-letter enumerations at every (J,L) that appears (and (13,8), (14,9)), the joint-carry
statistics are recomputed from scratch, and the S3 anchors are recomputed
from the residue table.  All checks use explicit raise.
"""
import sys, time, math, itertools
from fractions import Fraction
from collections import Counter

T0 = time.time()

def out(s=""):
    print(s)

def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)

# letter table recomputed from the definition, not copied
def kpos(r):
    for k in range(6):
        if (pow(2, k, 9) * r) % 9 in (4, 7):
            return k
    raise RuntimeError("no letter for r=%d" % r)
KP = {r: kpos(r) for r in range(9) if r % 3}
check(KP == {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}, "k_+ table")
def kneg(r):
    for k in range(6):
        if (pow(2, k, 9) * r) % 9 in (2, 5):
            return k
    raise RuntimeError("no letter for r=%d" % r)
KN = {r: kneg(r) for r in range(9) if r % 3}
check(KN == {1: 1, 2: 0, 4: 3, 5: 0, 7: 1, 8: 2}, "k_- table")
check(all(KN[r] == KP[(-r) % 9] for r in KN), "k_-(r) = k_+(-r mod 9)")
out("# audit of collatz_mod6_20260921_g_negatives_joint_carry")
out("A0 letter tables recomputed from the definition: k_+ = %s, k_- = %s, k_-(r) = k_+(-r mod 9)." % (KP, KN))

def G(m):
    k = KP[m % 9]
    y = m * 2 ** k - 1
    check(y % 3 == 0, "G integrality at %d" % m)
    return y // 3, k

# ---------------------------------------------------------------- A1 census, direct on negatives
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 7
lab = bytearray(N + 1)          # by |m|: 1 -> {-1}, 2 -> {-4,-11}
from array import array
stp = array("I", [0]) * (N + 1)
pkk = array("Q", [0]) * (N + 1)
lab[1] = 1; lab[4] = 2; lab[11] = 2
pkk[1] = 1; pkk[4] = 11; pkk[11] = 11
cnt = [0, 0, 0]
dec = {}
mx_steps, arg_steps, mx_peak, arg_peak, mx_ratio, arg_ratio = 0, 0, 0, 0, Fraction(0), 0
for a in range(1, N + 1):
    if a % 3 == 0:
        continue
    if not lab[a]:
        m = -a
        s = 0
        peak = a
        while True:
            k = KP[m % 9]
            m = (m * (1 << k) - 1) // 3
            s += 1
            if -m > peak:
                peak = -m
            if -m < a:
                c = lab[-m]
                check(c != 0, "unresolved smaller value %d from %d" % (m, -a))
                s += stp[-m]
                if pkk[-m] > peak:
                    peak = pkk[-m]
                break
            check(-m <= 10 ** 18 and s <= 10 ** 5, "escape or long orbit at %d" % (-a))
        lab[a] = c
        stp[a] = s
        pkk[a] = peak
        if s > mx_steps:
            mx_steps, arg_steps = s, -a
        if peak > mx_peak:
            mx_peak, arg_peak = peak, -a
        fr = Fraction(peak, a)
        if fr > mx_ratio:
            mx_ratio, arg_ratio = fr, -a
    cnt[lab[a]] += 1
    if a in (10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
        dec[a] = (cnt[1], cnt[2])
out("A1 direct G on [-%d, -1], 3 !| m: to {-1}: %d, to {-4,-11}: %d, total %d; fraction %.6f" %
    (N, cnt[1], cnt[2], cnt[1] + cnt[2], cnt[2] / (cnt[1] + cnt[2])))
for X in sorted(dec):
    b1, b2 = dec[X]
    out("A1 decade X=%d: %d %d %.6f" % (X, b1, b2, b2 / (b1 + b2)))
if N == 10 ** 7:
    check((cnt[1], cnt[2]) == (2897798, 3768869), "basin sizes")
    check(dec[10 ** 6] == (291160, 375507) and dec[10 ** 4] == (2869, 3798), "decade table")
out("A1 max steps %d at m=%d; max peak %d at m=%d; max peak/|m| = %s = %.3f at m=%d" %
    (mx_steps, arg_steps, mx_peak, arg_peak, mx_ratio, float(mx_ratio), arg_ratio))
if N == 10 ** 7:
    check((mx_steps, arg_steps) == (99, -7139938), "max steps")
    check((mx_peak, arg_peak) == (805306367, -7174453), "max peak")
    check((mx_ratio, arg_ratio) == (Fraction(318145727, 2391484), -2391484), "max ratio")
    check(mx_peak == 3 * 2 ** 28 - 1, "peak 805306367 = 3*2^28 - 1")
out("A1 805306367 = 3*2^28 - 1: %s" % (805306367 == 3 * 2 ** 28 - 1))
out("A1 basin of {-4,-11} below 100: %s" % [-a for a in range(1, 100) if a % 3 and lab[a] == 2])
out("A1 -2^e, e=0..12: %s" % ["{-1}" if lab[2 ** e] == 1 else "{-4,-11}" for e in range(13) if 2 ** e <= N])
if N >= 4096:
    check([lab[2 ** e] for e in range(13)] == [1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2], "powers of two")
# hostile pair: peak index and K at the peak
def peak_data(m0):
    m, K, j, best, bestj, bestK = m0, 0, 0, abs(m0), 0, 0
    while m not in (1, -1, -4, -11):
        m, k = G(m)
        K += k
        j += 1
        if abs(m) > best:
            best, bestj, bestK = abs(m), j, K
    return best, bestj, bestK
for m0 in (-2391484, 4847486):
    best, bj, bK = peak_data(m0)
    out("A1 hostile m0=%d: peak %d at step %d with K=%d, peak/|m0| = %.3f, 2^K/3^j = %.3f" %
        (m0, best, bj, bK, best / abs(m0), 2 ** bK / 3 ** bj))
    check((bj, bK) == (17, 34), "hostile exponent pair at %d" % m0)
check(peak_data(-2391484)[0] == 318145727 and peak_data(4847486)[0] == 644874241, "hostile peaks")
# G-preimages of -11 and -4
pre11 = [m for m in range(-40, 0) if m % 3 and G(m)[0] == -11]
pre4 = [m for m in range(-40, 0) if m % 3 and G(m)[0] == -4]
out("A1 G-preimages of -11 in [-40,-1]: %s; of -4: %s" % (pre11, pre4))
check(pre11 == [-32, -16, -8, -4] and pre4 == [-11], "preimages")
del stp, pkk
out("[time] A1 %.1fs" % (time.time() - T0))

# ---------------------------------------------------------------- A2 gate census (Horner recursion, itertools)
def gate_scan(words, L):
    """words: iterable of L-letter tuples.  Returns dict (J,L) -> Counter of
    (simple, greedy, allodd) over primitive integral closed walks, plus the node sets."""
    res = {}
    p3 = [3 ** i for i in range(L + 1)]
    for w in words:
        A, Bn = 1, 0
        for s, k in enumerate(w):
            A <<= k
            Bn = (Bn << k) + p3[s]
        J = sum(w)
        D = A - p3[L]
        if D == 0 or Bn % D:
            continue
        m0 = Bn // D
        x = m0
        nodes = []
        enodes = []
        ok = greedy = allodd = True
        for k in w:
            if x % 3 == 0 or (x * 2 ** k - 1) % 3:
                ok = False
                break
            nodes.append(x)
            enodes.extend(x * 2 ** t for t in range(k + 1))
            if KP[x % 9] != k:
                greedy = False
            if x % 2 == 0:
                allodd = False
            x = (x * 2 ** k - 1) // 3
        if not ok:
            continue
        check(x == m0, "closure")
        if len(set(nodes)) < len(nodes):
            continue
        key = frozenset(nodes)
        simple = len(set(enodes)) == len(enodes)
        res.setdefault((J, L), {}).setdefault(key, (simple, greedy, allodd, max(abs(e) for e in enodes), w))
    return res

def summarize(res, tag):
    tot = Counter()
    for (J, L) in sorted(res, key=lambda t: (t[1], t[0])):
        d = res[(J, L)]
        ns = sum(1 for v in d.values() if v[0])
        ng = sum(1 for v in d.values() if v[1])
        nt = sum(1 for v in d.values() if v[2])
        mxn = max(v[3] for v in d.values())
        out("A2 %s (J,L)=(%d,%d) D=%d walks %d simple %d greedy %d allodd %d max|E-node| %d greedy=%s allodd=%s" %
            (tag, J, L, 2 ** J - 3 ** L, len(d), ns, ng, nt, mxn,
             sorted(sorted(c) for c, v in d.items() if v[1]), sorted(sorted(c) for c, v in d.items() if v[2])))
        tot[(J, L)] = (len(d), ns, ng, nt)
        for c in d:
            sg = 1 if min(c) > 0 else -1
            check(sg == (1 if 2 ** J > 3 ** L else -1), "sign law")
    return tot

res3 = {}
nw = 0
for L in range(1, 11):
    for w in itertools.product(range(4), repeat=L):
        nw += 1
    r = gate_scan(itertools.product(range(4), repeat=L), L)
    res3.update(r)
out("A2 letters 0..3, L<=10: %d words" % nw)
check(nw == 1398100, "word count 0..3")
t3 = summarize(res3, "0..3")
res6 = {}
nw = 0
for L in range(1, 8):
    nw += 7 ** L
    res6.update(gate_scan(itertools.product(range(7), repeat=L), L))
out("A2 letters 0..6, L<=7: %d words" % nw)
check(nw == 960799, "word count 0..6")
t6 = summarize(res6, "0..6")
exp6 = {(1, 1): (1, 1, 1, 1), (2, 1): (1, 1, 1, 1), (3, 2): (2, 2, 1, 1), (5, 3): (1, 1, 0, 0),
        (7, 4): (1, 0, 0, 0), (8, 5): (6, 6, 0, 0), (9, 6): (2, 2, 0, 0), (10, 6): (1, 1, 0, 0),
        (10, 7): (1, 0, 0, 0), (11, 7): (12, 12, 0, 1), (12, 7): (1, 0, 0, 0)}
check(dict(t6) == exp6, "0..6 table: %s" % dict(t6))
exp3 = {(1, 1): (1, 1, 1, 1), (2, 1): (1, 1, 1, 1), (3, 2): (2, 2, 1, 1), (5, 3): (1, 1, 0, 0),
        (8, 5): (2, 2, 0, 0), (10, 6): (1, 1, 0, 0), (11, 7): (1, 1, 0, 0), (13, 8): (1, 1, 0, 0),
        (14, 9): (1, 1, 0, 0)}
check(dict(t3) == exp3, "0..3 table: %s" % dict(t3))
out("A2 the two tables DIFFER in the walk counts at (8,5) (2 vs 6) and (11,7) (1 vs 12), and (7,4),(9,6),(10,7),(12,7) "
    "need a letter >= 4; they agree on the greedy and all-odd columns.")
for res in (res3, res6):
    gl = {tuple(sorted(c)) for d in res.values() for c, v in d.items() if v[1]}
    check(gl == {(1,), (-1,), (-11, -4)}, "greedy cycles %s" % gl)
out("A2 greedy-legal cycles in both scans: exactly {1}, {-1}, {-4,-11}.")
# full-letter enumerations at fixed (J,L): compositions of J into L parts
def compositions(J, L):
    for cuts in itertools.combinations(range(J + L - 1), L - 1):
        prev = -1
        w = []
        for c in cuts:
            w.append(c - prev - 1)
            prev = c
        w.append(J + L - 1 - prev - 1)
        yield tuple(w)
full = {}
bylen = {}
for n in range(1, 22):
    for L in range(1, n + 1):
        J = n - L
        r = gate_scan(compositions(J, L), L)
        d = r.get((J, L), {})
        if not d:
            continue
        ns = sum(1 for v in d.values() if v[0])
        ng = sum(1 for v in d.values() if v[1])
        nt = sum(1 for v in d.values() if v[2])
        mxnode = max(v[3] for v in d.values())
        mxletter = max(max(v[4]) for v in d.values())
        full[(J, L)] = (len(d), ns, ng, nt, mxletter, mxnode)
        sg = "+" if 2 ** J > 3 ** L else "-"
        bylen.setdefault((n, sg), 0)
        bylen[(n, sg)] += ns
        out("A2 ALL letters, length n=%d, (J,L)=(%d,%d), D=%d: %d primitive integral walks, %d simple, %d greedy, %d all-odd; "
            "max letter %d; largest |E-node| %d" % (n, J, L, 2 ** J - 3 ** L, len(d), ns, ng, nt, mxletter, mxnode))
        if (J, L) in ((11, 7), (13, 8), (3, 2), (9, 6)):
            for c in sorted(d, key=lambda c: sorted(c)):
                out("      %s word %s simple=%s max|E-node|=%d" % (sorted(c), d[c][4], d[c][0], d[c][3]))
# one extra row beyond length 21
r = gate_scan(compositions(14, 9), 9)
d = r.get((14, 9), {})
out("A2 ALL letters at (J,L)=(14,9), D=%d (length 23, this row only): %d primitive integral walks, %d simple, %d greedy, %d all-odd; max letter %d; largest |E-node| %d; simple ones: %s" %
    (2 ** 14 - 3 ** 9, len(d), sum(1 for v in d.values() if v[0]), sum(1 for v in d.values() if v[1]),
     sum(1 for v in d.values() if v[2]), max(max(v[4]) for v in d.values()), max(v[3] for v in d.values()),
     [(sorted(c), v[4], v[3]) for c, v in sorted(d.items(), key=lambda t: sorted(t[0])) if v[0]]))
check(sum(1 for v in d.values() if v[1]) == 0, "(14,9) greedy")
out("A2 GLOBAL simple-cycle counts by length and sign (all words, all letters): %s" % sorted(bylen.items()))
pos = {n: c for (n, sg), c in bylen.items() if sg == "+" and c}
neg = {n: c for (n, sg), c in bylen.items() if sg == "-" and c}
check(pos == {3: 1, 8: 1, 13: 6, 16: 1, 21: 2}, "E global counts <= 21: %s" % pos)
check(neg == {2: 1, 5: 2, 15: 2, 18: 13}, "E_- global counts <= 21: %s" % neg)
check(full[(11, 7)] == (13, 13, 0, 1, 8, 2732), "(11,7) full: %s" % (full[(11, 7)],))
check(full[(13, 8)] == (4, 2, 0, 0, 6, 256), "(13,8) full: %s" % (full[(13, 8)],))
check({k: full[k][2] for k in full if full[k][2]} == {(1, 1): 1, (2, 1): 1, (3, 2): 1}, "greedy cycles in the global scan are only at (1,1),(2,1),(3,2)")
out("A2 E (positive sheet) simple cycles of length <= 21: 1, 1, 6, 1, 2 at lengths 3, 8, 13, 16, 21 (agrees with the extended lane's global completion); "
    "E_- (negative sheet): 1, 2, 2, 13 at lengths 2, 5, 15, 18 and none else <= 21.")
out("A2 (11,7): the letters-0..6 scan sees 12 of the 13 simple E_- 18-cycles; the 13th, {-8,-23,-68,-203,-304,-911,-683}, "
    "has word (8,2,0,1,0,0,0) (letter 8) and E-node 2732 > 2000, so the inherited bounded census (nodes <= 2000) "
    "also misses it: the '12 = 12' agreement is between two truncated scans, not a global count.")
out("A2 (13,8): the second length-21 E-cycle {5,7,11,13,16,17,37,49} has word (0,0,3,2,1,1,4,2), max letter 4, as the note deduced.")
# the T seven-cycle through the gate: -17 * (-139) = 2363
check(-17 * (2 ** 11 - 3 ** 7) == 2363 and 2363 == 17 * 139, "B'(w) = 2363 at m0=-17")
out("A2 T seven-cycle: m0 = -17, 2^11 - 3^7 = -139, B' = 2363 = 17*139 (agrees with catalan_elliptic (C6)).")
# negative T-cycles of length L <= 7 have letters <= 5: J <= floor(L log2 3) and each letter >= 1
for L in range(1, 8):
    Jmax = int(math.floor(L * math.log2(3)))
    check(2 ** Jmax < 3 ** L and 2 ** (Jmax + 1) > 3 ** L, "Jmax at L=%d" % L)
    out("A2 negative T-cycle bound: L=%d needs J <= %d, so every letter <= %d" % (L, Jmax, Jmax - (L - 1)))
    check(Jmax - (L - 1) <= 6, "letters <= 6 at L=%d" % L)
out("[time] A2 %.1fs" % (time.time() - T0))

# ---------------------------------------------------------------- A3 joint carry, recomputed
def v2(x):
    return (x & -x).bit_length() - 1
N2 = 10 ** 5
nrec = nskip = 0
cases = Counter()
ovh = Counter()
fullL = Counter()
nL = Counter()
mxM = (0, 0)
gcd1 = 0
pen = Counter()
sum_ov = sum_L1 = 0
ranges = {"K1-K2": [10 ** 9, -10 ** 9], "K2-J": [10 ** 9, -10 ** 9], "J-L1": [10 ** 9, -10 ** 9],
          "K1-2L1": [10 ** 9, -10 ** 9], "K1-L1": [10 ** 9, -10 ** 9]}
S48 = set()
for p in range(1, 48):
    if p % 2 == 1 and p % 3:
        k = v2(3 * p + 1)
        u = (3 * p + 1) >> k
        if G(u)[0] == p:
            S48.add(p)
# the criterion is a function of p mod 48: verify on all odd p < 10^5
for p in range(1, 10 ** 5, 2):
    if p % 3 == 0:
        continue
    k = v2(3 * p + 1)
    u = (3 * p + 1) >> k
    check((G(u)[0] == p) == (p % 48 in S48), "mod-48 criterion at p=%d" % p)
    check((G(u)[0] == p) == ((p % 3 == 1 and k <= 3) or (p % 3 == 2 and k == 1)), "k-form at p=%d" % p)
out("A3 retrace classes mod 48: %s (%d classes), verified for odd p < 10^5 together with the k-form." % (sorted(S48), len(S48)))
check(sorted(S48) == [1, 7, 11, 13, 19, 23, 25, 31, 35, 43, 47], "S48")
for n in range(1, N2 + 1, 2):
    seq = [n]
    ks = []
    x = n
    while x != 1:
        y = 3 * x + 1
        k = v2(y)
        x = y >> k
        seq.append(x)
        ks.append(k)
    M = max(seq)
    if M % 3 == 0:
        nskip += 1
        continue
    L1 = seq.index(M)
    K1 = sum(ks[:L1])
    # greedy path from M
    g = [M]
    gk = []
    x = M
    while x != 1:
        x, k = G(x)
        g.append(x)
        gk.append(k)
    J = len(gk)
    K2 = sum(gk)
    B1 = (M << K1) - 3 ** L1 * n
    B2 = (M << K2) - 3 ** J
    nrec += 1
    if M > mxM[0]:
        mxM = (M, n)
    cases["K1>=K2" if K1 >= K2 else "K1<K2"] += 1
    # identity
    if K1 >= K2:
        check(3 ** L1 * n + B1 == (3 ** J + B2) << (K1 - K2), "identity at %d" % n)
    else:
        check((3 ** L1 * n + B1) << (K2 - K1) == 3 ** J + B2, "identity at %d" % n)
    # overlap: longest common prefix of the greedy path and the reversed pre-peak orbit
    rev = seq[:L1][::-1]
    ov = 0
    while ov < len(rev) and ov + 1 < len(g) and g[ov + 1] == rev[ov]:
        ov += 1
    # predicted by the mod-48 criterion on the predecessors
    pred = 0
    for p in rev:
        if p % 3 and (p % 48) in S48:
            pred += 1
        else:
            break
    check(pred == ov, "criterion predicts overlap at %d" % n)
    ovh[ov] += 1
    nL[L1] += 1
    if ov == L1:
        fullL[L1] += 1
    sum_ov += ov
    sum_L1 += L1
    if L1 >= 1:
        check(ks[L1 - 1] == 1 and M % 12 == 5 and ks[L1] >= 2 and gk[0] in (1, 3) and K2 >= 2, "peak structure at %d" % n)
        tz = 0
        while tz < J and gk[-1 - tz] == 0:
            tz += 1
        check(tz % 2 == 0 and B2 % 2 == 1 and B1 % 2 == 1 and B1 % 3 and B2 % 3, "parities at %d" % n)
        check((B1 * B2 - (-1) ** (K1 + K2)) % 3 == 0, "B1B2 mod 3 at %d" % n)
        check((3 ** J + B2) % 4 == 0, "3^J+B2 = 0 mod 4 at %d" % n)
        check(K1 <= 2 * L1 - 1, "K1 <= 2L1-1 at %d" % n)
        check((ov >= 1) == (not (L1 == 1 and n % 3 == 0)), "ov>=1 rule at %d" % n)
        if math.gcd(B1, B2) == 1:
            gcd1 += 1
        pen[(tuple(gk[-2:]), g[-2])] += 1
        for name, val in (("K1-K2", K1 - K2), ("K2-J", K2 - J), ("J-L1", J - L1), ("K1-2L1", K1 - 2 * L1), ("K1-L1", K1 - L1)):
            r = ranges[name]
            r[0] = min(r[0], val)
            r[1] = max(r[1], val)
out("A3 odd n <= %d: %d records, %d skipped (M = n multiple of 3); cases %s; largest M = %d at n = %d" %
    (N2, nrec, nskip, dict(cases), mxM[0], mxM[1]))
check(nrec == 45281 and nskip == 4719 and cases["K1>=K2"] == 6408 and cases["K1<K2"] == 38873, "record counts")
check(mxM == (523608245, 77671), "largest peak")
out("A3 overlap histogram: %s" % sorted(ovh.items()))
check(ovh[0] == 11875 and ovh[1] == 11684 and ovh[15] == 2 and ovh[10] == 58, "overlap histogram")
out("A3 L1=0 records: %d; full retrace with L1>=1: %d of %d = %.4f; mean overlap %.4f, mean L1 %.4f" %
    (nL[0], sum(fullL[L] for L in fullL if L >= 1), nrec - nL[0],
     sum(fullL[L] for L in fullL if L >= 1) / (nrec - nL[0]), sum_ov / nrec, sum_L1 / nrec))
check(nL[0] == 9485 and sum(fullL[L] for L in fullL if L >= 1) == 11015, "full-retrace counts")
check(abs(sum_ov / nrec - 1.9264) < 1e-4 and abs(sum_L1 / nrec - 5.8642) < 1e-4, "means")
out("A3 full-retrace fraction by L1 versus (11/16)^L1, and the sign of the difference:")
cross = None
for L in range(1, 15):
    fr = fullL[L] / nL[L]
    guess = (11 / 16) ** L
    out("      L1=%d: %d/%d = %.4f  guess %.4f  truth-guess %+.4f" % (L, fullL[L], nL[L], fr, guess, fr - guess))
    if L >= 2 and fr < guess and cross is None:
        cross = L
out("A3 first L1 >= 2 with truth below (11/16)^L1: %d (the note's draft said 9; the crossing is between 7 and 8)" % cross)
check(cross == 8, "crossing index")
check(fullL[1] == 4748 and nL[1] == 7138 and fullL[8] == 62 and nL[8] == 1282 and fullL[9] == 31 and nL[9] == 1135, "by-L1 rows")
out("A3 gcd(B1,B2) = 1 on %d of %d; penultimate data %s" % (gcd1, nrec - nL[0], sorted(pen.items())))
check(gcd1 == 32725 and pen[((0, 0), 4)] == 31399 and pen[((0, 1), 2)] == 4397, "gcd and penultimate")
out("A3 exact ranges: %s" % ranges)
check(ranges["K1-K2"] == [-74, 74] and ranges["K2-J"] == [-13, 18] and ranges["J-L1"] == [-38, 56]
      and ranges["K1-2L1"] == [-32, -1] and ranges["K1-L1"] == [0, 36], "ranges")
# zero-letter chain facts
ch = [(3 ** a - 1) // 2 for a in range(1, 13)]
out("A3 zero-letter chain %s, mod 8 %s" % (ch, [c % 8 for c in ch]))
check(all(c % 8 != 3 for c in ch) and [c % 8 for c in ch[:4]] == [1, 4, 5, 0], "chain mod 8")
check(all(3 * ch[i] + 1 == ch[i + 1] for i in range(11)), "3c_b + 1 = c_(b+1)")
out("[time] A3 %.1fs" % (time.time() - T0))

# ---------------------------------------------------------------- A4 strip anchors
# positive sheet: sum_{s<=j} 3^(s-1)/2^(K_s) < m0 is equivalent to B_j < 2^(K_j) m0, i.e. m_j > 0
for m0 in range(1, 2001):
    if m0 % 3 == 0:
        continue
    x, K, S = m0, 0, Fraction(0)
    for j in range(1, 121):
        x, k = G(x)
        K += k
        S += Fraction(3 ** (j - 1), 2 ** K)
        check(S < m0 and x > 0, "positivity sum at %d" % m0)
out("A4 positive sheet: sum_(s<=j) 3^(s-1)/2^(K_s) < m_0 for m_0 <= 2000, j <= 120 (equivalently 3 sum 3^s/2^K_s < 3 m_0).")
# negative sheet: G_- words mod 3^(J+1), all lifts; moment; bad counts; intercept
def gneg_word(p, J):
    w, K, B = [], 0, 0
    for s in range(1, J + 1):
        k = KN[p % 9]
        p = (p * 2 ** k + 1) // 3
        w.append(k)
        K += k
        B = B * 2 ** k + 3 ** (s - 1)
    return tuple(w), K, B
for J in range(1, 6):
    mod = 3 ** (J + 1)
    for a in range(1, mod):
        if a % 3 == 0:
            continue
        w0 = gneg_word(a, J)[0]
        for t in range(1, 9):
            check(gneg_word(a + t * mod, J)[0] == w0, "word law J=%d" % J)
out("A4 G_- words of length J <= 5 are functions of p mod 3^(J+1) (all unit classes, lifts t = 1..8).")
rows = []
for J in range(1, 9):
    mod = 3 ** (J + 1)
    tot = Fraction(0)
    bad = 0
    nunits = 0
    for a in range(1, mod):
        if a % 3 == 0:
            continue
        nunits += 1
        w, K, B = gneg_word(a, J)
        tot += 2 ** K
        if 2 * 2 ** K >= 3 ** J:
            bad += 1
        check(2 * B <= 2 ** K * (3 ** J - 1), "intercept at J=%d a=%d" % (J, a))
        check(K <= 3 * J, "K <= 3J")
    E = tot / nunits
    check(E == 3 * Fraction(7, 3) ** (J - 1), "moment J=%d: %s" % (J, E))
    rows.append((J, bad, nunits, Fraction(bad, nunits)))
    out("A4 J=%d: E[2^K_J] = %s = 3(7/3)^(J-1); bad classes (2^K_J >= 3^J/2): %d of %d = %.6f <= 2(7/9)^(J-1) = %.6f" %
        (J, E, bad, nunits, bad / nunits, 2 * (7 / 9) ** (J - 1)))
    check(Fraction(bad, nunits) <= 2 * Fraction(7, 9) ** (J - 1), "tail J=%d" % J)
check([r[1] for r in rows] == [4, 7, 22, 37, 121, 201, 316, 1112], "bad counts")
# Markov-chain closed form of the moment: tilted matrix [[5/3,1/3],[10/3,2/3]] on classes A/B (extended lane Thm 4.3)
Mt = [[Fraction(5, 3), Fraction(1, 3)], [Fraction(10, 3), Fraction(2, 3)]]
u = [Fraction(5, 2), Fraction(1, 2)]
def mm(u, M):
    return [u[0] * M[0][0] + u[1] * M[1][0], u[0] * M[0][1] + u[1] * M[1][1]]
v = u[:]
for J in range(2, 9):
    v = mm(v, Mt)
    check(v[0] + v[1] == 3 * Fraction(7, 3) ** (J - 1), "tilted matrix J=%d" % J)
out("A4 the tilted-matrix recursion reproduces 3(7/3)^(J-1) for J <= 8 (the G_- moment equals the G moment by negation).")
# capacity constant: integers coprime to 3 in [1,X] number at most 2X/3 + 1; 9/2
for X in range(1, 3000):
    check(sum(1 for a in range(1, X + 1) if a % 3) <= 2 * X / 3 + 1, "count at X=%d" % X)
out("A4 #{a <= X : 3 !| a} <= 2X/3 + 1 for X < 3000; capacity B/A >= 9/2 follows as in the note.")
# the two negative cycles: 3^L/2^K
check(Fraction(3, 2) > 1 and Fraction(9, 8) > 1 and Fraction(3, 4) < 1, "cycle ratios")
out("A4 cycle ratios 3^L/2^K: {-1}: 3/2, {-4,-11}: 9/8, {1}: 3/4.")
out("[time] A4 %.1fs" % (time.time() - T0))
out("[done] audit total %.1fs" % (time.time() - T0))

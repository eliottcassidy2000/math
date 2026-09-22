#!/usr/bin/env python3
"""collatz_mod6_20260921_g_negatives_joint_carry.py

Lane g_negatives_joint_carry (session collatz-mod6-20260921).

S1  G(m) = (2^k m - 1)/3, k >= 0 minimal with 2^k m in {4,7} mod 9, on NEGATIVE
    integers coprime to 3: exact census of the cycles reached from every
    -10^7 <= m <= -1 (escape bound 10^18), basin sizes, worst step count and peak.
S1b Cycle-gate census m_0 (2^J - 3^L) = B'(w): every word w in {0..3}^L (L <= 10)
    and {0..6}^L (L <= 7) with q = 1 (integral chain), split into greedy-legal
    (G-cycles) and all-odd (reversed T-cycles), both signs; mirror table.
S2  Joint carry for odd n <= 10^5: forward T-word n -> orbit maximum M
    (2^{K1} M = 3^{L1} n + B1) and greedy G-word M -> 1 (2^{K2} M = 3^J + B2);
    conservation identity; overlap of the greedy path with the forward orbit;
    residue criterion for a full retrace; invariant search on (B1,B2,K1,K2,L1,J).
S3  G-analogue of the bounded-strip theorem: numerical anchors for the proof in
    the note (positivity forces 3^j/2^{K_j} -> 0 on the positive side; on the
    negative side the guards 2a transfer needs the word law mod 3^{J+1}, the
    moment 3(7/3)^{J-1}, the intercept bound and the (7/9)^{J-1} tail).

All checks use explicit raise.  python3 -O gives identical output modulo timing.
"""
import sys, time, hashlib, math
from array import array
from fractions import Fraction
from collections import Counter, defaultdict

T0 = time.time()

def out(s=""):
    print(s)

def tline(tag):
    out("[time] %s %.1fs" % (tag, time.time() - T0))

def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)

# greedy letters: G(m) = (2^k m - 1)/3 with 2^k m in {4,7} mod 9 (k minimal)
KPOS = {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}
# G_-(p) = (2^k p + 1)/3 with 2^k p in {2,5} mod 9 (k minimal); G(-p) = -G_-(p)
KNEG = {1: 1, 2: 0, 4: 3, 5: 0, 7: 1, 8: 2}

def G(m):
    r = m % 9
    if r % 3 == 0:
        raise ValueError("G undefined on multiples of 3: %d" % m)
    k = KPOS[r]
    y = (m << k) - 1
    if y % 3 != 0:
        raise RuntimeError("non-integral G step at %d" % m)
    return y // 3, k

def Gneg(p):
    r = p % 9
    if r % 3 == 0:
        raise ValueError("G_- undefined on multiples of 3: %d" % p)
    k = KNEG[r]
    y = (p << k) + 1
    if y % 3 != 0:
        raise RuntimeError("non-integral G_- step at %d" % p)
    return y // 3, k

src = open(__file__, "rb").read()
out("# collatz_mod6_20260921_g_negatives_joint_carry")
out("source sha256 %s" % hashlib.sha256(src).hexdigest())
out()

# ---------------------------------------------------------------- S1
out("## S1. G on negative integers coprime to 3: cycles and basins")
# S1.1 the residue table and the negation conjugacy
for r in (1, 2, 4, 5, 7, 8):
    check(KNEG[r] == KPOS[(-r) % 9], "conjugacy table r=%d" % r)
for p in range(1, 20001):
    if p % 3 == 0:
        continue
    g1, k1 = G(-p)
    g2, k2 = Gneg(p)
    check(g1 == -g2 and k1 == k2, "G(-p) = -G_-(p) at p=%d" % p)
out("S1.1 G(-p) = -G_-(p) with the same k for all p <= 20000, 3 !| p; k-tables "
    "k_+(r)=%s, k_-(r)=%s, k_-(r) = k_+(-r mod 9)." % (KPOS, KNEG))
# S1.2 the two cycles, by direct evaluation on negative m
c1 = G(-1)
c4 = G(-4)
c11 = G(-11)
out("S1.2 G(-1) = %d (k=%d); G(-4) = %d (k=%d); G(-11) = %d (k=%d)." %
    (c1[0], c1[1], c4[0], c4[1], c11[0], c11[1]))
check(c1 == (-1, 1) and c4 == (-11, 3) and c11 == (-4, 0), "the two negative cycles")
# gate check: m_0 (2^J - 3^L) = B'(w) with B' = sum 3^(i-1) 2^(J-J_i)
def Bprime(word):
    J = sum(word)
    Ji = 0
    B = 0
    for i, j in enumerate(word, 1):
        Ji += j
        B += 3 ** (i - 1) * 2 ** (J - Ji)
    return B
for m0, word in ((-1, (1,)), (-4, (3, 0)), (-11, (0, 3)), (1, (2,))):
    J, L = sum(word), len(word)
    check(m0 * (2 ** J - 3 ** L) == Bprime(word), "gate at m0=%d" % m0)
    out("S1.2 gate: m0=%d word=%s (J,L)=(%d,%d) 2^J-3^L=%d B'=%d" %
        (m0, word, J, L, 2 ** J - 3 ** L, Bprime(word)))

# S1.3 census on p in [1, 10^7] for G_- (= -G on [-10^7, -1]); resolve by drop-below
N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 7
STEPCAP = 10 ** 5
ESC = 10 ** 18
memo = bytearray(N1 + 1)        # 0 unresolved / multiple of 3, 1 -> {-1}, 2 -> {-4,-11}
tot = array("I", [0]) * (N1 + 1)  # total G_- steps to first hit of the cycle set
pk = array("Q", [0]) * (N1 + 1)   # orbit peak (uint64; raise if larger)
memo[1] = 1
memo[4] = 2
memo[11] = 2
tot[1] = 0
tot[4] = 0
tot[11] = 0
pk[1] = 1
pk[4] = 11
pk[11] = 11
escapes = []
newcycle = []
basin = [0, 0, 0]
maxsteps = 0
argsteps = []
maxpeak = 0
argpeak = 0
maxratio = Fraction(0)
argratio = 0
pow10_marks = {10 ** e: None for e in range(1, 8)}
basin_by_decade = []
kn = KNEG
for p in range(1, N1 + 1):
    if p % 3 == 0:
        continue
    if memo[p]:
        basin[memo[p]] += 1
    else:
        x = p
        steps = 0
        peak = p
        while True:
            k = kn[x % 9]
            x = ((x << k) + 1) // 3
            steps += 1
            if x > peak:
                peak = x
            if x < p:
                c = memo[x]
                steps += tot[x]
                if pk[x] > peak:
                    peak = pk[x]
                break
            if x == p:
                c = 3
                newcycle.append(p)
                break
            if x > ESC:
                c = 0
                escapes.append(p)
                break
            if steps > STEPCAP:
                c = 3
                newcycle.append(p)
                break
        if c == 3 or c == 0:
            memo[p] = 0
            continue
        memo[p] = c
        tot[p] = steps
        if peak >= (1 << 64):
            raise RuntimeError("peak overflow at p=%d" % p)
        pk[p] = peak
        basin[c] += 1
        if steps > maxsteps:
            maxsteps = steps
            argsteps = [p]
        elif steps == maxsteps:
            argsteps.append(p)
        if peak > maxpeak:
            maxpeak = peak
            argpeak = p
        fr = Fraction(peak, p)
        if fr > maxratio:
            maxratio = fr
            argratio = p
    if p in pow10_marks:
        basin_by_decade.append((p, basin[1], basin[2]))
tline("S1.3 census")
check(not escapes, "escapes: %s" % escapes[:5])
check(not newcycle, "new cycles: %s" % newcycle[:5])
out("S1.3 FINITE-EXACT: every m in [-10^7, -1] with 3 !| m reaches {-1} or {-4,-11} under G;"
    " escapes to |m| > 10^18: %d; new cycles: %d." % (len(escapes), len(newcycle)))
out("S1.3 basins: to {-1}: %d; to {-4,-11}: %d; total %d (= %d non-multiples of 3)." %
    (basin[1], basin[2], basin[1] + basin[2], N1 - N1 // 3))
out("S1.3 basin fraction to {-4,-11}: %.6f" % (basin[2] / (basin[1] + basin[2])))
out("S1.3 basin counts by decade (X, to {-1}, to {-4,-11}, fraction to {-4,-11}):")
for X, b1, b2 in basin_by_decade:
    out("      X=%d  %d  %d  %.6f" % (X, b1, b2, b2 / (b1 + b2)))
out("S1.3 max steps to the cycle set: %d at m = %s" % (maxsteps, [-p for p in argsteps]))
out("S1.3 max orbit peak: %d at m = -%d; max peak/|m| = %d/%d = %.3f at m = -%d" %
    (maxpeak, argpeak, maxratio.numerator, maxratio.denominator, float(maxratio), argratio))
# print the two worst words
def gneg_word(p):
    w = []
    x = p
    seq = [x]
    while x not in (1, 4, 11):
        x, k = Gneg(x)
        w.append(k)
        seq.append(x)
    return w, seq
w, seq = gneg_word(argsteps[0])
out("S1.3 word of m=-%d (%d letters): %s" % (argsteps[0], len(w), "".join(map(str, w))))
w, seq = gneg_word(argratio)
out("S1.3 orbit of m=-%d (negated): %s" % (argratio, [-v for v in seq]))
out("S1.3 word of m=-%d: %s" % (argratio, "".join(map(str, w))))
# compare with the positive-side worst ratio of the three-adic lane (m = 4847486)
xp = 4847486
x = xp
peakp = x
Kp = 0
wp = []
while x != 1:
    x, k = G(x)
    wp.append(k)
    if x > peakp:
        peakp = x
out("S1.3 positive control: G-orbit of %d peaks at %d, peak/m = %d/%d = %.3f, word %s" %
    (xp, peakp, Fraction(peakp, xp).numerator, Fraction(peakp, xp).denominator, peakp / xp, "".join(map(str, wp))))
out("S1.3 the two worst ratios: negative %s, positive %s, equal: %s; 2^34/3^17 = %.3f" %
    (maxratio, Fraction(peakp, xp), maxratio == Fraction(peakp, xp), 2 ** 34 / 3 ** 17))
# S1.4 direct check on the negative side for |m| <= 10^5 (no conjugation)
for p in range(1, 100001):
    if p % 3 == 0:
        continue
    x = -p
    steps = 0
    while x not in (-1, -4, -11):
        x, k = G(x)
        steps += 1
        if steps > 10 ** 6:
            raise RuntimeError("direct negative orbit too long at %d" % (-p))
    c = 1 if x == -1 else 2
    check(c == memo[p], "direct/conjugate basin mismatch at m=%d" % (-p))
out("S1.4 direct evaluation of G on m in [-10^5, -1] reproduces the conjugated basin labels.")
# small basin listing
out("S1.5 basin of {-4,-11} below 100: %s" %
    [-p for p in range(1, 100) if p % 3 and memo[p] == 2])
out("S1.5 powers of two: -2^e for e=0..12 go to %s" %
    [(e, "{-4,-11}" if memo[2 ** e] == 2 else "{-1}") for e in range(0, 13)])
del tot, pk
tline("S1 done")
out()

# ---------------------------------------------------------------- S1b
out("## S1b. Cycle-gate census: q = 1 words, greedy-legal (G) and all-odd (T), both signs")

def gate_census(maxL, maxletter):
    """(J,L) -> [set of primitive E-cycles (frozenset of nodes), greedy-legal subset,
    all-odd subset]; every word in {0..maxletter}^L, 1 <= L <= maxL."""
    res = defaultdict(lambda: [set(), set(), set(), set()])
    pow3 = [3 ** i for i in range(maxL + 2)]
    nwords = [0]
    def rec(word, J, Bp, L):
        if L > 0:
            nwords[0] += 1
            D = (1 << J) - pow3[L]
            if D != 0 and Bp % D == 0:
                m0 = Bp // D
                x = m0
                ok = True
                greedy = True
                allodd = True
                cyc = []
                for j in word:
                    if x % 3 == 0:
                        ok = False
                        break
                    cyc.append(x)
                    y = (x << j) - 1
                    if y % 3 != 0:
                        ok = False
                        break
                    if KPOS[x % 9] != j:
                        greedy = False
                    if x % 2 == 0:
                        allodd = False
                    x = y // 3
                if ok:
                    check(x == m0, "gate closure failed for word %s" % (word,))
                    if len(set(cyc)) == len(cyc):
                        key = frozenset(cyc)
                        # E-nodes of the closed walk: m_i, 3m_i+1, and its halvings
                        # E-nodes of the closed walk: m_i, 2m_i, ..., 2^{j_i} m_i = 3 m_{i+1} + 1
                        enodes = []
                        for i2, j2 in enumerate(word):
                            for t in range(j2 + 1):
                                enodes.append(cyc[i2] << t)
                            check(3 * cyc[(i2 + 1) % L] + 1 == cyc[i2] << j2, "E-walk arrow at word %s" % (word,))
                        check(len(enodes) == J + L, "E-walk length")
                        simple = (len(set(enodes)) == len(enodes))
                        res[(J, L)][0].add(key)
                        if simple:
                            res[(J, L)][3].add(key)
                        if greedy:
                            res[(J, L)][1].add(key)
                        if allodd:
                            res[(J, L)][2].add(key)
        if L == maxL:
            return
        for j in range(maxletter + 1):
            rec(word + (j,), J + j, (Bp << j) + pow3[L], L + 1)
    rec((), 0, 0, 0)
    return res, nwords[0]

for (maxL, maxletter) in ((10, 3), (7, 6)):
    res, nw = gate_census(maxL, maxletter)
    out("S1b census: letters 0..%d, 1 <= L <= %d, %d words scanned." % (maxletter, maxL, nw))
    out("  (J,L) | sign(2^J-3^L) | 2^J-3^L | #q=1 closed E-walks (simple E-cycles) | #greedy-legal (G-cycles) | "
        "#all-odd (reversed T-cycles) | greedy cycles | all-odd cycles")
    tot_g = Counter()
    tot_t = Counter()
    for (J, L) in sorted(res, key=lambda t: (t[1], t[0])):
        E, Gs, Ts, Ss = res[(J, L)]
        D = 2 ** J - 3 ** L
        sg = "+" if D > 0 else "-"
        gl = sorted(sorted(c) for c in Gs)
        tl = sorted(sorted(c) for c in Ts)
        out("  (%d,%d) | %s | %d | %d (simple %d) | %d | %d | %s | %s" %
            (J, L, sg, D, len(E), len(Ss), len(Gs), len(Ts), gl, tl))
        if len(E) <= 2:
            out("      q=1 closed walks at (%d,%d): %s" % (J, L, sorted(sorted(c) for c in E)))
        tot_g[sg] += len(Gs)
        tot_t[sg] += len(Ts)
    out("  totals: G-cycles + %d, - %d; reversed T-cycles + %d, - %d" %
        (tot_g["+"], tot_g["-"], tot_t["+"], tot_t["-"]))
    allG = set()
    for key in res:
        for c in res[key][1]:
            allG.add(tuple(sorted(c)))
    check(allG == {(1,), (-1,), (-11, -4)}, "greedy-legal cycles are not the three universal ones: %s" % allG)
    # every found cycle satisfies the sign law
    for (J, L) in res:
        D = 2 ** J - 3 ** L
        for c in res[(J, L)][0]:
            sgn = 1 if min(c) > 0 else -1
            check(sgn == (1 if D > 0 else -1), "sign law failed at (J,L)=(%d,%d)" % (J, L))
    out("  greedy-legal q=1 cycles are exactly {1}, {-1}, {-4,-11}; the sign law "
        "sign(m0) = sign(2^J-3^L) holds on every q=1 cycle.")
tline("S1b done")
out()

# ---------------------------------------------------------------- S1c (added at the audit of 2026-09-22)
out("## S1c. Global gate census by E-cycle length n = J + L <= 21, ALL letters (compositions of J into L parts)")
import itertools
def compositions(J, L):
    for cuts in itertools.combinations(range(J + L - 1), L - 1):
        prev = -1
        w = []
        for c in cuts:
            w.append(c - prev - 1)
            prev = c
        w.append(J + L - 1 - prev - 1)
        yield tuple(w)
bylen = Counter()
row_info = {}
extra = []
for n in range(1, 22):
    for L in range(1, n + 1):
        J = n - L
        D = (1 << J) - 3 ** L
        seen = {}
        for w in compositions(J, L):
            A, Bn = 1, 0
            for s_, k in enumerate(w):
                A <<= k
                Bn = (Bn << k) + 3 ** s_
            if Bn % D:
                continue
            x = m0 = Bn // D
            nodes, enodes, ok, greedy, allodd = [], [], True, True, True
            for k in w:
                if x % 3 == 0 or ((x << k) - 1) % 3:
                    ok = False
                    break
                nodes.append(x)
                enodes.extend(x << t for t in range(k + 1))
                if KPOS[x % 9] != k:
                    greedy = False
                if x % 2 == 0:
                    allodd = False
                x = ((x << k) - 1) // 3
            if not ok or len(set(nodes)) < len(nodes):
                continue
            check(x == m0, "S1c closure at %s" % (w,))
            key = frozenset(nodes)
            if key not in seen:
                seen[key] = (len(set(enodes)) == len(enodes), greedy, allodd, max(abs(e) for e in enodes), w)
        if seen:
            ns = sum(1 for v in seen.values() if v[0])
            ng = sum(1 for v in seen.values() if v[1])
            nt = sum(1 for v in seen.values() if v[2])
            sg = "+" if D > 0 else "-"
            bylen[(n, sg)] += ns
            row_info[(J, L)] = (len(seen), ns, ng, nt, max(max(v[4]) for v in seen.values()), max(v[3] for v in seen.values()))
            out("  n=%d (J,L)=(%d,%d) D=%d: %d primitive integral walks, %d simple, %d greedy, %d all-odd, max letter %d, max |E-node| %d" %
                ((n, J, L, D) + row_info[(J, L)]))
            for c, v in seen.items():
                if v[1]:
                    extra.append(tuple(sorted(c)))
                if (J, L) == (11, 7) and (max(v[4]) > 6 or v[3] > 2000):
                    out("      (11,7) cycle beyond the letters-0..6 scan and the nodes<=2000 census: %s word %s max|E-node| %d" %
                        (sorted(c), v[4], v[3]))
                if (J, L) == (13, 8) and v[0]:
                    out("      (13,8) simple cycle %s word %s max letter %d" % (sorted(c), v[4], max(v[4])))
out("S1c simple E-cycles (positive sheet) by length: %s" % sorted((n, c) for (n, sg), c in bylen.items() if sg == "+" and c))
out("S1c simple E_- cycles (negative sheet) by length: %s" % sorted((n, c) for (n, sg), c in bylen.items() if sg == "-" and c))
check({n: c for (n, sg), c in bylen.items() if sg == "+" and c} == {3: 1, 8: 1, 13: 6, 16: 1, 21: 2}, "S1c E counts")
check({n: c for (n, sg), c in bylen.items() if sg == "-" and c} == {2: 1, 5: 2, 15: 2, 18: 13}, "S1c E_- counts")
check(set(extra) == {(1,), (-1,), (-11, -4)}, "S1c greedy cycles")
check(row_info[(11, 7)] == (13, 13, 0, 1, 8, 2732) and row_info[(13, 8)] == (4, 2, 0, 0, 6, 256), "S1c rows")
out("S1c the letters-0..6 scan of S1b sees 12 of the 13 simple E_- cycles of length 18; the inherited bounded census "
    "(nodes <= 2000) also sees 12; the 13th needs the letter 8 and the E-node 2732. Greedy-legal cycles: still exactly "
    "{1}, {-1}, {-4,-11} (all letters, all lengths <= 21).")
tline("S1c done")
out()

# ---------------------------------------------------------------- S2
out("## S2. Joint carry: forward T-word to the peak, greedy G-word from the peak")
N2 = 10 ** 5
GCAP = 5000

def fwd_odd(n):
    seq = [n]
    ks = []
    while n != 1:
        y = 3 * n + 1
        k = (y & -y).bit_length() - 1
        n = y >> k
        seq.append(n)
        ks.append(k)
    return seq, ks

def greedy_path(M):
    seq = [M]
    ks = []
    x = M
    while x != 1:
        x, k = G(x)
        seq.append(x)
        ks.append(k)
        if len(ks) > GCAP:
            raise RuntimeError("greedy path from %d exceeds cap" % M)
    return seq, ks

def crit(p, k):
    # residue criterion for G(u) = p where u = (3p+1)/2^k
    r = p % 3
    return (r == 1 and k <= 3) or (r == 2 and k == 1)

records = []          # (n, M, L1, K1, B1, J, K2, B2, overlap)
skipped_mult3 = []
L1zero = 0
cnt_case = Counter()  # K1>=K2 vs K1<K2
overlap_hist = Counter()
gap_hist = Counter()
full = 0
full_by_L1 = Counter()
n_by_L1 = Counter()
crit_checks = 0
maxM = 0
argM = 0
peak_checks = 0
for n in range(1, N2 + 1, 2):
    seq, ks = fwd_odd(n)
    M = max(seq)
    L1 = seq.index(M)
    if M > maxM:
        maxM = M
        argM = n
    K1 = sum(ks[:L1])
    B1 = (M << K1) - 3 ** L1 * n
    # formula check
    Bf = 0
    Kp = 0
    for i in range(1, L1 + 1):
        Bf += 3 ** (L1 - i) * 2 ** Kp
        Kp += ks[i - 1]
    check(Bf == B1, "B1 formula at n=%d" % n)
    if M % 3 == 0:
        skipped_mult3.append(n)
        continue
    gseq, gks = greedy_path(M)
    J = len(gks)
    K2 = sum(gks)
    B2 = (M << K2) - 3 ** J
    Bg = 0
    Ks = 0
    for s in range(1, J + 1):
        Ks += gks[s - 1]
        Bg += 3 ** (s - 1) * 2 ** (K2 - Ks)
    check(Bg == B2, "B2 formula at n=%d" % n)
    # conservation identity
    if K1 >= K2:
        check(3 ** L1 * n + B1 == 2 ** (K1 - K2) * (3 ** J + B2), "identity (K1>=K2) at n=%d" % n)
        cnt_case["K1>=K2"] += 1
    else:
        check(2 ** (K2 - K1) * (3 ** L1 * n + B1) == 3 ** J + B2, "identity (K1<K2) at n=%d" % n)
        cnt_case["K1<K2"] += 1
    # overlap of the greedy path with the pre-peak forward orbit
    ov = 0
    while ov < L1 and ov < J and gseq[ov + 1] == seq[L1 - ov - 1]:
        ov += 1
    # residue criterion at every pre-peak step
    for i in range(L1, 0, -1):
        p = seq[i - 1]
        k = ks[i - 1]
        u = seq[i]
        actual = (G(u)[0] == p)
        check(actual == crit(p, k), "residue criterion failed at n=%d step %d" % (n, i))
        crit_checks += 1
    # the criterion predicts the overlap length
    pred = 0
    for i in range(L1, 0, -1):
        if crit(seq[i - 1], ks[i - 1]):
            pred += 1
        else:
            break
    check(pred == ov, "overlap prediction at n=%d" % n)
    if L1 == 0:
        L1zero += 1
    overlap_hist[ov] += 1
    gap_hist[L1 - ov] += 1
    n_by_L1[L1] += 1
    if ov == L1:
        full += 1
        full_by_L1[L1] += 1
    records.append((n, M, L1, K1, B1, J, K2, B2, ov, tuple(gks[-2:]), gseq[-2] if J >= 2 else None))
    if L1 >= 1:
        check(ks[L1 - 1] == 1, "k_L1 = 1 at n=%d" % n)
        check(M % 12 == 5, "M = 5 mod 12 at n=%d" % n)
        check(ks[L1] >= 2, "k_(L1+1) >= 2 at n=%d" % n)
        check(gks[0] in (1, 3), "first greedy letter at n=%d" % n)
        check(K2 >= 2, "K2 >= 2 at n=%d" % n)
        tz = 0
        while tz < J and gks[J - 1 - tz] == 0:
            tz += 1
        check(tz % 2 == 0, "even trailing zeros at n=%d" % n)
        check(B2 % 2 == 1, "B2 odd at n=%d" % n)
        check((ov >= 1) == (not (L1 == 1 and n % 3 == 0)), "overlap>=1 rule at n=%d" % n)
        peak_checks += 1
tline("S2 loop")
out("S2.1 odd n <= %d: %d records; %d skipped because M = n is a multiple of 3 (G undefined): first %s" %
    (N2, len(records), len(skipped_mult3), skipped_mult3[:8]))
out("S2.1 largest peak M = %d at n = %d; every greedy path from a peak reached 1 (cap %d)." % (maxM, argM, GCAP))
out("S2.1 carry formulas verified: B1 = sum_{i<=L1} 3^(L1-i) 2^(k_1+..+k_(i-1)), "
    "B2 = sum_{s<=J} 3^(s-1) 2^(K2-K_s), on every record.")
out("S2.2 conservation identity verified on every record; cases: %s" % dict(cnt_case))
out("S2.3 residue criterion G(u) = p  <=>  (p = 1 mod 3 and k <= 3) or (p = 2 mod 3 and k = 1), "
    "u = (3p+1)/2^k, checked at %d pre-peak steps (all agree)." % crit_checks)
# mod-48 form of the criterion
S48 = sorted(p for p in range(48) if p % 2 == 1 and p % 3 and
             ((p % 3 == 1 and p % 16 != 5) or (p % 12 == 11)))
for p in range(1, 20000, 2):
    if p % 3 == 0:
        continue
    y = 3 * p + 1
    k = (y & -y).bit_length() - 1
    check(crit(p, k) == ((p % 48) in S48), "mod-48 form at p=%d" % p)
out("S2.3 equivalently p mod 48 in %s (%d of the 16 odd classes coprime to 3), "
    "i.e. p = 1 mod 3 with p != 5 mod 16, or p = 11 mod 12; verified for odd p < 20000." % (S48, len(S48)))
out("S2.4 records with L1 = 0 (n is its own peak): %d; fully retraced (overlap = L1): %d of %d = %.4f; "
    "fully retraced with L1 >= 1: %d of %d = %.4f" %
    (L1zero, full, len(records), full / len(records), full - L1zero, len(records) - L1zero,
     (full - L1zero) / (len(records) - L1zero)))
out("S2.4 overlap histogram (overlap: count): %s" % sorted(overlap_hist.items()))
out("S2.4 gap L1 - overlap histogram (gap: count): %s" % sorted(gap_hist.items())[:25])
out("S2.4 full-retrace fraction by L1 (L1: full/count):")
for L1 in sorted(n_by_L1):
    if L1 <= 12 or full_by_L1[L1]:
        out("      L1=%d: %d/%d = %.4f  ((11/16)^L1 = %.4f)" %
            (L1, full_by_L1[L1], n_by_L1[L1], full_by_L1[L1] / n_by_L1[L1], (11 / 16) ** L1))
# mean overlap
mean_ov = sum(r[8] for r in records) / len(records)
mean_L1 = sum(r[2] for r in records) / len(records)
out("S2.4 mean overlap %.4f, mean L1 %.4f" % (mean_ov, mean_L1))
out("S2.7 peak structure verified on %d records with L1 >= 1: k_L1 = 1; M = 5 mod 12; k_(L1+1) >= 2; "
    "first greedy letter in {1,3}; K2 >= 2; number of trailing zero letters of the greedy word is even; "
    "B2 odd; overlap >= 1 iff not (L1 = 1 and 3 | n)." % peak_checks)
chain = [(3 ** a - 1) // 2 for a in range(1, 13)]
out("S2.7 zero-letter chain (3^a-1)/2, a=1..12: %s; residues mod 8: %s" % (chain, [c % 8 for c in chain]))
check(all(((3 ** a - 1) // 2) % 8 != 3 for a in range(1, 40)), "chain never 3 mod 8")
check(all(((3 ** (a + 1) - 1) // 2) % 2 == ((a + 1) % 2) for a in range(1, 40)), "chain parity law")
out("S2.7 (3^a-1)/2 is never 3 mod 8 (a < 40), and (3^(a+1)-1)/2 is even iff a is odd (a < 40).")

# S2.5 invariant search
out("S2.5 invariant search on (B1,B2,K1,K2,L1,J), records with L1 >= 1 unless stated")
rec1 = [r for r in records if r[2] >= 1]
def v2(x):
    return (x & -x).bit_length() - 1 if x else -1
def v3(x):
    if x == 0:
        return -1
    c = 0
    while x % 3 == 0:
        x //= 3
        c += 1
    return c
# (a) exact valuation identities
for (n, M, L1, K1, B1, J, K2, B2, ov, _l2, _pn) in rec1:
    check(v2(3 ** J + B2) == K2, "v2(3^J+B2)=K2 at n=%d" % n)
    check(v2(3 ** L1 * n + B1) == K1, "v2(3^L1 n + B1)=K1 at n=%d" % n)
    check(B1 % 2 == 1, "B1 odd at n=%d" % n)
    check(B1 % 3 != 0 and B2 % 3 != 0, "B1,B2 units mod 3 at n=%d" % n)
    check((B1 * B2 - (-1) ** (K1 + K2)) % 3 == 0, "B1 B2 = (-1)^(K1+K2) mod 3 at n=%d" % n)
    check(v3(2 ** K2 * M - B2) == J, "v3 at n=%d" % n)
out("S2.5a hold on all records: B1 odd; 3 !| B1 B2; v2(3^J + B2) = K2; v2(3^L1 n + B1) = K1; "
    "B1 B2 = (-1)^(K1+K2) mod 3; v3(2^K2 M - B2) = J.")
# (b) constancy tests of residues
exprs = {
    "B1": lambda r: r[4], "B2": lambda r: r[7], "B1+B2": lambda r: r[4] + r[7],
    "B1-B2": lambda r: r[4] - r[7], "B1*B2": lambda r: r[4] * r[7],
    "B1+3^L1": lambda r: r[4] + 3 ** r[2], "B2+3^J": lambda r: r[7] + 3 ** r[5],
    "B2-2^K2": lambda r: r[7] - 2 ** r[6], "B1-2^K1": lambda r: r[4] - 2 ** r[3],
}
consts = []
for name, f in exprs.items():
    for q in range(2, 49):
        vals = set()
        for r in rec1:
            vals.add(f(r) % q)
            if len(vals) > 1:
                break
        if len(vals) == 1:
            consts.append((name, q, vals.pop()))
out("S2.5b constant residues (expr, modulus, value) among moduli 2..48: %s" % consts)
# (c) residues functionally determined by exponent data
def functional(keyf, valf, label):
    table = {}
    for r in rec1:
        k = keyf(r)
        v = valf(r)
        if k in table and table[k] != v:
            return False
        table[k] = v
    return True
tests = [
    ("B1 mod 3 <- (K1-K2 mod 2, B2 mod 3)", lambda r: ((r[3] - r[6]) % 2, r[7] % 3), lambda r: r[4] % 3),
    ("B1 mod 9 <- (K1-K2 mod 6, B2 mod 9)", lambda r: ((r[3] - r[6]) % 6, r[7] % 9), lambda r: r[4] % 9),
    ("B1 mod 9 <- (K1,K2,L1,J mod 6, B2 mod 9)", lambda r: (r[3] % 6, r[6] % 6, r[2] % 6, r[5] % 6, r[7] % 9), lambda r: r[4] % 9),
    ("B2 mod 4 <- (K1,K2,L1,J mod 4, B1 mod 4)", lambda r: (r[3] % 4, r[6] % 4, r[2] % 4, r[5] % 4, r[4] % 4), lambda r: r[7] % 4),
    ("B2 mod 8 <- (K2 mod 2, J mod 2)", lambda r: (r[6] % 2, r[5] % 2), lambda r: r[7] % 8),
    ("B2 mod 2 <- (K2,J) mod 2", lambda r: (r[6] % 2, r[5] % 2), lambda r: r[7] % 2),
    ("B1 mod 4 <- (K1,L1) mod 2", lambda r: (r[3] % 2, r[2] % 2), lambda r: r[4] % 4),
    ("B1 mod 4 <- (K1,L1) mod 4", lambda r: (r[3] % 4, r[2] % 4), lambda r: r[4] % 4),
]
for label, kf, vf in tests:
    out("S2.5c %s : %s" % (label, "DETERMINED" if functional(kf, vf, label) else "not determined"))
# (d) linear forms in the exponents: exact ranges
forms = {}
for a in range(-2, 3):
    for b in range(-2, 3):
        for c in range(-2, 3):
            for d in range(-2, 3):
                if (a, b, c, d) == (0, 0, 0, 0):
                    continue
                lo = min(a * r[3] + b * r[6] + c * r[2] + d * r[5] for r in rec1)
                hi = max(a * r[3] + b * r[6] + c * r[2] + d * r[5] for r in rec1)
                forms[(a, b, c, d)] = (lo, hi)
zero_range = [k for k, v in forms.items() if v[0] == v[1]]
out("S2.5d linear forms a K1 + b K2 + c L1 + d J, |coeff| <= 2: constant ones: %s" % zero_range)
for key in ((1, -1, 0, 0), (0, 1, 0, -1), (0, 0, -1, 1), (1, 0, -2, 0), (1, 0, -1, 0), (0, 1, 0, -2), (1, -1, -1, 1)):
    out("S2.5d range of %s . (K1,K2,L1,J): [%d, %d]" % (key, forms[key][0], forms[key][1]))
# (e) gcd and sign structure
g_hist = Counter(math.gcd(r[4], r[7]) for r in rec1)
out("S2.5e gcd(B1,B2) = 1 on %d of %d records; other values (gcd: count, first 10): %s" %
    (g_hist[1], len(rec1), sorted((g, c) for g, c in g_hist.items() if g > 1)[:10]))
out("S2.5e last two greedy letters (k_(J-1), k_J) and penultimate node, counts: %s" %
    sorted(Counter((r[9], r[10]) for r in rec1 if r[5] >= 2).items()))
# sample rows
out("S2.6 sample rows (n, M, L1, K1, B1, J, K2, B2, overlap):")
for r in records:
    if r[0] in (1, 3, 7, 27, 97, 871, 6171, 77031):
        out("      %s" % (r,))
tline("S2 done")
out()

# ---------------------------------------------------------------- S3
out("## S3. G-analogue of the bounded-strip theorem: numerical anchors")
# S3.1 positive side: sum_{s>=1} 3^s/2^{K_s} <= 3 m_0 for every positive orbit (positivity)
worst = Fraction(0)
argw = 0
for m in range(1, 5001):
    if m % 3 == 0:
        continue
    x = m
    Ks = 0
    S = Fraction(0)
    steps = 0
    while steps < 300:
        x, k = G(x)
        Ks += k
        steps += 1
        S += Fraction(3 ** steps, 2 ** Ks)
    check(S < 3 * m, "positivity sum bound at m=%d" % m)
    fr = S / m
    if fr > worst:
        worst = fr
        argw = m
out("S3.1 sum_{s=1}^{300} 3^s/2^{K_s} < 3 m_0 verified for every positive m_0 <= 5000; "
    "largest ratio (sum)/m_0 = %.6f at m_0 = %d (bound 3)." % (float(worst), argw))
# S3.2 negative side (G_- on positives): word law mod 3^(J+1), moment, intercept, tail
def gneg_prefix(x, J):
    ks = []
    B = 0
    K = 0
    for s in range(J):
        x, k = Gneg(x)
        ks.append(k)
    return ks
def gneg_KB(x, J):
    K = 0
    B = 0
    for s in range(1, J + 1):
        x, k = Gneg(x)
        K += k
        B = (B << k) + 3 ** (s - 1)   # B_s = 2^k B_(s-1) + 3^(s-1)
    return K, B
for J in range(1, 7):
    mod = 3 ** (J + 1)
    for a in range(1, mod):
        if a % 3 == 0:
            continue
        w0 = gneg_prefix(a, J)
        for t in (1, 2, 5):
            check(gneg_prefix(a + t * mod, J) == w0, "word law mod 3^(J+1) fails J=%d a=%d" % (J, a))
out("S3.2a G_- words of length J are functions of m mod 3^(J+1) for J <= 6 (all unit classes, lifts t=1,2,5).")
rows = []
for J in range(1, 9):
    mod = 3 ** (J + 1)
    units = [a for a in range(1, mod) if a % 3]
    E = Fraction(0)
    bad = 0
    intercept_ok = True
    for a in units:
        K, B = gneg_KB(a, J)
        E += Fraction(2 ** K)
        if 2 ** K >= Fraction(3 ** J, 2):
            bad += 1
        if not (B <= 2 ** K * (3 ** J - 1) // 2):
            intercept_ok = False
    E /= len(units)
    check(E == 3 * Fraction(7, 3) ** (J - 1), "moment at J=%d" % J)
    check(intercept_ok, "intercept bound at J=%d" % J)
    tail = Fraction(bad, len(units))
    bound = 2 * Fraction(7, 9) ** (J - 1)
    check(tail <= bound, "tail bound at J=%d" % J)
    rows.append((J, str(E), bad, len(units), float(tail), float(bound)))
out("S3.2b G_-: E_Haar[2^{K_J}] = 3(7/3)^(J-1) exactly, B_J <= 2^{K_J}(3^J-1)/2 on every class, and "
    "Haar{2^{K_J} >= 3^J/2} <= 2(7/9)^(J-1), J <= 8:")
out("      J | E[2^K_J] | #classes with 2^K_J >= 3^J/2 | #units mod 3^(J+1) | tail | 2(7/9)^(J-1)")
for r in rows:
    out("      %d | %s | %d | %d | %.6f | %.6f" % r)
# S3.3 sign law for cycles (from S1b) and the strip-forcing ratio on the two known negative cycles
for m0, word in ((-1, (1,)), (-4, (3, 0))):
    J, L = sum(word), len(word)
    out("S3.3 negative cycle from %d: (K,L)=(%d,%d), 3^L/2^K = %s > 1, so 3^j/2^{K_j} -> infinity along it." %
        (m0, J, L, Fraction(3 ** L, 2 ** J)))
out("S3.3 positive cycle {1}: (K,L)=(2,1), 3^L/2^K = 3/4 < 1, ratio -> 0.")
# S3.4 capacity constant for G_- (integers coprime to 3): B/A >= 9/2
out("S3.4 capacity bound for a hypothetical G_- strip A <= 3^j/2^{K_j} <= B: images are coprime to 3, "
    "count <= 2X/3 + 1, giving B/A >= 9/2 (Collatz odd-and-coprime count gives 9).")
tline("S3 done")
out()
out("[done] total %.1fs" % (time.time() - T0))

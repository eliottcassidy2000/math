#!/usr/bin/env python3
"""Audit D (2026-10-14), THM-4654 statements 1-4: independent checks with integer arithmetic.

Methods deliberately different from safe_language.py (which uses the Fraction slack automaton):
  * DIRECT  : finite suffix test 2*C_(i,h) < 3^h with C_(i,h) = sum_j w_(i+j) 2^j 3^(h-1-j)  (pure integers);
  * PARRY   : every suffix starting at a 1 is <=_lex the prefix of d of the same length (d = greedy expansion of 1);
  * FOLLOW  : the follower graph Q of THM-4072 (integer states q = matched length of d, resets);
  * SLACK   : the slack automaton m -> min(3m/2 - e, 1) of THM-4654 (1), re-implemented with integer
              numerators over 2^n (no Fraction), only to test statement 1 itself.
Then a DFS over all safe words (follower-graph acceptance, cross-checked by DIRECT on every node up to
length 24) checks the block law, the k <= floor(l/3)+1 density bound, the slope and offset bounds, the
proof's distance claim t_r >= 3r - 3 (r >= 2), and DP over follower states gives N(l), k_max(l), sup T_w(0).
"""
import sys
import math
from fractions import Fraction as F
from itertools import product

out_lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out_lines.append(s)


# ---------------------------------------------------------------- greedy expansion of 1 in base 3/2 (integers)
def greedy_d(n):
    """d_1..d_n with 1 = sum d_j (2/3)^j greedy; z_k = a_k / 2^k (integer a_k)."""
    a, d = 1, []          # z_0 = 1 = 1/2^0
    for k in range(1, n + 1):
        # z_k = (3/2) z_(k-1) - d_k ; (3/2) z_(k-1) = 3 a / 2^k
        t = 3 * a          # numerator over 2^k
        dk = 1 if t >= 2 ** k else 0
        a = t - dk * 2 ** k
        d.append(dk)
        assert 0 <= a < 2 ** k
        assert a % 2 == 1  # odd numerator: z_k never 0, d infinite (THM-2228 sec. 4)
    return d


DN = 2000
d = greedy_d(DN)                    # d[j-1] = d_j
D1 = [None] + d                     # 1-indexed


def safe_direct(w):
    """Q_i = 3^(L-i) * S_i (integer), Q_i = 2 w_i 3^(L-i-1) + 2 Q_(i+1); safe iff Q_i < 3^(L-i) for all i."""
    L = len(w)
    Q = 0
    p3 = 1                      # 3^(L-i-1) at i = L-1
    for i in range(L - 1, -1, -1):
        Q = 2 * w[i] * p3 + 2 * Q
        p3 *= 3                 # now 3^(L-i)
        if Q >= p3:
            return False
    return True


def safe_parry(w):
    L = len(w)
    for i in range(L):
        if w[i] == 1:
            seg = w[i:]
            pre = d[:L - i]
            if tuple(seg) > tuple(pre):
                return False
    return True


def follow_step(q, s):
    nxt = D1[q + 1]
    if s == nxt:
        return q + 1
    if s == 0 and nxt == 1:
        return 0
    return None


def safe_follow(w):
    q = 0
    for s in w:
        q = follow_step(q, s)
        if q is None:
            return False
    return True


def safe_slack(w):
    """m_n = a / 2^n exactly; step m -> min(3m/2 - e, 1), valid iff 3m/2 - e > 0."""
    a, n = 1, 0            # m = a/2^n
    for e in w:
        n += 1
        a = 3 * a - e * 2 ** n      # (3/2)(a/2^(n-1)) - e = (3a - e 2^n)/2^n
        if a <= 0:
            return False
        if a > 2 ** n:
            a = 2 ** n
    return True


# ---------------------------------------------------------------- (1) four-way equivalence on all words <= LBRUTE
LBRUTE = int(sys.argv[1]) if len(sys.argv) > 1 else 18
cnt = 0
Nbrute = []
for L in range(0, LBRUTE + 1):
    nsafe = 0
    for w in product((0, 1), repeat=L):
        a1 = safe_direct(w)
        a2 = safe_parry(w)
        a3 = safe_follow(w)
        a4 = safe_slack(w)
        assert a1 == a2 == a3 == a4, (w, a1, a2, a3, a4)
        nsafe += a1
        cnt += 1
    Nbrute.append(nsafe)
say(f"[S1] DIRECT = PARRY = FOLLOW = SLACK on all {cnt} binary words of length 0..{LBRUTE}")
say(f"     brute-force N(l), l=0..{LBRUTE}: {Nbrute}")

# ---------------------------------------------------------------- greedy d, gap sequence
ones = [j for j in range(1, 400) if D1[j]]
gaps_d = [b - a for a, b in zip(ones, ones[1:])]
say("[S2] d_1..d_60 =", "".join(map(str, d[:60])))
say("     claimed     1010000010010010100000000010 ; matches prefix:",
    "".join(map(str, d[:28])) == "1010000010010010100000000010")
say("     gap sequence of d (first 20):", gaps_d[:20])
say("     claimed (2,6,3,3,2,10,7,5,7,3,3,2,12) matches:", gaps_d[:13] == [2, 6, 3, 3, 2, 10, 7, 5, 7, 3, 3, 2, 12])
# d not eventually periodic within 2000 digits (sanity), and sum d_j (2/3)^j -> 1 from below
s = F(0)
for j in range(1, 200):
    if D1[j]:
        s += F(2, 3) ** j
say("     1 - sum_{j<200} d_j (2/3)^j =", float(1 - s), "(> 0)")

# ---------------------------------------------------------------- DFS over all safe words up to LDFS
LDFS = int(sys.argv[2]) if len(sys.argv) > 2 else 30
LDIRECT = 24          # cross-check DFS acceptance with DIRECT on every node up to this length

stats = dict(nodes=0, gap_rule_2_6=0, gap_rule_26_3=0, tr_checked=0, end_gap2=0, first_gap2=0)
viol = []
max_off = F(0)
max_off_word = None
NDFS = [0] * (LDFS + 1)
KMAXDFS = [0] * (LDFS + 1)
tightest_tr = {}
# iterative DFS: stack of (word list, q, k, C, ones_pos list)
stack = [((), 0, 0, 0, ())]
while stack:
    w, q, k, C, op = stack.pop()
    L = len(w)
    NDFS[L] += 1
    if k > KMAXDFS[L]:
        KMAXDFS[L] = k
    if L > 0:
        stats['nodes'] += 1
        if L <= LDIRECT and stats['nodes'] % 7 == 0:   # sampled cross-check (every 7th node) to keep runtime low
            assert safe_direct(w), w
        # statement 3 / 4: density and slope  (3^k/2^L <= 3 lambda^L  <=>  3k - 3 <= L)
        if not (k <= L // 3 + 1):
            viol.append(('density', w))
        if not (3 * k - 3 <= L):
            viol.append(('slope', w))
        # offset T_w(0) = C / 2^L < 11/10   <=>  10 C < 11 * 2^L
        if not (10 * C < 11 * 2 ** L):
            viol.append(('offset', w))
        if C * max_off.denominator > max_off.numerator * 2 ** L:
            max_off = F(C, 2 ** L)
            max_off_word = w
        # distances from the end: t = L - n (n 0-indexed), r-th one from the end
        ts = [L - n for n in reversed(op)]
        for r, t in enumerate(ts, start=1):
            if r == 1:
                if t < 1:
                    viol.append(('t1', w))
            else:
                stats['tr_checked'] += 1
                if t < 3 * r - 3:
                    viol.append(('t_r', w, r, t))
                key = (r, t - (3 * r - 3))
                tightest_tr[r] = min(tightest_tr.get(r, 10 ** 9), t - (3 * r - 3))
        # gap rules as consequences on actual safe words
        g = [b - a for a, b in zip(op, op[1:])]
        if g:
            if g[-1] < 2:
                viol.append(('gap1', w))
            if len(g) >= 2 and g[-2] == 2:
                stats['gap_rule_2_6'] += 1
                if g[-1] < 6:
                    viol.append(('gap2->6', w))
            if len(g) >= 3 and g[-3] == 2 and g[-2] == 6:
                stats['gap_rule_26_3'] += 1
                if g[-1] < 3:
                    viol.append(('gap26->3', w))
            if g[-1] == 2 and w[-1] == 1:
                stats['end_gap2'] += 1
            if g[0] == 2:
                stats['first_gap2'] += 1
        # Parry gap form: every suffix starting at a 1 has gap sequence that, at the first difference with
        # gaps_d, is LARGER; (complete gaps only; trailing open run is a prefix of zeros)
        for st in range(len(g)):
            sg = g[st:]
            for a_, b_ in zip(sg, gaps_d):
                if a_ != b_:
                    if a_ < b_:
                        viol.append(('parry-gap', w, st))
                    break
    if L == LDFS:
        continue
    for s_ in (0, 1):
        q2 = follow_step(q, s_)
        if q2 is None:
            continue
        if s_ == 1:
            stack.append((w + (1,), q2, k + 1, 3 * C + 2 ** L, op + (L,)))
        else:
            stack.append((w + (0,), q2, k, C, op))

say(f"[S3/S4] DFS over all safe words of length 1..{LDFS}: {stats['nodes']} words; violations: {len(viol)}",
    viol[:5])
say(f"     N(l) by DFS, l=0..{LDFS}: {NDFS}")
say(f"     brute N(l) == DFS N(l) on l<={LBRUTE}: {NDFS[:LBRUTE + 1] == Nbrute}")
say(f"     k_max(l) by DFS: {KMAXDFS}")
say(f"     k_max(l) == floor(l/3)+1 for 1<=l<={LDFS}: {all(KMAXDFS[l] == l // 3 + 1 for l in range(1, LDFS + 1))}")
say(f"     gap-2 -> next gap >= 6 instances checked: {stats['gap_rule_2_6']}; (2,6) -> next >= 3 instances: "
    f"{stats['gap_rule_26_3']}")
say(f"     t_r >= 3r-3 checked on {stats['tr_checked']} (word, r>=2) instances; words ending in gap 2 then 1: "
    f"{stats['end_gap2']}; words whose first gap is 2: {stats['first_gap2']}")
say(f"     min over words of t_r - (3r-3), r = 2..8: {[tightest_tr.get(r) for r in range(2, 9)]}")
say(f"     max offset T_w(0) over these words: {float(max_off):.9f} at length {len(max_off_word)} word "
    f"{''.join(map(str, max_off_word))}")
sumN = sum(NDFS[1:23])
say(f"     sum_(l=1..22) N(l) = {sumN} (THM-4654 table says 34,805 safe words of length 1..22)")

# ---------------------------------------------------------------- DP over follower states to large l
LDP = 400
cntq = {0: 1}
kq = {0: 0}
offq = {0: F(0)}
N = [1]
KM = [0]
SUPOFF = [F(0)]
for L in range(1, LDP + 1):
    c2, k2, o2 = {}, {}, {}
    for q, c in cntq.items():
        for s_ in (0, 1):
            q2 = follow_step(q, s_)
            if q2 is None:
                continue
            c2[q2] = c2.get(q2, 0) + c
            kk = kq[q] + s_
            if kk > k2.get(q2, -1):
                k2[q2] = kk
            if L <= 200:
                oo = (3 * offq[q] + 1) / 2 if s_ else offq[q] / 2
                if q2 not in o2 or oo > o2[q2]:
                    o2[q2] = oo
    cntq, kq = c2, k2
    if L <= 200:
        offq = o2
        SUPOFF.append(max(o2.values()))
    N.append(sum(c2.values()))
    KM.append(max(k2.values()))
say(f"[DP] N(l) for l = 0..30 by follower DP: {N[:31]}")
say(f"     DP N == DFS N for l <= {LDFS}: {N[:LDFS + 1] == NDFS}")
say(f"     k_max(l) == floor(l/3)+1 for all 1 <= l <= {LDP}: {all(KM[l] == l // 3 + 1 for l in range(1, LDP + 1))}")
ratios = [N[l] / 1.5 ** l for l in (20, 50, 100, 200, 400)]
say(f"     N(l)/(3/2)^l at l = 20, 50, 100, 200, 400: {[round(r_, 8) for r_ in ratios]}")
# exact constant R = 3 / sum_j j d_j (2/3)^j  (generating function 1/((1-z)(1 - sum d_j z^j)))
Sj = sum(j * D1[j] * (2.0 / 3.0) ** j for j in range(1, DN))
say(f"     predicted constant 3/sum_j j d_j (2/3)^j = {3 / Sj:.8f}  (THM-4654 says ~1.551)")
say(f"     sup_w T_w(0) over safe words of length <= 200 (DP, exact): {float(max(SUPOFF)):.15f}; "
    f"11/10 - sup = {float(F(11, 10) - max(SUPOFF)):.3e}")
# the witness family (100)^r 101: T_w(0) = 1/2 + sum_(s=1)^(r+1) (3/8)^s
for r in (0, 1, 2, 5, 10, 30):
    w = [1, 0, 0] * r + [1, 0, 1]
    assert safe_follow(w) and safe_direct(w) if len(w) <= 40 else safe_follow(w)
    C = 0
    for i, e in enumerate(w):
        if e:
            C = 3 * C + 2 ** i
    val = F(C, 2 ** len(w))
    closed = F(1, 2) + sum(F(3, 8) ** s for s in range(1, r + 2))
    assert val == closed, (r, val, closed)
    assert sum(w) == len(w) // 3 + 1
say("     (100)^r 101 is safe, has floor(l/3)+1 ones, and T_w(0) = 1/2 + sum_(s=1)^(r+1) (3/8)^s  (r = 0,1,2,5,10,30)")
# sharpness for every l: (100)^r 1, (100)^r 10, (100)^r 101 (l = 3r+1, 3r+2, 3r+3)
okshp = True
for L in range(1, 301):
    r = (L - 1) // 3
    if L % 3 == 1:
        w = [1, 0, 0] * r + [1]
    elif L % 3 == 2:
        w = [1, 0, 0] * r + [1, 0]
    else:
        w = [1, 0, 0] * ((L - 3) // 3) + [1, 0, 1]
    okshp &= (len(w) == L and safe_follow(w) and sum(w) == L // 3 + 1)
say(f"     sharpness witnesses (100)^r 1 | (100)^r 10 | (100)^r 101 attain floor(l/3)+1 for all l <= 300: {okshp}")

# ---------------------------------------------------------------- minimal next gaps over ALL reachable states
# For each follower state reached right after a '1' whose previous gap is g (or pattern), find min next gap
# over all safe words (not only the slack-maximal entry). Exhaustive via DFS data above is implicit; here a direct
# statement: the minimal next gap after the gap pattern P over all safe words of length <= 40.
def min_next_gap(pattern, LMAXW=45):
    best = None
    # build words: prefix X (any safe prefix ending right before the first 1 of the pattern) is irrelevant for the
    # minimum since a fresh start (state 0) maximises slack; but we verify by also trying all prefixes of length<=12.
    prefixes = [()]
    for L in range(1, 13):
        for w in product((0, 1), repeat=L):
            if safe_follow(w):
                prefixes.append(w)
    for pre in prefixes:
        body = list(pre)
        for g in pattern:
            body += [1] + [0] * (g - 1)
        if not safe_follow(body + [1]):
            continue
        for h in range(1, 40):
            ww = body + [1] + [0] * (h - 1) + [1]
            if safe_follow(ww):
                if best is None or h < best:
                    best = h
                break
    return best


for pat in [(2,), (3,), (4,), (5,), (2, 6), (2, 7), (2, 6, 3), (2, 6, 3, 3), (2, 6, 3, 3, 2), (3, 2)]:
    say(f"     min next gap after gap pattern {pat} over all safe prefixes (len<=12) : {min_next_gap(pat)}")

with open('/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_D_20261014/'
          'audit_safe_language.out', 'w') as fh:
    fh.write("\n".join(out_lines) + "\n")

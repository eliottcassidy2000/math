#!/usr/bin/env python3
"""collatz_procgen_20260922_mirror_loops_family.py -- the splicing closure for loops through -1 (Q1 mirror
of loops_family.py and of Theorem 4.2 of the loops note).

Records.  eta(n) = ceil(n log_3 2) - n log_3 2 in (0,1); 3^{eta(n)} = 3^{a(n)}/2^n with a(n) = ceil(n log_3 2).
A record is an n with eta(n) < eta(n') for all n' < n (computed exactly: 3^{a(n)} 2^{n'} < 3^{a(n')} 2^n).
These are the LOWER best approximations n/a(n) < log_2 3 (3^a slightly above 2^n).
Theorem (mirror of 4.2): if every record n >= 11 has a base loop B_n (K = n-1 halvings, a0(n-1)
multiplications) and a cycle C_n with (p,q) = (n, a(n)) through a hub visited by every B_{n*}, n* >= n,
then every K with n*(K+1) >= 11 has an a0-loop.  Greedy: n* = largest record <= K+1, m = K+1-n*,
m = sum of records (largest record < m first), splice those cycles into B_{n*}.
Input: a file with 'RECON h=-H P=p q=q word: ...' lines (base loops have h=-1).
Output: one line per K (K, a, word length, verified) plus a summary; with --words prints the words.
Every word is checked independently: E-simulation from -1 (halvings only at even values), return to -1,
K halvings, a = a0(K) exactly (3^a > 2^(K+1) >= 3^(a-1)), identity 2^K + B = 3^a, cost formula.
Usage: python3 ..._mirror_loops_family.py RECONFILE KMAX [--words]
"""
import sys, re
from fractions import Fraction as F

import math
def aceil(n):                       # ceil(n log_3 2): least a with 3^a > 2^n (n >= 1), exact
    a = max(0, int(n * math.log(2) / math.log(3)) - 2)
    while 3 ** a <= 2 ** n: a += 1
    return a

def a0(K): return aceil(K + 1)       # least a with 3^a > 2^(K+1)

def records_upto(nmax):
    recs = []; best = None          # best = (a, n) with minimal 3^a/2^n
    for n in range(1, nmax + 1):
        a = aceil(n)
        if best is None or 3 ** a * 2 ** best[1] < 3 ** best[0] * 2 ** n:
            recs.append(n); best = (a, n)
    return recs

def parse_word(s):
    w = []
    for tok in s.split():
        m = re.fullmatch(r'([MH])(\d+)', tok)
        if m: w += [m.group(1)] * int(m.group(2))
    return w

def traj(word, start=-1):
    v = start; vals = [v]
    for c in word:
        if c == 'M': v = 3 * v + 1
        else:
            if v % 2: return None
            v //= 2
        vals.append(v)
    return vals

def verify(word, K):
    """E-simulation from -1 (halvings only at even values), return to -1, K halvings, a = a0(K) exactly,
    identity 2^K + B = 3^a (B by Horner: B <- 3B + 2^(K_i)), and the cost formula
    3^a/2^K = prod 3v/(3v+1), checked as prod(v) * 2^K == prod(3v+1)."""
    vals = traj(word)
    if vals is None or vals[-1] != -1: return False
    a = word.count('M')
    if word.count('H') != K or a != a0(K): return False
    B = 0; b = 0
    for c in word:
        if c == 'H': b += 1
        else: B = 3 * B + (1 << b)
    if (1 << K) + B != 3 ** a: return False
    pv = 1; pw = 1; v = -1
    for c in word:
        if c == 'M': pv *= v; pw *= 3 * v + 1; v = 3 * v + 1
        else: v //= 2
    return pv * (1 << K) == pw

def splice(word, hub, cyc):
    vals = traj(word)
    for i, v in enumerate(vals):
        if v == hub: return word[:i] + cyc + word[i:]
    return None

def compress(word):
    out = []; last = None; run = 0
    for c in word:
        if c == last: run += 1
        else:
            if last: out.append(f"{last}{run}")
            last = c; run = 1
    if last: out.append(f"{last}{run}")
    return " ".join(out)

def main():
    fn = sys.argv[1]; KMAX = int(sys.argv[2]); show = '--words' in sys.argv
    base = {}; cyc = {}
    for line in open(fn):
        m = re.match(r'RECON h=-(\d+) P=(\d+) q=(\d+) word:(.*)$', line)
        if not m: continue
        h, p, q, w = int(m.group(1)), int(m.group(2)), int(m.group(3)), parse_word(m.group(4))
        vals = traj(w, -h)
        assert vals is not None and vals[-1] == -h and w.count('H') == p and w.count('M') == q
        if h == 1 and p >= 2: base[p + 1] = w          # base loop B_n, n = K+1
        elif q == aceil(p): cyc.setdefault(p, (-h, w))   # cycle C_p with (p, ceil(p log_3 2)) through hub -h
    recs = records_upto(KMAX + 1)
    print(f"records n <= {KMAX+1} (lower best approximations n/a < log_2 3): {recs}")
    print(f"base loops available for n = {sorted(base)}; cycles for p = {sorted(cyc)} (hubs {[cyc[p][0] for p in sorted(cyc)]})")
    base[3] = ['M', 'H', 'M', 'H']
    ok_K = []; fail_K = []
    for K in range(2, KMAX + 1):
        n = K + 1
        nstar = max(r for r in recs if r <= n)
        if nstar not in base: fail_K.append((K, 'no base')); continue
        word = list(base[nstar]); m = n - nstar; parts = []
        while m > 0:
            if m in recs: parts.append(m); break
            r = max(r for r in recs if r < m); parts.append(r); m -= r
        good = True
        for r in parts:
            if r not in cyc: good = False; break
            hub, cw = cyc[r]
            nw = splice(word, hub, cw)
            if nw is None: good = False; break
            word = nw
        if good and verify(word, K):
            ok_K.append(K)
            if show: print(f"K={K} a={word.count('M')} parts={parts} base=B_{nstar} word: {compress(word)}")
        else:
            fail_K.append((K, f"parts {parts} base B_{nstar}"))
    def ranges(L):
        out = []; s = None; prev = None
        for x in L:
            if s is None: s = prev = x
            elif x == prev + 1: prev = x
            else: out.append((s, prev)); s = prev = x
        if s is not None: out.append((s, prev))
        return ", ".join(f"{a}" if a == b else f"{a}..{b}" for a, b in out)
    print(f"verified a0-loops (independent exact check) for K in: {ranges(ok_K)}  ({len(ok_K)} values)")
    print(f"not produced by the splicing family: {ranges([k for k, _ in fail_K])}")
    for K, why in fail_K[:12]: print(f"   K={K}: {why}")

if __name__ == '__main__':
    main()

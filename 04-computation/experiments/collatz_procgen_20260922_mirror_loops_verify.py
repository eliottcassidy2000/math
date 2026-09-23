#!/usr/bin/env python3
"""collatz_procgen_20260922_mirror_loops_verify.py -- exact checks for loops through -1 (Q1 mirror).

A loop through -1 is a word in {M: v->3v+1, H: v->v/2 (v even)} from -1 back to -1.
Checks (all exact integer arithmetic, independent of the C DP):
 (1) Loop equation and automatic legality (mirror of Prop 1.1): for every K <= KEX and every a <= AEX,
     enumerate ALL nondecreasing (K_0..K_{a-1}) in [0,K]; the identity 2^K + B = 3^a holds iff the
     word is a legal E-loop through -1 (simulated); also K_0 = 0, the initial M-run is odd and the
     final H-run is odd.
 (2) Ratio bound (mirror of Prop 1.2): every solution has 2^(K+1) <= 3^a + 1, equality only (a,K)=(1,1);
     so a >= 2 gives 3^a > 2^(K+1), i.e. a >= a0(K) = ceil((K+1) log_3 2).
 (3) Exact amin(K) for K <= KEX (no height cap): the exceptions K = 1, 5..9.
 (4) Words from a file ('RECON K=.. word: M3 H1 ...' lines, or 'K=.. word: ...'): legality, endpoint,
     a = a0(K) (exact), identity B = 3^a - 2^K, cost formula 3^a/2^K = prod 3|v|/(3|v|-1), height.
Usage: python3 ..._mirror_loops_verify.py [WORDFILE ...]
"""
import sys, re, itertools
from fractions import Fraction as F

def a0(K):                       # least a with 3^a > 2^(K+1)
    a = 0
    while 3 ** a <= 2 ** (K + 1): a += 1
    return a

def word_from_Ks(Ks, K):
    """Ks nondecreasing halving counts before each multiplication; total K halvings."""
    w = []; b = 0
    for Ki in Ks:
        w += ['H'] * (Ki - b); b = Ki; w.append('M')
    w += ['H'] * (K - b)
    return w

def simulate(word, start=-1):
    v = start; pts = []; ok = True
    for c in word:
        if c == 'M': pts.append(v); v = 3 * v + 1
        else:
            if v % 2: ok = False; break
            v //= 2
    return ok, v, pts

def check_equation(KEX=12, AEX=12):
    print(f"(1)-(3) exhaustive loop equation, K <= {KEX}, a <= {AEX} (all nondecreasing K_i in [0,K]):")
    sols = {}; bad = 0; nwords = 0
    for K in range(1, KEX + 1):
        for a in range(1, AEX + 1):
            for Ks in itertools.combinations_with_replacement(range(K + 1), a):
                nwords += 1
                B = sum(3 ** (a - 1 - i) * 2 ** Ks[i] for i in range(a))
                eq = (2 ** K + B == 3 ** a)
                w = word_from_Ks(Ks, K)
                ok, v, pts = simulate(w)
                legal_loop = ok and v == -1
                if eq != legal_loop: bad += 1
                if eq:
                    sols.setdefault(K, []).append((a, Ks))
                    # structure: K_0=0, initial M-run odd, final H-run odd, ratio bound
                    run_m = sum(1 for x in Ks if x == 0)
                    assert Ks[0] == 0 and run_m % 2 == 1 and (K - Ks[-1]) % 2 == 1
                    assert 2 ** (K + 1) <= 3 ** a + 1 and (2 ** (K + 1) != 3 ** a + 1 or (a, K) == (1, 1))
    print(f"   words tested {nwords}; identity <=> legal loop through -1 on all of them: {bad == 0}")
    print("   K: a0(K), exact amin(K) (a <= %d), all (a, #solutions):" % AEX)
    for K in range(1, KEX + 1):
        s = sols.get(K, [])
        amin = min(a for a, _ in s) if s else None
        cnt = {}
        for a, _ in s: cnt[a] = cnt.get(a, 0) + 1
        flag = '=a0' if amin == a0(K) else ('<a0 (a=1 loop)' if amin is not None and amin < a0(K) else '>a0')
        print(f"   K={K:2d} a0={a0(K):2d} amin={amin} {flag}  solutions by a: {dict(sorted(cnt.items()))}")
    return sols

def parse_word(s):
    w = []
    for tok in s.split():
        m = re.fullmatch(r'([MH])(\d+)', tok)
        if not m: continue
        w += [m.group(1)] * int(m.group(2))
    return w

def verify_word(word, K_expect=None, start=-1, cycle=False, quiet=False):
    ok, v, pts = simulate(word, start)
    a = word.count('M'); K = word.count('H')
    res = {'legal': ok, 'closed': v == start, 'a': a, 'K': K}
    if not (ok and v == start): return False, res
    # identity: B = 3^a*start... for start=-1: 2^K + B = 3^a, B = sum 3^(a-1-i) 2^(K_i)
    Ks = []; b = 0
    for c in word:
        if c == 'H': b += 1
        else: Ks.append(b)
    B = sum(3 ** (a - 1 - i) * 2 ** Ks[i] for i in range(a))
    ident = (3 ** a * start + B == start * 2 ** K)
    cost = F(1)
    for p in pts: cost *= F(3 * p, 3 * p + 1)
    costok = (cost == F(3 ** a, 2 ** K))
    H = max((-p for p in pts), default=0)
    res.update(ident=ident, cost=costok, H=H)
    if cycle:
        good = ident and costok
    else:
        good = ident and costok and (K_expect is None or K == K_expect) and a == a0(K)
    res['a0'] = (a == a0(K))
    return good, res

def check_files(files):
    for fn in files:
        n = 0; good = 0; bad = []
        for line in open(fn):
            m = re.search(r'K=(\d+).*word:\s*(.*)$', line)
            if not m: continue
            K = int(m.group(1)); w = parse_word(m.group(2)); n += 1
            g, res = verify_word(w, K)
            if g: good += 1
            else: bad.append((K, res))
        print(f"(4) {fn}: {good}/{n} words verified (legal E-loop through -1, identity, cost formula, a = a0(K))")
        for K, res in bad[:10]: print("    FAILED", K, res)

if __name__ == '__main__':
    check_equation(12, 12)
    if len(sys.argv) > 1: check_files(sys.argv[1:])

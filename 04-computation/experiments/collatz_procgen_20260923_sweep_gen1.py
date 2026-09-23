#!/usr/bin/env python3
"""collatz_procgen_20260923_sweep_gen1.py  (HYP-9126: generation-1 points above 3/2)

Backward (3-adic, Q2 side): T_i^{-1}(2) = 1/2 + 3^i/2^a, a = K0(i-1) = floor(i log2 3), value 1/2 + 3/c(i-1)
  in (3/2, 5/2), where c(i-1) = 2^a/3^(i-1).
  * if c(i-1) < 9/4: explicit certificate  (k_1+1, k_2, ..., k_{i-1}, 0, 1)  built from a K0-loop
    (k_1..k_{i-1}) of length i-1 (spliced from the record objects), multiplier 2^(a+2)/3^(i+1) = 4c/9 < 1;
  * if c(i-1) > 9/4: the Psi-landing obstruction (proved in the note) forbids every descent through the
    landing value 2; we search "first move k1, then the greedy map G" certificates.
Forward (2-adic, Q1 side): X_i = -1 - 2^i/3^(A-1), A = a0(i) = ceil((i+1) log3 2), |X_i| = 1 + 3*2^i/3^A.
  * if 3^A < 2^(i+2): explicit certificate: the forward word of an a0-loop through -1 with i halvings,
    then H H (lands -4 -> -2 -> -1), multiplier 3^A/2^(i+2) < 1;
  * else: obstruction through the landing -2 (Lemma S); route search M^a + negative Collatz.
Every certificate is replayed with exact rational arithmetic.

usage: python3 collatz_procgen_20260923_sweep_gen1.py IMAX [K1MAX]
"""
import sys, os, math
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collatz_procgen_20260923_sweep_splice as sp
from collatz_procgen_20260923_sweep_records import eps_rec, eta_rec, ceil_ntheta

def is_unit3(x):
    return x.numerator % 3 != 0 and x.denominator % 3 != 0

def replay_backward(y0, ks):
    """reverse moves x -> (2^k x - 1)/3 on 3-adic units (dyadic rationals); returns min prefix multiplier info"""
    x = y0; K = 0; best = None
    for s, k in enumerate(ks, 1):
        t = x * 2**k - 1
        if t.numerator % 3: return False, None
        x = t / 3; K += k
        if not is_unit3(x): return False, None
        if 2**K < 3**s and best is None: best = (s, K)
    return best is not None, best

def greedy_k(m):
    for k in range(0, 7):
        if (2**k * m) % 9 in (4, 7): return k
    raise ValueError

L2_3 = math.log2(3)
def _desc(K, s):
    """exact test 2^K < 3^s, with a float pre-filter"""
    d = s * L2_3 - K
    if d > 1e-6: return True
    if d < -1e-6: return False
    return 2**K < 3**s

def slow_descent(y0, K1MAX, SMAX=30000, K1MIN=1):
    for k1 in range(K1MIN, K1MAX + 1):
        t = y0 * 2**k1 - 1
        if t.numerator % 3: continue
        x = t / 3
        if not is_unit3(x) or x.denominator != 1: continue
        m = x.numerator; s = 1; K = k1; ks = [k1]
        if _desc(K, s): return ks
        while m != 1 and s < SMAX:
            k = greedy_k(m); m = ((m << k) - 1) // 3; s += 1; K += k; ks.append(k)
            if _desc(K, s): return ks
    return None

def forward_replay(x, word):
    """word: list of 'M'/'H' applied to the 2-adic rational x; returns True iff legal and some halving
    has multiplier 3^a/2^b < 1"""
    a = b = 0; v = x; desc = False
    for c in word:
        if c == 'M': v = 3*v + 1; a += 1
        else:
            if v.numerator % 2: return False, None
            v = v / 2; b += 1
            if 3**a < 2**b: desc = True
    return desc, (a, b, v)

def fwd_route(x, AEXTRA=600, STEPS=40000):
    j = 0; d = x.denominator
    while d % 3 == 0: d //= 3; j += 1
    for a in range(j, j + AEXTRA):
        v = x * 3**a + F(3**a - 1, 2)
        v = v.numerator; at = a; b = 0
        for _ in range(STEPS):
            if v in (-1, -5, -7, -17): break
            if v % 2 == 0:
                v //= 2; b += 1
                if 3**at < 2**b: return (a, at, b)
            else:
                v = 3*v + 1; at += 1
    return None

def run_collatz_extra(v, at, b, STEPS, allow_extra):
    """negative Collatz map from the integer v (halve even, 3v+1 odd), optionally with ONE extra x3-move at an
    even value of absolute value < 10^7 (an E-only arrow); returns (a_tot, b) at the first halving with
    3^a_tot < 2^b, or None"""
    for _ in range(STEPS):
        if v in (-1, -5, -7, -17): return None          # entered a negative Collatz cycle
        if allow_extra and -v <= 200: return None       # main trajectory reached the small region
        if v % 2 == 0:
            if allow_extra and -v < 10**7:
                r = run_collatz_extra(3*v + 1, at + 1, b, STEPS, False)
                if r: return r
            v //= 2; b += 1
            if 3**at < 2**b: return (at, b)
        else:
            v = 3*v + 1; at += 1
    return None

def deep_fwd(i, AMAX, STEPS):
    A = ceil_ntheta(i + 1); x = -1 - F(2**i, 3**(A - 1)); j = A - 1
    for a in range(j, j + AMAX):
        v = (x * 3**a + F(3**a - 1, 2)).numerator
        r = run_collatz_extra(v, a, 0, STEPS, True)
        if r: return (a,) + r
    return None

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] in ('deep-back', 'deep-fwd'):
    for i in [int(t) for t in sys.argv[2].split(',')]:
        if sys.argv[1] == 'deep-back':
            a = (3**i).bit_length() - 1; x = F(1, 2) + F(3**i, 2**a)
            ks = slow_descent(x, int(sys.argv[3]), K1MIN=int(sys.argv[4]) if len(sys.argv) > 4 else 1)
            if ks is None: print("deep-back i=%d: no certificate k1<=%s" % (i, sys.argv[3]))
            else:
                ok, best = replay_backward(x, ks)
                print("deep-back i=%d: DESCENDS k1=%d s=%d K=%d replay=%s" % (i, ks[0], best[0], best[1], ok))
        else:
            r = deep_fwd(i, int(sys.argv[3]), int(sys.argv[4]))
            print("deep-fwd i=%d: %s" % (i, ("DESCENDS a=%d 3^%d<2^%d" % r) if r else "none (a<=j+%s, one extra x3-move)" % sys.argv[3]))
        sys.stdout.flush()
    sys.exit(0)

if __name__ == '__main__':
    IMAX = int(sys.argv[1]); K1MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    Bp, Cp = sp.load('pos'); Bn, Cn = sp.load('neg')
    recp = eps_rec(max(max(Bp), max(Cp)) + 1); recn = eta_rec(max(max(Bn), max(Cn)) + 1)
    L = math.log2(3)
    nb_easy = nb_hard = nb_hard_ok = 0; hard_fail = []
    print("# backward T_i^{-1}(2) = 1/2 + 3^i/2^floor(i log2 3)")
    for i in range(3, IMAX + 1):
        a = (3**i).bit_length() - 1                      # floor(i log2 3)
        c = F(2**a, 3**(i-1)); x = F(1, 2) + F(3**i, 2**a)
        if c < F(9, 4):
            # K0-loop of length i-1 from the record objects
            if i - 1 == 2: loop = [2, 2]
            else:
                ok, info = sp.build('pos', i - 1, Bp, Cp, recp)
                if not ok: print("  i=%d: no loop" % i); continue
                qs, parts = sp.decompose('pos', i, recp); loop = list(Bp[qs])
                for q in parts: h, cyc = Cp[q]; loop = sp.splice(1, loop, h, cyc)
            ks = [loop[0] + 1] + loop[1:] + [0, 1]
            ok, best = replay_backward(x, ks)
            nb_easy += 1
            if not ok: print("  i=%d easy route FAILED" % i)
        else:
            nb_hard += 1
            ks = slow_descent(x, K1MAX)
            if ks is None: hard_fail.append(i); st = "no certificate k1<=%d" % K1MAX
            else:
                ok, best = replay_backward(x, ks); nb_hard_ok += ok
                st = "DESCENDS k1=%d s=%d K=%d replay=%s" % (ks[0], best[0], best[1], ok)
            print("  i=%d c=%.4f value=%.5f %s" % (i, float(c), float(x), st), flush=True)
    print("backward i<=%d: easy clocks (c<9/4) certified %d; hard clocks %d, certified %d, open %s"
          % (IMAX, nb_easy, nb_hard, nb_hard_ok, hard_fail))
    print("# forward X_i = -1 - 2^i/3^(a0(i)-1)")
    nf_easy = nf_hard = nf_hard_ok = 0; ffail = []
    for i in range(10, IMAX + 1):
        A = ceil_ntheta(i + 1); x = -1 - F(2**i, 3**(A - 1))
        if 3**A < 2**(i + 2):
            ok, info = sp.build('neg', i, Bn, Cn, recn)
            if not ok: print("  i=%d: no mirror loop" % i); continue
            ns, parts = sp.decompose('neg', i + 1, recn); loop = list(Bn[ns])
            for q in parts: h, cyc = Cn[q]; loop = sp.splice(-1, loop, h, cyc)
            word = []
            for idx in range(len(loop), 0, -1): word += ['M'] + ['H'] * loop[idx-1]
            word += ['H', 'H']
            desc, st = forward_replay(x, word)
            nf_easy += 1
            if not desc or st[2] != -1: print("  i=%d forward easy route FAILED %s" % (i, st))
        else:
            nf_hard += 1
            r = fwd_route(x)
            if r is None: ffail.append(i); st = "no route certificate"
            else: nf_hard_ok += 1; st = "DESCENDS a=%d 3^%d<2^%d" % r
            print("  i=%d |x|=%.5f %s" % (i, float(-x), st), flush=True)
    print("forward 10<=i<=%d: easy clocks certified %d; hard clocks %d, certified %d, open %s"
          % (IMAX, nf_easy, nf_hard, nf_hard_ok, ffail))

#!/usr/bin/env python3
"""collatz_procgen_20260922_loops_verify.py -- independent exact check of E-loops through 1.

Input: lines containing 's=<s> ... k:<k_1,...,k_s>' (reverse moves from x_0 = 1:
x_i = (2^{k_i} x_{i-1} - 1)/3), e.g. the output of loops_dp ... recon or of loops_family.py.
For every loop it checks, with exact integers and two independent code paths:
  (1) reverse legality: every x_i is a positive integer not divisible by 3, x_s = 1;
  (2) forward E-simulation: from 1 apply, for i = s..1, the arrow n -> 3n+1 then k_i arrows n -> n/2,
      each halving applied to an even number, and the walk ends at 1 (an E-cycle through 1);
  (3) the representation identity 2^K = sum_{i=0}^{s} 3^i 2^{e_i} with e nonincreasing, e_s = 0,
      where e_i = K - (k_1 + ... + k_{s-i})  (e_i = number of halvings after the (s-i)-th multiplication
      counted from the end);
  (4) K = sum k_i equals K0(s) = max{K : 2^K < 3^(s+1)} (exact), so 2^K < 3^(s+1) (1-escape) and
      2^(K+1) > 3^(s+1) (the lower bound: no loop with smaller K exists for s >= 2).
Also reports: climb length (number of leading pure 3n+1 steps from 1), whether the loop passes through 13,
primitive (returns to 1 only at the end), height (max multiplication point).
Usage: python3 ..._loops_verify.py FILE [--quiet]
"""
import re, sys

def K0(s):
    K = 0
    t = 3 ** (s + 1)
    while 2 ** (K + 1) < t:
        K += 1
    return K

def check(s, ks):
    assert len(ks) == s, "length"
    # (1) reverse legality
    xs = [1]
    for k in ks:
        assert k >= 0
        t = (1 << k) * xs[-1] - 1
        assert t % 3 == 0, "non-integer"
        y = t // 3
        assert y >= 1 and y % 3 != 0, "illegal value %d" % y
        xs.append(y)
    assert xs[-1] == 1, "does not return to 1"
    K = sum(ks)
    # (2) forward E-simulation (independent of (1))
    n = 1; visited = [1]; mults = []
    for i in range(s, 0, -1):
        mults.append(n)
        n = 3 * n + 1
        for _ in range(ks[i - 1]):
            assert n % 2 == 0, "halving an odd number"
            n //= 2
        visited.append(n)
    assert n == 1, "forward walk does not end at 1"
    # (3) representation identity.  With K_i = k_1+...+k_i the reverse recursion gives
    #     3^s x_s = 2^K x_0 - B,  B = sum_{i=1}^s 3^(i-1) 2^(K-K_i);  with x_0 = x_s = 1:
    #     2^K = sum_{i=0}^{s} 3^i 2^(e_i),  e_i = K - K_(i+1) (i < s), e_s = 0, e nonincreasing.
    pref = [0]
    for k in ks:
        pref.append(pref[-1] + k)
    B = sum(3 ** (i - 1) * 2 ** (K - pref[i]) for i in range(1, s + 1))
    assert 2 ** K == 3 ** s + B, "carry identity"
    ee = [K - pref[i + 1] for i in range(s)] + [0]
    assert sum(3 ** i * 2 ** ee[i] for i in range(s + 1)) == 2 ** K, "representation"
    assert all(ee[i] >= ee[i + 1] for i in range(s)) and ee[s - 1] == 0, "exponents"
    # (4) K = K0(s)
    k0 = K0(s)
    climb = 0
    for i in range(s, 0, -1):
        if ks[i - 1] == 0 and i > 1:
            climb += 1
        else:
            break
    info = dict(K=K, K0=k0, eqK0=(K == k0), ratio=2 ** K / 3 ** s, height=max(xs),
                through13=(13 in xs), primitive=all(x != 1 for x in xs[1:-1]),
                climb=climb + 1)
    return info

def main():
    fn = sys.argv[1]; quiet = '--quiet' in sys.argv
    n = bad = 0; notK0 = []
    for line in open(fn):
        m = re.search(r's=\s*(\d+).*?k:(\S+)', line)
        if not m:
            continue
        s = int(m.group(1)); ks = [int(t) for t in m.group(2).split(',')]
        try:
            info = check(s, ks)
        except AssertionError as ex:
            bad += 1; print("FAIL s=%d: %s" % (s, ex)); continue
        n += 1
        if s >= 2 and not info['eqK0']:
            notK0.append(s)
        if not quiet:
            print("s=%4d K=%4d K0=%4d %s ratio=%.6f height=%d climb=%d prim=%d thru13=%d" % (
                s, info['K'], info['K0'], 'OK' if info['eqK0'] else '--', info['ratio'], info['height'],
                info['climb'], info['primitive'], info['through13']))
    print("verified loops: %d, failures: %d, s>=2 with K != K0(s): %s" % (n, bad, notK0 if notK0 else 'none'))

if __name__ == '__main__':
    main()

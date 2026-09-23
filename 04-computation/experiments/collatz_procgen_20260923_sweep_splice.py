#!/usr/bin/env python3
"""collatz_procgen_20260923_sweep_splice.py  (HYP-9122 / HYP-9125)

Explicit splicing closure (loops note Theorem 4.2, mirror note Theorem 2.2) from the record
objects written by collatz_procgen_20260923_sweep_records.py, with an independent exact check
of every produced loop.

  python3 collatz_procgen_20260923_sweep_splice.py SIDE s1,s2,...     (SIDE = pos | neg)
     pos: s = number of multiplications of a loop through 1 with K0(s)=floor((s+1)log2 3) halvings
     neg: s = K = number of halvings of a loop through -1 with a0(K)=ceil((K+1)log3 2) multiplications
  python3 ... SIDE hubcheck      -> checks the hub hypothesis of the reduction theorem for all objects

Splicing rule (proved in the note): if a reverse-move object visits h as the value 2^j x_{t-1}
(0 <= j <= k_t, forward reading x_t -> 3x_t+-1 = 2^{k_t} x_{t-1} -> ... -> x_{t-1}), then replacing
k_t by (j+c_1, c_2, ..., c_q, k_t-j), where (c_i) are the reverse moves of a cycle at h, splices the
cycle in.  (For j = 0 this inserts the cycle right after the M-point x_{t-1}.)
"""
import os, sys, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collatz_procgen_20260923_sweep_records import (verify, ceil_qL, ceil_ntheta, eps_rec, eta_rec,
                                                     OBJDIR)
sys.set_int_max_str_digits(0) if hasattr(sys, 'set_int_max_str_digits') else None

def load(side):
    B = {}; C = {}
    for f in glob.glob(os.path.join(OBJDIR, '%s_B_*.txt' % side)):
        q = int(re.search(r'_B_(\d+)\.txt', f).group(1)); B[q] = [int(t) for t in open(f).read().split(',')]
    for f in glob.glob(os.path.join(OBJDIR, '%s_C_*_*.txt' % side)):
        m = re.search(r'_C_(\d+)_(\d+)\.txt', f); q, h = int(m.group(1)), int(m.group(2))
        C[q] = (h, [int(t) for t in open(f).read().split(',')])
    # the smallest records are fixed by hand (and re-verified here):
    #   pos: C_1 = trivial loop 1 -> 4 -> 2 -> 1 (reverse move [2]); B_3 = C_1 twice
    #   neg: C_1 = basic loop -1 -> -2 -> -1 ([1]); C_3 = negative Collatz 2-cycle {-5,-7} at u=5 ([2,1])
    if side == 'pos':
        C.setdefault(1, (1, [2])); B.setdefault(3, [2, 2])
        assert verify(1, 1, [2], 2)[0] and verify(1, 1, [2, 2], 4)[0]
    else:
        C.setdefault(1, (1, [1])); C.setdefault(3, (5, [2, 1]))
        assert verify(-1, 1, [1], 1)[0] and verify(-1, 5, [2, 1], 3)[0]
    return B, C

def values(sign, h, moves):
    x = h; xs = [h]
    for k in moves:
        x = ((x << k) + (-1 if sign > 0 else 1)) // 3; xs.append(x)
    return xs

def find_visit(sign, moves, hub):
    """return (t, j) with 2^j x_{t-1} == hub, 0<=j<=k_t, or None"""
    xs = values(sign, 1, moves)
    for t in range(1, len(moves)+1):
        xp = xs[t-1]
        if xp > hub: continue
        if hub % xp: continue
        r = hub // xp
        if r & (r-1): continue
        j = r.bit_length()-1
        if j <= moves[t-1]: return (t, j)
    return None

def splice(sign, moves, hub, cyc):
    tj = find_visit(sign, moves, hub)
    if tj is None: return None
    t, j = tj
    k = moves[t-1]
    new = moves[:t-1] + [j + cyc[0]] + cyc[1:] + [k - j] + moves[t:]
    return new

def decompose(side, s_plus_1, recs):
    """greedy decomposition of Theorem 4.2 / 2.2: returns (q*, [q_2..q_t])"""
    rs = [r for r, _ in recs]
    qs = max(r for r in rs if r <= s_plus_1)
    n = s_plus_1 - qs; parts = []
    while n > 0:
        if n in rs: parts.append(n); break
        qq = max(r for r in rs if r < n); parts.append(qq); n -= qq
    return qs, parts

def build(side, s, B, C, recs):
    sign = 1 if side == 'pos' else -1
    qs, parts = decompose(side, s+1, recs)
    if qs not in B: return None, 'no base %d' % qs
    moves = list(B[qs])
    for q in parts:
        if q not in C: return None, 'no cycle %d' % q
        h, cyc = C[q]
        new = splice(sign, moves, h, cyc)
        if new is None: return None, 'hub %d of C_%d not visited by B_%d' % (h, q, qs)
        moves = new
    if side == 'pos':
        K = (3**(s+1)).bit_length()-1          # K0(s) = floor((s+1) log2 3)
        if len(moves) != s: return None, 'length'
    else:
        K = s; a = ceil_ntheta(s+1)
        if len(moves) != a: return None, 'length'
    ok, info = verify(sign, 1, moves, K)
    return ok, (qs, parts, info['height'] if ok else info)

if __name__ == '__main__':
    side = sys.argv[1]
    B, C = load(side)
    recmax = max(max(B), max(C))
    recs = eps_rec(recmax+1) if side == 'pos' else eta_rec(recmax+1)
    if sys.argv[2] == 'reverify':
        # independent exact re-verification of every stored record object
        sign = 1 if side == 'pos' else -1
        for q in sorted(B):
            K = ((3**q).bit_length() - 1) if side == 'pos' else q - 1      # pos: K0(q-1)=floor(q log2 3); neg: n-1
            ok, info = verify(sign, 1, B[q], K)
            print('%s B_%d: moves=%d K=%d verified=%s height=%s' % (side, q, len(B[q]), K, ok, info['height'] if ok else info), flush=True)
        for q in sorted(C):
            h, cyc = C[q]
            P = (3**q).bit_length() if side == 'pos' else q                 # pos: ceil(q log2 3); neg: n
            ok, info = verify(sign, h, cyc, P)
            print('%s C_%d at %d: moves=%d K=%d verified=%s height=%s' % (side, q, h, len(cyc), P, ok, info['height'] if ok else info), flush=True)
        sys.exit(0)
    if sys.argv[2] == 'tightcheck':
        # tightness lemmas T1-T3: every prefix carries the extremal K (exact integers), objects with record <= 6000
        sign = 1 if side == 'pos' else -1
        def kap(i, x, h):
            t = 3**i*x - (3**i-1)//2 if sign < 0 else 3**i*x + (3**i-1)//2
            if sign > 0:                      # least m with 2^m h >= t
                m = max(0, (t//h).bit_length()-1)
                while (h << m) < t: m += 1
                while m > 0 and (h << (m-1)) >= t: m -= 1
            else:                             # largest m with 2^m h <= t
                m = max(0, (t//h).bit_length()-1)
                while (h << (m+1)) <= t: m += 1
                while m > 0 and (h << m) > t: m -= 1
            return m
        objs = [(1, B[q], 'B_%d' % q) for q in sorted(B) if q <= 6000] + [(C[q][0], C[q][1], 'C_%d' % q) for q in sorted(C) if q <= 6000]
        bad = 0
        for h, mv, name in objs:
            x = h; K = 0
            for i, k in enumerate(mv, 1):
                x = ((x << k) + (-1 if sign > 0 else 1)) // 3; K += k
                if kap(i, x, h) != K: bad += 1; print('not tight:', name, i); break
        print('tightcheck %s: %d objects, violations %d' % (side, len(objs), bad))
        sys.exit(0)
    if sys.argv[2] == 'hubcheck':
        sign = 1 if side == 'pos' else -1
        rs = [r for r, _ in recs if r in B or r in C]
        bad = 0
        for q in sorted(B):
            for q2 in sorted(C):
                if q2 > q: continue
                h = C[q2][0]
                if find_visit(sign, B[q], h) is None:
                    bad += 1; print('B_%d does not visit hub %d of C_%d' % (q, h, q2))
        print('hubcheck %s: bases %s cycles %s failures %d' % (side, sorted(B), sorted(C), bad))
        sys.exit(0)
    for s in [int(t) for t in sys.argv[2].split(',')]:
        ok, info = build(side, s, B, C, recs)
        print('%s s=%d ok=%s base=%s parts=%s height=%s' % (side, s, ok, info[0] if ok else '-',
              (len(info[1]), sorted(set(info[1]))) if ok else info, info[2] if ok else '-'), flush=True)

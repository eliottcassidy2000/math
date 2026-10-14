#!/usr/bin/env python3
"""The repair algebra at the -1 anchor: collisions of Collatz words at -1, primitive collisions as ladder pairs,
the Gersonides carries, and an exact census.

Conventions (THM-4555): f_c(x) = (3x+1)/2^c; a word u = (a_1..a_k) acts by f_u = f_{a_k} o ... o f_{a_1}; u is reduced if
a_1 >= 2; reduced u, u' collide if f_u(-1) = f_u'(-1).

Numerator form (proved in THM-4640): for a reduced word u = (a_1..a_k) with partial sums S_j, total A = S_k,
    f_u(-1) = K(s) / 2^(A-1),   s = (e_1 < ... < e_(k-1)),  e_j = S_j - 1,
    K(s) = -3^m + sum_{j=1..m} 3^(m-j) 2^(e_j)   (m = k-1),  i.e. K(()) = -1, K(t.e) = 3 K(t) + 2^e.
So s is the increasing sequence of a word's HEAD (all letters but the last); the last letter only fixes the total.

Parts (each prints PASS/FAIL lines and counts; all checks exact):
  A  word-level validation with Fractions for every reduced word of total <= TA
  B  census of primitive collisions (K with representations of >= 2 distinct last exponents) for max exponent <= EB
  C  two-position (4-term) relations 3^D (2^x - 1) = 2^y - 1: only (D,x,y) = (1,1,2), (2,3,6) in a large box
  D  no primitive collision of admissible sequences differs in <= 2 padded positions (census up to ED)
  E  the uniform trailing-ones switches (THM-4555 (iv)) induced by the smallest primitive collisions
Usage: python3 repair_algebra.py [TA] [EB] [ED]"""
import sys, time, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np

TA = int(sys.argv[1]) if len(sys.argv) > 1 else 16
EB = int(sys.argv[2]) if len(sys.argv) > 2 else 24
ED = int(sys.argv[3]) if len(sys.argv) > 3 else 16
CHECKS = 0
def check(cond, msg):
    global CHECKS
    CHECKS += 1
    if not cond:
        print("FAIL:", msg); sys.exit(1)

def K_of(s):
    K = -1
    for e in s: K = 3 * K + (1 << e)
    return K

def head_seq(u):
    """increasing sequence of the head of reduced word u (all letters but the last)"""
    S = 0; out = []
    for a in u[:-1]:
        S += a; out.append(S - 1)
    return tuple(out)

def word_from(s, c):
    """reduced word with head sequence s and last letter c"""
    if not s: return (c,)
    letters = [s[0] + 1] + [s[j] - s[j - 1] for j in range(1, len(s))] + [c]
    return tuple(letters)

# ---------------------------------------------------------------- Part A
def part_A():
    t0 = time.time()
    vals = defaultdict(list)            # value -> words
    stack = [((a,), Fr(3 * -1 + 1, 2 ** a), a) for a in range(2, TA + 1)]
    while stack:
        u, y, A = stack.pop()
        vals[(y, A)].append(u)
        for c in range(1, TA - A + 1):
            stack.append((u + (c,), (3 * y + 1) / 2 ** c, A + c))
    nwords = sum(len(v) for v in vals.values())
    expect = sum(2 ** (A - 2) for A in range(2, TA + 1))
    check(nwords == expect, f"reduced word count {nwords} != {expect}")
    # numerator formula and collision criterion
    for (y, A), words in vals.items():
        for u in words:
            check(sum(u) == A, "total")
            check(y == Fr(K_of(head_seq(u)), 2 ** (A - 1)), f"numerator formula for {u}")
    # distinct values at equal total <-> distinct K (denominator 2^(A-1) exactly: K odd)
    for (y, A), words in vals.items():
        check(K_of(head_seq(words[0])) % 2 == 1, "K odd")
    # collisions across different totals are impossible (THM-4555 (ii))
    byval = defaultdict(set)
    for (y, A) in vals: byval[y].add(A)
    check(all(len(s) == 1 for s in byval.values()), "colliding words with different totals")
    # primitive collisions <-> >= 2 distinct last letters <-> ladder identity on heads
    nprim = 0; nroot = 0
    for (y, A), words in vals.items():
        lasts = {u[-1] for u in words}
        if len(lasts) < 2: continue
        nprim += 1
        heads = {}
        for u in words: heads.setdefault(u[-1], u[:-1])
        L = sorted(heads)
        for c1, c2 in itertools.combinations(L, 2):
            h1, h2 = heads[c1], heads[c2]
            y1 = Fr(-1) if not h1 else f_word(h1)
            y2 = Fr(-1) if not h2 else f_word(h2)
            d = c2 - c1                       # (3 y1 + 1) / 2^c1 = (3 y2 + 1) / 2^c2  ->  3 y2 + 1 = 2^d (3 y1 + 1)
            check(3 * y2 + 1 == 2 ** d * (3 * y1 + 1), "ladder identity")
            check(d % 2 == 0, "ladder exponent even")
            i = d // 2
            check(y2 == 4 ** i * y1 + Fr(4 ** i - 1, 3), "ladder in compiler form")
        if any(not u[:-1] for u in words): nroot += 1
    print(f"A: {nwords} reduced words of total <= {TA}; {len(vals)} distinct values; {nprim} values carry a primitive "
          f"collision ({nroot} of them involve a one-letter word, i.e. the root collision); all identities exact  "
          f"[{time.time() - t0:.1f}s]", flush=True)

def f_word(u):
    y = Fr(-1)
    for a in u: y = (3 * y + 1) / 2 ** a
    return y

# ---------------------------------------------------------------- Part B
def census(E):
    """all increasing sequences with max <= E: arrays K, last, length"""
    K = np.array([-1], dtype=np.int64); L = np.array([0], dtype=np.int8); M = np.array([0], dtype=np.int8)
    for e in range(1, E + 1):
        K = np.concatenate([K, 3 * K + (1 << e)])
        L = np.concatenate([L, np.full(len(L), e, dtype=np.int8)])
        M = np.concatenate([M, M + 1])
    return K, L, M

def primitive_table(E):
    K, L, M = census(E)
    o = np.lexsort((L, K)); Ks, Ls, Ms = K[o], L[o], M[o]
    same = np.concatenate([[False], Ks[1:] == Ks[:-1]])
    newlast = same & np.concatenate([[False], Ls[1:] != Ls[:-1]])
    # birth of each primitive K = last exponent at its first 'newlast' mark
    idx = np.nonzero(newlast)[0]
    firstmark = idx[np.concatenate([[True], Ks[idx[1:]] != Ks[idx[:-1]]])] if len(idx) else idx
    return Ks, Ls, Ms, firstmark

def part_B():
    t0 = time.time()
    # cross-check the vectorized census against brute force for small E
    for E in range(1, 13):
        Ks, Ls, Ms, fm = primitive_table(E)
        reps = defaultdict(set)
        for n in range(E + 1):
            for s in itertools.combinations(range(1, E + 1), n):
                reps[K_of(s)].add(s[-1] if s else 0)
        brute = sorted(k for k, ls in reps.items() if len(ls) > 1)
        check(sorted(Ks[fm].tolist()) == brute, f"census vs brute force at E = {E}")
    Ks, Ls, Ms, fm = primitive_table(EB)
    births = Ls[fm].astype(int)
    print(f"B: census of primitive collisions at -1 (max exponent <= {EB}; {len(Ks)} representations): "
          f"{len(fm)} primitive K  [{time.time() - t0:.1f}s]", flush=True)
    cum = 0; prev = None
    print("   E   new   cumulative   new/prev", flush=True)
    for E in range(1, EB + 1):
        new = int((births == E).sum()); cum += new
        ratio = f"{new / prev:.4f}" if prev else "-"
        print(f"   {E:2d} {new:7d} {cum:10d}   {ratio}", flush=True)
        prev = new if new else prev
    # smallest primitive K with their canonical representations (strip the root twin: drop a leading exponent 1)
    small = sorted(set(Ks[fm].tolist()), key=abs)[:30]
    print("   smallest |K| primitive collisions (canonical representations; head words; ladder exponent i):")
    for k in small:
        # recover representations by backward search (exact)
        rs = sorted({r[1:] if r and r[0] == 1 else r for r in backward(k, EB + 1)}, key=lambda r: (len(r), r))
        lasts = sorted({r[-1] if r else 0 for r in rs})
        heads = [word_from(r, 1)[:-1] if r else () for r in rs]
        print(f"     K = {k:6d}: " + "  ~  ".join(str(r) for r in rs) + f"   heads {heads}   last exponents {lasts}", flush=True)
    return Ks, Ls, Ms, fm

def backward(K, emax):
    """all increasing representations of K with exponents < emax (exact backward search)"""
    out = []
    if K == -1: out.append(())
    for e in range(1, emax):
        if (K - (1 << e)) % 3 == 0:
            P = (K - (1 << e)) // 3
            if P < -(3 ** e) * 2: continue
            for r in backward(P, e): out.append(r + (e,))
    return out

# ---------------------------------------------------------------- Part C
def part_C(DMAX=60, XMAX=400):
    t0 = time.time()
    sols = []
    for D in range(1, DMAX + 1):
        p3 = 3 ** D
        for x in range(1, XMAX + 1):
            v = p3 * ((1 << x) - 1) + 1
            if v & (v - 1) == 0: sols.append((D, x, v.bit_length() - 1))
    check(sols == [(1, 1, 2), (2, 3, 6)], f"two-position relations {sols}")
    print(f"C: 3^D (2^x - 1) = 2^y - 1 with 1 <= D <= {DMAX}, 1 <= x <= {XMAX}: solutions {sols} "
          f"(unit carry 3*1 = 4 - 1 and Catalan carry 9*7 = 64 - 1)  [{time.time() - t0:.1f}s]", flush=True)

# ---------------------------------------------------------------- Part D
def part_D():
    t0 = time.time()
    reps = defaultdict(list)
    for n in range(ED + 1):
        for s in itertools.combinations(range(1, ED + 1), n):
            reps[K_of(s)].append(s)
    hist = defaultdict(int); npairs = 0
    for k, rs in reps.items():
        for s, t in itertools.combinations(rs, 2):
            if (s[-1] if s else 0) == (t[-1] if t else 0): continue
            m = max(len(s), len(t))
            sp = (1,) * (m - len(s)) + s; tp = (1,) * (m - len(t)) + t
            r = sum(1 for a, b in zip(sp, tp) if a != b)
            hist[r] += 1; npairs += 1
            # r = 0 only for padding twins (the root collision () ~ (1)); otherwise at least three positions differ
            check(r >= 3 or (r == 0 and k == -1), f"primitive collision differing in {r} padded positions: {s} ~ {t}")
    print(f"D: {npairs} primitive pairs (distinct last exponents) up to exponent {ED}; padded differing positions "
          f"histogram {dict(sorted(hist.items()))}  [{time.time() - t0:.1f}s]", flush=True)

# ---------------------------------------------------------------- Part E
def part_E(n=12):
    """uniform switches: a primitive pair of heads (h, h') with words (h, c) ~ (h', c'), lengths differing by D >= 1,
    gives (THM-4555 (iv)) U^j(N) = U^j((N+1)/2^D - 1) for every odd N whose word starts 1^r (shorter word), r >= D."""
    t0 = time.time()
    shown = 0
    for k in sorted({x for x in primitive_K_small()}, key=abs):
        rs = sorted({r[1:] if r and r[0] == 1 else r for r in backward(k, 24)}, key=lambda r: (len(r), r))
        for s, t in itertools.combinations(rs, 2):
            es, et = (s[-1] if s else 0), (t[-1] if t else 0)
            if es == et: continue
            # word lengths: |head| + 1;  last letters c_s, c_t with e_s + c_s = e_t + c_t  (choose the smallest valid)
            c_s = max(1, 1 + et - es); c_t = c_s + es - et
            ws, wt = word_from(s, c_s), word_from(t, c_t)
            check(f_word(ws) == f_word(wt) and sum(ws) == sum(wt), "induced word collision")
            short, long_ = (ws, wt) if len(ws) <= len(wt) else (wt, ws)
            D = len(long_) - len(short)
            if D < 1: continue
            # verify the switch on actual integers: odd N with word starting 1^r short (r = D), a few instances
            r = D; ok = 0
            for N in odd_sources_with_prefix((1,) * r + short, 4):
                Mrt = (N + 1) // 2 ** D - 1
                j = r + len(short)
                check(Uj(N, j) == Uj(Mrt, j), f"switch {short} ~ {long_} at N = {N}")
                ok += 1
            print(f"   {short} ~ {long_}: D = {D}; verified on {ok} sources N (merge with (N+1)/2^{D} - 1 at step {r + len(short)})")
            shown += 1
            if shown >= n: break
        if shown >= n: break
    print(f"E: {shown} switch families verified on actual integers  [{time.time() - t0:.1f}s]", flush=True)

def primitive_K_small():
    reps = defaultdict(set)
    for n in range(15):
        for s in itertools.combinations(range(1, 15), n):
            reps[K_of(s)].add(s[-1] if s else 0)
    return [k for k, ls in reps.items() if len(ls) > 1 and k != -1]

def Uj(n, j):
    for _ in range(j):
        n = 3 * n + 1
        while n % 2 == 0: n //= 2
    return n

def odd_sources_with_prefix(w, count):
    """odd N >= 3 whose U-word starts with w: N = I_w(z) for endpoints z = 2^A... (THM-4555 (i)); found by inverse map"""
    out = []
    A = sum(w); k = len(w)
    # sources: N = (2^A z - B)/3^k with z odd, N odd positive and word exactly w (checked forward)
    z = 1
    while len(out) < count and z < 10 ** 7:
        B = 0
        # inverse of the composition along w applied to z
        N = Fr(z)
        for a in reversed(w): N = (2 ** a * N - 1) / 3
        if N.denominator == 1 and N > 1 and N.numerator % 2 == 1:
            Ni = N.numerator
            if word_prefix(Ni, k) == tuple(w): out.append(Ni)
        z += 2
    return out

def word_prefix(n, k):
    out = []
    for _ in range(k):
        n = 3 * n + 1; a = 0
        while n % 2 == 0: n //= 2; a += 1
        out.append(a)
    return tuple(out)

if __name__ == '__main__':
    part_A()
    part_B()
    part_C()
    part_D()
    part_E()
    print(f"TOTAL exact checks: {CHECKS}")

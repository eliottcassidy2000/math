#!/usr/bin/env python3
"""hankel_bf.py -- brute-force Hankel matrices of the Hamiltonian-path series on G(k,w), independent of the automaton:
entries H(u r) for prefix index-words u and suffix index-words r, computed by decoding u r with the auditor's decoder
and counting with the subset DP (hp.exe batch).  Rank mod two primes (lower bounds on the rational rank) and, when the
matrix is small enough, the exact rational rank.
Usage: python3 -I hankel_bf.py k w Lp Ls [seam_c]
  seam_c: restrict rows to prefixes whose cut has exactly seam_c pending slots (the seam block).
"""
import sys, os, itertools
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from decode import letters, decode, adjacency, batch_counts

def partial_decode(k, w, word, pending0=None):
    """returns pending list after reading word from the given pending state, or None if invalid"""
    pending = list(pending0 or []); groups = 100
    for t, li in enumerate(word):
        c = len(pending); L = letters(c, w, k)
        if li < 0 or li >= len(L): return None
        T, b, A = L[li]
        gs = [pending[i][1] for i in range(c) if T >> i & 1]
        if len(gs) != len(set(gs)): return None
        newp = [pending[i] for i in range(c) if not (T >> i & 1)]
        for j in range(b): newp.append((t, groups))
        if b: groups += 1
        pending = newp
    return pending

def enum_words(k, w, L, pending0, end_c=None, require_complete=False):
    out = []
    def rec(word, pending):
        if len(word) > 0 or not require_complete:
            if (end_c is None or len(pending) == end_c) and (not require_complete or len(pending) == 0):
                out.append(tuple(word))
        if len(word) == L: return
        c = len(pending); Ls = letters(c, w, k)
        for li, (T, b, A) in enumerate(Ls):
            gs = [pending[i][1] for i in range(c) if T >> i & 1]
            if len(gs) != len(set(gs)): continue
            newp = [pending[i] for i in range(c) if not (T >> i & 1)]
            for j in range(b): newp.append((len(word), 500 + len(word)))
            rec(word + [li], newp)
    rec([], list(pending0))
    return out

def rank_mod(M, p):
    A = [[x % p for x in row] for row in M]
    n = len(A); m = len(A[0]) if A else 0; r = 0; row = 0
    for col in range(m):
        if row >= n: break
        piv = None
        for i in range(row, n):
            if A[i][col]: piv = i; break
        if piv is None: continue
        A[row], A[piv] = A[piv], A[row]
        inv = pow(A[row][col], p - 2, p)
        A[row] = [(x * inv) % p for x in A[row]]
        prow = A[row]
        for i in range(n):
            if i != row and A[i][col]:
                f = A[i][col]; Ai = A[i]
                A[i] = [(a - f * b) % p for a, b in zip(Ai, prow)]
        row += 1; r += 1
    return r

def rank_exact(M):
    rows = [{j: Fraction(v) for j, v in enumerate(r) if v} for r in M]
    basis = []
    for v in rows:
        for piv, row in basis:
            if piv in v:
                f = v[piv] / row[piv]
                for k2, val in row.items():
                    nv = v.get(k2, Fraction(0)) - f * val
                    if nv == 0: v.pop(k2, None)
                    else: v[k2] = nv
        if v: basis.append((min(v), v))
    return len(basis)

def main():
    empty0 = '--empty0' in sys.argv; sys.argv = [a for a in sys.argv if a != '--empty0']
    k, w, Lp, Ls = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    seam = int(sys.argv[5]) if len(sys.argv) > 5 else None
    prefixes = enum_words(k, w, Lp, [], end_c=seam)
    if seam is not None: prefixes = [u for u in prefixes if len(u) > 0]
    # suffixes: index words complete from some start cut (all group patterns of the start cut)
    starts = []
    cvals = [seam] if seam is not None else list(range(w + 1))
    for c in cvals:
        for pat in set(itertools.product(range(c), repeat=c)) if c > 0 else [()]:
            # group pattern: slots i with group pat[i]; only patterns where groups appear as intervals in order
            pend = [(-1 - i, 1000 + pat[i]) for i in range(c)]
            starts.append(pend)
    suf = set()
    for pend in starts:
        for r in enum_words(k, w, Ls, pend, require_complete=True): suf.add(r)
    suffixes = sorted(suf, key=lambda r: (len(r), r))
    # entries
    cache = {}; items = []; keys = []
    for u in prefixes:
        for r in suffixes:
            word = list(u) + list(r)
            if tuple(word) in cache: continue
            dec = decode(k, w, word) if word else None
            if word == []: cache[()] = 0 if empty0 else 1; continue
            if dec is None: cache[tuple(word)] = 0; continue
            arcs, apex = dec
            if len(word) + k > 16: cache[tuple(word)] = None; continue
            keys.append(tuple(word)); items.append(adjacency(k, len(word), arcs, apex))
    hs = batch_counts(items)
    for key, h in zip(keys, hs): cache[key] = h
    M = []; dropped = 0
    for u in prefixes:
        row = []
        for r in suffixes:
            v = cache[tuple(list(u) + list(r))]
            if v is None: v = 0; dropped += 1
            row.append(v)
        M.append(row)
    M = [row for row in M if any(row)]
    cols = [j for j in range(len(suffixes)) if any(row[j] for row in M)]
    M = [[row[j] for j in cols] for row in M]
    r1 = rank_mod(M, 2147483647); r2 = rank_mod(M, 2147483629)
    ex = rank_exact(M) if len(M) <= 400 and len(M[0]) <= 400 else None
    print(f"G({k},{w}) brute-force Hankel{' seam c=' + str(seam) if seam is not None else ''}: {len(prefixes)} prefixes (len<={Lp}), {len(suffixes)} suffixes (len<={Ls}), nonzero matrix {len(M)}x{len(M[0])}, {len(keys)} DPs, {dropped} entries dropped (too large); rank mod p1 = {r1}, mod p2 = {r2}, exact = {ex}")

main()

"""For each order-12 escaper: is there a vertex v such that T - v is a substitution OF FIRST CREATORS, i.e. T - v is
strong and non-prime, its quotient Q over the maximal modules and every module M are first creators of their orders
(a single vertex counts as the first creator of 1)?  Cross-tabulated with the re-connected-sum description
(two first creators joined by an ordinal sum and re-connected by 3 arcs).  (genesis lane; exact per class)"""
import sys
from collections import Counter
from verify_escapers import adj_from, ham
from spectra import in_S

def sub_out(out, verts):
    vs = [v for v in range(len(out)) if verts >> v & 1]
    idx = {v: i for i, v in enumerate(vs)}
    o = [0] * len(vs)
    for v in vs:
        for w in vs:
            if out[v] >> w & 1: o[idx[v]] |= 1 << idx[w]
    return o

def is_fc(o):
    k = len(o)
    if k == 1: return True
    return in_S(k - 1, ham(o, (1 << k) - 1)) is False

def strong(o):
    k = len(o); full = (1 << k) - 1
    r = 1; prev = 0
    while r != prev:
        prev = r
        for v in range(k):
            if r >> v & 1: r |= o[v]
    if r != full: return False
    inn = [0] * k
    for v in range(k):
        for w in range(k):
            if o[v] >> w & 1: inn[w] |= 1 << v
    r = 1; prev = 0
    while r != prev:
        prev = r
        for v in range(k):
            if r >> v & 1: r |= inn[v]
    return r == full

def maximal_modules(o):
    k = len(o); full = (1 << k) - 1
    mods = []
    for M in range(1, full):
        if bin(M).count('1') < 2: continue
        ok = True
        for x in range(k):
            if M >> x & 1: continue
            t = o[x] & M
            if t and t != M: ok = False; break
        if ok: mods.append(M)
    maxi = [M for M in mods if not any(M != N and (M & N) == M for N in mods)]
    covered = 0
    for M in maxi: covered |= M
    parts = maxi + [1 << v for v in range(k) if not covered >> v & 1]
    return parts

def substitution_of_fcs(o):
    if not strong(o): return False
    parts = maximal_modules(o)
    if len(parts) == len(o): return False          # prime
    # parts must partition V for a strong tournament
    if sum(bin(p).count('1') for p in parts) != len(o): return None
    reps = [(p & -p).bit_length() - 1 for p in parts]
    Q = [0] * len(parts)
    for i, a in enumerate(reps):
        for j, b in enumerate(reps):
            if o[a] >> b & 1: Q[i] |= 1 << j
    if not is_fc(Q): return False
    for p in parts:
        if not is_fc(sub_out(o, p)): return False
    return True

def reconnected(out, n):
    full = (1 << n) - 1
    for S in range(1, full):
        k = bin(S).count('1')
        if k < 2 or k > n - 2: continue
        c = sum(bin(out[v] & ~S & full).count('1') for v in range(n) if S >> v & 1)
        if c != 3: continue
        if is_fc(sub_out(out, full & ~S)) and is_fc(sub_out(out, S)): return True
    return False

tab = Counter(); examples = {}
for line in open(sys.argv[1]):
    n, code = line.split()[:2]; n = int(n)
    out = adj_from(code, n); full = (1 << n) - 1
    a = reconnected(out, n)
    b = any(substitution_of_fcs(sub_out(out, full & ~(1 << v))) for v in range(n))
    tab[(a, b)] += 1
    if (a, b) not in examples: examples[(a, b)] = ham(out, full)
print('(re-connected sum of two first creators, substitution of first creators plus a vertex) -> classes, example H')
for k, v in sorted(tab.items()): print('  ', k, v, examples[k])

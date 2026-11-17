"""Wave 57: toothpick caps, the orphan 63 as a simultaneous doubling, and a first-creator reproduction.
Needs w57_reps.txt (w57_classes.py), w57_pairs.txt + w57_first8.txt (w57_capdata.c), w57_first9.txt + w57_spec9.txt
(w57_first9.c).  Prior art: THM-3121 (substitution kernel), THM-5172 / wave 25 (cyclic join), THM-4876 (6) (one-centre
counts), THM-4748 (Fibonacci fan), THM-4889 (census, least orders), THM-4922 (first creators, genealogy, orphan 63)."""
import itertools, collections
from w45_hk import H_of, sub

def F2(n, out):
    """ordered two-path covers sum_S H(X[S]) H(X[V-S]) (THM-3166's chi_X(2))."""
    tot = 0
    for S in range(1 << n):
        a = [v for v in range(n) if S >> v & 1]; b = [v for v in range(n) if not S >> v & 1]
        ha = H_of(*sub(out, a)) if a else 1; hb = H_of(*sub(out, b)) if b else 1
        tot += ha * hb
    return tot

def substitute(qn, qout, blocks):
    """Q[S_1..S_q]; blocks = list of (n_i, out_i)."""
    off = [0]
    for (m, _) in blocks: off.append(off[-1] + m)
    N = off[-1]; out = [0] * N
    for i, (m, o) in enumerate(blocks):
        for a in range(m):
            v = off[i] + a
            for b in range(m):
                if o[a] >> b & 1: out[v] |= 1 << (off[i] + b)
            for j in range(len(blocks)):
                if j != i and qout[i] >> j & 1:
                    for b in range(blocks[j][0]): out[v] |= 1 << (off[j] + b)
    return N, out

TT = lambda k: (k, [sum(1 << j for j in range(i + 1, k)) for i in range(k)])
ONE = (1, [0])
C3 = (3, [0b010, 0b100, 0b001])

def one_centre(w):
    """TT_m + u, u -> c_i iff w_i = 1 (chain c_0 -> c_1 -> ...); u is the last vertex."""
    m = len(w); out = [sum(1 << j for j in range(i + 1, m)) for i in range(m)] + [0]
    for i, b in enumerate(w):
        if b: out[m] |= 1 << i
        else: out[i] |= 1 << m
    return m + 1, out

def converse(n, out):
    return [sum(1 << w for w in range(n) if out[w] >> v & 1) for v in range(n)]

def canon(n, out):
    return min(tuple(out[p[i]] >> p[j] & 1 for i in range(n) for j in range(i + 1, n)) for p in itertools.permutations(range(n)))

def modules(n, out):
    res = []
    for m in range(1, (1 << n) - 1):
        if bin(m).count('1') < 2: continue
        if all(len({out[x] >> y & 1 for y in range(n) if m >> y & 1}) == 1 for x in range(n) if not m >> x & 1):
            res.append([v for v in range(n) if m >> v & 1])
    return res

reps = collections.defaultdict(list)
for line in open('w57_reps.txt'):
    t = list(map(int, line.split())); reps[t[0]].append(t[1:])

print("== A. Toothpick caps C3[TT_k, X, 1]:  H = M_k (F2(X) - H(X)) + 2 H(X),  M_k = 2^k - 1  (k = 1: H + F2) ==")
bad = cnt = 0; capcanon = collections.defaultdict(set)
for n in range(1, 6):
    for out in reps[n]:
        h = H_of(n, out); f = F2(n, out)
        for k in (1, 2, 3):
            N, o = substitute(3, C3[1], [TT(k), (n, out), ONE])
            cnt += 1
            if H_of(N, o) != ((1 << k) - 1) * (f - h) + 2 * h: bad += 1
        if n <= 4:
            N, o = substitute(3, C3[1], [ONE, (n, out), ONE]); capcanon[n].add(canon(N, o))
print(f"  checked {cnt} (class, k) cases with X of order <= 5, k <= 3: mismatches {bad}")
print("  cap X -> C3[1, X, 1] is injective on classes (orders 1-4):", {n: (len(capcanon[n]), len(reps[n])) for n in capcanon})

print("\n== B. The orphan 63: the Fibonacci fan F_6 is a simultaneous TT2-doubling of a 9-path quotient ==")
Q = one_centre([1, 0, 1, 0])
print(f"  Q = TT4 + u with contact word 1010: H = {H_of(*Q)}, self-converse: {canon(Q[0], converse(*Q)) == canon(*Q)}, "
      f"nontrivial modules: {modules(*Q)}")
# Q's vertices: c0..c3 (chain) and u = 4.  Double c0, c1, c2 (blocks in vertex order c0,c1,c2,c3,u)
blocks = [TT(2), TT(2), TT(2), ONE, ONE]
F6 = substitute(5, Q[1], blocks)
W = one_centre([1, 1, 0, 0, 1, 1, 0])
print(f"  Q[TT2, TT2, TT2, 1, u] has H = {H_of(*F6)}; equal to the one-centre word 11.00.11.0: {canon(*F6) == canon(*W)}")
Wc = one_centre([1, 0, 0, 1, 1, 0, 0])
print(f"  converse = doubling the last three chain vertices = word 1.00.11.00: {canon(Wc[0], converse(*Wc)) == canon(*W)}, H = {H_of(*Wc)}")
print(f"  modules of F_6: {modules(*W)}")
dels = sorted({H_of(*sub(W[1], [u for u in range(8) if u != v])) for v in range(8)})
print(f"  deletion values of F_6: {dels} (all first created by order 6; THM-4922 orphan)")
print("  doubling subsets of Q's chain (one TT2 per chosen chain vertex):")
for S in range(16):
    bl = [TT(2) if S >> i & 1 else ONE for i in range(4)] + [ONE]
    T = substitute(5, Q[1], bl)
    print(f"    double {[i for i in range(4) if S >> i & 1]}: order {T[0]}, H = {H_of(*T)}")
def onecentre_count(w):
    """THM-4876 (6): H = 1 + #10 + sum over 1 before 0 of 2^(gap-1)... in word form."""
    acc = tot = 0
    for b in w:
        if b == 0: tot += acc
        acc = 2 * acc + b
    return 1 + sum(1 for i in range(len(w) - 1) if w[i] == 1 and w[i + 1] == 0) + tot
fan = []
for k in range(4, 13):
    a = k // 2; b = k + 1 - a; w = [0 if i in (a, b, k + 1) else 1 for i in range(1, k + 2)]
    fan.append((k, onecentre_count(w), (1 << k) - 1))
print("  Fibonacci fans (THM-4748: apex loses to a, b, k+1): (k, H, 2^k - 1):", fan)

print("\n== C. How sparse the toothpick caps are (cores of order <= 8; exact below 513) ==")
pairs = set(); spec = collections.defaultdict(set); sspec = collections.defaultdict(set)
for line in open('w57_pairs.txt'):
    n, h, f, s = map(int, line.split()); pairs.add((h, f)); spec[n].add(h)
    if s: sspec[n].add(h)
B = 20000
def fam(kmax, smax):
    got = set()
    for h, f in pairs:
        for s in range(smax + 1):
            Fs = f << s
            if Fs + h > 3 * B: break
            for k in range(1, kmax + 1):
                N = ((1 << k) - 1) * (Fs - h) + 2 * h
                if N > B: break
                got.add(N)
    return got
for label, km, sm in (("caps H + F2", 1, 0), ("caps with doubling rays H + 2^s F2", 1, 30), ("toothpick caps M_k(2^s F2 - H) + 2H", 30, 30)):
    g = fam(km, sm); miss = [N for N in range(1, B + 1, 2) if N not in g]
    print(f"  {label}: misses {sum(1 for N in miss if N < 513)} of 256 odd N < 513 (exact); {len(miss)} of {B // 2} odd N <= {B}; first {miss[:12]}")

print("\n== D. First creators, reproduced (THM-4889, THM-4922) ==")
for line in open('w57_spec9.txt'):
    h, s = map(int, line.split()); spec[9].add(h)
    if s: sspec[9].add(h)
first = {}
for n in sorted(spec):
    for h in sorted(spec[n]): first.setdefault(h, n)
M = {n: max(spec[n]) for n in spec}
print("  max H(n):", [M[n] for n in sorted(M)], "| strong floor f(n), n >= 3:", [min(sspec[n]) for n in sorted(sspec) if n >= 3])
print("  new values per order:", [sum(1 for h in spec[n] if first[h] == n) for n in sorted(spec)])
for n in range(6, 10):
    late = [h for h in range(1, M[n] + 1, 2) if h not in spec[n] and h not in (7, 21)]
    print(f"  order {n}: odd values below max H({n}) not attained: {len(late)}" + (f" {late}" if len(late) <= 16 else f" in [{min(late)}, {max(late)}]")
          + f"; all attained at order {n + 1}: {all(h in spec[n + 1] for h in late) if n + 1 in spec else 'n/a'}")
for fn, n in (("w57_first8.txt", 8), ("w57_first9.txt", 9)):
    orph = []
    for line in open(fn):
        t = [x for x in line.split() if x != 'big']; h = int(t[0]); par = list(map(int, t[2:]))
        if int(t[1]) != 1: print("   non-strong creator of a new value at", n, h)
        if not any(first.get(p) == n - 1 for p in par): orph.append((h, par))
    print(f"  order {n}: new values with no deletion parent new at order {n - 1}: {orph}")
W = {h: [n for n in sorted(sspec) if n >= 3 and min(sspec[n]) <= h <= M[n]] for h in range(3, 31, 2)}
print("  windows {n : f(n) <= h <= max H(n)}:", W)

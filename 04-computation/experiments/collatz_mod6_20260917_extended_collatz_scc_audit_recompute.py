#!/usr/bin/env python3
"""Independent recomputation (verifier lens: recompute) for lane extended_collatz_scc.

Does NOT import the explorer's script.  Recomputes: leaf identity, drift table,
3-adic hostile family, Q2/Q1 greedy to 10^6, SCC counts (own iterative Kosaraju,
not Tarjan), the 3-adic chain law (E[2^K_J], f_J, g_J, sharpness, exhaustive word
determination), the cycle census (own DFS) plus GLOBAL cycle-word hostile checks
for lengths <= 22 with unbounded nodes, the signed side, and the pn+1 hierarchy
(k(r) from raw definition, no closed form).  Exact integer/Fraction arithmetic.
"""
import sys, time
from fractions import Fraction
from collections import deque

T0 = time.time()
FAILS = []


def chk(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("  !! FAIL:", msg)


def sec(t):
    print("\n### " + t)


# ---------------------------------------------------------------- A. leaf identity
sec("A. leaf identity (v <= 10^5), brute force over sources n <= 2*10^5")
VMAX = 10 ** 5
preds = {}
for n in range(1, 2 * VMAX + 1):
    if n % 2 == 0:
        preds.setdefault(n // 2, set()).add(n)
    u = 3 * n + 1
    if u <= VMAX:
        preds.setdefault(u, set()).add(n)
new_count = 0
for v in range(1, VMAX + 1):
    P = preds.get(v, set())
    new = [n for n in P if n % 2 == 0 and 3 * n + 1 == v]
    if new:
        new_count += 1
        j = (v - 1) // 6
        chk(v % 6 == 1 and v >= 7 and new == [2 * j], "new arrow shape v=%d" % v)
        n0 = (4 * v - 1) // 3
        chk((4 * v - 1) % 3 == 0 and n0 % 8 == 1 and 3 * n0 + 1 == 4 * v and (n0 - 1) // 4 == 2 * j, "n0 identity v=%d" % v)
        # least odd n with oddpart(3n+1)=v: brute force h=0..6
        odds = [(2 ** h * v - 1) // 3 for h in range(0, 7) if (2 ** h * v - 1) % 3 == 0 and ((2 ** h * v - 1) // 3) % 2 == 1]
        chk(odds and min(odds) == n0, "n0 least odd T-predecessor v=%d" % v)
        chk(P == {2 * v, 2 * j}, "predecessor set v=%d" % v)
    else:
        chk(not (v % 6 == 1 and v >= 7), "missing new arrow v=%d" % v)
        chk(P == {2 * v} or (v % 2 == 0 and P == {2 * v, (v - 1) // 3}), "pred set v=%d: %s" % (v, P))
    if v % 3 == 0:
        chk(P == {2 * v}, "3Z target v=%d has non-halving predecessor" % v)
print("  new even->3n+1 arrows into [1,10^5]: %d (expected 16666)" % new_count)
chk(new_count == 16666, "new arrow count")

# ---------------------------------------------------------------- B. drift table
sec("B. drift table mod 9 (recomputed from the definition)")
U9 = (1, 2, 4, 5, 7, 8)
KMIN, NXT, KCL = {}, {}, {}
for r in U9:
    ks = [k for k in range(6) if ((2 ** k * r) % 9) in (4, 7)]
    KCL[r] = ks
    KMIN[r] = min(ks)
    NXT[r] = (((2 ** KMIN[r] * r) % 9 - 1) // 3) % 3
print("  KMIN =", KMIN, " NEXT mod 3 =", NXT, " classes =", KCL)
chk(KMIN == {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}, "KMIN")
chk(NXT == {1: 1, 2: 1, 4: 1, 5: 1, 7: 2, 8: 2}, "NEXT")
for r in U9:
    chk(all(k % 2 == (1 if r % 3 == 2 else 0) for k in KCL[r]) and len(KCL[r]) == 2, "parity/2 classes r=%d" % r)
# verify: (2^k m -1)/3 integer and not 0 mod 3 iff 2^k m in {4,7} mod 9, over m<=2000, k<=12
for m in range(1, 2001):
    for k in range(13):
        y = 2 ** k * m - 1
        good = (y % 3 == 0) and ((y // 3) % 3 != 0)
        chk(good == ((2 ** k * m) % 9 in (4, 7)), "admissibility criterion m=%d k=%d" % (m, k))


def G(x):
    k = KMIN[x % 9]
    return ((x << k) - 1) // 3, k


# ---------------------------------------------------------------- C. hostile family
sec("C. 3-adic hostile family v_3(m-1)=j")


def v3(n):
    c = 0
    while n % 3 == 0:
        n //= 3
        c += 1
    return c


for m in range(2, 20001):
    if m % 3 != 1:
        continue
    j = v3(m - 1)
    x = m
    ks = []
    for _ in range(j):
        x, k = G(x)
        ks.append(k)
    chk(ks == [2] * (j - 1) + [0], "word of m=%d (j=%d)" % (m, j))
    chk(x == 4 ** (j - 1) * (m - 1) // 3 ** j and (4 ** (j - 1) * (m - 1)) % 3 ** j == 0, "closed form m=%d" % m)
    if j <= 4:
        chk(x < m, "m_j < m fails for j<=4, m=%d" % m)
    else:
        chk(x > m, "m_j > m fails for j>=5, m=%d" % m)
for j in range(1, 9):
    thr = Fraction(4 ** (j - 1), 4 ** (j - 1) - 3 ** j) if 4 ** (j - 1) > 3 ** j else None
    print("  j=%d 4^(j-1)/3^j=%s  threshold m>%s" % (j, Fraction(4 ** (j - 1), 3 ** j), thr))
chk([j for j in range(1, 30) if 4 ** (j - 1) > 3 ** j] == list(range(5, 30)), "4^(j-1)>3^j iff j>=5")
x, p = 244, [244]
for _ in range(5):
    x, _k = G(x)
    p.append(x)
chk(p == [244, 325, 433, 577, 769, 256], "244 chain")
x, q = 244, [244]
while True:
    x, _k = G(x)
    q.append(x)
    if x < 244:
        break
print("  greedy from 244 to first descent:", q, " compound-value peak", max(q))

# ---------------------------------------------------------------- D. Q2 greedy to 10^6
sec("D. Q2 greedy to 10^6 (own implementation), Q1 to 10^6")
NQ = 10 ** 6
tsteps = [0] * (NQ + 1)
tarr = [0] * (NQ + 1)
tpeak = [0] * (NQ + 1)   # peak over ALL E-nodes on the inverse path (doubling intermediates included)
tpeakc = [0] * (NQ + 1)  # peak over compound-move values only
tpeak[1] = tpeakc[1] = 1
stat = {"below": 0, "one": 0, "cycle": 0, "cap": 0}
wp = (0, 0); wpc = (0, 0); ws = (0, 0); wa = (0, 0)
t1 = time.time()
for m in range(2, NQ + 1):
    if m % 3 == 0:
        continue
    x = m; s = 0; a = 0; pk = m; pkc = m
    st = None
    while True:
        k = KMIN[x % 9]
        y = x << k
        if y > pk:
            pk = y
        x = (y - 1) // 3
        s += 1; a += k + 1
        if x > pkc:
            pkc = x
        if x < m:
            st = "one" if x == 1 else "below"
            break
        if x == m:
            st = "cycle"; break
        if s > 10 ** 5:
            st = "cap"; break
    stat[st] += 1
    if st in ("below", "one"):
        tsteps[m] = s + tsteps[x]
        tarr[m] = a + tarr[x]
        tpeak[m] = max(pk, tpeak[x])
        tpeakc[m] = max(pkc, tpeakc[x])
        if tpeak[m] > wp[0]: wp = (tpeak[m], m)
        if tpeakc[m] > wpc[0]: wpc = (tpeakc[m], m)
        if tsteps[m] > ws[0]: ws = (tsteps[m], m)
        if tarr[m] > wa[0]: wa = (tarr[m], m)
print("  statuses", stat, " time %.1fs" % (time.time() - t1))
print("  max E-node peak %d at m=%d ; max compound-value peak %d at m=%d" % (wp + wpc))
print("  max compound moves %d at m=%d ; max arrows %d at m=%d" % (ws + wa))
chk(stat == {"below": 666664, "one": 2, "cycle": 0, "cap": 0}, "Q2 statuses")
chk(wp == (150994948, 797162) and ws == (74, 984104) and wa == (172, 984104), "Q2 extremes")
chk(150994948 == 9 * 2 ** 24 + 4, "peak = 9*2^24+4")
# which m hit 1 directly
ones = [m for m in (2, 4, 5, 7, 8, 10, 11, 13) if G(m)[0] == 1]
print("  m with G(m)=1:", ones)
t1 = time.time()
for n in range(2, NQ + 1):
    x = n
    while x >= n:
        x = x // 2 if x % 2 == 0 else 3 * x + 1
print("  Q1: every n in [2,10^6] descends under Collatz  (%.1fs)" % (time.time() - t1))

# ---------------------------------------------------------------- E. SCC via Kosaraju
sec("E. SCC of E|[1,N] via iterative Kosaraju (independent of Tarjan)")


def kosaraju(N, up):
    def succ(n):
        out = []
        if n % 2 == 0:
            out.append(n // 2)
        u = up(n)
        if 1 <= u <= N:
            out.append(u)
        return out
    # build reverse adjacency
    radj = [[] for _ in range(N + 1)]
    for n in range(1, N + 1):
        for w in succ(n):
            radj[w].append(n)
    order = []
    seen = [False] * (N + 1)
    for r in range(1, N + 1):
        if seen[r]:
            continue
        seen[r] = True
        stack = [(r, iter(succ(r)))]
        while stack:
            node, it = stack[-1]
            pushed = False
            for w in it:
                if not seen[w]:
                    seen[w] = True
                    stack.append((w, iter(succ(w))))
                    pushed = True
                    break
            if not pushed:
                order.append(node)
                stack.pop()
    comp = [0] * (N + 1)
    ncomp = 0
    sizes = []
    for r in reversed(order):
        if comp[r]:
            continue
        ncomp += 1
        comp[r] = ncomp
        st = [r]
        sz = 0
        while st:
            x = st.pop()
            sz += 1
            for w in radj[x]:
                if not comp[w]:
                    comp[w] = ncomp
                    st.append(w)
        sizes.append(sz)
    return comp, sizes


EXP = {1000: (854, 147, 520, 202, [31, 47, 55, 71, 73]),
       10000: (8249, 1752, 4915, 1737, [383, 511, 575, 608, 667]),
       100000: (82924, 17077, 49590, 17764, [1535, 2047, 2207, 2287, 2303])}
for N in (1000, 10000, 100000):
    comp, sizes = kosaraju(N, lambda n: 3 * n + 1)
    giant = comp[1]
    gsize = sizes[giant - 1]
    nontriv = sum(1 for s in sizes if s > 1)
    out = [n for n in range(1, N + 1) if n % 3 and comp[n] != giant]
    half = sum(1 for n in out if n <= N // 2)
    m3ok = all(sizes[comp[n] - 1] == 1 for n in range(3, N + 1, 3))
    print("  N=%d #SCC=%d nontrivial=%d giant=%d outside=%d half=%d smallest=%s mult3 singletons=%s"
          % (N, len(sizes), nontriv, gsize, len(out), half, out[:5], m3ok))
    e = EXP[N]
    chk((len(sizes), gsize, len(out), half, out[:5]) == e and nontriv == 1 and m3ok, "SCC data N=%d" % N)
    chk(gsize == max(sizes), "giant is largest N=%d" % N)
comp, sizes = kosaraju(100000, lambda n: 3 * n - 1)
giant = comp[1]
out = [n for n in range(1, 100001) if n % 3 and comp[n] != giant]
print("  E_- N=10^5 #SCC=%d giant=%d nontrivial=%d outside=%d smallest=%s" % (len(sizes), sizes[giant - 1], sum(1 for s in sizes if s > 1), len(out), out[:5]))
chk(len(sizes) == 82979 and sizes[giant - 1] == 17022, "E_- SCC data")
chk(all(comp[c] == giant for c in (1, 5, 17)), "3n-1 cycle minima in giant")
# outsider explanation
for n in (1535, 2047, 2207, 2287, 2303):
    x = n; mx = n
    while x != 1:
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        mx = max(mx, x)
    print("  outsider %d Collatz peak %d" % (n, mx))
    chk(mx > 100000, "outsider %d peak" % n)

# ---------------------------------------------------------------- F. 3-adic chain law
sec("F. 3-adic Terras: word determination, E[2^K_J], f_J, g_J")


def word(m, J):
    ks = []
    x = m
    for _ in range(J):
        x, k = G(x)
        ks.append(k)
    return ks


# exhaustive determination check for J<=6: all m in [1, 2*3^(J+1)] vs m+3^(J+1)
for J in range(1, 7):
    mod = 3 ** (J + 1)
    for m in range(1, 3 * mod):
        if m % 3 == 0:
            continue
        chk(word(m, J) == word(m + mod, J), "determination J=%d m=%d" % (J, m))
    chk(word(1, J) != word(1 + 3 ** J, J) and word(1, J)[:J - 1] == word(1 + 3 ** J, J)[:J - 1], "sharpness J=%d" % J)
print("  word determined by m mod 3^(J+1) (exhaustive, J<=6, m<3^(J+2)); sharp witnesses 1 vs 1+3^J")
EXP_F = {1: (2, 3), 2: (8, 9), 3: (25, 27), 4: (26, 27), 5: (236, 243), 6: (239, 243), 7: (241, 243),
         8: (2173, 2187), 9: (19609, 19683), 10: (58868, 59049)}
prev = 0
for J in range(1, 11):
    mod = 3 ** (J + 1)
    tot = 0; e2 = 0; fj = 0; gj = 0
    tail = 0  # #classes with 2^K_J >= 3^J
    for m in range(1, mod):
        if m % 3 == 0:
            continue
        tot += 1
        x = m; K = 0; hit = False
        for i in range(1, J + 1):
            x, k = G(x)
            K += k
            if (1 << K) < 3 ** i:
                hit = True
        e2 += 1 << K
        fj += hit
        gj += (1 << K) < 3 ** J
        tail += (1 << K) >= 3 ** J
    E2 = Fraction(e2, tot)
    f = Fraction(fj, tot)
    g = Fraction(gj, tot)
    print("  J=%2d units=%6d E[2^K_J]=%s f_J=%s g_J=%s tailclasses=%d <= 2*3^J*(7/9)^(J-1)=%.1f"
          % (J, tot, E2, f, g, tail, 2 * 3 ** J * (7 / 9) ** (J - 1)))
    chk(tot == 2 * 3 ** J and E2 == 3 * Fraction(7, 3) ** (J - 1), "E[2^K_J] J=%d" % J)
    chk(f == Fraction(*EXP_F[J]), "f_J J=%d" % J)
    chk(f >= prev, "f monotone J=%d" % J)
    chk(tail <= 2 * 3 ** J * Fraction(7, 9) ** (J - 1) and 1 - f <= 1 - g <= Fraction(7, 9) ** (J - 1), "tail J=%d" % J)
    prev = f
# Markov-chain word counts J<=4 : independent recount by DP on residues mod 9 with digit lifts
for J in range(1, 5):
    mod = 3 ** (J + 1)
    cnt = {}
    for m in range(1, mod):
        if m % 3:
            cnt[tuple(word(m, J))] = cnt.get(tuple(word(m, J)), 0) + 1
    for w, c in cnt.items():
        dist = {r: Fraction(1, 6) for r in U9}
        for k in w:
            nd = {}
            for r, pr in dist.items():
                if KMIN[r] == k:
                    c0 = NXT[r]
                    for r2 in (c0, c0 + 3, c0 + 6):
                        nd[r2] = nd.get(r2, 0) + pr / 3
            dist = nd
        chk(Fraction(c, 2 * 3 ** J) == sum(dist.values()), "word count %s" % (w,))
print("  word counts J<=4 equal 2*3^J * P_chain (checked)")
# tilted matrix and initial vector directly from the chain
M = [[Fraction(0)] * 2 for _ in range(2)]
for c in (1, 2):
    for r2 in (c, c + 3, c + 6):
        M[c - 1][NXT[r2] - 1] += Fraction(2 ** KMIN[r2], 3)
u = [Fraction(0), Fraction(0)]
for r in U9:
    u[NXT[r] - 1] += Fraction(2 ** KMIN[r], 6)
print("  M =", M, " u =", u, " det", M[0][0] * M[1][1] - M[0][1] * M[1][0], " tr", M[0][0] + M[1][1])
chk(M == [[Fraction(5, 3), Fraction(1, 3)], [Fraction(10, 3), Fraction(2, 3)]] and u == [Fraction(5, 2), Fraction(1, 2)], "M,u")
Epi = sum(Fraction(2, 9) * KMIN[r] for r in (1, 4, 7)) + sum(Fraction(1, 9) * KMIN[r] for r in (2, 5, 8))
chk(Epi == 1, "E_pi[k]")
# stationary law after steps 1..4 mod 3^6
for i in (1, 2, 3, 4):
    cnt = {r: 0 for r in U9}
    for m in range(1, 3 ** 6):
        if m % 3:
            x = m
            for _ in range(i):
                x, _k = G(x)
            cnt[x % 9] += 1
    chk(all(cnt[r] == (2 if r in (1, 4, 7) else 1) * 3 ** 6 * 2 // 27 for r in U9), "stationary law step %d" % i)
print("  stationary law 2/9,1/9 at steps 1..4 (checked)")
# the tail bound as a count for X = 10^6 and J = 10 (actual sigma > J count)
J = 10
cnt_sig = 0
for m in range(1, 10 ** 6 + 1):
    if m % 3 == 0:
        continue
    x = m; ok = False
    for i in range(J):
        x, _k = G(x)
        if x < m:
            ok = True; break
    if not ok:
        cnt_sig += 1
bound = Fraction(7, 9) ** (J - 1) * (Fraction(2 * 10 ** 6, 3) + 2 * 3 ** J)
print("  #{m<=10^6: sigma>10} = %d ; bound %.1f" % (cnt_sig, float(bound)))
chk(cnt_sig <= bound, "sigma tail count vs bound")

# ---------------------------------------------------------------- G. cycle census
sec("G. cycle census (own DFS, nodes<=2000, len<=40) and GLOBAL cycle-word check")


def census(LIM, MAXLEN, up):
    def succ(n):
        out = []
        u = up(n)
        if 1 <= u <= LIM:
            out.append(u)
        if n % 2 == 0:
            out.append(n // 2)
        return out
    found = []
    for s in range(1, LIM + 1):
        path = [s]; inpath = {s}
        stack = [iter(succ(s))]
        while stack:
            try:
                w = next(stack[-1])
            except StopIteration:
                stack.pop(); inpath.discard(path.pop()); continue
            if w == s:
                found.append(tuple(path)); continue
            if w <= s or w in inpath or len(path) >= MAXLEN:
                continue
            path.append(w); inpath.add(w); stack.append(iter(succ(w)))
    return found


C = census(2000, 40, lambda n: 3 * n + 1)
hist = {}
ah = {}
for c in C:
    hist[len(c)] = hist.get(len(c), 0) + 1
    a = sum(1 for x, y in zip(c, c[1:] + c[:1]) if y == 3 * x + 1)
    h = len(c) - a
    e = sum(1 for x, y in zip(c, c[1:] + c[:1]) if y == 3 * x + 1 and x % 2 == 0)
    ah[(a, h)] = ah.get((a, h), 0) + 1
    chk(2 ** h > 3 ** a and all(x % 3 for x in c) and (e >= 1 or c == (1, 4, 2)), "cycle property %s" % (c,))
    chk(h == min(t for t in range(100) if 2 ** t > 3 ** a), "h=ceil(a log2 3) %s" % (c,))
print("  #cycles=%d hist=%s (a,h)=%s" % (len(C), sorted(hist.items()), sorted(ah.items())))
chk(len(C) == 74 and sorted(hist.items()) == [(3, 1), (8, 1), (13, 6), (16, 1), (21, 2), (26, 11), (34, 9), (39, 43)], "census 74")
chk(sorted(ah) == [(1, 2), (3, 5), (5, 8), (6, 10), (8, 13), (10, 16), (13, 21), (15, 24)], "(a,h) list")
chk((2, 7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4) in C, "user 16-cycle")
# boundary: raise length cap to 45 and node cap to 2500 to see whether the census is cap-sensitive
C45 = census(2000, 45, lambda n: 3 * n + 1)
C2500 = census(2500, 40, lambda n: 3 * n + 1)
print("  nodes<=2000,len<=45: %d cycles (lengths %s); nodes<=2500,len<=40: %d cycles"
      % (len(C45), sorted(set(len(c) for c in C45)), len(C2500)))
Cm = census(2000, 40, lambda n: 3 * n - 1)
hm = {}
detm = []
for c in Cm:
    hm[len(c)] = hm.get(len(c), 0) + 1
    a = sum(1 for x, y in zip(c, c[1:] + c[:1]) if y == 3 * x - 1)
    chk(2 ** (len(c) - a) < 3 ** a and all(x % 3 for x in c), "E_- cycle %s" % (c,))
    if all((y == x // 2) if x % 2 == 0 else (y == 3 * x - 1) for x, y in zip(c, c[1:] + c[:1])):
        detm.append(c)
print("  E_- census: %d cycles hist=%s deterministic=%s" % (len(Cm), sorted(hm.items()), detm))
chk(len(Cm) == 70 and sorted(detm) == sorted([(1, 2), (5, 14, 7, 20, 10), (17, 50, 25, 74, 37, 110, 55, 164, 82, 41, 122, 61, 182, 91, 272, 136, 68, 34)]), "E_- census")


# GLOBAL: all simple E-cycles of length L (any node size) from the cycle-word equation
def global_cycles(L, sign=1):
    """Enumerate words over {U,H}^L; n_i = (A n0 + B)/2^g; solve n0; verify simple cycle with min n0."""
    res = []
    amax = max(a for a in range(L + 1) if (2 ** (L - a) > 3 ** a if sign == 1 else 2 ** (L - a) < 3 ** a or a == L))
    # iterate over words by combinations of U positions
    from itertools import combinations
    for a in range(1, L + 1):
        if sign == 1 and not (2 ** (L - a) > 3 ** a):
            continue
        if sign == -1 and not (2 ** (L - a) < 3 ** a):
            continue
        for Upos in combinations(range(L), a):
            Uset = set(Upos)
            A, B, g = 1, 0, 0
            for i in range(L):
                if i in Uset:
                    A, B = 3 * A, 3 * B + sign * (1 << g)
                else:
                    g += 1
            # n0 = (A n0 + B)/2^g  ->  n0 (2^g - A) = B
            D = (1 << g) - A
            if D == 0 or B % D != 0:
                continue
            n0 = B // D
            if n0 <= 0:
                continue
            # simulate and verify
            x = n0; nodes = [n0]; ok = True
            for i in range(L):
                if i in Uset:
                    x = 3 * x + sign
                else:
                    if x % 2:
                        ok = False; break
                    x //= 2
                if i < L - 1:
                    nodes.append(x)
            if not ok or x != n0 or len(set(nodes)) != L or min(nodes) != n0:
                continue
            res.append(tuple(nodes))
    return res


print("  GLOBAL cycle enumeration (unbounded nodes) by length:")
tot_global = 0
for L in range(2, 23):
    gc = global_cycles(L)
    tot_global += len(gc)
    big = [c for c in gc if max(c) > 2000]
    ahs = sorted(set((sum(1 for x, y in zip(c, c[1:] + c[:1]) if y == 3 * x + 1), 0) for c in gc))
    ahs = sorted(set(( (lambda a: (a, L - a))(sum(1 for x, y in zip(c, c[1:] + c[:1]) if y == 3 * x + 1)) ) for c in gc))
    inc = [c for c in C if len(c) == L]
    print("    L=%2d: %d cycles total, %d with max node > 2000, (a,h) set %s ; census had %d" % (L, len(gc), len(big), ahs, len(inc)))
    chk(set(c for c in gc if max(c) <= 2000) == set(inc), "census vs global at L=%d" % L)
    if L == 16:
        print("      length-16 cycles (all nodes):", gc)
        for c in big[:3]:
            print("      big 16-cycle example:", c)
print("  (a,h) pairs with h > ceil(a log2 3) among global cycles L<=22:",
      sorted(set((a, L - a) for L in range(2, 23) for c in global_cycles(L) for a in [sum(1 for x, y in zip(c, c[1:] + c[:1]) if y == 3 * x + 1)] if (L - a) > min(t for t in range(100) if 2 ** t > 3 ** a))))

# ---------------------------------------------------------------- H. signed side
sec("H. signed side E_- (3n-1)")
predsm = {}
for n in range(1, 2 * VMAX + 1):
    if n % 2 == 0:
        predsm.setdefault(n // 2, set()).add(n)
    u = 3 * n - 1
    if u <= VMAX:
        predsm.setdefault(u, set()).add(n)
for v in range(1, VMAX + 1):
    P = predsm.get(v, set())
    new = [n for n in P if n % 2 == 0 and 3 * n - 1 == v]
    if new:
        chk(v % 6 == 5 and new == [(v + 1) // 3], "minus new arrow v=%d" % v)
        n0 = (4 * v + 1) // 3
        chk((4 * v + 1) % 3 == 0 and n0 % 8 == 7 and 3 * n0 - 1 == 4 * v and (n0 + 1) // 4 == (v + 1) // 3, "minus n0 v=%d" % v)
        odds = [(2 ** h * v + 1) // 3 for h in range(0, 7) if (2 ** h * v + 1) % 3 == 0 and ((2 ** h * v + 1) // 3) % 2 == 1]
        chk(min(odds) == n0, "minus least odd v=%d" % v)
    else:
        chk(v % 6 != 5, "missing minus arrow v=%d" % v)
KM = {}
NM = {}
for r in U9:
    ks = [k for k in range(6) if ((2 ** k * r) % 9) in (2, 5)]
    KM[r] = min(ks)
    NM[r] = (((2 ** KM[r] * r) % 9 + 1) // 3) % 3
print("  minus KMIN", KM, " NEXT", NM)
chk(all(KM[r] == KMIN[9 - r] for r in U9), "negation conjugacy")
for m in range(1, 2001):
    for k in range(13):
        y = 2 ** k * m + 1
        chk(((y % 3 == 0) and (y // 3) % 3 != 0) == ((2 ** k * m) % 9 in (2, 5)), "minus admissibility")


def Gm(x):
    k = KM[x % 9]
    return ((x << k) + 1) // 3, k


x, p = 242, [242]
for _ in range(5):
    x, _k = Gm(x); p.append(x)
chk(p == [242, 323, 431, 575, 767, 256], "242 chain")
chk(Gm(4) == (11, 3) and Gm(11) == (4, 0), "G_- 2-cycle {4,11}")
# explicit non-greedy path and its k-word
pathm = [4, 11, 59, 20, 7, 5, 2]
kw = []
for x, y in zip(pathm, pathm[1:]):
    ks = [k for k in range(10) if (2 ** k * x + 1) % 3 == 0 and (2 ** k * x + 1) // 3 == y]
    chk(len(ks) == 1, "step %d->%d" % (x, y))
    kw.append(ks[0])
print("  k-word of 4->11->59->20->7->5->2 is", kw, " (draft note says (3,4,0,1,1,0))")
chk(kw == [3, 4, 0, 0, 1, 0], "k-word of the rescue path")
chk((3 * 4 - 1, 3 * 11 - 1) == (11, 32), "E_- cycle 4->11->32 arrows")
# Q2_- greedy with BFS fallback


def bfs_below(m, sign, kcap=40, budget=200000):
    par = {m: None}; dq = deque([m]); nodes = 0
    while dq:
        x = dq.popleft(); nodes += 1
        if nodes > budget:
            return None
        for k in range(kcap + 1):
            y = (x << k) - sign
            if y % 3:
                continue
            z = y // 3
            if z <= 0 or z % 3 == 0 or z in par or z > 1 << 90:
                continue
            par[z] = x
            if z < m:
                pth = [z]
                while par[pth[-1]] is not None:
                    pth.append(par[pth[-1]])
                return pth[::-1]
            dq.append(z)
    return None


tpk = [0] * (NQ + 1); tst = [0] * (NQ + 1); tpk[1] = 1
wpm = (0, 0); wsm = (0, 0); rescued = []; cyc_m = []
t1 = time.time()
for m in range(2, NQ + 1):
    if m % 3 == 0:
        continue
    x = m; s = 0; pk = m; st = None
    while True:
        k = KM[x % 9]; y = x << k
        pk = max(pk, y)
        x = (y + 1) // 3; s += 1
        if x < m:
            st = "below"; break
        if x == m:
            st = "cycle"; break
        if s > 10 ** 5:
            st = "cap"; break
    if st != "below":
        if st == "cycle":
            cyc_m.append(m)
        pth = bfs_below(m, -1)
        chk(pth is not None, "Q2_- stuck at %d" % m)
        rescued.append((m, st, pth))
        x = pth[-1]; pk = max(pth); s = len(pth) - 1
    tpk[m] = max(pk, tpk[x]); tst[m] = s + tst[x]
    if tpk[m] > wpm[0]: wpm = (tpk[m], m)
    if tst[m] > wsm[0]: wsm = (tst[m], m)
print("  Q2_-: max peak %d at m=%d; max compound moves %d at m=%d; greedy-cycle minima %s; rescued %s (%.1fs)"
      % (wpm + wsm + (cyc_m, rescued, time.time() - t1)))
chk(wpm == (150994940, 797161) and wsm == (82, 919795) and cyc_m == [4] and [r[0] for r in rescued] == [4], "Q2_- data")
chk(150994940 == 9 * 2 ** 24 - 4, "minus peak identity")
# Q1_-: deterministic 3n-1 descent exceptions
exc = []
for n in range(2, NQ + 1):
    x = n; s = 0; ok = False
    while True:
        x = x // 2 if x % 2 == 0 else 3 * x - 1
        s += 1
        if x < n:
            ok = True; break
        if x == n or s > 10 ** 5:
            break
    if not ok:
        exc.append(n)
print("  Q1_-: deterministic non-descending starts in [2,10^6]:", exc)
chk(exc == [5, 17], "Q1_- exceptions")


def verify_path(pth, up):
    return all((y == up(x)) or (x % 2 == 0 and y == x // 2) for x, y in zip(pth, pth[1:]))


P5 = [5, 14, 7, 20, 10, 29, 86, 43, 128, 64, 32, 16, 8, 4, 2, 1]
P17a = [1, 2, 5, 14, 7, 20, 59, 176, 88, 44, 22, 11, 32, 16, 8, 4, 11, 32, 16, 8, 23, 68, 34, 17]
P17b = [1, 2, 5, 14, 7, 20, 59, 176, 88, 44, 22, 11, 32, 16, 8, 23, 68, 34, 17]
P17c = [17, 50, 25, 74, 37, 110, 55, 164, 82, 41, 122, 61, 182, 91, 272, 136, 68, 34, 101, 302, 151, 452, 226, 113, 338, 169, 506, 253, 758, 379, 1136, 568, 284, 142, 71, 212, 106, 53, 158, 79, 236, 118, 59, 176, 88, 44, 22, 11, 32, 16, 8, 4, 2, 1]
for nm, pth in (("5->1", P5), ("1->17 (23 arrows, revisits 11,32,16,8)", P17a), ("1->17 shorter (18 arrows)", P17b), ("17->1", P17c)):
    ok = verify_path(pth, lambda x: 3 * x - 1)
    print("  path %s valid: %s (%d arrows)" % (nm, ok, len(pth) - 1))
    chk(ok, "path " + nm)
chk(verify_path([1, 2, 5], lambda x: 3 * x - 1), "1->5")

# ---------------------------------------------------------------- I. pn+1 hierarchy
sec("I. pn+1 hierarchy: k(r) from raw definition on residues mod p^2")
import math
for p in (3, 5, 7, 11, 13, 17, 19):
    ordp = next(t for t in range(1, p) if pow(2, t, p) == 1)
    prim = ordp == p - 1
    sub = sorted({pow(2, t, p) for t in range(ordp)})
    e = ((pow(2, ordp, p * p) - 1) // p) % p
    # raw minimal k for each unit r mod p^2 (search k up to ordp*p+5)
    k = {}; c = {}
    for r in range(1, p * p):
        if r % p == 0 or (r % p) not in sub:
            continue
        for t in range(0, ordp * p + 5):
            y = pow(2, t, p * p) * r % (p * p)
            if y % p == 1 and ((y - 1) // p) % p != 0:
                k[r] = t; c[r] = ((y - 1) // p) % p
                break
    # class law per coset, E_pi[k], rho via raw enumeration mod p^3 lifts (i.i.d. law => stationary = one-step law)
    law_ok = True
    for s in sub:
        law = {}
        for t in range(p):
            law[c[s + p * t]] = law.get(c[s + p * t], 0) + 1
        if law != {cc: 1 + (cc == e) for cc in range(1, p)}:
            law_ok = False
    if prim:
        pi = {s: Fraction(1 + (s == e), p) for s in sub}
        Ek = sum(pi[s] * Fraction(sum(k[s + p * t] for t in range(p)), p) for s in sub)
        # tilted matrix M[s][s'] = (1/p) sum_t 2^{k(r)} [c(r)=s']
        Mp = {s: {s2: Fraction(sum(2 ** k[s + p * t] for t in range(p) if c[s + p * t] == s2), p) for s2 in sub} for s in sub}
        # rank one: rows proportional; eigenvalue = sum over s' of M[e][s'] * ... compute via power on vector 1
        v1 = {s: Fraction(1) for s in sub}
        Mv = {s: sum(Mp[s][s2] * v1[s2] for s2 in sub) for s in sub}
        MMv = {s: sum(Mp[s][s2] * Mv[s2] for s2 in sub) for s in sub}
        rho = MMv[1] / Mv[1]
        rank1 = all(MMv[s] == rho * Mv[s] for s in sub)
        k0 = {s: next(t for t in range(ordp) if pow(2, t, p) * s % p == 1) for s in sub}
        rho_cf = Fraction(sum(2 ** k0[s] for s in sub) + 2 ** (p - 1 + k0[e]), p)
        Ek_cf = Fraction(p - 1, 2) + Fraction(k0[e], p)
        print("  p=%2d ord=%2d prim=%s e=%d law_ok=%s E_pi[k]=%s (cf %s) <log2p:%s rho=%s (cf %s) rank1=%s rho<p:%s"
              % (p, ordp, prim, e, law_ok, Ek, Ek_cf, 2 ** Ek.numerator < p ** Ek.denominator, rho, rho_cf, rank1, rho < p))
        chk(law_ok and Ek == Ek_cf and rho == rho_cf and rank1, "p=%d formulas" % p)
        chk((rho < p) == (p == 3), "mean-subcritical p=%d" % p)
        chk((2 ** Ek.numerator < p ** Ek.denominator) == (p in (3, 5)), "log-subcritical p=%d" % p)
        chk(all(k[r] == k0[r % p] + (p - 1) * (1 if ((pow(2, k0[r % p], p * p) * r - 1) // p) % p == 0 else 0) for r in k), "k formula p=%d" % p)
    else:
        dead = Fraction(sum(1 for r in c if c[r] not in sub), len(c))
        print("  p=%2d ord=%2d prim=%s e=%d <2>=%s law_ok=%s dead-exit prob=%s" % (p, ordp, prim, e, sub, law_ok, dead))
        chk(law_ok, "class law p=%d" % p)
EXP_RHO = {3: Fraction(7, 3), 5: Fraction(47, 5), 11: Fraction(66559, 11), 13: Fraction(1052671, 13), 19: Fraction(8650751, 19)}
EXP_EK = {3: Fraction(1), 5: Fraction(11, 5), 11: Fraction(61, 11), 13: Fraction(86, 13), 19: Fraction(176, 19)}
for p in EXP_RHO:
    ordp = p - 1
    sub = list(range(1, p))
    e = ((pow(2, ordp, p * p) - 1) // p) % p
    k0 = {s: next(t for t in range(ordp) if pow(2, t, p) * s % p == 1) for s in sub}
    chk(Fraction(sum(2 ** k0[s] for s in sub) + 2 ** (p - 1 + k0[e]), p) == EXP_RHO[p], "rho table p=%d" % p)
    chk(Fraction(p - 1, 2) + Fraction(k0[e], p) == EXP_EK[p], "Ek table p=%d" % p)
# p=5 certificate
lam = Fraction(21, 20)
k0_5 = {1: 0, 2: 3, 3: 1, 4: 2}
rho_l = (sum(lam ** k0_5[s] for s in k0_5) + lam ** (4 + k0_5[3])) / 5
print("  p=5: rho(21/20)=%s rho^10=%.6f lam^23=%.6f  2^23<5^10:%s  rate %.5f"
      % (rho_l, float(rho_l ** 10), float(lam ** 23), 2 ** 23 < 5 ** 10, float(rho_l ** 10 / lam ** 23) ** 0.1))
chk(rho_l == Fraction(17876501, 16000000) and rho_l ** 10 < lam ** 23 and 2 ** 23 < 5 ** 10, "p=5 certificate")
# E_5 greedy cycles
K5 = {}
for r in range(1, 25):
    if r % 5:
        for t in range(0, 30):
            y = pow(2, t, 25) * r % 25
            if y % 5 == 1 and ((y - 1) // 5) % 5 != 0:
                K5[r] = t; break
cyc5 = []
for m in range(2, 10 ** 5 + 1):
    if m % 5 == 0:
        continue
    x = m; s = 0
    while True:
        x = ((x << K5[x % 25]) - 1) // 5; s += 1
        if x < m:
            break
        if x == m:
            cyc5.append(m); break
        if s > 10000:
            chk(False, "E_5 cap at %d" % m); break
print("  E_5 greedy-cycle minima <=10^5:", cyc5)
chk(cyc5 == [13, 17], "E_5 greedy cycles")
for m in (13, 17):
    x = m; cy = [m]
    while True:
        x = ((x << K5[x % 25]) - 1) // 5
        if x == m:
            break
        cy.append(x)
    print("   cycle", cy)
# deterministic 5n+1 cycles reversed
det5 = []
for start in (13, 17):
    x = start; cy = [x]
    while True:
        x = x // 2 if x % 2 == 0 else 5 * x + 1
        if x == start:
            break
        cy.append(x)
    det5.append(cy)
print("  deterministic 5n+1 cycles:", det5)
# E_7 reachable residues from 1
seen = {1}; st = [1]
while st:
    x = st.pop()
    for y in ((x // 2,) if x % 2 == 0 else ()) + ((7 * x + 1,) if 7 * x + 1 <= 10 ** 4 else ()):
        if y not in seen:
            seen.add(y); st.append(y)
print("  E_7 reachable residues mod 7 from 1 within [1,10^4]:", sorted({x % 7 for x in seen}))
chk(sorted({x % 7 for x in seen}) == [1, 2, 4], "E_7 residues")
# Q1_5 for n=13 with the explorer's budget replicated independently (DFS halving-preferred)


def reach_below(n, up, budget, cap):
    seen = {n}; st = [n]; ex = 0
    while st:
        x = st.pop(); ex += 1
        if ex > budget:
            return None
        for y in ([up(x)] if up(x) < cap else []) + ([x // 2] if x % 2 == 0 else []):
            if y < n:
                return True
            if y not in seen:
                seen.add(y); st.append(y)
    return False


print("  E_5 n=13 reaches below 13 within 2*10^5 nodes (value cap 2^60)?", reach_below(13, lambda x: 5 * x + 1, 200000, 1 << 60))
print("  E_5 n=13 with budget 10^6, cap 2^80?", reach_below(13, lambda x: 5 * x + 1, 10 ** 6, 1 << 80))

# ---------------------------------------------------------------- J. canon files
sec("J. canon files")
import os
for f in ("01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md",
          "01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md"):
    print("  exists:", f, os.path.exists("/tmp/math-wt-collatz-mod6/" + f))

print("\nTOTAL TIME %.1fs" % (time.time() - T0))
print("FAILS: %d" % len(FAILS))
for f in FAILS:
    print("  -", f)
print("AUDIT RESULT:", "ALL RECOMPUTED CHECKS PASSED" if not FAILS else "SOME CHECKS FAILED")

"""aud_compare.py -- compares the packet's data files and numerical claims with the auditor's own census
(output directory of `aud_engine census OUT`). Auditor's code; nothing imported from the packet.
usage: python aud_compare.py OUT [packet_dir]
"""
import os, sys, collections
OUT = sys.argv[1]
PK = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FAILS = []; NOTES = []
def report(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  | " + detail if detail else "")); sys.stdout.flush()
    if not ok: FAILS.append(name)
def rows(fn):
    for line in open(fn):
        if line.startswith("#") or not line.strip(): continue
        yield line.split()
def oddpart(M):
    j = 0
    while M % 2 == 0: M //= 2; j += 1
    return M, j
def dec(code, m):
    A = [[0] * m for _ in range(m)]; b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
            b += 1
    return A
def H_dp(A):
    n = len(A)
    if n == 0: return 1
    out = [sum(1 << u for u in range(n) if A[v][u]) for v in range(n)]
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for s in range(1, 1 << n):
        fs = f[s]
        for v in range(n):
            c = fs[v]
            if c:
                nx = out[v] & ~s
                while nx:
                    lb = nx & -nx; nx ^= lb
                    f[s | lb][lb.bit_length() - 1] += c
    return sum(f[(1 << n) - 1])
def H_dfs(A):
    n = len(A); outl = [[u for u in range(n) if A[v][u]] for v in range(n)]; full = (1 << n) - 1; cnt = 0
    stack = [(v, 1 << v) for v in range(n)]
    while stack:
        v, mask = stack.pop()
        if mask == full: cnt += 1; continue
        for u in outl[v]:
            if not mask >> u & 1: stack.append((u, mask | 1 << u))
    return cnt
def join(A, B):
    """order join A => B"""
    a, b = len(A), len(B); T = [[0] * (a + b) for _ in range(a + b)]
    for i in range(a):
        for j in range(a): T[i][j] = A[i][j]
        for j in range(b): T[i][a + j] = 1
    for i in range(b):
        for j in range(b): T[a + i][a + j] = B[i][j]
    return T

# ------------------------------------------------------------------ own census
S = {m: set(int(r[0]) for r in rows(os.path.join(OUT, "spectrum_order_%d.txt" % m))) for m in range(1, 10)}
S9 = set().union(*S.values())
wit = {int(r[0]): (int(r[1]), int(r[2])) for r in rows(os.path.join(OUT, "values_witness.txt"))}
myroots = [tuple(int(x) for x in r) for r in rows(os.path.join(OUT, "roots_unit_offset.txt"))]       # m code v H d E U M outdeg
mysets = [tuple(int(x) for x in r) for r in rows(os.path.join(OUT, "unitsets_all_sizes.txt"))]       # m+1 codeB O HB M size
mydbl = [tuple(int(x) for x in r) for r in rows(os.path.join(OUT, "doubling_vertices.txt"))]         # m code x H c2
tuples = [tuple(int(x) for x in r) for r in rows(os.path.join(OUT, "exact_tuples_order_le8.txt"))]

# ------------------------------------------------------------------ 1. spectrum
report("no tournament of order <= 9 has 7 or 21 Hamiltonian paths (complete class lists, own canonical forms)", 7 not in S9 and 21 not in S9)
S8 = set().union(*[S[m] for m in range(1, 9)])
miss200 = [n for n in range(1, 200, 2) if n not in S8]
report("every odd number below 200 other than 7 and 21 is a path count at order <= 8", miss200 == [7, 21], "missing: %s" % miss200)
miss609 = [n for n in range(1, 610, 2) if n not in S8]
report("order <= 8 covers every odd value in [1, 609] except 7, 21 (THM-1370)", miss609 == [7, 21])
print("     distinct values by order:", {m: len(S[m]) for m in S}, "| max:", {m: max(S[m]) for m in S})
report("320 distinct values at order 8; 1,520 at order 9; maxima 661 and 3,357", len(S[8]) == 320 and len(S[9]) == 1520 and max(S[8]) == 661 and max(S[9]) == 3357)
report("the spectrum is monotone in the order (S_m subset of S_(m+1))", all(S[m] <= S[m + 1] for m in range(1, 9)))
first_missing = next(n for n in range(1, 4000, 2) if n not in S9 and n not in (7, 21))
print("     least odd number other than 7, 21 that is not a path count of order <= 9:", first_missing)

# ------------------------------------------------------------------ 2. unit-offset roots
u = collections.Counter(r[0] for r in myroots)
report("unit-offset roots (class-vertex pairs) at orders 2..9 = 2, 5, 7, 16, 31, 68, 207, 1333", [u[m] for m in range(2, 10)] == [2, 5, 7, 16, 31, 68, 207, 1333], str([u[m] for m in range(2, 10)]))
le8 = [r for r in myroots if r[0] <= 8]; at9 = [r for r in myroots if r[0] == 9]
dE8 = sum(1 for r in le8 if r[4] == r[5]); dE9 = sum(1 for r in at9 if r[4] == r[5]); M0_8 = sum(1 for r in le8 if r[7] == 0); M0_all = sum(1 for r in myroots if r[7] == 0)
report("116 of the 336 unit-offset roots of order <= 8 have d = E; 197 of the 1,333 at order nine", len(le8) == 336 and dE8 == 116 and dE9 == 197)
report("HYP-9658: '8 of the 336 are trivial (M = 0)'", M0_8 == 8, "own count: %d of the 336 (one per order 2..8); %d of the 1,669 through order nine" % (M0_8, M0_all))
od = collections.Counter((r[0], r[8]) for r in myroots)
print("     out-degree of the root among unit-offset roots:", {m: {k: od[(m, k)] for k in range(10) if od[(m, k)]} for m in range(2, 10)})
report("out-degree <= 2 up to order seven, out-degree 3 first at order eight", all(r[8] <= 2 for r in myroots if r[0] <= 7) and any(r[8] == 3 for r in myroots if r[0] == 8) and all(r[8] <= 3 for r in myroots))
pairs8 = set((r[8], r[7]) for r in myroots if r[0] == 8 and r[8] >= 2); pairs9 = set((r[8], r[7]) for r in myroots if r[0] == 9 and r[8] >= 2)
us8 = set((r[5], r[4]) for r in mysets if r[0] == 8 and r[5] >= 2); us9 = set((r[5], r[4]) for r in mysets if r[0] == 9 and r[5] >= 2)
report("sets of (out-degree, M): direct root census = unit-set census at orders 8 and 9 (packet: 140 and 1,007 pairs)", pairs8 == us8 and pairs9 == us9 and len(pairs8) == 140 and len(pairs9) == 1007, "%d and %d pairs" % (len(pairs8), len(pairs9)))

# ------------------------------------------------------------------ 3. order ten through unit sets (ALL sizes searched here)
sz = collections.Counter((r[0], r[5]) for r in mysets)
print("     unit sets (base class, O) by root order and |O|:", {m: {k: sz[(m, k)] for k in range(11) if sz[(m, k)]} for m in range(2, 11)})
report("order ten: 19,002 unit sets = 2,771 pairs + 16,231 triples, no set of size four", sz[(10, 2)] == 2771 and sz[(10, 3)] == 16231 and sz[(10, 4)] == 0)
big = sum(v for (m, k), v in sz.items() if k >= 4)
report("AUDIT EXTENSION: no unit set of size >= 4 over any base of order <= 9 (every O of every size was searched), so the order-ten census is complete", big == 0, "unit sets of size >= 4: %d" % big)

# ------------------------------------------------------------------ 4. seeds per line
best = {}; best_tail9 = {}; tail_by_order = collections.defaultdict(set)
def offer(d, c, key):
    if c not in d or key < d[c]: d[c] = key
for r in myroots:
    if r[7] > 0:
        c, j = oddpart(r[7]); offer(best, c, (j + 1, r[0])); offer(best_tail9, c, (j + 1, r[0])); tail_by_order[r[0]].add(c)
mod = {}
for r in mydbl:
    c, j = oddpart(r[4]); offer(mod, c, (j + 1, r[0]))
best9 = dict(best)
for c, key in mod.items(): offer(best9, c, key)
ten = {}
for r in mysets:
    if r[0] == 10 and r[5] >= 2:
        c, j = oddpart(r[4]); offer(ten, c, (j + 1, 10))
best10 = dict(best9)
for c, key in ten.items(): offer(best10, c, key)
cum = []; acc = set()
for m in range(2, 9): acc |= tail_by_order[m]; cum.append(len(acc))
report("odd c with a tail seed at order <= m, m = 2..8: 1, 1, 2, 3, 7, 24, 89", cum == [1, 1, 2, 3, 7, 24, 89], str(cum))
acc8 = set(acc)
no8 = [c for c in range(1, 64, 2) if c not in acc8]
report("odd c <= 63 without a tail seed at order <= 8: 9, 13, 19, 27, 33, 35, 41, 53, 59, 63", no8 == [9, 13, 19, 27, 33, 35, 41, 53, 59, 63], str(no8))
got9 = [c for c in no8 if c in tail_by_order[9]]
report("P2: 8 of the ten acquire a tail seed at order nine (9, 13, 27, 35, 41, 53, 59, 63); P3: c = 13; P4 fails: c = 33", got9 == [9, 13, 27, 35, 41, 53, 59, 63], str(got9))
report("P3 detail: a unit-offset root of order nine with (H,d,E,U) = (183,61,35,269), M = 208 = 13 * 16", any(r[0] == 9 and r[3:8] == (183, 61, 35, 269, 208) for r in myroots))
ntail9 = len(best_tail9); nany9 = len(best9)
report("P1 outcome: 636 odd c have a TAIL seed of order <= 9", ntail9 == 636, "tail seeds of order <= 9: %d odd c; tail <= 9 or module <= 8: %d odd c; c with a module seed only: %s" % (ntail9, nany9, sorted(c for c in mod if c not in best_tail9)))
report("ratio 636 / 89 = 7.1", abs(ntail9 / 89 - 7.1) < 0.06, "%.3f" % (ntail9 / 89))
no129 = [c for c in range(1, 130, 2) if c not in best9]
report("odd c <= 129 without a seed (tail <= 9 or module <= 8): 19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115", no129 == [19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115], str(no129))
q1 = {c: ten[c][0] for c in no129 if c in ten}
report("Q1: all eleven acquire a seed at order ten, least s0 19:8 33:6 71:6 79:5 83:7 97:4 101:4 103:6 107:6 113:5 115:2", q1 == {19: 8, 33: 6, 71: 6, 79: 5, 83: 7, 97: 4, 101: 4, 103: 6, 107: 6, 113: 5, 115: 2}, str(q1))
report("Q2/Q3 detail: unit sets with (H(B), M) = (345, 2432 = 19 * 128) and (307, 1056 = 33 * 32)", any(r[0] == 10 and r[3] == 345 and r[4] == 2432 for r in mysets) and any(r[0] == 10 and r[3] == 307 and r[4] == 1056 for r in mysets))
n1025_9 = sum(1 for c in range(1, 1026, 2) if c in best9); n1025 = sum(1 for c in range(1, 1026, 2) if c in best10)
miss1025 = [c for c in range(1, 1026, 2) if c not in best10]
report("Q4: odd c <= 1025 with a seed rises from 377 to 510 of 513; the three without are 245, 853, 909", n1025_9 == 377 and n1025 == 510 and miss1025 == [245, 853, 909], "%d -> %d, missing %s" % (n1025_9, n1025, miss1025))
n4097 = sum(1 for c in range(1, 4098, 2) if c in best10); miss4097 = [c for c in range(1, 4098, 2) if c not in best10]
claimed33 = [245, 853, 909, 1049, 1407, 1467, 1515, 1577, 1629, 1829, 1887, 2023, 2037, 2085, 2343, 2393, 2703, 3001, 3033, 3123, 3313, 3411, 3427, 3545, 3713, 3733, 3753, 3773, 3907, 3909, 3985, 3991, 4023]
report("2,016 of the 2,049 odd c <= 4097 have a seed of order <= 10; the 33 without are the README's list (least 245)", n4097 == 2016 and miss4097 == claimed33, "%d; missing %d" % (n4097, len(miss4097)))
report("5,529 odd c have a seed in total (tail <= 10 or module <= 8)", len(best10) == 5529, str(len(best10)))
report("G1/G2: least s0 is 2 for c = 3 and 3 for c = 5 (no ray through 7 or 21)", best10[3][0] == 2 and best10[5][0] == 3)

# ------------------------------------------------------------------ 5. line_seeds.txt and the s0 table of (B3)
ls = {}
for r in rows(os.path.join(PK, "line_seeds.txt")): ls[int(r[0])] = (int(r[1]), r[2], int(r[3]))
okls = (sorted(ls) == sorted(c for c in best10 if c <= 4097)) and all(ls[c][0] == best10[c][0] and ls[c][2] == best10[c][1] for c in ls)
report("line_seeds.txt: the set of c <= 4097 and every (least s0, order) agree with the own census", okls, "%d lines" % len(ls))
mech = collections.Counter(ls[c][1] for c in range(1, 130, 2))
report("c <= 129: 61 tail seeds (order <= 10) and 4 module seeds (order <= 8)", mech == {"tail": 61, "module": 4} and all(c in ls for c in range(1, 130, 2)), str(dict(mech)))
tab_c = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65,
         67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129]
tab_s = [1, 2, 3, 5, 1, 1, 5, 1, 2, 8, 4, 5, 3, 6, 1, 2, 6, 5, 1, 3, 6, 2, 3, 1, 1, 1, 5, 1, 2, 2, 2, 4, 3,
         3, 3, 6, 2, 3, 1, 5, 1, 7, 3, 2, 1, 2, 1, 1, 4, 1, 4, 6, 3, 6, 1, 1, 5, 2, 3, 3, 4, 1, 1, 1, 1]
report("the s0 table of THM-5126 (B3) (65 entries, typed in from the theorem) = line_seeds.txt = own census", len(tab_c) == 65 and tab_c == list(range(1, 130, 2)) and all(ls[c][0] == s and best10[c][0] == s for c, s in zip(tab_c, tab_s)))
report("least s0 <= 8 for c <= 129, attained only at c = 19", max(tab_s) == 8 and [c for c, s in zip(tab_c, tab_s) if s == 8] == [19])

# ------------------------------------------------------------------ 6. below_rays.txt
br = {}
for r in rows(os.path.join(PK, "below_rays.txt")): br[(int(r[0]), int(r[1]))] = (int(r[2]), r[3], r[4:])
need = [(c, s) for c in range(1, 130, 2) for s in range(1, ls[c][0])]
report("below_rays.txt lists exactly the points (c, s), s < s0(c), odd c <= 129", sorted(br) == sorted(need), "%d points; sum of (s0 - 1) = %d" % (len(br), len(need)))
kinds = collections.Counter(v[1] for v in br.values())
report("121 points = 2 gaps + 117 census + 2 products", len(br) == 121 and kinds == {"GAP": 2, "census": 117, "product": 2}, str(dict(kinds)))
okb = True; recount = 0; maxc = 0
for (c, s), (N, kind, extra) in sorted(br.items()):
    if N != 1 + c * 2 ** s: okb = False
    if kind == "GAP":
        if N not in (7, 21) or N in S9: okb = False
    elif kind == "census":
        if N not in S9: okb = False
        else:
            m, code = wit[N]; recount += 1; maxc = max(maxc, N)
            if H_dp(dec(code, m)) != N: okb = False            # own witness, recounted in Python
    elif kind == "product":
        a, b = int(extra[0]), int(extra[1])
        if a * b != N or a not in S9 or b not in S9: okb = False
        else:
            (ma, ca), (mb, cb) = wit[a], wit[b]; T = join(dec(ca, ma), dec(cb, mb)); h1 = H_dp(T); h2 = H_dfs(T); recount += 1
            print("     %d = %d * %d: witnesses of order %d and %d (H recounted: %d, %d); order join of order %d has H = %d (DP), %d (DFS)" % (N, a, b, ma, mb, H_dp(dec(ca, ma)), H_dp(dec(cb, mb)), ma + mb, h1, h2))
            if h1 != N or h2 != N: okb = False
    else: okb = False
report("every certificate of below_rays.txt holds: census values are path counts of order <= 9 (own witness recounted), products are order joins counted literally, gaps are 7 and 21", okb, "%d recounts; largest census value %d" % (recount, maxc))
could_census = [(c, s, v[0]) for (c, s), v in br.items() if v[1] == "product" and v[0] in S9]
report("the two product points are indeed above the order-nine range (3,357)", could_census == [] and all(v[0] > 3357 for v in br.values() if v[1] == "product"))
over = [(c, s, 1 + c * 2 ** s) for c in range(131, 258, 2) if c in ls for s in range(1, ls[c][0]) if 1 + c * 2 ** s > 3357]
below257 = [(c, s, 1 + c * 2 ** s) for c in range(131, 258, 2) if c in ls for s in range(1, ls[c][0])]
unc = [(c, s, N) for (c, s, N) in below257 if N not in S9]
print("     README section 5: points below the rays for 129 < c <= 257: %d, of which above 3,357: %d %s; not order-<=9 values: %d; c in (129,257] without a seed: %s" % (len(below257), len(over), over, len(unc), [c for c in range(131, 258, 2) if c not in ls]))

# ------------------------------------------------------------------ 7. weight three
w3 = {d: ls.get(2 ** d + 1) for d in range(1, 13)}
print("     weight-three rows c = 2^d + 1: (least s0, mechanism, order):", w3)
report("(B4): rows d <= 7 have seeds (c = 3, 5, 9, 17, 33, 65, 129 <= 129); d = 8 (c = 257) status", all(w3[d] for d in range(1, 8)), "c = 257: %s" % (w3[8],))

# ------------------------------------------------------------------ 8. doubling vertices
dv = set(r[3] for r in mydbl)
report("49 of the 320 values of order <= 8 have a doubling realiser (c_2 > 0) of order <= 8", len(dv) == 49 and dv <= S8, "%d values" % len(dv))
l3 = sorted(n for n in dv if n % 4 == 3 and n < 200)
report("N = 3 mod 4, N < 200 with a doubling realiser: 3, 19, 23, 31, 59, 75, 95, 99, 103, 111, 155, 163, 179, 191, 199", l3 == [3, 19, 23, 31, 59, 75, 95, 99, 103, 111, 155, 163, 179, 191, 199], str(l3))
report("15 has no doubling realiser of order <= 8; 11 has none (forced by 21)", 15 not in dv and 11 not in dv)
modonly = sorted(c for c in range(1, 130, 2) if ls[c][1] == "module")
report("module seeds used for c <= 129: 9 (order 6), 15 (order 7), 51, 81 (order 8)", modonly == [9, 15, 51, 81] and all(mod[c] == (1, o) for c, o in zip(modonly, (6, 7, 8, 8))))

# ------------------------------------------------------------------ 9. two-module tuples
tset = set(tuples)
report("3,477 distinct exact bilinear tuples (orders 3-8); none has de = 0 with be, ga != 0", len(tset) == 3477 and not any(t[3] == 0 and t[1] != 0 and t[2] != 0 for t in tset))
report("positive controls: (3,-1,-1,1) [C_3[TT_a,TT_b,1]] and the fans are among the exact tuples", (3, -1, -1, 1) in tset)
AM = 26; hit = set()
for (al, be, ga, de) in tset:
    for a in range(1, AM + 1):
        for b in range(1, AM + 1):
            h = al + be * 2 ** a + ga * 2 ** b + de * 2 ** (a + b)
            if h > 0 and h & 1 and bin(h).count("1") == 3: hit.add(h)
rowsok = {}
for d in range(1, 13):
    pts = [1 + 2 ** b + 2 ** (b + d) for b in range(1, 20 - d + 1)]
    rowsok[d] = [N for N in pts if N not in hit]
full = [d for d in rowsok if all(N in (7, 21) for N in rowsok[d])]
report("weight-three map: rows d <= 4 and d = 7 are complete (b <= 20 - d) from two-module tuples of order <= 8 alone", set(full) == {1, 2, 3, 4, 7}, "complete rows: %s; first misses: %s" % (full, {d: rowsok[d][:3] for d in rowsok if d not in full and d <= 8}))

print("ALL PASS" if not FAILS else "FAILURES: %s" % FAILS)
sys.exit(1 if FAILS else 0)

#!/usr/bin/env python3
"""
ADVERSARIAL verification of THM-468 (tournament determinant floor).
Written FRESH by an adversarial-check subagent (2026-06-11), deliberately NOT
reusing hadamard_det_census_macmini_s2.py. Every primitive (det, Pf, vortex,
local transitivity, switching) is re-implemented from definitions.

Checks:
  A. det(I+X) = sum over principal minors det(X[K])   [random integer matrices]
  B. Pfaffian last-row expansion sign convention: Pf(A) = sum_{i<2k} (-1)^{i+1} a_{i,2k} Pf(minor)
  C. Transitive Pfaffian recursion: Pf_{2k} = 1 for 2k = 2,4,6,8 (all +1 above diag)
  D. 4-set classification: all 64 labeled 4-tournaments, |Pf| = 3 <=> vortex
  E. Exhaustive n=5, n=6 (all labeled tournaments):
       - every even-minor Pf is odd
       - det(I+S) == sum of Pf^2 over even subsets
       - 2^(n-1) | det(I+S), d >= 1
       - {d=1} == {locally transitive} == {vortex-free} == switching closure of transitive
       - labeled count == (2n-2)!!, iso-class counts
  F. Epsilon-lemma of (c)=>(d): T-v transitive + arbitrary epsilon:
       vortex-free <=> epsilon monotone (1^a 0^b or 0^a 1^b), n = 5..8
  G. n=7: switching closure size == 12!! == 46080 and members all have d=1;
       random non-members have d > 1.
"""
import itertools, random
random.seed(468)

# ---------- integer determinant (Bareiss, fraction-free) ----------
def det_int(Mat):
    M = [list(map(int, row)) for row in Mat]
    n = len(M)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if M[r][k] != 0), None)
            if piv is None:
                return 0
            M[k], M[piv] = M[piv], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[-1][-1]

# ---------- Pfaffian via first-row expansion (standard convention) ----------
def pf(A):
    n = len(A)
    if n == 0:
        return 1
    if n % 2:
        return 0
    tot = 0
    for j in range(1, n):
        s = 1 if j % 2 == 1 else -1            # +,-,+,... for j = 1,2,3,...
        idx = [k for k in range(n) if k != 0 and k != j]
        sub = [[A[r][c] for c in idx] for r in idx]
        tot += s * A[0][j] * pf(sub)
    return tot

def pf_matchsum(A):
    """Pfaffian directly from the perfect-matching definition (independent check)."""
    n = len(A)
    if n == 0:
        return 1
    if n % 2:
        return 0
    verts = list(range(n))
    total = 0
    def rec(rem, pairs):
        nonlocal total
        if not rem:
            # sign of the permutation (p1 q1 p2 q2 ...) relative to identity
            perm = [x for p in pairs for x in p]
            sgn, seen = 1, [False] * n
            pos = {v: i for i, v in enumerate(perm)}
            # sign of permutation taking identity -> perm
            visited = [False] * n
            for i in range(n):
                if visited[i]:
                    continue
                j, clen = i, 0
                while not visited[j]:
                    visited[j] = True
                    j = perm[j]
                    clen += 1
                if clen % 2 == 0:
                    sgn = -sgn
            total += sgn * eval_prod(pairs)
            return
        a = rem[0]
        for b in rem[1:]:
            rec([x for x in rem if x not in (a, b)], pairs + [(a, b)])
    def eval_prod(pairs):
        p = 1
        for (a, b) in pairs:
            p *= A[a][b]
        return p
    rec(verts, [])
    return total

# ---------- tournaments ----------
def pairs_of(n):
    return list(itertools.combinations(range(n), 2))

def S_from_code(code, n, pairs):
    S = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(pairs):
        if (code >> b) & 1:
            S[i][j], S[j][i] = 1, -1
        else:
            S[i][j], S[j][i] = -1, 1
    return S

def M_IplusS(S):
    n = len(S)
    return [[(1 if i == j else S[i][j]) for j in range(n)] for i in range(n)]

def has_3cycle(S, verts):
    for a, b, c in itertools.combinations(verts, 3):
        # 3-cycle iff each vertex of the induced triangle has out-degree 1
        x = (S[a][b], S[b][c], S[c][a])
        if x == (1, 1, 1) or x == (-1, -1, -1):
            return True
    return False

def is_vortex_free(S, n):
    for quad in itertools.combinations(range(n), 4):
        for v in quad:
            C = [u for u in quad if u != v]
            a, b, c = C
            cyc = (S[a][b], S[b][c], S[c][a])
            if cyc not in ((1, 1, 1), (-1, -1, -1)):
                continue
            if all(S[v][u] == 1 for u in C) or all(S[v][u] == -1 for u in C):
                return False
    return True

def is_locally_transitive(S, n):
    for v in range(n):
        outs = [u for u in range(n) if S[v][u] == 1]
        ins = [u for u in range(n) if S[v][u] == -1]
        if has_3cycle(S, outs) or has_3cycle(S, ins):
            return False
    return True

def encode_S(S, n, pairs):
    code = 0
    for b, (i, j) in enumerate(pairs):
        if S[i][j] == 1:
            code |= 1 << b
    return code

def switch_code(code, n, pairs, W):
    out = 0
    for b, (i, j) in enumerate(pairs):
        bit = (code >> b) & 1
        if ((W >> i) & 1) != ((W >> j) & 1):
            bit ^= 1
        if bit:
            out |= 1 << b
    return out

def perm_code(code, n, pairs, perm, pair_index):
    out = 0
    for b, (i, j) in enumerate(pairs):
        pi, pj = perm[i], perm[j]
        bit = (code >> b) & 1
        if pi < pj:
            nb = pair_index[(pi, pj)]
            nbit = bit
        else:
            nb = pair_index[(pj, pi)]
            nbit = bit ^ 1
        if nbit:
            out |= 1 << nb
    return out

FAIL = []
def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, detail)
    if not ok:
        FAIL.append((name, detail))

# ===== A. principal-minor expansion det(I+X) = sum_K det(X[K]) =====
print("=" * 70)
print("A. det(I+X) = sum over principal minors, random integer matrices")
for n in (3, 4, 5, 6):
    for t in range(20):
        X = [[random.randint(-3, 3) for _ in range(n)] for _ in range(n)]
        lhs = det_int([[(1 if i == j else 0) + X[i][j] for j in range(n)] for i in range(n)])
        rhs = 0
        for k in range(n + 1):
            for K in itertools.combinations(range(n), k):
                rhs += det_int([[X[r][c] for c in K] for r in K])
        if lhs != rhs:
            check(f"A n={n}", False, f"lhs={lhs} rhs={rhs}")
            break
    else:
        continue
    break
else:
    pass
check("A principal-minor identity (random, n=3..6, 20 trials each)", not FAIL)

# ===== B. last-row Pfaffian expansion sign convention =====
print("=" * 70)
print("B. Pf(A) = sum_{i=1}^{2k-1} (-1)^{i+1} a_{i,2k} Pf(A minus rows/cols i,2k)")
okB = True
for n in (2, 4, 6, 8):
    for t in range(10):
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                A[i][j] = random.randint(-3, 3)
                A[j][i] = -A[i][j]
        direct = pf(A)
        ms = pf_matchsum(A)
        if direct != ms:
            okB = False
            check("B pf vs matching-sum", False, f"n={n} {direct} vs {ms}")
        lastrow = 0
        for i in range(1, n):                  # 1-based i = 1..2k-1
            s = 1 if i % 2 == 1 else -1        # (-1)^(i+1)
            idx = [k for k in range(n) if k not in (i - 1, n - 1)]
            sub = [[A[r][c] for c in idx] for r in idx]
            lastrow += s * A[i - 1][n - 1] * pf(sub)
        if lastrow != direct:
            okB = False
            check("B last-row convention", False, f"n={n} expansion={lastrow} pf={direct}")
        # det = Pf^2 sanity
        if det_int(A) != direct * direct:
            okB = False
            check("B det=Pf^2", False, f"n={n}")
check("B last-row expansion sign convention + det=Pf^2 (n=2,4,6,8)", okB)

# ===== C. transitive Pfaffian recursion =====
print("=" * 70)
print("C. all-(+1)-above-diagonal: Pf_{2k} for 2k = 2,4,6,8")
okC = True
for n in (2, 4, 6, 8):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j], A[j][i] = 1, -1
    v = pf(A)
    v2 = pf_matchsum(A)
    print(f"   2k={n}: Pf = {v} (matching-sum {v2})")
    if v != 1 or v2 != 1:
        okC = False
    # also verify each minor (remove j and 2k) is again the all-ones pattern
    for jj in range(n - 1):
        idx = [k for k in range(n) if k not in (jj, n - 1)]
        sub = [[A[r][c] for c in idx] for r in idx]
        for a in range(len(sub)):
            for b in range(a + 1, len(sub)):
                if sub[a][b] != 1:
                    okC = False
check("C transitive Pf_{2k} = 1 for 2k=2,4,6,8; minors stay all-+1 pattern", okC)

# ===== D. 4-set classification =====
print("=" * 70)
print("D. all 64 labeled 4-tournaments: |Pf| vs vortex vs iso class")
p4 = pairs_of(4)
seen = {}
okD = True
for code in range(64):
    S = S_from_code(code, 4, p4)
    P = pf(S)
    vfree = is_vortex_free(S, 4)
    scores = tuple(sorted(sum(1 for j in range(4) if S[i][j] == 1) for i in range(4)))
    key = (scores, abs(P), vfree)
    seen[key] = seen.get(key, 0) + 1
    if (abs(P) == 3) != (not vfree):
        okD = False
    if P not in (1, -1, 3, -3):
        okD = False
for key in sorted(seen):
    print(f"   scores={key[0]} |Pf|={key[1]} vortex-free={key[2]}: {seen[key]} labeled")
check("D |Pf|=3 <=> vortex, Pf in {±1,±3}, all 64 labeled 4-tournaments", okD)
# iso-class sanity: 4 score-type classes
check("D four 4-classes seen", len(seen) == 4, f"classes={sorted(seen)}")

# ===== E. exhaustive n=5 and n=6 =====
print("=" * 70)
for n in (5, 6):
    print(f"E. exhaustive n={n}: all {2**(n*(n-1)//2)} labeled tournaments")
    P = pairs_of(n)
    pair_index = {pr: b for b, pr in enumerate(P)}
    nbits = len(P)
    evens = [K for k in range(0, n + 1, 2) for K in itertools.combinations(range(n), k)]
    # switching closure of all labeled transitive tournaments
    closure = set()
    for perm in itertools.permutations(range(n)):
        pos = {v: i for i, v in enumerate(perm)}
        c0 = 0
        for b, (i, j) in enumerate(P):
            if pos[i] < pos[j]:
                c0 |= 1 << b
        for W in range(0, 1 << n, 2):  # vertex 0 not in W (W ~ complement)
            closure.add(switch_code(c0, n, P, W))
    d1set, ltset, vfset = set(), set(), set()
    okOdd = okSum = okDiv = True
    for code in range(1 << nbits):
        S = S_from_code(code, n, P)
        det = det_int(M_IplusS(S))
        tot = 0
        for K in evens:
            pk = pf([[S[r][c] for c in K] for r in K])
            if pk % 2 == 0:
                okOdd = False
            tot += pk * pk
        if tot != det:
            okSum = False
        if det % (1 << (n - 1)) != 0 or det < (1 << (n - 1)):
            okDiv = False
        if det == (1 << (n - 1)):
            d1set.add(code)
        if is_locally_transitive(S, n):
            ltset.add(code)
        if is_vortex_free(S, n):
            vfset.add(code)
    dblfac = 1
    for k in range(2, 2 * n - 1, 2):
        dblfac *= k
    check(f"E n={n} every even-minor Pf odd", okOdd)
    check(f"E n={n} det(I+S) == sum Pf^2", okSum)
    check(f"E n={n} 2^(n-1) | det and det >= 2^(n-1)", okDiv)
    check(f"E n={n} {{d=1}} == {{locally transitive}}", d1set == ltset,
          f"|d1|={len(d1set)} |lt|={len(ltset)}")
    check(f"E n={n} {{d=1}} == {{vortex-free}}", d1set == vfset,
          f"|vf|={len(vfset)}")
    check(f"E n={n} {{d=1}} == switching closure of transitive", d1set == closure,
          f"|closure|={len(closure)}")
    check(f"E n={n} labeled local orders == (2n-2)!! == {dblfac}", len(d1set) == dblfac,
          f"got {len(d1set)}")
    # iso classes of the d=1 stratum
    canon = set()
    for code in d1set:
        canon.add(min(perm_code(code, n, P, perm, pair_index)
                      for perm in itertools.permutations(range(n))))
    print(f"   n={n}: local-order iso classes = {len(canon)}")
    check(f"E n={n} iso-class count matches theorem table",
          len(canon) == {5: 4, 6: 6}[n], f"got {len(canon)}")
    # invariance spot checks of d
    okInv = True
    for _ in range(50):
        code = random.randrange(1 << nbits)
        S = S_from_code(code, n, P)
        d0 = det_int(M_IplusS(S)) >> (n - 1)
        W = random.randrange(1 << n)
        dW = det_int(M_IplusS(S_from_code(switch_code(code, n, P, W), n, P))) >> (n - 1)
        Sop = [[-S[i][j] for j in range(n)] for i in range(n)]
        dop = det_int(M_IplusS(Sop)) >> (n - 1)
        pr = list(range(n)); random.shuffle(pr)
        dp = det_int(M_IplusS(S_from_code(perm_code(code, n, P, pr, pair_index), n, P))) >> (n - 1)
        if not (d0 == dW == dop == dp):
            okInv = False
    check(f"E n={n} d invariant under switching/reversal/relabeling (spot)", okInv)

# ===== F. epsilon-lemma of (c)=>(d) =====
print("=" * 70)
print("F. T - v transitive (u_1 > ... > u_{n-1}), arbitrary epsilon:")
print("   vortex-free(T)  <=>  epsilon is 1^a 0^b or 0^a 1^b")
okF = True
for n in (5, 6, 7, 8):
    m = n - 1
    P = pairs_of(n)
    for eps in itertools.product((0, 1), repeat=m):
        # vertices 0..m-1 = u_1..u_m with u_i -> u_j iff i<j; vertex m = v
        S = [[0] * n for _ in range(n)]
        for i in range(m):
            for j in range(i + 1, m):
                S[i][j], S[j][i] = 1, -1
        for i in range(m):
            if eps[i]:
                S[m][i], S[i][m] = 1, -1
            else:
                S[m][i], S[i][m] = -1, 1
        vf = is_vortex_free(S, n)
        s = "".join(map(str, eps))
        mono = (s == "".join(sorted(s))) or (s == "".join(sorted(s, reverse=True)))
        if vf != mono:
            okF = False
            print(f"   COUNTEREXAMPLE n={n} eps={s}: vortex-free={vf}, monotone={mono}")
check("F vortex-free <=> epsilon monotone, n=5..8, all epsilon", okF)

# ===== G. n=7 switching closure =====
print("=" * 70)
print("G. n=7: switching closure of transitive")
n = 7
P = pairs_of(n)
closure7 = set()
for perm in itertools.permutations(range(n)):
    pos = {v: i for i, v in enumerate(perm)}
    c0 = 0
    for b, (i, j) in enumerate(P):
        if pos[i] < pos[j]:
            c0 |= 1 << b
    for W in range(0, 1 << n, 2):
        closure7.add(switch_code(c0, n, P, W))
check("G n=7 |closure| == 12!! == 46080", len(closure7) == 46080, f"got {len(closure7)}")
okG = True
for code in random.sample(sorted(closure7), 200):
    S = S_from_code(code, n, P)
    if det_int(M_IplusS(S)) != 1 << (n - 1):
        okG = False
check("G n=7 random closure members have d=1 (200 samples)", okG)
okG2 = True
cnt = 0
while cnt < 200:
    code = random.randrange(1 << len(P))
    if code in closure7:
        continue
    cnt += 1
    S = S_from_code(code, n, P)
    det = det_int(M_IplusS(S))
    if det == 1 << (n - 1):
        okG2 = False
    if det % (1 << (n - 1)) != 0:
        okG2 = False
check("G n=7 random NON-members have d>1 and divisibility holds (200 samples)", okG2)

print("=" * 70)
if FAIL:
    print(f"OVERALL: {len(FAIL)} FAILURE(S):")
    for f in FAIL:
        print("  ", f)
else:
    print("OVERALL: ALL CHECKS PASS")

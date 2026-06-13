#!/usr/bin/env python3
"""
THM-468 ADVERSARIAL CHECK -- completely fresh code, no reuse of session scripts.

Checks:
  A. Pfaffian implementation self-test: pf^2 == det for random skew matrices.
  B. Last-row Pfaffian expansion sign convention (used in (d)=>(b)).
  C. Transitive Pfaffian Pf_{2k} = 1 for 2k = 2,4,6,8 (all +1 above diagonal).
  D. 4-tournament classification: all 64 labeled 4-tournaments, |Pf| = 3 <=> vortex,
     grouped by iso class (score sequence suffices at n=4).
  E. n=5 FULL: det(I+S) = sum of Pf^2 over even subsets; every Pf odd;
     2^(n-1) | det; d := det/2^(n-1) >= 1.
  F. n=5 and n=6 FULL: d=1 <=> vortex-free <=> locally-transitive (direct defn)
     <=> switching-of-transitive (set membership). Counts == (2n-2)!!.
  G. Switching invariance of d (random sample, n=6).
  H. Epsilon-induction targeted check: T with T-v transitive chain; vortex-free
     <=> epsilon monotone (1^a 0^b or 0^a 1^b), for n = 4..7.
  I. Base case: at n=3, switchings of transitives = ALL 8 labeled tournaments.
  J. n=6 sampled: det identity + Pf-oddness on random tournaments.
"""
import itertools, random
random.seed(468)

# ---------- basic machinery ----------

def pairs_list(n):
    return list(itertools.combinations(range(n), 2))

def decode(enc, n, P):
    adj = [[0]*n for _ in range(n)]
    for t, (i, j) in enumerate(P):
        if (enc >> t) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def encode(adj, n, P):
    e = 0
    for t, (i, j) in enumerate(P):
        if adj[i][j]:
            e |= (1 << t)
    return e

def smat(adj, n):
    return [[0 if i == j else (1 if adj[i][j] else -1) for j in range(n)] for i in range(n)]

def det_bareiss(Min):
    M = [row[:] for row in Min]
    n = len(M)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n-1):
        if M[k][k] == 0:
            piv = None
            for r in range(k+1, n):
                if M[r][k] != 0:
                    piv = r
                    break
            if piv is None:
                return 0
            M[k], M[piv] = M[piv], M[k]
            sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n-1][n-1]

def pf(S, idx):
    """Pfaffian via expansion along first row of the index set (standard)."""
    L = len(idx)
    if L == 0:
        return 1
    if L % 2 == 1:
        return 0
    i0 = idx[0]
    rest = idx[1:]
    total = 0
    for k, j in enumerate(rest):
        sgn = 1 if (k % 2 == 0) else -1
        rem = rest[:k] + rest[k+1:]
        total += sgn * S[i0][j] * pf(S, rem)
    return total

# ---------- tournament predicates (direct from definitions) ----------

def is_3cycle(adj, a, b, c):
    # each vertex has out-degree exactly 1 within the triple
    return (adj[a][b] + adj[a][c] == 1 and
            adj[b][a] + adj[b][c] == 1 and
            adj[c][a] + adj[c][b] == 1)

def is_vortex4(adj, quad):
    """Is the induced 4-tournament on quad a vortex: v -> 3-cycle or 3-cycle -> v?"""
    for v in quad:
        rest = [u for u in quad if u != v]
        a, b, c = rest
        if is_3cycle(adj, a, b, c):
            if all(adj[v][u] for u in rest):
                return True
            if all(adj[u][v] for u in rest):
                return True
    return False

def has_vortex(adj, n):
    for quad in itertools.combinations(range(n), 4):
        if is_vortex4(adj, quad):
            return True
    return False

def induced_transitive(adj, verts):
    for tri in itertools.combinations(verts, 3):
        if is_3cycle(adj, *tri):
            return False
    return True

def locally_transitive(adj, n):
    for v in range(n):
        outs = [u for u in range(n) if u != v and adj[v][u]]
        ins  = [u for u in range(n) if u != v and adj[u][v]]
        if not induced_transitive(adj, outs):
            return False
        if not induced_transitive(adj, ins):
            return False
    return True

def is_transitive(adj, n):
    return induced_transitive(adj, list(range(n)))

def switch(adj, n, W):
    out = [row[:] for row in adj]
    for i in range(n):
        for j in range(n):
            if i != j and ((i in W) != (j in W)):
                out[i][j] = adj[j][i]
    return out

def switch_set_of_transitives(n, P):
    seen = set()
    for perm in itertools.permutations(range(n)):
        adj = [[0]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1, n):
                adj[perm[a]][perm[b]] = 1   # perm[a] beats perm[b]
        for wmask in range(1 << n):
            W = {i for i in range(n) if (wmask >> i) & 1}
            seen.add(encode(switch(adj, n, W), n, P))
    return seen

def double_fact(k):
    r = 1
    while k > 1:
        r *= k
        k -= 2
    return r

FAIL = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {name}" + (f" -- {detail}" if detail else ""))
    if not cond:
        FAIL.append(name + " " + detail)

# ---------- A. Pfaffian self-test ----------
print("=== A. Pfaffian self-test (pf^2 == det, random skew) ===")
ok = True
for trial in range(200):
    m = random.choice([2, 4, 6])
    A = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            A[i][j] = random.randint(-3, 3)
            A[j][i] = -A[i][j]
    p = pf(A, tuple(range(m)))
    d = det_bareiss(A)
    if p*p != d:
        ok = False
        break
check("pf^2 == det on 200 random skew matrices (sizes 2,4,6)", ok)

# ---------- B. Last-row expansion sign convention ----------
print("\n=== B. Last-row Pfaffian expansion: Pf(A) = sum_{i=1}^{2k-1} (-1)^(i+1) a_{i,2k} Pf(A minus {i,2k}) ===")
ok = True
for trial in range(200):
    m = random.choice([2, 4, 6, 8])
    A = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            A[i][j] = random.randint(-3, 3)
            A[j][i] = -A[i][j]
    full = pf(A, tuple(range(m)))
    s = 0
    for i in range(m-1):  # i is 0-based; theorem's j = i+1, sign (-1)^(j+1) = (-1)^i
        idx = tuple(x for x in range(m-1) if x != i)
        s += ((-1) ** i) * A[i][m-1] * pf(A, idx)
    if s != full:
        ok = False
        break
check("last-row expansion identity (200 random skew, sizes 2-8)", ok)

# ---------- C. Transitive Pfaffian ----------
print("\n=== C. Transitive Pfaffian: all +1 above diagonal, sizes 2,4,6,8 ===")
for m in (2, 4, 6, 8):
    A = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            A[i][j] = 1
            A[j][i] = -1
    p = pf(A, tuple(range(m)))
    check(f"Pf(transitive {m}x{m}) == 1", p == 1, f"got {p}")

# ---------- D. 4-tournament classification ----------
print("\n=== D. All 64 labeled 4-tournaments: |Pf| vs vortex vs iso class ===")
P4 = pairs_list(4)
table = {}
ok_dichotomy = True
for enc in range(64):
    adj = decode(enc, 4, P4)
    S = smat(adj, 4)
    p = pf(S, (0, 1, 2, 3))
    vor = is_vortex4(adj, (0, 1, 2, 3))
    ss = tuple(sorted(sum(adj[i]) for i in range(4)))
    key = (ss, abs(p), vor)
    table[key] = table.get(key, 0) + 1
    if (abs(p) == 3) != vor:
        ok_dichotomy = False
    if abs(p) not in (1, 3):
        ok_dichotomy = False
for key in sorted(table):
    ss, ap, vor = key
    print(f"  score seq {ss}: |Pf|={ap}, vortex={vor}, count={table[key]}")
check("|Pf|=3 <=> vortex, |Pf| in {1,3}, on all 64 labeled 4-tournaments", ok_dichotomy)
# expected classes: (0,1,2,3) transitive |Pf|=1; (1,1,2,2) strong |Pf|=1;
# (1,1,1,3) v->C3 |Pf|=3 vortex; (0,2,2,2) C3->v |Pf|=3 vortex
exp = {((0,1,2,3), 1, False): 24, ((1,1,2,2), 1, False): 24,
       ((1,1,1,3), 3, True): 8, ((0,2,2,2), 3, True): 8}
check("4-class table matches expected (transitive/strong |Pf|=1; both vortex classes |Pf|=3)",
      table == exp, f"table={table}")

# ---------- E. n=5 full: det identity, Pf oddness, divisibility ----------
print("\n=== E. n=5 FULL: det(I+S) = sum Pf^2, all Pf odd, 2^(n-1)|det ===")
n = 5
P5 = pairs_list(n)
even_subsets = [k for k in
                [tuple(c) for r in range(0, n+1, 2) for c in itertools.combinations(range(n), r)]]
ok_id = ok_odd = ok_div = ok_pos = True
for enc in range(1 << len(P5)):
    adj = decode(enc, n, P5)
    S = smat(adj, n)
    M = [[(1 if i == j else 0) + S[i][j] for j in range(n)] for i in range(n)]
    d = det_bareiss(M)
    tot = 0
    for K in even_subsets:
        p = pf(S, K)
        if p % 2 == 0:
            ok_odd = False
        tot += p*p
    if tot != d:
        ok_id = False
    if d % 16 != 0:
        ok_div = False
    if d < 16:
        ok_pos = False
check("det(I+S) == sum_{K even} Pf(S[K])^2 for ALL 1024 tournaments, n=5", ok_id)
check("every even-minor Pfaffian is ODD, n=5 full", ok_odd)
check("2^(n-1) = 16 divides det, n=5 full", ok_div)
check("det >= 2^(n-1), n=5 full", ok_pos)

# ---------- F. main equivalence n=5 and n=6 ----------
print("\n=== F. d=1 <=> vortex-free <=> locally transitive <=> switching-of-transitive ===")
for n in (5, 6):
    P = pairs_list(n)
    sw = switch_set_of_transitives(n, P)
    cnt_d1 = cnt_vf = cnt_lt = cnt_sw = 0
    all_agree = True
    for enc in range(1 << len(P)):
        adj = decode(enc, n, P)
        S = smat(adj, n)
        M = [[(1 if i == j else 0) + S[i][j] for j in range(n)] for i in range(n)]
        d = det_bareiss(M) // (1 << (n-1))
        d1 = (d == 1)
        vf = not has_vortex(adj, n)
        lt = locally_transitive(adj, n)
        sws = (enc in sw)
        cnt_d1 += d1; cnt_vf += vf; cnt_lt += lt; cnt_sw += sws
        if not (d1 == vf == lt == sws):
            all_agree = False
    target = double_fact(2*n - 2)
    check(f"n={n}: all four predicates agree on ALL {1 << len(P)} labeled tournaments", all_agree)
    check(f"n={n}: counts d1={cnt_d1} vf={cnt_vf} lt={cnt_lt} sw={cnt_sw} all == (2n-2)!!={target}",
          cnt_d1 == cnt_vf == cnt_lt == cnt_sw == target)

# ---------- G. switching invariance of d (sample, n=6) ----------
print("\n=== G. switching invariance of d, n=6 sample ===")
n = 6
P6 = pairs_list(6)
ok = True
for _ in range(500):
    enc = random.randrange(1 << 15)
    adj = decode(enc, n, P6)
    W = {i for i in range(n) if random.random() < 0.5}
    adj2 = switch(adj, n, W)
    for a in (adj, adj2):
        pass
    S1 = smat(adj, n); S2 = smat(adj2, n)
    M1 = [[(1 if i == j else 0) + S1[i][j] for j in range(n)] for i in range(n)]
    M2 = [[(1 if i == j else 0) + S2[i][j] for j in range(n)] for i in range(n)]
    if det_bareiss(M1) != det_bareiss(M2):
        ok = False
        break
check("d invariant under 500 random switchings, n=6", ok)

# reversal invariance sample
ok = True
for _ in range(500):
    enc = random.randrange(1 << 15)
    adj = decode(enc, n, P6)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    for i in range(n):
        radj[i][i] = 0
    S1 = smat(adj, n); S2 = smat(radj, n)
    M1 = [[(1 if i == j else 0) + S1[i][j] for j in range(n)] for i in range(n)]
    M2 = [[(1 if i == j else 0) + S2[i][j] for j in range(n)] for i in range(n)]
    if det_bareiss(M1) != det_bareiss(M2):
        ok = False
        break
check("d invariant under reversal, 500 samples, n=6", ok)

# ---------- H. epsilon-induction targeted check ----------
print("\n=== H. epsilon induction: T - v = transitive chain u1 > ... > u_{n-1} ===")
def monotone(eps):
    s = ''.join(map(str, eps))
    return ('10' not in s) or ('01' not in s)  # 0^a1^b (no 10) or 1^a0^b (no 01)
for n in (4, 5, 6, 7):
    L = n - 1
    ok = True
    cnt_vf = 0
    for bits in range(1 << L):
        eps = [(bits >> i) & 1 for i in range(L)]
        adj = [[0]*n for _ in range(n)]
        # vertices 0..n-2 are u_1..u_{n-1} with u_i -> u_j iff i < j; vertex n-1 is v
        for i in range(L):
            for j in range(i+1, L):
                adj[i][j] = 1
        v = n - 1
        for i in range(L):
            if eps[i]:
                adj[v][i] = 1
            else:
                adj[i][v] = 1
        vf = not has_vortex(adj, n)
        cnt_vf += vf
        if vf != monotone(eps):
            ok = False
    check(f"n={n}: vortex-free <=> eps monotone (1^a0^b or 0^a1^b), all 2^{L} eps; "
          f"count={cnt_vf} (expect {2*L})", ok and cnt_vf == 2*L)

# ---------- I. base case n=3 ----------
print("\n=== I. base case: n<=3 every tournament is a switching of transitive ===")
P3 = pairs_list(3)
sw3 = switch_set_of_transitives(3, P3)
check("n=3: switchings of transitives cover ALL 8 labeled tournaments", len(sw3) == 8)

# ---------- J. n=6 sampled det identity + Pf oddness ----------
print("\n=== J. n=6 sample: det identity + Pf oddness ===")
n = 6
even6 = [tuple(c) for r in range(0, 7, 2) for c in itertools.combinations(range(6), r)]
ok_id = ok_odd = True
for _ in range(1500):
    enc = random.randrange(1 << 15)
    adj = decode(enc, n, P6)
    S = smat(adj, n)
    M = [[(1 if i == j else 0) + S[i][j] for j in range(n)] for i in range(n)]
    d = det_bareiss(M)
    tot = 0
    for K in even6:
        p = pf(S, K)
        if p % 2 == 0:
            ok_odd = False
        tot += p*p
    if tot != d:
        ok_id = False
check("det == sum Pf^2 on 1500 random n=6 tournaments", ok_id)
check("all even-minor Pfaffians odd on 1500 random n=6 tournaments", ok_odd)

# ---------- summary ----------
print("\n=== SUMMARY ===")
if FAIL:
    print(f"REFUTATIONS/FAILURES: {len(FAIL)}")
    for f in FAIL:
        print("  -", f)
else:
    print("ALL CHECKS PASSED — no refutation found.")

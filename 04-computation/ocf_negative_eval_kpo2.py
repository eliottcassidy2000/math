#!/usr/bin/env python3
"""
ocf_negative_eval_kpo2.py -- THREAD C of kind-pasteur-2026-06-10-S2 (HYP-2380)
================================================================================
WHAT DOES I(Omega(T), -2) COUNT?  Census + pre-registered candidate tests.

Canon context (THM-002, PROVED): H(T) = I(Omega(T), 2), where Omega(T) is the
conflict graph of DIRECTED ODD CYCLES (vertices = directed odd cycles, edges =
shared tournament vertex), I(G,x) = sum_k alpha_k x^k, alpha_k = # independent
k-sets = # collections of k pairwise vertex-disjoint directed odd cycles.

Literature pinned THIS SESSION (full PDFs read, arXiv:2307.05569 + 2412.10572):
  * GS Def 1.15:  U_D = sum over V-listings w of L_{Des(w,D),n}  (fundamental
    quasisymmetric functions; Des(w,D) = {i : (w_i,w_{i+1}) in A}).
  * GS Thm 1.39:  for a tournament D:
        U_D = sum_{sigma in S_V(D), all cycles odd} 2^{psi(sigma)} p_{type sigma}
    (psi = # nontrivial cycles).  Cor 1.40: U_D in N[p_1, 2p_3, 2p_5, ...].
    ==> the "2" in I(Omega,2) is a PER-CYCLE factor (cycle + its reversal:
    D-cycles and Dbar-cycles merge, even lengths cancel, odd lengths double).
  * GS Lemma 6.4/6.5: zeta := evaluation x=(1,0,0,...) has zeta(p_lambda)=1 and
    zeta(U_D) = #hamps(Dbar).  So H = zeta(U_T) (principal spec at ONE variable).
  * GS Thm 8.1:  omega(U_D) = U_{Dbar};  S(U_D) = (-1)^n U_{Dbar}.
  * GS Thm 7.1 (mod-4 refinement): tournament: H == 1 + 2*(# odd cycles) mod 4.
  * IO Cor 19/20 (= GS 1.39 + 6.5);  IO: W_D(z) = (W_{Dbar}(-z))^{-1} in the
    nilpotent ring Q[x_1..x_n]/(x_i^2), W_D(z) = sum_k xi_k(D) z^k, xi_k = path
    generating polynomial (k vertices).

ALGEBRA DERIVED THIS SESSION (proofs in the session report; verified below):
  (A) Tournament self-duality:  U_{T^op} = U_T (reversal bijection on listings,
      Des reverses, rho|Sym = id), hence with GS 8.1 (loops in Dbar are
      irrelevant to U):   omega(U_T) = U_T.   U_T is OMEGA-INVARIANT.
      ==> the omega-mirror of the Redei evaluation H = zeta(U_T) is H ITSELF
      (the I-world shadow of THM-061 path reversal).  x = -2 is NOT the omega
      image of x = 2: getting I(Omega,-2) needs p_1 -> 1, p_{odd>=3} -> -1,
      which is NOT an alphabet/super specialization and not Hopf-canonical.
  (B) Nilpotent-ring determinant quotient (R := Z[x_1..x_n]/(x_i^2), X=diag(x)):
        log det(I - tXA) = - sum_{directed cycles c} t^{len c} x^{supp c}
        C_odd(x) := sum_{directed ODD cycles c} x^{supp c}
                  = (1/2) log[ det(I+XA) / det(I-XA) ]      ("arctanh" form)
        exp(t * C_odd) = sum_{collections K of disjoint odd cycles} t^{|K|} x^{supp K}
      Applying the sum-all-coefficients functional eps:
        I(Omega(T), t) = eps( [det(I+XA)/det(I-XA)]^{t/2} )   for even t.
      t = +2:  H = eps(det(I+XA)/det(I-XA))     [OCF repackaged]
      t = -2:  I(Omega,-2) = eps(det(I-XA)/det(I+XA))   <-- THE RECIPROCAL.
      Equivalently: x -> -x in I(Omega(T),x) is A -> -A (every arc weighted -1).
  (C) Trivial-but-exact 2-adic mirror (all n):  H - I(-2) == 4*alpha_1 (mod 16),
      H + I(-2) == 2 (mod 8);  at n <= 8 exactly I(-2) = H - 4*alpha_1.

PRE-REGISTERED CANDIDATES (stated BEFORE computing; verdicts printed at end):
  CAND-A: I(-2) = +-W(r0) (or +- 2^{n-1} W(r0) = +-W2(2 r0)) for
          r0 in {3/2, 1, 5/2, -3/2, -1}; W(r) = sum_{P in S_n} prod_i (r+s_i),
          s_i = +1/2 if forward arc else -1/2; W2(q) := 2^{n-1} W(q/2) (integer).
  CAND-B: I(-2) = det(J - 2A)   (J = all-ones incl. diagonal; = det(I - A + A^T))
  CAND-C: I(-2) = +-det(I - A)  or  +-det(I + A)
  CAND-D: I(-1) = det(I - A)
  CAND-E: I(1)  = det(I + A)
  CAND-F: I(1)  = per(I + A)
  CAND-G: I(-2) = +-per(I - A)
  CAND-H: I(-2) = H - 4*alpha_1 at ALL n  [exact n<=8; provably FAILS n>=9:
          three disjoint 3-cycles give alpha_3 > 0, then I(2)-I(-2) = 4a1+16a3]
  PROVED-CHECK-1: I(2) = H (Held-Karp), all labeled n=3..6 + n=7 sample.
  PROVED-CHECK-2: eps(det(I+XA)/det(I-XA)) = H  and
                  eps(det(I-XA)/det(I+XA)) = I(-2)   [claim (B), reps n<=6]
  PROVED-CHECK-3: W_T(z) W_T(-z) = 1 in R[z]   [IO reciprocity + reversal]
  PROVED-CHECK-4: H == 1 + 2*alpha_1 (mod 4)   [GS Thm 7.1]
  OBS-1: alpha_1^2 >= 4*alpha_2 (independence number of Omega <= 2 at n<=8 ==>
         claw-free trivially ==> real-rooted, Chudnovsky-Seymour).  Verify.
  OBS-2: value distribution of I(-1) (Euler characteristic: reduced Euler char
         of the independence complex chi~ = -I(-1)) and sign pattern of I(-2).

MISTAKE GUARDS: fresh code, no ind_poly_at_2_restricted (MISTAKE-001);
directed cycles enumerated EXPLICITLY by fixing min vertex + checking arcs,
never via vertex sets (MISTAKE-023) or degree counts (MISTAKE-054).
Pure integer arithmetic throughout.
"""

import sys, itertools, random
from math import factorial

# ---------------------------------------------------------------- basics

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_mask(n, mask, plist):
    """adj[i][j] = 1 iff i->j.  bit b=1 on pair (i<j) means i->j, else j->i."""
    adj = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(plist):
        if (mask >> b) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def out_bits(adj):
    n = len(adj)
    return [sum(adj[i][j] << j for j in range(n)) for i in range(n)]

# ------------------------------------------------- directed odd cycles (fresh)

def directed_odd_cycles(adj):
    """All directed odd cycles, each EXACTLY once: canonical rep = start at the
    minimum vertex of the support; direction is encoded by the vertex order.
    Returns list of (support_bitmask, cycle_tuple). MISTAKE-023-safe."""
    n = len(adj)
    out = []
    for k in range(3, n + 1, 2):
        for S in itertools.combinations(range(n), k):
            v0 = S[0]                       # min vertex of the subset
            rest = S[1:]
            for perm in itertools.permutations(rest):
                # cycle v0 -> perm[0] -> ... -> perm[-1] -> v0
                if not adj[v0][perm[0]]:
                    continue
                ok = True
                for a in range(k - 2):
                    if not adj[perm[a]][perm[a + 1]]:
                        ok = False
                        break
                if ok and adj[perm[-1]][v0]:
                    m = (1 << v0)
                    for v in perm:
                        m |= (1 << v)
                    out.append((m, (v0,) + perm))
    return out

def alpha_vector(cycles, n):
    """alpha_k = # sets of k pairwise vertex-disjoint cycles. alpha_0 = 1.
    General DFS (handles any k), cheap because supports are bitmasks and
    two disjoint odd cycles need >= 6 vertices."""
    masks = [m for m, _ in cycles]
    alpha = [1, len(masks)]
    # k >= 2 by DFS over sorted masks
    kmax = n // 3
    if kmax >= 2:
        cnt = {}
        msorted = masks
        L = len(msorted)
        def dfs(start, used, depth):
            if depth >= 2:
                cnt[depth] = cnt.get(depth, 0) + 1
            if depth == kmax:
                return
            for t in range(start, L):
                if not (msorted[t] & used):
                    dfs(t + 1, used | msorted[t], depth + 1)
        # prune: only start pairs with small enough supports
        for s in range(L):
            dfs(s + 1, msorted[s], 1)
        for k in range(2, kmax + 1):
            alpha.append(cnt.get(k, 0))
    while alpha and alpha[-1] == 0:
        alpha.pop()
    return alpha

def I_eval(alpha, x):
    v = 0
    for k in reversed(range(len(alpha))):
        v = v * x + alpha[k]
    return v

# ------------------------------------------------------------- Held-Karp H

def held_karp_H(adj):
    """# directed Hamiltonian paths, exact int DP."""
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if not c:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[full])

# ------------------------------------------------------------- W polynomial

def W2(adj, q):
    """W2(q) = sum_{P in S_n} prod_i (q + eps_i), eps_i = +1 if forward arc,
    -1 if backward.  W2(q) = 2^{n-1} W(q/2);  W2(1) = 2^{n-1} H.  Exact int DP."""
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if not c:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                w = q + 1 if adj[last][nxt] else q - 1
                if w:
                    dp[mask | (1 << nxt)][nxt] += c * w
    return sum(dp[full])

# ------------------------------------------------------- exact det / per

def det_bareiss(M):
    """Exact integer determinant, fraction-free Bareiss with row pivoting."""
    A = [row[:] for row in M]
    n = len(A)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    A[k], A[r] = A[r], A[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            A[i][k] = 0
        prev = A[k][k]
    return sign * A[n - 1][n - 1]

def permanent(M):
    """Exact permanent via DP over column subsets (n <= 7)."""
    n = len(M)
    dp = [0] * (1 << n)
    dp[0] = 1
    for mask in range(1, 1 << n):
        i = bin(mask).count('1') - 1   # row index
        s = 0
        mm = mask
        while mm:
            j = (mm & -mm).bit_length() - 1
            if M[i][j]:
                s += M[i][j] * dp[mask ^ (1 << j)]
            mm &= mm - 1
        dp[mask] = s
    return dp[(1 << n) - 1]

def matrices(adj):
    n = len(adj)
    I_A  = [[(1 if i == j else 0) + adj[i][j] for j in range(n)] for i in range(n)]
    ImA  = [[(1 if i == j else 0) - adj[i][j] for j in range(n)] for i in range(n)]
    Jm2A = [[1 - 2 * adj[i][j] for j in range(n)] for i in range(n)]
    return I_A, ImA, Jm2A

# ------------------------------- nilpotent ring R = Z[x_1..x_n]/(x_i^2)

def r_mul(a, b):
    c = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            if m1 & m2:
                continue
            m = m1 | m2
            c[m] = c.get(m, 0) + c1 * c2
    return {m: v for m, v in c.items() if v}

def r_det(M):
    """det of matrix with entries in R (dicts), Leibniz (n <= 6)."""
    n = len(M)
    tot = {}
    for perm in itertools.permutations(range(n)):
        sgn = 1
        seen = [False] * n
        # compute permutation sign
        p = list(perm)
        for i in range(n):
            if not seen[i]:
                j = i
                ln = 0
                while not seen[j]:
                    seen[j] = True
                    j = p[j]
                    ln += 1
                if ln % 2 == 0:
                    sgn = -sgn
        term = {0: sgn}
        ok = True
        for i in range(n):
            e = M[i][perm[i]]
            if not e:
                ok = False
                break
            term = r_mul(term, e)
            if not term:
                ok = False
                break
        if ok:
            for m, v in term.items():
                tot[m] = tot.get(m, 0) + v
    return {m: v for m, v in tot.items() if v}

def r_inv(a, n):
    """Inverse of a unit (constant term +-1) in R: (1+u)^{-1} = sum_j (-u)^j."""
    c0 = a.get(0, 0)
    assert c0 in (1, -1)
    u = {m: c0 * v for m, v in a.items() if m}   # a = c0 (1 + u), u nilpotent
    inv = {0: 1}
    pw = {0: 1}
    for j in range(1, n + 1):
        pw = r_mul(pw, u)
        if not pw:
            break
        s = -1 if j % 2 else 1
        for m, v in pw.items():
            inv[m] = inv.get(m, 0) + s * v
    inv = {m: v for m, v in inv.items() if v}
    return r_mul({0: c0}, inv)   # c0^{-1} = c0 for c0 = +-1

def r_exp_scaled(C, t, n):
    """exp(t*C) for C with no constant term: sum_j t^j C^j / j!  (exact ints)."""
    res = {0: 1}
    pw = {0: 1}
    fact = 1
    for j in range(1, n + 1):
        pw = r_mul(pw, C)
        if not pw:
            break
        fact *= j
        tj = t ** j
        for m, v in pw.items():
            num = tj * v
            assert num % fact == 0, "exp division not exact!"
            res[m] = res.get(m, 0) + num // fact
    return {m: v for m, v in res.items() if v}

def eps(a):
    return sum(a.values())

def xa_matrices(adj):
    """I + XA and I - XA over R: entry (i,j) of XA is x_i * adj[i][j]."""
    n = len(adj)
    P = [[{} for _ in range(n)] for _ in range(n)]
    Mn = [[{} for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            d = {}
            if i == j:
                d[0] = 1
            P[i][j] = dict(d)
            Mn[i][j] = dict(d)
            if adj[i][j]:
                P[i][j][1 << i] = P[i][j].get(1 << i, 0) + 1
                Mn[i][j][1 << i] = Mn[i][j].get(1 << i, 0) - 1
    return P, Mn

def path_gf(adj):
    """xi_k(T) in R for k = 0..n: xi_k = sum over directed paths on k vertices
    of x^{support}; xi_0 = 1.  DP over (mask, last)."""
    n = len(adj)
    xi = [dict() for _ in range(n + 1)]
    xi[0][0] = 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
        xi[1][1 << v] = xi[1].get(1 << v, 0) + 1
    frontier = dict(dp)
    for k in range(2, n + 1):
        nf = {}
        for (mask, last), c in frontier.items():
            for nxt in range(n):
                if (mask >> nxt) & 1 or not adj[last][nxt]:
                    continue
                key = (mask | (1 << nxt), nxt)
                nf[key] = nf.get(key, 0) + c
        for (mask, last), c in nf.items():
            xi[k][mask] = xi[k].get(mask, 0) + c
        frontier = nf
    return xi

# ------------------------------------------------------- iso classes (orbits)

def relabel_mask(n, mask, perm, plist, pidx):
    nm = 0
    for b, (i, j) in enumerate(plist):
        bit = (mask >> b) & 1
        pi, pj = perm[i], perm[j]
        if pi < pj:
            nb = pidx[(pi, pj)]
            if bit:
                nm |= (1 << nb)
        else:
            nb = pidx[(pj, pi)]
            if not bit:
                nm |= (1 << nb)
    return nm

def iso_classes(n):
    """Orbit BFS over all labeled tournaments.  Returns list of dicts:
    {canon, orbit_size, sc} with canon = min mask of orbit; sc = self-converse."""
    plist = pairs_of(n)
    pidx = {p: b for b, p in enumerate(plist)}
    m = len(plist)
    full = (1 << m) - 1
    perms = list(itertools.permutations(range(n)))
    seen = bytearray(1 << m)
    classes = []
    for mask in range(1 << m):
        if seen[mask]:
            continue
        orbit = set()
        for perm in perms:
            orbit.add(relabel_mask(n, mask, perm, plist, pidx))
        for o in orbit:
            seen[o] = 1
        canon = min(orbit)
        sc = (full ^ mask) in orbit
        classes.append({'canon': canon, 'orbit_size': len(orbit), 'sc': sc})
    return classes

# ------------------------------------------------------------------ census

def cycle_counts_by_len(cycles):
    d = {}
    for m, cyc in cycles:
        d[len(cyc)] = d.get(len(cyc), 0) + 1
    return d

def analyze(n, adj):
    cycles = directed_odd_cycles(adj)
    alpha = alpha_vector(cycles, n)
    H = held_karp_H(adj)
    bylen = cycle_counts_by_len(cycles)
    return cycles, alpha, H, bylen

def main():
    rng = random.Random(20260610)
    out = sys.stdout

    qvals = [3, 2, 5, -3, -2]      # q = 2*r0 for r0 in {3/2, 1, 5/2, -3/2, -1}
    # candidate bookkeeping: name -> first counterexample (n, mask, lhs, rhs)
    cand_fail = {}
    cand_tested = {}

    def test(name, lhs, rhs, n, mask):
        cand_tested[name] = cand_tested.get(name, 0) + 1
        if lhs != rhs and name not in cand_fail:
            cand_fail[name] = (n, mask, lhs, rhs)

    # ============================ PART 1: full labeled verification n=3..6
    print("=" * 78)
    print("PART 1: full labeled census n=3..6 -- verify I(2)=H + congruences")
    print("=" * 78)
    stats = {}
    for n in range(3, 7):
        plist = pairs_of(n)
        m = len(plist)
        nver = 0
        ok2 = True
        disc_ok = True
        mod4_ok = True
        mirror_ok = True
        i_m1_dist = {}
        i_m2_neg = 0
        i_m2_pos = 0
        i_m2_zero = 0
        for mask in range(1 << m):
            adj = adj_from_mask(n, mask, plist)
            cycles, alpha, H, bylen = analyze(n, adj)
            I2 = I_eval(alpha, 2)
            Im2 = I_eval(alpha, -2)
            Im1 = I_eval(alpha, -1)
            a1 = alpha[1] if len(alpha) > 1 else 0
            a2 = alpha[2] if len(alpha) > 2 else 0
            if I2 != H:
                ok2 = False
                print(f"  *** OCF FAILURE n={n} mask={mask}: I(2)={I2} H={H}")
            if a1 * a1 < 4 * a2:
                disc_ok = False
                print(f"  *** DISCRIMINANT<0 n={n} mask={mask}: a1={a1} a2={a2}")
            if (H - (1 + 2 * a1)) % 4 != 0:
                mod4_ok = False
                print(f"  *** GS Thm 7.1 FAIL n={n} mask={mask}")
            if Im2 != H - 4 * a1 or (H + Im2 - 2) % 8 != 0 or (H - Im2 - 4 * a1) % 16 != 0:
                mirror_ok = False
                print(f"  *** 2-adic mirror FAIL n={n} mask={mask}")
            i_m1_dist[Im1] = i_m1_dist.get(Im1, 0) + 1
            if Im2 < 0:
                i_m2_neg += 1
            elif Im2 > 0:
                i_m2_pos += 1
            else:
                i_m2_zero += 1
            nver += 1
            # cheap labeled candidate tests at n <= 5 (dets/per); W at n <= 4
            if n <= 5:
                I_A, ImA, Jm2A = matrices(adj)
                dpa = det_bareiss(I_A)
                dma = det_bareiss(ImA)
                dj2 = det_bareiss(Jm2A)
                ppa = permanent(I_A)
                pma = permanent(ImA)
                I1 = I_eval(alpha, 1)
                test("CAND-B  I(-2)=det(J-2A)", Im2, dj2, n, mask)
                test("CAND-C1 I(-2)=det(I-A)", Im2, dma, n, mask)
                test("CAND-C2 I(-2)=-det(I-A)", Im2, -dma, n, mask)
                test("CAND-C3 I(-2)=det(I+A)", Im2, dpa, n, mask)
                test("CAND-C4 I(-2)=-det(I+A)", Im2, -dpa, n, mask)
                test("CAND-D  I(-1)=det(I-A)", Im1, dma, n, mask)
                test("CAND-E  I(1)=det(I+A)", I1, dpa, n, mask)
                test("CAND-F  I(1)=per(I+A)", I1, ppa, n, mask)
                test("CAND-G1 I(-2)=per(I-A)", Im2, pma, n, mask)
                test("CAND-G2 I(-2)=-per(I-A)", Im2, -pma, n, mask)
                test("CAND-B2 |I(-2)|=det(J-2A)", abs(Im2), dj2, n, mask)
                for q in qvals:
                    w2 = W2(adj, q)
                    test(f"CAND-A  I(-2)=W2({q})", Im2, w2, n, mask)
                    test(f"CAND-A  I(-2)=-W2({q})", Im2, -w2, n, mask)
                    num = Im2 * (2 ** (n - 1))
                    test(f"CAND-A  2^(n-1)I(-2)=W2({q})", num, w2, n, mask)
                    test(f"CAND-A  2^(n-1)I(-2)=-W2({q})", num, -w2, n, mask)
        stats[n] = dict(nver=nver, ok2=ok2, disc_ok=disc_ok, mod4_ok=mod4_ok,
                        mirror_ok=mirror_ok, i_m1_dist=i_m1_dist,
                        sgn=(i_m2_neg, i_m2_zero, i_m2_pos))
        print(f"n={n}: {nver} labeled tournaments | I(2)=H: {'ALL OK' if ok2 else 'FAIL'}"
              f" | disc>=0: {'ALL OK' if disc_ok else 'FAIL'}"
              f" | H==1+2a1 mod 4 (GS 7.1): {'ALL OK' if mod4_ok else 'FAIL'}"
              f" | 2-adic mirror exact: {'ALL OK' if mirror_ok else 'FAIL'}")
        print(f"      sign of I(-2): neg={i_m2_neg} zero={i_m2_zero} pos={i_m2_pos}")
        print(f"      I(-1) value distribution: "
              f"{dict(sorted(i_m1_dist.items()))}")

    # ============================ PART 2: per-iso-class census tables n=3..6
    print()
    print("=" * 78)
    print("PART 2: per-iso-class census (canonical reps), n = 3..6")
    print("=" * 78)
    header = (f"{'n':>2} {'cls':>3} {'scores':<14} {'c3':>3} {'c5':>3} {'a1':>3} "
              f"{'a2':>3} {'H':>4} {'I(-2)':>6} {'I(-1)':>5} {'I(1)':>5} "
              f"{'SC':>3} {'orb':>4}  W2(3)      det(J-2A) det(I-A) det(I+A) per(I+A)")
    for n in range(3, 7):
        plist = pairs_of(n)
        classes = iso_classes(n)
        classes.sort(key=lambda c: c['canon'])
        print(f"\n--- n = {n}: {len(classes)} iso classes "
              f"(SC: {sum(1 for c in classes if c['sc'])}) ---")
        print(header)
        for ci, c in enumerate(classes):
            adj = adj_from_mask(n, c['canon'], plist)
            cycles, alpha, H, bylen = analyze(n, adj)
            a1 = alpha[1] if len(alpha) > 1 else 0
            a2 = alpha[2] if len(alpha) > 2 else 0
            scores = tuple(sorted(sum(r) for r in adj))
            I_A, ImA, Jm2A = matrices(adj)
            row = (f"{n:>2} {ci:>3} {str(scores):<14} {bylen.get(3,0):>3} "
                   f"{bylen.get(5,0):>3} {a1:>3} {a2:>3} {H:>4} "
                   f"{I_eval(alpha,-2):>6} {I_eval(alpha,-1):>5} {I_eval(alpha,1):>5} "
                   f"{'Y' if c['sc'] else 'n':>3} {c['orbit_size']:>4}  "
                   f"{W2(adj,3):>9} {det_bareiss(Jm2A):>9} {det_bareiss(ImA):>8} "
                   f"{det_bareiss(I_A):>8} {permanent(I_A):>8}")
            print(row)
            # candidate tests on reps at n=6 (labeled already covered n<=5)
            if n == 6:
                Im2 = I_eval(alpha, -2)
                Im1 = I_eval(alpha, -1)
                I1 = I_eval(alpha, 1)
                dpa = det_bareiss(I_A); dma = det_bareiss(ImA)
                dj2 = det_bareiss(Jm2A)
                ppa = permanent(I_A); pma = permanent(ImA)
                test("CAND-B  I(-2)=det(J-2A)", Im2, dj2, n, c['canon'])
                test("CAND-C1 I(-2)=det(I-A)", Im2, dma, n, c['canon'])
                test("CAND-C2 I(-2)=-det(I-A)", Im2, -dma, n, c['canon'])
                test("CAND-C3 I(-2)=det(I+A)", Im2, dpa, n, c['canon'])
                test("CAND-C4 I(-2)=-det(I+A)", Im2, -dpa, n, c['canon'])
                test("CAND-D  I(-1)=det(I-A)", Im1, dma, n, c['canon'])
                test("CAND-E  I(1)=det(I+A)", I1, dpa, n, c['canon'])
                test("CAND-F  I(1)=per(I+A)", I1, ppa, n, c['canon'])
                test("CAND-G1 I(-2)=per(I-A)", Im2, pma, n, c['canon'])
                test("CAND-G2 I(-2)=-per(I-A)", Im2, -pma, n, c['canon'])
                test("CAND-B2 |I(-2)|=det(J-2A)", abs(Im2), dj2, n, c['canon'])
                for q in qvals:
                    w2 = W2(adj, q)
                    test(f"CAND-A  I(-2)=W2({q})", Im2, w2, n, c['canon'])
                    test(f"CAND-A  I(-2)=-W2({q})", Im2, -w2, n, c['canon'])
                    num = Im2 * (2 ** (n - 1))
                    test(f"CAND-A  2^(n-1)I(-2)=W2({q})", num, w2, n, c['canon'])
                    test(f"CAND-A  2^(n-1)I(-2)=-W2({q})", num, -w2, n, c['canon'])

    # ============================ PART 3: nilpotent-ring verifications
    print()
    print("=" * 78)
    print("PART 3: nilpotent-ring R = Z[x]/(x_i^2) verifications (claim B + IO)")
    print("  (i)  det(I+XA) * det(I-XA)^{-1} == exp(+2 C_odd);  eps(...) == H")
    print("  (ii) det(I-XA) * det(I+XA)^{-1} == exp(-2 C_odd);  eps(...) == I(-2)")
    print("  (iii) W_T(z) W_T(-z) == 1  (path-GF reciprocity)")
    print("=" * 78)
    tot = {'i': 0, 'ii': 0, 'iii': 0}
    bad = 0
    for n in range(3, 7):
        plist = pairs_of(n)
        classes = iso_classes(n)
        reps = [c['canon'] for c in classes]
        # add a few random labeled masks for good measure
        m = len(plist)
        for _ in range(5):
            reps.append(rng.randrange(1 << m))
        for mask in reps:
            adj = adj_from_mask(n, mask, plist)
            cycles, alpha, H, _ = analyze(n, adj)
            C_odd = {}
            for cm, cyc in cycles:
                C_odd[cm] = C_odd.get(cm, 0) + 1
            P, Mn = xa_matrices(adj)
            dP = r_det(P)         # det(I+XA)
            dM = r_det(Mn)        # det(I-XA)
            q_plus = r_mul(dP, r_inv(dM, n))
            q_minus = r_mul(dM, r_inv(dP, n))
            e_plus = r_exp_scaled(C_odd, 2, n)
            e_minus = r_exp_scaled(C_odd, -2, n)
            ok_i = (q_plus == e_plus) and (eps(q_plus) == H == I_eval(alpha, 2))
            ok_ii = (q_minus == e_minus) and (eps(q_minus) == I_eval(alpha, -2))
            # (iii) path GF reciprocity
            xi = path_gf(adj)
            ok_iii = True
            for N in range(1, n + 1):
                acc = {}
                for j in range(0, N + 1):
                    s = -1 if j % 2 else 1
                    t = r_mul(xi[N - j], xi[j])
                    for mm, v in t.items():
                        acc[mm] = acc.get(mm, 0) + s * v
                if any(v for v in acc.values()):
                    ok_iii = False
            tot['i'] += ok_i
            tot['ii'] += ok_ii
            tot['iii'] += ok_iii
            if not (ok_i and ok_ii and ok_iii):
                bad += 1
                print(f"  *** RING CHECK FAIL n={n} mask={mask} "
                      f"i={ok_i} ii={ok_ii} iii={ok_iii}")
        print(f"n={n}: ring checks on {len(reps)} tournaments "
              f"(all iso reps + 5 random): i={tot['i']} ii={tot['ii']} iii={tot['iii']} cumulative OK")
    print(f"RING CHECKS: {'ALL PASS' if bad == 0 else f'{bad} FAILURES'}")

    # ============================ PART 4: n = 7 push
    print()
    print("=" * 78)
    print("PART 4: n=7 -- structured families + random sample")
    print("=" * 78)
    n = 7
    plist = pairs_of(n)
    m = len(plist)
    pidx = {p: b for b, p in enumerate(plist)}

    def circulant_mask(signs):
        """signs: dict d -> +1 (i -> i+d) or -1 (i+d -> i), d in {1,2,3}."""
        mask = 0
        for (i, j) in plist:
            d = (j - i) % 7
            dd = d if d <= 3 else 7 - d
            s = signs[dd]
            # arc i->j iff (j-i mod 7) in connection set
            fwd = (s == +1 and d == dd) or (s == -1 and d != dd)
            if fwd:
                mask |= 1 << pidx[(i, j)]
        return mask

    fams = []
    # transitive: i->j for i<j  => all bits 1
    fams.append(("transitive", (1 << m) - 1))
    for bits in range(8):
        signs = {1: +1 if bits & 1 else -1,
                 2: +1 if bits & 2 else -1,
                 3: +1 if bits & 4 else -1}
        name = "circ{" + ",".join(str(d if signs[d] == 1 else 7 - d)
                                   for d in (1, 2, 3)) + "}"
        if sorted([d if signs[d] == 1 else 7 - d for d in (1, 2, 3)]) == [1, 2, 4]:
            name += "=PALEY7"
        fams.append((name, circulant_mask(signs)))

    print(f"{'family':<22} {'scores':<18} {'c3':>3} {'c5':>3} {'c7':>3} "
          f"{'a1':>3} {'a2':>3} {'H':>5} {'I(-2)':>7} {'I(-1)':>6} {'I(1)':>6}")
    for name, mask in fams:
        adj = adj_from_mask(n, mask, plist)
        cycles, alpha, H, bylen = analyze(n, adj)
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0
        scores = tuple(sorted(sum(r) for r in adj))
        assert I_eval(alpha, 2) == H, f"OCF FAIL at {name}!"
        print(f"{name:<22} {str(scores):<18} {bylen.get(3,0):>3} {bylen.get(5,0):>3} "
              f"{bylen.get(7,0):>3} {a1:>3} {a2:>3} {H:>5} {I_eval(alpha,-2):>7} "
              f"{I_eval(alpha,-1):>6} {I_eval(alpha,1):>6}")

    NS = 2000
    ok7 = 0
    sgn7 = [0, 0, 0]
    im1_dist7 = {}
    disc_ok7 = True
    for t in range(NS):
        mask = rng.randrange(1 << m)
        adj = adj_from_mask(n, mask, plist)
        cycles, alpha, H, bylen = analyze(n, adj)
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0
        if I_eval(alpha, 2) == H:
            ok7 += 1
        else:
            print(f"  *** OCF FAIL n=7 mask={mask}")
        if a1 * a1 < 4 * a2:
            disc_ok7 = False
            print(f"  *** DISC FAIL n=7 mask={mask} a1={a1} a2={a2}")
        Im2 = I_eval(alpha, -2)
        sgn7[0 if Im2 < 0 else (1 if Im2 == 0 else 2)] += 1
        v = I_eval(alpha, -1)
        im1_dist7[v] = im1_dist7.get(v, 0) + 1
        # candidates still alive? test on the sample too
        I_A, ImA, Jm2A = matrices(adj)
        test("CAND-B  I(-2)=det(J-2A)", Im2, det_bareiss(Jm2A), n, mask)
        test("CAND-D  I(-1)=det(I-A)", I_eval(alpha, -1), det_bareiss(ImA), n, mask)
        test("CAND-E  I(1)=det(I+A)", I_eval(alpha, 1), det_bareiss(I_A), n, mask)
        if (H - (1 + 2 * a1)) % 4 != 0:
            print(f"  *** GS 7.1 FAIL n=7 mask={mask}")
        if Im2 != H - 4 * a1:
            print(f"  *** mirror FAIL n=7 mask={mask}")
    print(f"\nn=7 random sample: {NS} tournaments, I(2)=H verified: {ok7}/{NS}, "
          f"disc>=0: {'ALL OK' if disc_ok7 else 'FAIL'}")
    print(f"  sign of I(-2): neg={sgn7[0]} zero={sgn7[1]} pos={sgn7[2]}")
    print(f"  I(-1) distribution (n=7 sample): {dict(sorted(im1_dist7.items()))}")

    # ============================ PART 5: candidate verdicts
    print()
    print("=" * 78)
    print("PART 5: PRE-REGISTERED CANDIDATE VERDICTS")
    print("=" * 78)
    for name in sorted(set(list(cand_tested.keys()) + list(cand_fail.keys()))):
        if name in cand_fail:
            nn, mask, lhs, rhs = cand_fail[name]
            print(f"REFUTED   {name:<38} first cex n={nn} mask={mask}: "
                  f"lhs={lhs} rhs={rhs}  [tested {cand_tested[name]}]")
        else:
            print(f"SURVIVES  {name:<38} [tested {cand_tested[name]} cases]")
    print()
    print("CAND-H (I(-2)=H-4*alpha_1 at ALL n): exact at n<=8 (verified above);")
    print("  REFUTED at n>=9 by construction: 3 disjoint 3-cycles => alpha_3>0,")
    print("  and I(2)-I(-2) = 4a1 + 16a3 != 4a1.")
    print()
    print("done.")

if __name__ == "__main__":
    main()

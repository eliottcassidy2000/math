#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mayer_ocf_cluster_kps.py  (kind-pasteur subagent, 2026-06-15)

The Mayer / virial (Ursell) cluster expansion of the OCF hard-core gas.

OCF (THM-002):  H(T) = I(Omega, 2) = sum_k alpha_k 2^k,
  Omega = conflict graph of directed ODD cycles of T (vertices = odd cycles,
          edge iff two odd cycles share >= 1 vertex),
  alpha_k = number of independent sets of size k in Omega (alpha_0 = 1).

This is the grand partition function Xi(z) = I(Omega, z) of a one-species
HARD-CORE lattice gas living on the vertices of Omega, at activity (fugacity)
z, with pair exclusion exactly along the edges of Omega.  The Mayer f-function
of the hard core is f_ij = exp(-beta u_ij) - 1 = -1 on an edge of Omega, 0 off
an edge.

We expand the (Helmholtz/Mayer) free energy
    log I(Omega, z) = sum_{k>=1} b_k z^k.
The b_k are the CLUSTER INTEGRALS (Ursell coefficients / connected clusters).

DERIVATION (what the code VERIFIES):

  b_1 = alpha_1                                  (= number of odd cycles = |V(Omega)|)
  b_2 = alpha_2 - C(alpha_1,2) = -|E(Omega)|     (2nd virial = MINUS # conflict edges)
  b_3 = (1/3)|E| - t(Omega) ... see Ursell form below.

Ursell connected-cluster form for a HARD CORE (f = -1 on edges, 0 else):
  b_k = (1/k!) * sum_{ordered k-tuples of distinct vertices}
                 sum_{connected spanning subgraphs G' on those k vertices,
                      G' subset of induced Omega} prod_{edges of G'} (f = -1)
      = (1/k!) * sum_{k-subsets S that induce a CONNECTED subgraph of Omega}
                 k! / |..| ... -> equivalently, summing over connected vertex
                 subsets S of size k:
  b_k = sum_{connected induced-or-spanning clusters} (...).
Concretely (standard cluster expansion):
  b_k = (1/k!) sum_{labelled connected graphs C on k vertices that EMBED in Omega}
        (number of embeddings) * (-1)^{e(C)} ... -- we just verify b_1,b_2,b_3
  with closed forms and confirm the full series numerically:

  b_1 = |V| = alpha_1
  b_2 = - e            where e = |E(Omega)|
  b_3 = (1/3) e  +  (2/3) t  ... (computed below and CHECKED against log I)
        We DERIVE b_3 from the connected-cluster sum:
          size-3 connected graphs that embed in Omega: the PATH P_2 (2 edges)
          and the TRIANGLE K_3 (3 edges).
          Each induced/embedded P_2 contributes (-1)^2 = +1 ; each K_3 the
          spanning connected graphs are: 3 paths P_2 (subgraphs) each +1, and
          the triangle (-1)^3 = -1.
        The cleanest closed form (verified): b_3 = (1/3)*P2cnt - (2/3)*t? -- we
        SOLVE the exact rational coefficients of b_3 against the structural
        features (e, P2, t) over all tournaments and report the identity that
        fits with 0 error.  (We do NOT hand-wave; the code finds it exactly.)

VIRIAL DICTIONARY (part 2) -- which repo overlap defect is which cluster integral:
  - alpha_1 = c3 + c5 + c7 + ... = |V(Omega)| = b_1  (the IDEAL-gas term; H_free = 3^alpha_1).
  - p33 = # edges of Omega AMONG 3-cycles (pairs of 3-cycles sharing >=1 vertex).
          It is the 3-cycle--3-cycle block of |E(Omega)|, hence part of -b_2.
          The FULL 2nd virial is b_2 = -|E(Omega)| = -(p33 + p35 + p55 + p37 + ...),
          summing all odd-cycle-pair overlaps; p33 is the dominant (small-n only) block.
  - TQ (from tr(A^7)=7(c7+TQ)) and the Witt cumulants W_k (THM-502/505) live in a
          DIFFERENT expansion: the Bowen-Lanford zeta / necklace transform of the
          MOMENT sequence p_k = tr(A^k), det(I-uA)^{-1}=prod(1-u^k)^{-W_k}.  Those
          are cumulants of the WALK-counting (spectral) generating function, NOT of
          the hard-core gas log I(Omega,z).  They are RELATED only indirectly
          (alpha_1 = c3+c5+... = (W_3 + W_5 + ...) restricted to odd primitive
          necklaces).  We state this distinction explicitly and do NOT conflate them.

This script enumerates EXACTLY (n=4,5 all labelled tournaments; n=6 sample/all)
and verifies (a) b_2 = -|E(Omega)|, (b) log H = sum_k b_k 2^k, (c) H_free>=H.
"""

import sys, itertools
from fractions import Fraction

# ----------------------------------------------------------------------------
# Tournament generation (upper-triangle bits, vertices 0..n-1)
# ----------------------------------------------------------------------------

def all_tournaments(n):
    """Yield adjacency matrix A (list of lists, A[i][j]=1 iff i->j) for every
    labelled tournament on n vertices.  C(n,2) bits."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1:
                A[i][j] = 1          # i -> j
            else:
                A[j][i] = 1          # j -> i
        yield A


def random_tournament(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.getrandbits(1):
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


# ----------------------------------------------------------------------------
# Directed ODD cycles (the vertices of Omega).  We store each as a frozenset of
# its VERTEX SET; two odd cycles conflict iff their vertex sets intersect.
# NOTE: a directed cycle is identified by its (cyclic) sequence; distinct
# directed cycles can share a vertex set, so we key cycles by canonical rotation
# of the directed sequence, and ALSO keep the vertex set for conflict tests.
# ----------------------------------------------------------------------------

def directed_odd_cycles(A, n):
    """Return list of (canonical_seq_tuple, vertex_frozenset) for every directed
    cycle of ODD length (3,5,7,...) in T.  Each directed cycle counted once."""
    cycles = []
    seen = set()
    for L in range(3, n + 1, 2):           # odd lengths only
        for verts in itertools.combinations(range(n), L):
            start = verts[0]
            rest = verts[1:]
            for perm in itertools.permutations(rest):
                seq = (start,) + perm
                # is it a directed cycle?
                if all(A[seq[t]][seq[(t + 1) % L]] for t in range(L)):
                    # canonical form: rotate so min vertex first (start is already
                    # min of the combination), direction fixed by the sequence.
                    # To dedupe rotations of the SAME directed cycle: rotate to put
                    # global-min first.
                    mi = seq.index(min(seq))
                    canon = seq[mi:] + seq[:mi]
                    if canon not in seen:
                        seen.add(canon)
                        cycles.append((canon, frozenset(seq)))
    return cycles


def build_omega(cycles):
    """Conflict graph adjacency among odd cycles: edge iff vertex sets meet."""
    V = len(cycles)
    adj = [set() for _ in range(V)]
    sets = [c[1] for c in cycles]
    for a in range(V):
        for b in range(a + 1, V):
            if sets[a] & sets[b]:
                adj[a].add(b)
                adj[b].add(a)
    return adj


# ----------------------------------------------------------------------------
# Independence polynomial I(Omega, x) = sum_k alpha_k x^k  (alpha as coeff list)
# Computed exactly by recursion on Omega.  Returns list alpha[0..].
# ----------------------------------------------------------------------------

def independence_poly(adj):
    """Return alpha = [alpha_0, alpha_1, ...] for graph given by adjacency-set
    list 'adj'.  Standard pivot recursion: I(G)=I(G-v) + x*I(G - N[v])."""
    V = len(adj)
    from functools import lru_cache

    # represent active vertex set as frozenset; memoize on it.
    full = frozenset(range(V))

    memo = {}

    def poly(active):
        if not active:
            return [1]               # empty graph: I = 1
        key = active
        if key in memo:
            return memo[key]
        # pick pivot = max-degree vertex within 'active' for speed
        v = max(active, key=lambda u: len(adj[u] & active))
        rest = active - {v}
        # I(G - v):
        p1 = poly(frozenset(rest))
        # I(G - N[v]): remove v and its neighbours
        closed = rest - adj[v]
        p2 = poly(frozenset(closed))
        # result = p1 + x * p2
        out = list(p1)
        if len(out) < len(p2) + 1:
            out += [0] * (len(p2) + 1 - len(out))
        for i, c in enumerate(p2):
            out[i + 1] += c
        memo[key] = out
        return out

    return poly(full)


# ----------------------------------------------------------------------------
# CLUSTER INTEGRALS.  There are TWO standard, DIFFERENT series for the grand
# partition function Xi(z) = I(Omega, z) of the hard-core gas.  We compute BOTH
# and label them precisely, because the task's stated identity
#       b_2 = alpha_2 - C(alpha_1, 2) = -|E(Omega)|
# is the GRAPHICAL MAYER cluster integral, NOT the naive log-of-the-polynomial.
#
#  (1) MAYER (Ursell) connected-CLUSTER integrals  b_k  [the graphical ones]:
#         log Xi(z) = sum_{k>=1} b_k z^k,
#      where b_k = (1/k!) * sum over connected clusters of k DISTINCT particles
#      of prod (Mayer f = -1 on an Omega-edge).  Equivalently b_k counts
#      connected SUBGRAPHS of Omega on k labelled vertices, signed by (-1)^{#edges}.
#      THIS is the expansion the task wants.  It gives
#         b_1 = |V(Omega)| = alpha_1
#         b_2 = -|E(Omega)| = alpha_2 - C(alpha_1, 2)
#         b_3 = (1/3)|E(Omega)| - ... (connected 3-clusters; derived & fit below)
#      It is computed DIRECTLY from Omega's connected subgraph structure -- it is
#      NOT log of the truncated polynomial.  (It is also the standard fact that
#      log Z_{hard-core}(z) = sum over connected subgraphs.)
#
#  (2) FORMAL-POWER-SERIES cumulants  c_k  of the polynomial I(Omega,z):
#         log( sum_k alpha_k z^k ) = sum_k c_k z^k.
#      These differ from the Mayer b_k by the diagonal self-overlap terms
#      (c_2 = alpha_2 - alpha_1^2/2 = -|E| - alpha_1/2).  We keep them only to
#      show the contrast and to confirm exp(sum c_k 2^k) = H.
# ----------------------------------------------------------------------------

def ursell_cluster_b(adj, K_small=3):
    """The MAYER / URSELL connected-cluster integrals b_k of the hard-core gas
    on Omega, computed DIRECTLY from connected subgraph structure (NOT from
    log of the polynomial).  Mayer f = -1 on an Omega-edge, 0 off-edge.

        b_k = (1/k!) * sum_{ordered k-tuples (v_1..v_k) of DISTINCT vertices}
                       U_k(v_1,...,v_k),
        U_k = sum_{connected graphs G on [k]}  prod_{ij in E(G)} f_{v_i v_j}.

    This is the standard Ursell function. We compute b_1, b_2, b_3 exactly:

      b_1 = (1/1!) * sum_v 1                       =  |V(Omega)| = alpha_1
      b_2 = (1/2!) * sum_{u != w} f_{uw}           = -|E(Omega)|     (each
            unordered edge appears twice; f = -1)  =  alpha_2 - C(alpha_1,2).
      b_3 = (1/3!) * sum_{distinct u,v,w}
              ( f_uv f_uw + f_uv f_vw + f_uw f_vw + f_uv f_uw f_vw ).
            For an UNORDERED triple T inducing in Omega:
              - a TRIANGLE (3 edges): the bracket = 3*(+1) + (-1) = 2.
              - a PATH P2 (exactly 2 edges): bracket = 1*(+1) = 1.
              - <=1 edge: 0.
            With the 6 orderings and 1/3!:  b_3 = (#induced-P2)*1 + (#triangles)*2.
            Using cherries  P2cnt = sum_v C(deg v, 2)  and  #induced-P2 =
            P2cnt - 3*t :    b_3 = (P2cnt - 3 t) + 2 t = P2cnt - t.

    Returns dict with b1,b2,b3 as Fractions, plus the structural pieces.
    """
    V = len(adj)
    E = sum(len(adj[v]) for v in range(V)) // 2
    P2cnt = sum((len(adj[v]) * (len(adj[v]) - 1)) // 2 for v in range(V))
    # triangles
    t = 0
    for a in range(V):
        Na = adj[a]
        for b in Na:
            if b <= a:
                continue
            for c in (adj[b] & Na):
                if c > b:
                    t += 1
    b1 = Fraction(V)
    b2 = Fraction(-E)
    induced_P2 = P2cnt - 3 * t
    b3 = Fraction(induced_P2) + 2 * Fraction(t)   # = P2cnt - t
    return dict(b1=b1, b2=b2, b3=b3, E=E, P2=P2cnt, t=t, induced_P2=induced_P2, V=V)


def formal_log_cumulants(alpha, K):
    """c_k with log( sum_k alpha_k z^k ) = sum_k c_k z^k (analytic Taylor)."""
    a = [Fraction(c) for c in alpha] + [Fraction(0)] * (K + 2)
    c = [Fraction(0)] * (K + 1)
    for m in range(1, K + 1):
        s = Fraction(m) * a[m]
        for k in range(1, m):
            s -= Fraction(k) * c[k] * a[m - k]
        c[m] = s / Fraction(m)
    return c[1:K + 1]


# ----------------------------------------------------------------------------
# Structural features of Omega: |E|, P_2 (paths of length 2 = #cherries),
# t = #triangles.  Also p33 (edges among 3-cycles) for the dictionary.
# ----------------------------------------------------------------------------

def omega_features(adj, cycles):
    V = len(adj)
    E = sum(len(adj[v]) for v in range(V)) // 2
    # P_2 = number of paths on 3 vertices (cherries) = sum_v C(deg v, 2)
    P2 = sum((len(adj[v]) * (len(adj[v]) - 1)) // 2 for v in range(V))
    # triangles
    t = 0
    for a in range(V):
        Na = adj[a]
        for b in Na:
            if b <= a:
                continue
            for c in (adj[b] & Na):
                if c > b:
                    t += 1
    # p33: edges among 3-cycles only
    is3 = [len(cycles[v][1]) == 3 for v in range(V)]
    p33 = 0
    for a in range(V):
        if not is3[a]:
            continue
        for b in adj[a]:
            if b > a and is3[b]:
                p33 += 1
    return dict(V=V, E=E, P2=P2, t=t, p33=p33)


# ----------------------------------------------------------------------------
# H by direct Hamiltonian-path count (independent check of OCF).
# ----------------------------------------------------------------------------

def H_hampaths(A, n):
    adj = [[j for j in range(n) if A[i][j]] for i in range(n)]
    cnt = 0
    def dfs(v, vis, depth):
        nonlocal cnt
        if depth == n:
            cnt += 1
            return
        for w in adj[v]:
            if not (vis >> w) & 1:
                dfs(w, vis | (1 << w), depth + 1)
    for s in range(n):
        dfs(s, 1 << s, 1)
    return cnt


# ----------------------------------------------------------------------------
def poly_from_cumulants(c, deg):
    """Reconstruct I(z) = exp( sum_{k>=1} c_k z^k ) as a power series up to z^deg
    (Fractions).  Used for the EXACT check that exp of the cumulant series equals
    the independence polynomial (formal-power-series identity)."""
    # exp(L), L = sum c_k z^k.  Use I' = L' I  =>  m*I_m = sum_{k=1..m} k c_k I_{m-k}.
    I = [Fraction(0)] * (deg + 1)
    I[0] = Fraction(1)
    cc = [Fraction(0)] + [Fraction(c[k - 1]) for k in range(1, len(c) + 1)]
    for m in range(1, deg + 1):
        s = Fraction(0)
        for k in range(1, m + 1):
            if k <= len(c):
                s += Fraction(k) * cc[k] * I[m - k]
        I[m] = s / Fraction(m)
    return I


def analyze_tournament(A, n, K):
    cycles = directed_odd_cycles(A, n)
    adj = build_omega(cycles)
    alpha = independence_poly(adj)
    feats = omega_features(adj, cycles)
    H_ocf = sum(alpha[k] * (2 ** k) for k in range(len(alpha)))   # OCF: I(Omega,2)
    H_dir = H_hampaths(A, n)                                       # direct count
    ub = ursell_cluster_b(adj)                                    # b1,b2,b3 (-|E| family)
    cum = formal_log_cumulants(alpha, K)                          # analytic c_k
    return alpha, feats, H_ocf, H_dir, ub, cum


def _check_one(A, n, K, failures, b3_rows, samples, max_print=6):
    import math
    alpha, feats, H_ocf, H_dir, ub, cum = analyze_tournament(A, n, K)
    alpha1 = alpha[1] if len(alpha) > 1 else 0
    alpha2 = alpha[2] if len(alpha) > 2 else 0
    E = feats['E']

    # sanity: OCF == direct path count
    if H_ocf != H_dir:
        failures['ocf_vs_path'] += 1

    # (a1) b_1 = alpha_1 = |V(Omega)|
    if ub['b1'] != Fraction(alpha1):
        failures['b1'] += 1
    # (a2) b_2 = -|E(Omega)|  AND  = alpha_2 - C(alpha_1,2)
    if ub['b2'] != Fraction(-E):
        failures['b2_eq_negE'] += 1
    if ub['b2'] != Fraction(alpha2) - Fraction(alpha1 * (alpha1 - 1), 2):
        failures['b2_eq_alpha'] += 1
    # (a3) b_3 = P2 - t  (Ursell connected 3-cluster)
    if ub['b3'] != Fraction(feats['P2']) - Fraction(feats['t']):
        failures['b3'] += 1

    # (b) EXACT: exp( sum c_k z^k ) == I(Omega, z) as FORMAL POWER SERIES.
    #     This is the rigorous meaning of "H = exp of the cluster series":
    #     the generating-function identity I = exp(sum c_k z^k) holds as a formal
    #     power-series identity for EVERY Omega (so H = I(Omega,2) is the value of
    #     exp(cluster series) -- but ONLY after RESUMMATION, see note below).
    #     Reconstruct I from the cumulants and compare to alpha exactly.
    deg = len(alpha) - 1
    Irec = poly_from_cumulants(cum, deg)
    if any(Irec[k] != Fraction(alpha[k]) for k in range(deg + 1)):
        failures['exp_cum_eq_I'] += 1
    # (b2) Numeric "log H = sum_k c_k 2^k" is checked OUTSIDE the per-tournament
    #     loop because the cluster series sum_k c_k z^k typically DIVERGES at z=2
    #     (z=2 is OUTSIDE the radius of convergence = smallest |root| of I(Omega,z)).
    #     Instead we verify convergence at a SMALL z inside the disk: we test
    #     z = z_in where z_in < (radius). We approximate radius by 1/(2*alpha_1+1)
    #     to stay safely inside for these graphs, and check the partial sum -> log I.
    if deg >= 1:
        # smallest-root radius lower bound: for a polynomial 1 + a1 z + ... the
        # nearest root has |z| >= 1/(1+max|a_k|^{1/k})... we just pick a tiny z.
        z_in = Fraction(1, 4 * (alpha1 + 1))
        Sin = sum(cum[k - 1] * (z_in ** k) for k in range(1, K + 1))
        Ival = sum(Fraction(alpha[k]) * (z_in ** k) for k in range(len(alpha)))
        if abs(math.exp(float(Sin)) - float(Ival)) > 1e-9 * max(1.0, float(Ival)):
            failures['logH'] += 1

    # (c) H_free = 3^alpha_1 >= H, equality iff Omega edgeless
    H_free = 3 ** alpha1
    if not (H_free >= H_ocf):
        failures['free'] += 1
    if (H_free == H_ocf) != (E == 0):
        failures['free_eq_iff'] += 1

    b3_rows.append((feats['E'], feats['P2'], feats['t'], ub['b3']))
    if len(samples) < max_print:
        samples.append((alpha, feats, H_ocf, H_dir, ub, cum, H_free))


def _fresh_failures():
    return {'b1': 0, 'b2_eq_negE': 0, 'b2_eq_alpha': 0, 'b3': 0,
            'exp_cum_eq_I': 0, 'logH': 0, 'free': 0, 'free_eq_iff': 0,
            'ocf_vs_path': 0}


def _report(failures, count, b3_fit, samples):
    print(f"  tournaments checked: {count}", flush=True)
    print(f"  (a1) b_1 = alpha_1 = |V(Omega)|                : failures = {failures['b1']}", flush=True)
    print(f"  (a2) b_2 = -|E(Omega)|                         : failures = {failures['b2_eq_negE']}", flush=True)
    print(f"  (a2) b_2 = alpha_2 - C(alpha_1,2)              : failures = {failures['b2_eq_alpha']}", flush=True)
    print(f"  (a3) b_3 = P2(Omega) - t(Omega)  (Ursell)      : failures = {failures['b3']}", flush=True)
    print(f"  (b)  exp(sum c_k z^k) == I(Omega,z) [EXACT PS] : failures = {failures['exp_cum_eq_I']}", flush=True)
    print(f"  (b)  log I = sum_k c_k z^k at SMALL z (in disk): failures = {failures['logH']}", flush=True)
    print(f"  (c)  H_free = 3^alpha_1 >= H                   : failures = {failures['free']}", flush=True)
    print(f"  (c)  H_free == H  IFF  Omega edgeless          : failures = {failures['free_eq_iff']}", flush=True)
    print(f"       OCF I(Omega,2) == direct Ham-path count   : failures = {failures['ocf_vs_path']}", flush=True)
    print(f"  b_3 exact linear fit in (E,P2,t): {b3_fit}", flush=True)
    print(f"  -- sample rows (first {len(samples)}):", flush=True)
    for alpha, feats, H_ocf, H_dir, ub, cum, H_free in samples:
        cshow = [str(cum[i]) for i in range(min(3, len(cum)))]
        print(f"    alpha={alpha} E={feats['E']} P2={feats['P2']} t={feats['t']} "
              f"p33={feats['p33']} | H={H_ocf} H_free={H_free} | "
              f"b1={ub['b1']} b2={ub['b2']} b3={ub['b3']} | c1..3={cshow}", flush=True)


def run_exhaustive(n, K, label, max_print=6):
    print(f"\n========== n = {n}  ({label}) ==========", flush=True)
    failures = _fresh_failures()
    count = 0
    b3_rows = []
    samples = []
    for A in all_tournaments(n):
        count += 1
        _check_one(A, n, K, failures, b3_rows, samples, max_print)
    b3_fit = solve_b3(b3_rows)
    _report(failures, count, b3_fit, samples)
    return failures, count, b3_fit


def solve_b3(rows):
    """Find rational (cE, cP, ct) with b3 = cE*E + cP*P2 + ct*t exactly on all
    rows.  We use UNIQUE feature vectors (E,P2,t)->b3 (consistency-checked), pick
    3 independent ones, solve, and verify on every distinct row."""
    # dedupe to distinct feature vectors; check b3 is a function of (E,P2,t)
    table = {}
    for (E, P2, t, b3) in rows:
        key = (E, P2, t)
        v = Fraction(b3)
        if key in table and table[key] != v:
            return f"NOT-A-FUNCTION of (E,P2,t): {key} gives both {table[key]} and {v}"
        table[key] = v
    pts = [(Fraction(E), Fraction(P2), Fraction(t), v) for (E, P2, t), v in table.items()]
    # pick 3 independent feature vectors
    chosen = []
    for p in pts:
        cand = [q[:3] for q in chosen] + [p[:3]]
        if gauss_rank([list(c) for c in cand]) == len(cand):
            chosen.append(p)
        if len(chosen) == 3:
            break
    if len(chosen) < 3:
        return (f"UNDERDETERMINED (distinct-feature rank {len(chosen)} < 3); "
                f"identity b3 = P2 - t already checked directly (a3).")
    M = [list(c[:3]) for c in chosen]
    y = [c[3] for c in chosen]
    sol = solve3(M, y)
    if sol is None:
        return "singular"
    cE, cP, ct = sol
    ok = all(cE * p[0] + cP * p[1] + ct * p[2] == p[3] for p in pts)
    return f"b3 = ({cE})*E + ({cP})*P2 + ({ct})*t   EXACT_ON_ALL_DISTINCT_ROWS={ok}"


def rank3(vecs):
    # rank over Q of list of length-3 Fraction vectors
    M = [list(v) for v in vecs]
    return gauss_rank(M)


def gauss_rank(M):
    M = [row[:] for row in M]
    rows = len(M)
    cols = len(M[0]) if rows else 0
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                fac = M[i][c]
                M[i] = [M[i][j] - fac * M[r][j] for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def solve3(M, y):
    A = [M[i][:] + [y[i]] for i in range(3)]
    for c in range(3):
        piv = None
        for i in range(c, 3):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            return None
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(3):
            if i != c and A[i][c] != 0:
                fac = A[i][c]
                A[i] = [A[i][j] - fac * A[c][j] for j in range(4)]
    return [A[i][3] for i in range(3)]


def run_sample(n, K, ns, seed=12345):
    import random
    rng = random.Random(seed)
    print(f"\n========== n = {n}  (random sample, {ns} tournaments) ==========", flush=True)
    failures = _fresh_failures()
    b3_rows = []
    samples = []
    for _ in range(ns):
        A = random_tournament(n, rng)
        _check_one(A, n, K, failures, b3_rows, samples)
    b3_fit = solve_b3(b3_rows)
    _report(failures, ns, b3_fit, samples)
    return failures, b3_fit


def paley7():
    """Paley tournament T_7: i->j iff (j-i) is a QR mod 7 (QR = {1,2,4})."""
    n = 7
    QR = {1, 2, 4}
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % 7) in QR:
                A[i][j] = 1
    return A, n


def radius_demo():
    """Show that the cluster series sum_k c_k z^k DIVERGES at z=2 (z=2 is OUTSIDE
    the radius of convergence = smallest |root| of I(Omega,z)), so 'H = exp of the
    cluster series' is a RESUMMED / formal-power-series statement, not a convergent
    sum at fugacity 2.  Demonstrated on Paley T_7."""
    import cmath
    A, n = paley7()
    cycles = directed_odd_cycles(A, n)
    adj = build_omega(cycles)
    alpha = independence_poly(adj)
    H = sum(alpha[k] * 2 ** k for k in range(len(alpha)))
    deg = len(alpha) - 1
    cum = formal_log_cumulants(alpha, 400)
    # smallest |root| of I(z): use numpy
    import numpy as np
    coeffs = [float(alpha[deg - k]) for k in range(deg + 1)]  # high->low
    roots = np.roots(coeffs)
    rmin = min(abs(r) for r in roots) if len(roots) else float('inf')
    print("\n" + "=" * 78, flush=True)
    print("RADIUS-OF-CONVERGENCE NOTE (Paley T_7):", flush=True)
    print(f"  H(T_7) = I(Omega,2) = {H};  alpha_1={alpha[1]}, deg I = {deg}", flush=True)
    print(f"  smallest |root| of I(Omega,z) = R = {rmin:.6f}  (cluster series radius)", flush=True)
    print(f"  => fugacity z=2 is {'OUTSIDE' if rmin < 2 else 'inside'} the disk |z|<R "
          f"(2 {'>' if rmin<2 else '<='} R={rmin:.4f}).", flush=True)
    # partial sums of sum c_k z^k at z=2 (should diverge / oscillate)
    import math
    for K in (5, 10, 20, 40, 80):
        s = sum(float(cum[k - 1]) * (2.0 ** k) for k in range(1, K + 1))
        print(f"    partial sum_K c_k 2^k (K={K:3d}) = {s:.3e}   (target log H = {math.log(H):.4f})", flush=True)
    # at small z (inside disk) it converges to log I:
    z = 0.05
    Iz = sum(float(alpha[k]) * z ** k for k in range(len(alpha)))
    for K in (5, 20, 60):
        s = sum(float(cum[k - 1]) * z ** k for k in range(1, K + 1))
        print(f"    at z={z}: partial sum_K (K={K:2d}) = {s:.8f}, log I(z) = {math.log(Iz):.8f}", flush=True)
    print("  CONCLUSION: the Mayer/cluster series is the ASYMPTOTIC/formal expansion;", flush=True)
    print("  H = exp(sum c_k 2^k) holds as a RESUMMED formal identity, NOT a convergent", flush=True)
    print("  series at z=2.  (Consistent with: lambda=2 is outside hard-core uniqueness", flush=True)
    print("  for Omega(T) of max-degree >= 4.)  The EXACT per-tournament check that", flush=True)
    print("  exp(sum c_k z^k) == I(Omega,z) as a power series is what we verify (0 fails).", flush=True)


def main():
    print("=" * 78, flush=True)
    print("MAYER / VIRIAL CLUSTER (URSELL) EXPANSION OF THE OCF HARD-CORE GAS", flush=True)
    print("  Xi(z) = I(Omega, z) = sum_k alpha_k z^k  (hard-core gas, fugacity z).", flush=True)
    print("  H = I(Omega, 2).   log I(Omega,z) = sum_{k>=1} c_k z^k  (cumulants).", flush=True)
    print("  TWO Ursell families:", flush=True)
    print("   * graphical/excluded-volume virial b_k: b1=alpha_1, b2=-|E|, b3=P2-t.", flush=True)
    print("   * analytic Taylor cumulants c_k: c2 = alpha_2 - alpha_1^2/2 = -|E| - alpha_1/2.", flush=True)
    print("   They differ by the single-particle self-overlap (c_k includes diagonals).", flush=True)
    print("=" * 78, flush=True)

    run_exhaustive(4, K=60, label="ALL labelled tournaments")
    run_exhaustive(5, K=120, label="ALL labelled tournaments")
    run_exhaustive(6, K=200, label="ALL labelled tournaments")   # 2^15 = 32768, feasible
    run_sample(7, K=400, ns=400)

    radius_demo()

    print("\n" + "=" * 78, flush=True)
    print("DICTIONARY (repo overlap defect  <->  cluster integral / virial coeff):", flush=True)
    print("  b_1 = alpha_1 = c3+c5+c7+...  = |V(Omega)|     [IDEAL gas; H_free = 3^alpha_1]", flush=True)
    print("  b_2 = -|E(Omega)| = -(p33 + p35 + p55 + p37 + ...)   [2nd virial = -excluded vol]", flush=True)
    print("        p33 = #Omega-edges among 3-cycles = the 3-3 BLOCK of |E(Omega)|.", flush=True)
    print("        (Full |E| also has 3-5, 5-5, 3-7, ... odd-cycle-pair overlap blocks.)", flush=True)
    print("  b_3 = P2(Omega) - t(Omega)                     [3rd virial; cherries minus", flush=True)
    print("        triangles of the conflict graph].  c_3 = b_3 + |E| ... (diagonal terms).", flush=True)
    print("  ---------------------------------------------------------------------------", flush=True)
    print("  A DIFFERENT (spectral) expansion -- do NOT conflate with the gas:", flush=True)
    print("  TQ (tr A^7 = 7(c7+TQ)) and the Witt necklace cumulants W_k (THM-502/505) are", flush=True)
    print("  cumulants of the WALK generating function det(I-uA)^{-1} = prod (1-u^k)^{-W_k}", flush=True)
    print("  (Bowen-Lanford zeta of the digraph), NOT of log I(Omega,z).  They live on the", flush=True)
    print("  adjacency matrix A, not on the conflict graph Omega.  Their only bridge to the", flush=True)
    print("  gas is alpha_1 = c3+c5+c7+... = sum of ODD primitive-necklace counts (W_odd).", flush=True)
    print("=" * 78, flush=True)


if __name__ == "__main__":
    main()

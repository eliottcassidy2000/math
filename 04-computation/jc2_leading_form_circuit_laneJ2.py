#!/usr/bin/env python3
"""Lane J2 -- the Newton circuit of the JC(2) leading form as a gate.

Context (established upstream this session, restated, not re-derived here):

  L0  For a Jacobian pair (P,Q) in C[x,y] with Jac(P,Q)=const, n=deg P,
      m=deg Q, n+m>=3: the degree n+m-2 part of the Jacobian gives
      Jac(P_n,Q_m)=0, hence P_n = c H^a, Q_m = c' H^b with H a binary form,
      g = deg H = gcd(n,m), a=n/g, b=m/g coprime.
  L1  The degree n+m-3 part gives Jac(P_n,Q_(m-1)) + Jac(P_(n-1),Q_m) = 0.
  L2  (Jung--van der Kulk, classical)  a polynomial AUTOMORPHISM of C^2 with
      deg > 1 has n | m or m | n, i.e. a=1 or b=1.

This file studies H itself.  h(t)=H(1,t) has the DIRECTIONS AT INFINITY as its
roots, so H carries a repo Newton circuit (THM-3000 sec 1 convention):

      H(x,y) = sum_k A_k x^(g-k) y^k,   h_k = A_k / binom(g,k),
      R_k = h_k^2/(h_(k-1) h_(k+1))  (k=1..g-1),   c_k = log(R_k/R_(k-1)).

Everything below is exact (Fraction / sympy).  Nothing here is a reduction of
JC(2) to anything; MISTAKE-237 quarantines that move.  Every statement is a
stratification or a criterion with explicit scope.

PARTS
  A  automorphism census: K = #distinct directions is ALWAYS 1 (H = c L^g).
  B  exact proposition:  R_k = 1 for all k  <=>  K = 1.  (the gate)
  C  THM-3003 rigidity extended to complex roots, and the identity
     R_k(H o swap) = R_(g-k)(H): the P<->Q swap IS the repo reversal.
  D  coordinate-free form of the antipalindromic locus: a PGL_2 involution
     preserving the direction divisor.  K<=2 always has one; K>=3 with
     pairwise-distinct multiplicities never does.
  E  the exact order law at each direction extracted from L1, and its
     corollary rad(H) | P_(n-1) when n > m.
  F  leading-form anatomy of the only known JC counterexample (THM-1300,
     dim 3) and of the repo's homogeneous-Keller shear family.
"""
from fractions import Fraction as Fr
from math import gcd, comb
import itertools
import random
import sys

# sympy perturbs the GLOBAL random stream, so use a private generator
RNG = random.Random(20260731)

# ---------------------------------------------------------------- univariate
# poly = list of Fraction, index = degree, no trailing zeros


def u_trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def u_add(p, q):
    n = max(len(p), len(q))
    r = [Fr(0)] * n
    for i, c in enumerate(p):
        r[i] += c
    for i, c in enumerate(q):
        r[i] += c
    return u_trim(r)


def u_scal(p, s):
    return u_trim([c * s for c in p])


def u_mul(p, q):
    if not p or not q:
        return []
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            r[i + j] += a * b
    return u_trim(r)


def u_divmod(p, q):
    p = list(p)
    assert q, "division by zero polynomial"
    dq = len(q) - 1
    out = [Fr(0)] * max(0, len(p) - dq)
    while len(p) - 1 >= dq and p:
        k = len(p) - 1 - dq
        c = p[-1] / q[-1]
        out[k] = c
        for j, b in enumerate(q):
            p[k + j] -= c * b
        u_trim(p)
    return u_trim(out), u_trim(p)


def u_gcd(p, q):
    p, q = list(p), list(q)
    while q:
        p, q = q, u_divmod(p, q)[1]
    if p:
        p = u_scal(p, Fr(1) / p[-1])
    return p


def u_diff(p):
    return u_trim([p[i] * i for i in range(1, len(p))])


def u_deg(p):
    return len(p) - 1 if p else -1


def sqfree_decomp(f):
    """Yun.  Returns {i: a_i} with f = unit * prod_i a_i^i, a_i squarefree."""
    f = list(f)
    if u_deg(f) <= 0:
        return {}
    a0 = u_gcd(f, u_diff(f))
    b = u_divmod(f, a0)[0]
    c = u_divmod(u_diff(f), a0)[0]
    d = u_add(c, u_scal(u_diff(b), Fr(-1)))
    out = {}
    i = 1
    while u_deg(b) > 0:
        a = u_gcd(b, d)
        if u_deg(a) > 0:
            out[i] = a
        b = u_divmod(b, a)[0]
        c = u_divmod(d, a)[0]
        d = u_add(c, u_scal(u_diff(b), Fr(-1)))
        i += 1
        if i > 200:
            raise RuntimeError("sqfree runaway")
    return out


# ---------------------------------------------------------------- bivariate
# poly = dict {(i,j): Fraction} for x^i y^j


def b_clean(P):
    return {k: v for k, v in P.items() if v != 0}


def b_add(*Ps):
    R = {}
    for P in Ps:
        for k, v in P.items():
            R[k] = R.get(k, Fr(0)) + v
    return b_clean(R)


def b_scal(P, s):
    return b_clean({k: v * s for k, v in P.items()})


def b_mul(P, Q):
    R = {}
    for (i, j), a in P.items():
        for (k, l), b in Q.items():
            key = (i + k, j + l)
            R[key] = R.get(key, Fr(0)) + a * b
    return b_clean(R)


def b_pow(P, e):
    R = {(0, 0): Fr(1)}
    for _ in range(e):
        R = b_mul(R, P)
    return R


def b_dx(P):
    return b_clean({(i - 1, j): a * i for (i, j), a in P.items() if i >= 1})


def b_dy(P):
    return b_clean({(i, j - 1): a * j for (i, j), a in P.items() if j >= 1})


def b_jac(P, Q):
    return b_add(b_mul(b_dx(P), b_dy(Q)), b_scal(b_mul(b_dy(P), b_dx(Q)), Fr(-1)))


def b_deg(P):
    return max((i + j for (i, j) in P), default=-1)


def b_top(P, d=None):
    if d is None:
        d = b_deg(P)
    return b_clean({(i, j): a for (i, j), a in P.items() if i + j == d})


def b_sub(F, G1, G2):
    """F(G1, G2)."""
    di = max((i for (i, j) in F), default=0)
    dj = max((j for (i, j) in F), default=0)
    p1 = [{(0, 0): Fr(1)}]
    for _ in range(di):
        p1.append(b_mul(p1[-1], G1))
    p2 = [{(0, 0): Fr(1)}]
    for _ in range(dj):
        p2.append(b_mul(p2[-1], G2))
    R = {}
    for (i, j), a in F.items():
        R = b_add(R, b_scal(b_mul(p1[i], p2[j]), a))
    return R


def compose(F, G):
    """F o G as pairs."""
    return (b_sub(F[0], G[0], G[1]), b_sub(F[1], G[0], G[1]))


X = {(1, 0): Fr(1)}
Y = {(0, 1): Fr(1)}


# ------------------------------------------------- binary form -> direction data


def form_to_uni(H):
    """H = sum A_k x^(g-k) y^k  ->  (g, [A_0..A_g])."""
    g = b_deg(H)
    A = [Fr(0)] * (g + 1)
    for (i, j), a in H.items():
        assert i + j == g, "not homogeneous"
        A[j] = a
    return g, A


def direction_data(H):
    """Distinct directions of the binary form H in P^1, with multiplicities.

    Returns (K, sorted multiplicity vector, has_point_at_[0:1]).
    Exact: uses squarefree decomposition over Q, so 'distinct' means distinct
    in the algebraic closure.  No factorisation into irreducibles is needed.
    """
    g, A = form_to_uni(H)
    f = u_trim(list(A))                       # h(t) = H(1,t)
    inf_mult = g - u_deg(f)                   # multiplicity of the point [0:1]
    mults = []
    for i, a in sqfree_decomp(f).items():
        mults += [i] * u_deg(a)
    if inf_mult > 0:
        mults.append(inf_mult)
    mults.sort(reverse=True)
    return len(mults), mults, inf_mult > 0


def newton_R(H):
    """R_k of a binary form, k=1..g-1.  None if some A_k = 0 (undefined)."""
    g, A = form_to_uni(H)
    if g < 2 or any(a == 0 for a in A):
        return None
    h = [A[k] / Fr(comb(g, k)) for k in range(g + 1)]
    return [h[k] * h[k] / (h[k - 1] * h[k + 1]) for k in range(1, g)]


def roots_to_form(roots):
    """H = prod (x + r_i y)."""
    H = {(0, 0): Fr(1)} if not isinstance(roots[0], Fr) else {(0, 0): Fr(1)}
    H = {(0, 0): Fr(1)}
    for r in roots:
        H = b_mul(H, {(1, 0): Fr(1), (0, 1): r})
    return H


def R_from_roots(roots):
    """R_k for the multiset {r_i} directly (works over any exact field)."""
    d = len(roots)
    e = [Fr(0)] * (d + 1)
    e[0] = Fr(1)
    for r in roots:
        for k in range(d, 0, -1):
            e[k] = e[k] + e[k - 1] * r
    if any(x == 0 for x in e):
        return None
    h = [e[k] / Fr(comb(d, k)) for k in range(d + 1)]
    return [h[k] * h[k] / (h[k - 1] * h[k + 1]) for k in range(1, d)]


# ======================================================================= PART A
def random_affine():
    """SL2(Z) linear part (det 1) with a translation."""
    M = [[1, 0], [0, 1]]
    for _ in range(RNG.randint(1, 3)):
        k = RNG.choice([-2, -1, 1, 2])
        if RNG.random() < 0.5:
            N = [[1, k], [0, 1]]
        else:
            N = [[1, 0], [k, 1]]
        M = [[M[0][0] * N[0][0] + M[0][1] * N[1][0], M[0][0] * N[0][1] + M[0][1] * N[1][1]],
             [M[1][0] * N[0][0] + M[1][1] * N[1][0], M[1][0] * N[0][1] + M[1][1] * N[1][1]]]
    t1, t2 = RNG.randint(-2, 2), RNG.randint(-2, 2)
    F1 = {(1, 0): Fr(M[0][0]), (0, 1): Fr(M[0][1]), (0, 0): Fr(t1)}
    F2 = {(1, 0): Fr(M[1][0]), (0, 1): Fr(M[1][1]), (0, 0): Fr(t2)}
    return (b_clean(F1), b_clean(F2))


def random_triangular(dmax=3):
    d = RNG.randint(2, dmax)
    p = {(k, 0): Fr(RNG.randint(-3, 3)) for k in range(d)}
    p[(d, 0)] = Fr(RNG.choice([-3, -2, -1, 1, 2, 3]))
    if RNG.random() < 0.5:
        return (X, b_add(Y, b_clean(p)))                       # (x, y + p(x))
    q = {(0, k): v for (k, _), v in [((k, 0), p[(k, 0)]) for k in range(d + 1)]}
    return (b_add(X, b_clean(q)), Y)                            # (x + q(y), y)


def random_automorphism(nfac, dmax=3):
    F = (X, Y)
    for _ in range(nfac):
        F = compose(F, random_affine())
        F = compose(F, random_triangular(dmax))
    F = compose(F, random_affine())
    return F


def leading_form_data(P, Q):
    """L0 data for a Jacobian pair: (n,m,g,a,b,H,ok_L0)."""
    n, m = b_deg(P), b_deg(Q)
    Pn, Qm = b_top(P, n), b_top(Q, m)
    g = gcd(n, m)
    a, b = n // g, m // g
    # Jac(Pn,Qm) = 0 is L0's hypothesis; recover H.
    jac0 = (b_jac(Pn, Qm) == {})
    H = None
    if a == 1:
        c = Pn[max(Pn)] if Pn else None
        # normalise H monic-ish in the lexicographically largest monomial
        key = max(Pn)
        H = b_scal(Pn, Fr(1) / Pn[key])
    elif b == 1:
        key = max(Qm)
        H = b_scal(Qm, Fr(1) / Qm[key])
    return n, m, g, a, b, H, jac0


def partA(N=120, verbose_examples=4):
    print("=" * 74)
    print("PART A -- automorphism census: K = #distinct directions at infinity")
    print("=" * 74)
    bad = []
    rows = []
    seen = 0
    nL1 = 0
    nR = 0
    for trial in range(N):
        F = random_automorphism(RNG.randint(1, 3), dmax=RNG.choice([2, 2, 3]))
        P, Q = F
        n, m = b_deg(P), b_deg(Q)
        if n + m < 3:
            continue
        J = b_jac(P, Q)
        assert J == {(0, 0): Fr(1)} or J == {(0, 0): Fr(-1)} or len(J) == 1, ("Jac not constant", J)
        jc = list(J.values())[0]
        n_, m_, g, a, b, H, jac0 = leading_form_data(P, Q)
        divis = (m % n == 0) or (n % m == 0)                    # L2
        K, mults, at_inf = direction_data(H)
        # verify L0 fully: P_n = c H^a and Q_m = c' H^b
        Pn, Qm = b_top(P, n), b_top(Q, m)
        Ha, Hb = b_pow(H, a), b_pow(H, b)
        cP = Pn[max(Pn)] / Ha[max(Ha)]
        cQ = Qm[max(Qm)] / Hb[max(Hb)]
        okL0 = (b_add(Pn, b_scal(Ha, -cP)) == {}) and (b_add(Qm, b_scal(Hb, -cQ)) == {})
        # verify L1 -- valid only when n+m-3 >= 1 (else that graded piece is the
        # constant Jacobian itself, not zero)
        if n + m >= 4:
            Pn1, Qm1 = b_top(P, n - 1), b_top(Q, m - 1)
            okL1 = (b_add(b_jac(Pn, Qm1), b_jac(Pn1, Qm)) == {})
            nL1 += 1
        else:
            okL1 = True
        # circuit of H in a general-position coordinate frame
        Hgp = b_sub(H, {(1, 0): Fr(1), (0, 1): Fr(1)}, {(1, 0): Fr(1), (0, 1): Fr(2)})
        Rgp = newton_R(Hgp)
        okR = (Rgp is None) or all(r == 1 for r in Rgp)
        if Rgp is not None:
            nR += 1
        seen += 1
        rows.append((n, m, g, a, b, K, tuple(mults)))
        if not (jac0 and divis and okL0 and okL1 and K == 1 and okR):
            bad.append((n, m, g, a, b, K, mults, jac0, divis, okL0, okL1, okR))
    print(f"  automorphisms tested          : {seen}")
    print(f"  Jac(P_n,Q_m)=0 (L0)           : {seen - sum(1 for r in bad if not r[7])}/{seen}")
    print(f"  n|m or m|n     (L2, JvdK)     : {seen - sum(1 for r in bad if not r[8])}/{seen}")
    print(f"  P_n=cH^a, Q_m=c'H^b           : {seen - sum(1 for r in bad if not r[9])}/{seen}")
    print(f"  L1 identity (n+m>=4 only)     : {nL1 - sum(1 for r in bad if not r[10])}/{nL1}")
    print(f"  K = 1  (H = c L^g)            : {seen - sum(1 for r in bad if r[5] != 1)}/{seen}")
    print(f"  R(H) == 1 (general position)  : {nR - sum(1 for r in bad if not r[11])}/{nR}")
    print(f"  VIOLATIONS                    : {len(bad)}")
    degpairs = sorted(set((r[0], r[1]) for r in rows))
    print(f"  distinct degree pairs (n,m)   : {len(degpairs)}")
    print(f"    sample: {degpairs[:14]}")
    print(f"  observed (a,b) multiset       : {sorted(set((r[3], r[4]) for r in rows))}")
    print(f"  observed multiplicity vectors : {sorted(set(r[6] for r in rows))[:10]}")
    print("  --> every automorphism tested has a SINGLE direction at infinity,")
    print("      of full multiplicity g; multiplicity vector = (g).")
    return len(bad) == 0


def partA_proof_check():
    """The induction step used in the K=1 proof, checked mechanically:
    if (P,Q) is an automorphism with n<=m, n|m, then Q'=Q-lam P^(m/n) has
    smaller degree and (P,Q') is again an automorphism with the same P."""
    print()
    print("  induction-step check (Q' = Q - lam P^(m/n) drops the degree):")
    ok = True
    for trial in range(24):
        F = random_automorphism(RNG.randint(1, 2), dmax=3)
        P, Q = F
        n, m = b_deg(P), b_deg(Q)
        if n + m < 3:
            continue
        if n > m:
            P, Q, n, m = Q, P, m, n
        if m % n != 0:
            ok = False
            continue
        k = m // n
        Pn, Qm = b_top(P, n), b_top(Q, m)
        Pk = b_pow(Pn, k)
        lam = Qm[max(Qm)] / Pk[max(Pk)]
        Qp = b_add(Q, b_scal(b_pow(P, k), -lam))
        if b_deg(Qp) >= m:
            ok = False
        # (P,Q') is still a Jacobian pair
        J = b_jac(P, Qp)
        if len(J) != 1 or list(J.keys())[0] != (0, 0):
            ok = False
    print(f"    degree strictly drops and Jac preserved in every trial: {ok}")
    return ok


# ======================================================================= PART B
def partB(N=400):
    print()
    print("=" * 74)
    print("PART B -- the gate:  R_k = 1 for all k  <=>  K = 1  (H = c L^g)")
    print("=" * 74)
    bad = 0
    tested = 0
    hist = {}
    for trial in range(N):
        g = RNG.randint(2, 7)
        # bias towards small root pools so repeats (hence small K) occur often
        pool = [Fr(r) for r in [1, 2, 3, -1, -2, 5, Fr(1, 2), Fr(3, 2), -3, 7]]
        roots = [RNG.choice(pool) for _ in range(g)]
        H = roots_to_form(roots)
        K, mults, _ = direction_data(H)
        R = newton_R(H)
        if R is None:
            continue
        tested += 1
        trivial = all(r == 1 for r in R)
        hist[K] = hist.get(K, 0) + 1
        if trivial != (K == 1):
            bad += 1
    print(f"  random binary forms tested    : {tested}")
    print(f"  (R == all ones)  ==  (K == 1) : {tested - bad}/{tested}   violations: {bad}")
    print(f"  K histogram                   : {dict(sorted(hist.items()))}")
    print("  PROOF (exact, all fields, A_k != 0): R_k=1 for all k <=> h_k geometric")
    print("        <=> A_k = binom(g,k) r^k <=> H = A_0 (x + r y)^g <=> K=1.")
    print("  COROLLARY (with PART A): an automorphism of C^2 of degree > 1 has")
    print("        R(H) == 1 identically, i.e. leading-form circuit c(H) == 0.")
    print("        Contrapositive: c(H) != 0  ==>  (P,Q) is NOT an automorphism.")
    return bad == 0


# ======================================================================= PART C
def partC():
    print()
    print("=" * 74)
    print("PART C -- THM-3003 over C, and: the P<->Q swap IS the repo reversal")
    print("=" * 74)

    # C1. swap identity  R_k(H o sigma) = R_(g-k)(H), exactly, over Q
    ok1 = True
    for trial in range(60):
        g = RNG.randint(3, 7)
        roots = [Fr(RNG.randint(1, 9), RNG.randint(1, 5)) for _ in range(g)]
        R = R_from_roots(roots)
        Rsw = R_from_roots([1 / r for r in roots])
        if R is None or Rsw is None:
            continue
        if Rsw != list(reversed(R)):
            ok1 = False
    print(f"  C1  R_k(H o (x<->y)) = R_(g-k)(H) exactly            : {ok1}")
    print("      (H o sigma reverses the coefficient sequence A_k -> A_(g-k),")
    print("       i.e. r_i -> 1/r_i.  So the involution tau(P,Q)=(Q o sigma,")
    print("       P o sigma) -- a Jacobian pair with (n,m) swapped -- acts on the")
    print("       leading-form circuit as THM-3001's reversal.)")

    # C2. rigidity, complex roots: reciprocal-closed  =>  palindromic R
    import sympy as sp
    ok2 = True
    cases = []
    for trial in range(40):
        j = RNG.randint(1, 3)
        mu = sp.Rational(RNG.randint(1, 6), RNG.randint(1, 4))
        if RNG.random() < 0.5:
            mu = mu * (1 + sp.I) ** 2 / 2      # a genuinely complex mu
        rs = []
        for _ in range(j):
            r = sp.Rational(RNG.randint(1, 5), RNG.randint(1, 4))
            if RNG.random() < 0.5:
                r = r + sp.I * sp.Rational(RNG.randint(1, 3), 2)
            rs += [r, mu / r]
        if RNG.random() < 0.5:
            rs += [sp.sqrt(mu)]                 # a fixed point of r -> mu/r
        d = len(rs)
        e = [sp.Integer(0)] * (d + 1)
        e[0] = sp.Integer(1)
        for r in rs:
            for k in range(d, 0, -1):
                e[k] = sp.expand(e[k] + e[k - 1] * r)
        if any(sp.simplify(x) == 0 for x in e):
            continue
        h = [sp.simplify(e[k] / sp.binomial(d, k)) for k in range(d + 1)]
        R = [sp.simplify(h[k] ** 2 / (h[k - 1] * h[k + 1])) for k in range(1, d)]
        pal = all(sp.simplify(R[k - 1] - R[d - k - 1]) == 0 for k in range(1, d))
        cases.append((d, pal))
        if not pal:
            ok2 = False
    print(f"  C2  reciprocal-closed (complex mu) => R palindromic  : {ok2}"
          f"  ({len(cases)} multisets, degrees {sorted(set(c[0] for c in cases))})")

    # C3. converse over an exhaustive rational pool
    pool = [Fr(1), Fr(2), Fr(3), Fr(1, 2), Fr(1, 3), Fr(-1), Fr(-2), Fr(3, 2), Fr(2, 3)]
    ok3 = True
    npal = 0
    ntot = 0
    for g in (4, 5):
        for roots in itertools.combinations_with_replacement(pool, g):
            R = R_from_roots(list(roots))
            if R is None:
                continue
            ntot += 1
            if R == list(reversed(R)):
                npal += 1
                # must be reciprocal-closed: {r} = {mu/r} with mu^g = e_g^2
                e_g = Fr(1)
                for r in roots:
                    e_g *= r
                # mu is determined up to a g-th root of unity; over Q test all
                # rational mu with mu^g = e_g^2 by checking the multiset identity
                found = False
                for mu in set([e_g ** 2] if g == 1 else []):
                    pass
                # direct search: mu must be r_i * r_j for some pair (the
                # involution pairs roots), so test those finitely many values
                for r1 in roots:
                    for r2 in roots:
                        mu = r1 * r2
                        if mu == 0:
                            continue
                        if sorted([mu / r for r in roots]) == sorted(roots):
                            found = True
                            break
                    if found:
                        break
                if not found:
                    ok3 = False
                    print("    COUNTEREXAMPLE to converse:", roots)
    print(f"  C3  palindromic => reciprocal-closed, exhaustive     : {ok3}"
          f"  ({npal} palindromic of {ntot} multisets, g=4,5)")
    return ok1 and ok2 and ok3


# ======================================================================= PART D
def mobius_from_three(src, dst):
    """Unique Mobius map sending src[i] -> dst[i] (3 distinct points of P^1,
    given as Fractions or the string 'inf').  Returns a 2x2 matrix or None."""
    import sympy as sp

    def cross_matrix(p):
        # matrix sending p0,p1,p2 -> 0,1,inf
        a, b, c = p
        INF = 'inf'
        if a == INF:
            return sp.Matrix([[0, b - c], [1, -c]])
        if b == INF:
            return sp.Matrix([[1, -a], [1, -c]])
        if c == INF:
            return sp.Matrix([[1, -a], [0, b - a]])
        return sp.Matrix([[b - c, -a * (b - c)], [b - a, -c * (b - a)]])

    A = cross_matrix(src)
    B = cross_matrix(dst)
    if A.det() == 0 or B.det() == 0:
        return None
    return sp.simplify(B.inv() * A)


def apply_mobius(M, p):
    import sympy as sp
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    if p == 'inf':
        return 'inf' if c == 0 else sp.simplify(a / c)
    den = c * p + d
    if sp.simplify(den) == 0:
        return 'inf'
    return sp.simplify((a * p + b) / den)


def has_involutive_symmetry(points, mults):
    """Does some NON-IDENTITY involution of P^1 preserve the divisor
    sum m_i [p_i]?   Exact.  points: list of Fraction / 'inf'."""
    import sympy as sp
    K = len(points)
    idx = range(K)
    if K == 1:
        return True, "K=1: e.g. any involution fixing the single point"
    if K == 2:
        # the unique involution fixing both p,q: conjugate t -> -t by the map
        # sending p,q -> 0,inf.  Build it exactly and verify.
        p, q = points
        if p == 'inf':
            M = sp.Matrix([[-1, 2 * q], [0, 1]])          # t -> 2q - t : fixes q, inf
        elif q == 'inf':
            M = sp.Matrix([[-1, 2 * p], [0, 1]])
        else:
            M = sp.Matrix([[p + q, -2 * p * q], [2, -(p + q)]])   # fixes p and q
        good = (apply_mobius(M, p) == p and apply_mobius(M, q) == q
                and sp.simplify(M[0, 0] + M[1, 1]) == 0 and sp.simplify(M.det()) != 0)
        return good, f"K=2: involution fixing both points, matrix {list(M)}"
    for perm in itertools.permutations(idx):
        if any(perm[perm[i]] != i for i in idx):
            continue                       # must be an involution as a permutation
        if any(mults[perm[i]] != mults[i] for i in idx):
            continue                       # must preserve multiplicities
        M = mobius_from_three([points[i] for i in (0, 1, 2)],
                              [points[perm[i]] for i in (0, 1, 2)])
        if M is None:
            continue
        if sp.simplify(M[0, 0] + M[1, 1]) != 0:
            continue                       # trace 0 <=> involution (non-identity)
        if all(apply_mobius(M, points[i]) == points[perm[i]] for i in idx):
            return True, f"permutation {perm}, matrix {list(M)}"
    return False, "none"


def partD():
    print()
    print("=" * 74)
    print("PART D -- the antipalindromic locus, coordinate-free")
    print("=" * 74)
    print("  THM-3003 (extended in PART C): the circuit of H is antipalindromic")
    print("  <=> the direction multiset is invariant under t -> mu/t.  Every")
    print("  non-identity involution of P^1 is PGL_2-conjugate to such a map, and")
    print("  a linear change of (x,y) acts on directions by the full PGL_2.  So:")
    print("    'SOME coordinate system makes the leading-form circuit")
    print("     antipalindromic'  <=>  'the direction divisor has a nontrivial")
    print("     involutive PGL_2 symmetry'.")
    print()
    tests = [
        ("K=2, mults (1,3)", [Fr(2), Fr(5)], [1, 3]),
        ("K=2, mults (2,2)", [Fr(2), Fr(5)], [2, 2]),
        ("K=3, mults (1,1,3)", [Fr(0), 'inf', Fr(1)], [1, 1, 3]),
        ("K=3, mults (1,2,3)", [Fr(0), 'inf', Fr(1)], [1, 2, 3]),
        ("K=4, mults (1,1,2,2) harmonic", [Fr(0), 'inf', Fr(1), Fr(-1)], [1, 1, 2, 2]),
        ("K=4, mults (1,2,3,4)", [Fr(0), 'inf', Fr(1), Fr(2)], [1, 2, 3, 4]),
        ("K=5, mults (1,1,1,1,1) sym", [Fr(0), 'inf', Fr(1), Fr(-1), Fr(2)], [1] * 5),
    ]
    for name, pts, mm in tests:
        ok, why = has_involutive_symmetry(pts, mm)
        print(f"  {name:38s} involutive symmetry: {str(ok):5s}  {why[:44]}")
    print()
    print("  PROVED: K>=3 with pairwise-DISTINCT multiplicities has NO nontrivial")
    print("  involutive symmetry (the only multiplicity-preserving involutive")
    print("  permutation is the identity, and a Mobius map fixing 3 points is the")
    print("  identity).  Hence for such H the leading-form circuit is")
    print("  antipalindromic in NO coordinate system whatsoever.")
    print()
    # direct exact confirmation: K=3 multiplicities 1,2,3, many random SL2 changes
    roots = [Fr(1)] + [Fr(2)] * 2 + [Fr(5)] * 3      # K=3, mults (1,2,3), g=6
    H0 = roots_to_form(roots)
    K, mults, _ = direction_data(H0)
    npal = 0
    ntried = 0
    for _ in range(150):
        M = [[1, 0], [0, 1]]
        for _ in range(RNG.randint(1, 4)):
            k = RNG.choice([-2, -1, 1, 2])
            N = [[1, k], [0, 1]] if RNG.random() < 0.5 else [[1, 0], [k, 1]]
            M = [[M[0][0] * N[0][0] + M[0][1] * N[1][0], M[0][0] * N[0][1] + M[0][1] * N[1][1]],
                 [M[1][0] * N[0][0] + M[1][1] * N[1][0], M[1][0] * N[0][1] + M[1][1] * N[1][1]]]
        Hm = b_sub(H0, {(1, 0): Fr(M[0][0]), (0, 1): Fr(M[0][1])},
                       {(1, 0): Fr(M[1][0]), (0, 1): Fr(M[1][1])})
        R = newton_R(Hm)
        if R is None:
            continue
        ntried += 1
        if R == list(reversed(R)):
            npal += 1
    print(f"  exact control: H with K=3, mults {mults}, {ntried} random SL2(Z)")
    print(f"  coordinate changes -> palindromic circuits found: {npal}  (expected 0)")
    # and the mirror control: mults (1,1,3) IS symmetric
    roots2 = [Fr(2), Fr(1, 2)] + [Fr(1)] * 3          # {r}={1/r}: reciprocal-closed
    R2 = R_from_roots(roots2)
    print(f"  mirror control: roots {[str(r) for r in roots2]} (mults 1,1,3, K=3)")
    print(f"    R = {[str(r) for r in R2]}")
    print(f"    palindromic: {R2 == list(reversed(R2))}   (predicted True)")

    # how big is the symmetric locus?  census over reduced (multiplicity-free)
    # configurations of K points, K = 3..6
    print()
    print("  census: fraction of REDUCED K-point direction divisors (all m_i=1)")
    print("  admitting a nontrivial involutive PGL_2 symmetry")
    print("  (= fraction for which SOME coordinate frame gives an antipalindromic")
    print("   leading-form circuit):")
    pool = [Fr(k) for k in range(-6, 7)] + [Fr(1, 2), Fr(3, 2), Fr(1, 3), Fr(5, 2)]
    for K in (3, 4, 5, 6):
        hit = 0
        tot = 0
        for _ in range(60):
            pts = RNG.sample(pool, K)
            ok, _why = has_involutive_symmetry(pts, [1] * K)
            tot += 1
            hit += int(ok)
        dimM = K - 3
        print(f"    K={K}: {hit:3d}/{tot} symmetric   (moduli dim {dimM};"
              f" symmetric locus dim ~ {max(0, (K // 2) - 1)})")
    return npal == 0 and R2 == list(reversed(R2))


# ======================================================================= PART E
def order_at(P, L):
    """Order of vanishing of the form P along the linear form L (both binary,
    L = (l0,l1) meaning l0*x + l1*y).  Exact, by repeated division."""
    Lp = {(1, 0): L[0], (0, 1): L[1]}
    k = 0
    cur = dict(P)
    while cur:
        # divide cur by Lp if possible: substitute the line and test
        q = form_divide(cur, Lp)
        if q is None:
            break
        cur = q
        k += 1
        if k > 200:
            break
    return k


def form_divide(P, L):
    """Exact division of a binary form by a linear form, or None."""
    d = b_deg(P)
    if d < 1:
        return None
    # write in the univariate variable: use l0*x + l1*y ; if l0 != 0 divide as
    # polynomials in x with y as parameter
    l0, l1 = L[(1, 0)], L[(0, 1)]
    # convert P to coefficient list in y (degree j) : P = sum p_j x^(d-j) y^j
    p = [Fr(0)] * (d + 1)
    for (i, j), a in P.items():
        p[j] = a
    # L = l0 x + l1 y -> as univariate in t=y/x : l0 + l1 t
    if l1 != 0:
        # divide sum p_j t^j by (l0 + l1 t) with a possible root at t=inf
        # handle by synthetic division at t = -l0/l1
        root = -l0 / l1
        # p as poly in t of degree <= d
        q = [Fr(0)] * d
        rem = Fr(0)
        acc = Fr(0)
        # Horner from top
        coeffs = list(reversed(p))          # highest degree first (degree d)
        acc = coeffs[0]
        out = [acc]
        for c in coeffs[1:]:
            acc = acc * root + c
            out.append(acc)
        rem = out[-1]
        if rem != 0:
            return None
        qc = list(reversed(out[:-1]))       # quotient in t, degree d-1
        qc = [c / l1 for c in qc]
        return b_clean({(d - 1 - j, j): qc[j] for j in range(d)})
    else:
        # L = l0 x : divide by x -> every monomial needs i >= 1
        if any(i == 0 for (i, j) in P):
            return None
        return b_clean({(i - 1, j): a / l0 for (i, j), a in P.items()})


def partE(N=200):
    print()
    print("=" * 74)
    print("PART E -- the exact ORDER LAW at each direction, extracted from L1")
    print("=" * 74)
    print("  Write H = prod_i L_i^(m_i), e_i = ord_(L_i) P_(n-1), f_i = ord_(L_i) Q_(m-1).")
    print("  Lemma (exact):  ord_(L_i) Jac(F, H) = ord_(L_i)F + m_i - 1 + (extra),")
    print("  with extra >= 1 exactly when  ord_(L_i)(F) * g = m_i * deg F.")
    print("  Feeding this into L1 (Jac(P_n,Q_(m-1)) + Jac(P_(n-1),Q_m) = 0) gives")
    print()
    print("        e_i - f_i  =  m_i (a - b)  +  eps_i - delta_i            (J2-LAW)")
    print()
    print("  with delta_i = 0 unless e_i g = m_i (n-1), eps_i = 0 unless f_i g =")
    print("  m_i (m-1).  COROLLARY (n > m, both Jacobians nonzero): e_i >= 1 for")
    print("  every i, i.e. rad(H) divides P_(n-1) -- every direction at infinity")
    print("  reappears in the SUB-leading form of the higher-degree component.")
    print()
    ok = 0
    tot = 0
    nonres = 0
    nonres_ok = 0
    lemma_ok = 0
    lemma_tot = 0
    shown = 0
    for trial in range(N):
        F = random_automorphism(RNG.randint(1, 2), dmax=3)
        P, Q = F
        n, m = b_deg(P), b_deg(Q)
        if n < 2 or m < 2:
            continue
        n_, m_, g, a, b, H, jac0 = leading_form_data(P, Q)
        Pn1, Qm1 = b_top(P, n - 1), b_top(Q, m - 1)
        if not Pn1 or not Qm1:
            continue
        # K=1 for automorphisms: the single direction is the linear factor of H
        # recover L from H = c L^g by taking the (g-1)-st derivative direction
        Hx, Hy = b_dx(H), b_dy(H)
        # L is proportional to (H_y, -H_x)? no: H=cL^g -> H_x = c g l0 L^(g-1),
        # H_y = c g l1 L^(g-1); so (l0:l1) = (leading coeff ratio)
        # take L = l0 x + l1 y from the ratio of the top monomials
        g_, A = form_to_uni(H)
        # H = A_0 (x + r y)^g  =>  r = A_1/(g A_0)  (A_0 != 0 case)
        if A[0] != 0:
            r = A[1] / (g_ * A[0])
            L = (Fr(1), r)
        else:
            L = (Fr(0), Fr(1))
        e = order_at(Pn1, L)
        f = order_at(Qm1, L)
        mi = g
        JP = b_jac(Pn1, H)
        JQ = b_jac(H, Qm1)
        if JP == {} or JQ == {}:
            continue                          # degenerate branch of L1
        delta = order_at(JP, L) - (e + mi - 1)
        eps = order_at(JQ, L) - (f + mi - 1)
        # the LEMMA: extra >= 1 exactly when ord*g = m_i*deg
        lemma_tot += 2
        lemma_ok += int((delta >= 1) == (e * g == mi * (n - 1)))
        lemma_ok += int((eps >= 1) == (f * g == mi * (m - 1)))
        lhs = e - f
        rhs = mi * (a - b) + eps - delta
        tot += 1
        good = (lhs == rhs)
        if good:
            ok += 1
        resonant = (delta > 0) or (eps > 0)
        if not resonant:
            nonres += 1
            nonres_ok += int(lhs == mi * (a - b))
        if shown < 7:
            tag = "  (resonant)" if resonant else ""
            print(f"   (n,m)=({n:3d},{m:3d}) a,b=({a},{b}) g={g:2d}  e={e:2d} f={f:2d}"
                  f"  d={delta} eps={eps}  e-f={lhs:3d}  law_rhs={rhs:3d}"
                  f"  match={good}{tag}")
            shown += 1
    print(f"  (J2-LAW) exact, ALL instances                      : {ok}/{tot}")
    print(f"  non-resonant instances with e-f = m_i(a-b)         : {nonres_ok}/{nonres}")
    print(f"  order lemma (extra>=1 <=> ord*g = m_i*deg)         : {lemma_ok}/{lemma_tot}")

    # The order lemma itself is NOT special to K=1 -- it is a statement about
    # arbitrary binary forms.  Verify it directly on forms with K >= 2, which is
    # the regime a JC(2) counterexample would live in.
    print()
    print("  order lemma on forms with K >= 2 (the counterexample regime):")
    lok = 0
    ltot = 0
    kseen = set()
    for _ in range(300):
        pool = [Fr(1), Fr(2), Fr(3), Fr(-1), Fr(1, 2), Fr(5)]
        g_ = RNG.randint(3, 6)
        rts = [RNG.choice(pool) for _ in range(g_)]
        H = roots_to_form(rts)
        K, mults, _ = direction_data(H)
        if K < 2:
            continue
        dF = RNG.randint(2, 6)
        F = roots_to_form([RNG.choice(pool) for _ in range(dF)])
        JJ = b_jac(F, H)
        if JJ == {}:
            continue
        for r in set(rts):
            L = (Fr(1), r)
            mi = rts.count(r)
            eF = order_at(F, L)
            extra = order_at(JJ, L) - (eF + mi - 1)
            pred = (eF * g_ == mi * dF)
            ltot += 1
            lok += int((extra >= 1) == pred)
            kseen.add(K)
    print(f"    exact instances (K in {sorted(kseen)}) : {lok}/{ltot}")
    return (ok == tot and nonres_ok == nonres and lemma_ok == lemma_tot
            and lok == ltot)


# ======================================================================= PART F
def partF():
    print()
    print("=" * 74)
    print("PART F -- leading-form anatomy of actual repo objects")
    print("=" * 74)
    import sympy as sp
    x, y, z = sp.symbols('x y z')

    # F1: the only known JC counterexample (THM-1300, dim 3)
    u = 1 + x * y
    F1 = sp.expand(u ** 3 * z + y ** 2 * u * (4 + 3 * x * y))
    F2 = sp.expand(y + 3 * x * u ** 2 * z + 3 * x * y ** 2 * (4 + 3 * x * y))
    F3 = sp.expand(2 * x - 3 * x ** 2 * y - x ** 3 * z)
    J = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in (F1, F2, F3)])
    det = sp.expand(J.det())
    print(f"  THM-1300 map: det JF = {det}   (expected -2)")

    def top_form(F):
        P = sp.Poly(F, x, y, z)
        d = P.total_degree()
        return sp.expand(sum(c * x ** m[0] * y ** m[1] * z ** m[2]
                             for m, c in zip(P.monoms(), P.coeffs())
                             if sum(m) == d)), d

    T1, d1 = top_form(F1)
    T2, d2 = top_form(F2)
    T3, d3 = top_form(F3)
    print(f"    deg F1 = {d1}, top form = {T1}")
    print(f"    deg F2 = {d2}, top form = {T2}")
    print(f"    deg F3 = {d3}, top form = {T3}")
    Jt = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in (T1, T2, T3)])
    print(f"    det of the top-form Jacobian = {sp.expand(Jt.det())}   (L0 analogue: 0)")
    G = sp.gcd(sp.gcd(T1, T2), T3)
    print(f"    gcd of the three top forms   = {sp.factor(G)}")
    print(f"    cofactors                    = "
          f"{sp.simplify(T1/G)}, {sp.simplify(T2/G)}, {sp.simplify(T3/G)}")
    print("    --> all three leading forms are C[y]-multiples of the SINGLE form")
    print("        x^3 z: the dim-3 analogue of P_n = cH^a, Q_m = c'H^b.")
    print("    --> the direction divisor at infinity (the gcd) is  x^3 * z^1  :")
    print("        K = 2 distinct components, multiplicity vector (3,1).")
    # circuit of the generic-line restriction of the direction divisor x^3 z:
    # a line meets {x=0} in one point (multiplicity 3) and {z=0} in one point
    # (multiplicity 1), so the binary form is (t-p)^3 (t-q) -- a two-cluster
    # profile.  Scan the ratio q/p; the alternation is the invariant.
    print("    generic-line restriction = two-cluster profile (mults 3,1);")
    print("    circuit sign words over a scan of the ratio q/p:")
    words = {}
    for num, den in [(2, 1), (3, 1), (5, 1), (10, 1), (100, 1), (1, 2), (1, 5),
                     (1, 100), (7, 3), (-2, 1), (-1, 3)]:
        roots = [Fr(1)] * 3 + [Fr(num, den)]
        R = R_from_roots(roots)
        if R is None:
            continue
        sgn = ''.join('+' if R[k] > R[k - 1] else ('-' if R[k] < R[k - 1] else '0')
                      for k in range(1, len(R)))
        ch = sum(1 for i in range(1, len(sgn))
                 if sgn[i] != sgn[i - 1] and '0' not in (sgn[i], sgn[i - 1]))
        words[str(Fr(num, den))] = (sgn, ch)
    for k, v in words.items():
        print(f"      q/p = {k:>7s} : word {v[0]:>4s}  sign changes {v[1]}")
    mx = max(v[1] for v in words.values())
    print(f"    observed max alternation = {mx};  THM-3004 cap 2K-3 = {2*2-3}")

    # F2: the repo's homogeneous-Keller family (jc2_resonance_dictionary_...)
    print()
    print("  homogeneous Keller family (jc2_resonance_dictionary_leadingform_S107):")
    print("    F = I + H, H homogeneous of degree d, Keller <=> H = c(bx-ay)^d (a,b)")
    for (aa, bb, cc, dd) in [(1, 2, 1, 3), (2, -1, 3, 4), (1, 1, 1, 5)]:
        ell = {(1, 0): Fr(bb), (0, 1): Fr(-aa)}
        Hd = b_pow(ell, dd)
        P = b_add(X, b_scal(Hd, Fr(cc * aa)))
        Q = b_add(Y, b_scal(Hd, Fr(cc * bb)))
        J = b_jac(P, Q)
        n, m = b_deg(P), b_deg(Q)
        n_, m_, g, a, b, H, jac0 = leading_form_data(P, Q)
        K, mults, _ = direction_data(H)
        print(f"    (a,b,c,d)=({aa},{bb},{cc},{dd}): Jac={dict(J)}  (n,m)=({n},{m})"
              f"  g={g} (a,b)=({a},{b})  K={K}  mults={mults}")
    print("    --> K=1 throughout, as PART A predicts (shears are automorphisms).")

    # F3: the degree lattice of the repo's own JC(2) frontier
    print()
    print("  repo JC(2) degree-pair lattice (jc2_degree_lattice_hyp9070 / THM-2823):")
    print("    a counterexample needs a,b >= 2 coprime (L2), so g = gcd(n,m) >= 2")
    print("    and deg H = g >= 2, hence 1 <= K <= g.  Automorphisms sit in K=1.")
    for (n, m) in [(178, 288), (16, 24), (26, 39), (55, 89), (64, 96)]:
        g = gcd(n, m)
        a, b = n // g, m // g
        cfq = []
        A, B = a, b
        while B:
            cfq.append(A // B)
            A, B = B, A % B
        print(f"    (n,m)=({n},{m}): g={g:3d} (a,b)=({a:3d},{b:3d}) coprime={gcd(a,b)==1}"
              f"  cf(a/b)={cfq}  1<=K<={g}")
    return True


# ======================================================================= main
def main():
    okA = partA()
    okA2 = partA_proof_check()
    okB = partB()
    okC = partC()
    okD = partD()
    okE = partE()
    okF = partF()
    print()
    print("=" * 74)
    print("SUMMARY")
    print("=" * 74)
    for name, v in [("A automorphism census K=1", okA),
                    ("A induction step", okA2),
                    ("B R==1 <=> K==1", okB),
                    ("C THM-3003 over C + swap=reversal", okC),
                    ("D involutive-symmetry criterion", okD),
                    ("E order law from L1", okE),
                    ("F repo object anatomy", okF)]:
        print(f"  {name:38s} {v}")
    print()
    print("  ALL:", all([okA, okA2, okB, okC, okD, okE, okF]))


if __name__ == "__main__":
    main()

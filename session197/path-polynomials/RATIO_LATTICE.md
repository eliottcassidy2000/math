# Constructive cancellation through prime-valuation lattices

Session197. **PROVED exact all-order inverse criterion for a specified module language.** The problem solved here is cancellation of the binary exponential mode, not production of an arbitrary prime count or a specified ternary amplitude.

## 1. Physical inputs and multiplicative responses

Fix the seven-vertex quotient `[14,32,2,22,103,13,47]` from `UNIT_OFFSET_RIGIDITY.md`, with nonempty modules `(T,P,U,K1,R,K1,V)`. Fix the physical leaf tournaments P,R. Let a finite nonempty bank consist of nonempty physical tournaments G_1,...,G_k. Allow U and V to be any nonempty ordered joins of copies of these bank graphs, in a fixed chosen order.

For every nonempty graph X write

    h_X=H(X), b_X=B_X(2), rho_X=b_X/h_X.

These are positive integers h_X,b_X and a positive rational rho_X. The two choices of an occupied bin give b_X>=2h_X, so b_X-h_X>0.

**Multiplicativity.** If A->B is an ordered join, then

    B_(A->B)(d)=B_A(d)B_B(d),
    h_(A->B)=h_A h_B, b_(A->B)=b_A b_B,
    rho_(A->B)=rho_A rho_B.                                (1)

In each labelled bin, a directed path can visit A only before B. Restricting the path to the two parts gives two possibly empty directed paths. Conversely concatenate those two paths, using the universal A-to-B arc if both are nonempty. This is a bijection independently in every labelled bin, proving (1). It does not assert graph commutativity; the construction chooses an actual ordered join.

## 2. Exact cancellation, including zero cases

Set

    d_P=b_P-h_P>0, d_R=b_R-h_R>0,
    e_P=4h_P-b_P, e_R=4h_R-b_R.

With T=Tr_m, the audited arbitrary-module theorem gives

    H_m=A 3^m+B 2^m+C, A>0,

and

    B=b_V d_R e_P h_U+b_U d_P e_R h_V.

Dividing only by the strictly positive h_U h_V gives

    B/(h_U h_V)=d_R e_P rho_V+d_P e_R rho_U.                (2)

Therefore:

- If e_P=e_R=0, every allowed U,V cancels B.
- If exactly one e is zero, cancellation is impossible.
- If both e values are nonzero with the same sign, cancellation is impossible.
- If e_P e_R<0, cancellation is equivalent to

      rho_U/rho_V=kappa,
      kappa=-(d_R e_P)/(d_P e_R)>0.                        (3)

No zero or negative quantity is sent to a prime-valuation calculation. Since A>0, distinct exponential bases show that cancellation B=0 is exactly the condition for an eventual pure ternary ray A3^m+C; it then holds for every positive m. The constant C is not required to be a unit or a prime here.

## 3. A finite integer-lattice decision

Assume the opposite-sign case. Factor the positive reduced rationals rho_1,...,rho_k and kappa. Let p_1,...,p_s be every prime appearing in any numerator or denominator, including primes appearing only in kappa. Define the integer matrix and target

    M_(a,i)=v_(p_a)(rho_i), t_a=v_(p_a)(kappa).

Here denominator valuations are negative. If U uses u_i copies and V uses v_i copies of G_i, then (1) and unique prime factorization give

    rho_U/rho_V=kappa  iff  M(u-v)=t.                      (4)

Thus an allowed physical cancellation exists if and only if the integer system

    Mz=t, z in Z^k                                       (5)

is solvable. This is lattice membership, not rational-span membership: integrality is essential.

For any solution z, put u_i=max(z_i,0) and v_i=max(-z_i,0). If either sum is zero, add one copy of any fixed bank graph to **both** U and V. This common padding leaves their ratio unchanged and makes both modules nonempty. Construct U,V as their specified ordered joins. Equation(2) now proves an actual pure ternary family in the fixed quotient. There is no virtual-state realization step.

The quotient has a directed spanning cycle, so every resulting graph is strong, regardless of whether an individual bank graph is weak. Common padding changes A and C and is not claimed to preserve a prescribed count or amplitude.

## 4. Smith normal form gives witnesses or exact obstructions

Compute integer unimodular matrices E,F with

    D=E M F=diag(d_1,...,d_r,0,...),
    d_i>0 and d_i divides d_(i+1).

Put t'=E t. Then (5) is solvable precisely when

    d_i divides t'_i for i<=r,
    t'_i=0 for all remaining rows.

If these conditions hold, choose y_i=t'_i/d_i for i<=r and every free coordinate zero. The integer vector z=F y solves (5), and Section3 turns it directly into graph modules.

If a divisibility condition fails, the corresponding row of E gives an integer linear functional whose value on every generator column is divisible by d_i but whose value on t is not. If a zero-row condition fails, that row annihilates every generator column but not t. Either is an exact all-order certificate that no number of ordered-join copies can cancel B. The checker validates each displayed obstruction directly, rather than relying only on a solver's return status.

The included standard-library routine computes E,D,F by integer row and column operations. At each pivot it uses Euclidean reduction to clear that row and column. If the pivot does not divide a remaining entry, adding that entry's row to the pivot row exposes a nondivisible remainder, and the next Euclidean step strictly decreases the pivot's absolute value. Hence the pivot stage terminates. After it divides the entire lower block, later integer operations preserve that divisibility, giving the Smith divisibility chain. Swaps, integer row/column additions, and sign changes are all unimodular and are recorded in E,F. Thus the algorithm is finite and constructive, not a bounded search for exponents.

## 5. Physical positive and impossible examples

Let K=K1 and let W have rows `[8,1,3,6]`. Then

    (h_K,b_K,rho_K)=(1,2,2),
    (h_W,b_W,rho_W)=(5,32,32/5).

### Positive witness

Fix P=K1,R=Tr4. Then d_P=1,e_P=2,d_R=15,e_R=-12 and kappa=5/2. For the bank(K,W), in the prime order(2,5),

    M=[[1,5],[0,-1]], t=[-1,1].

The integer solution z=(4,-1) gives U=Tr4 and V=W. Substituting these actual graphs gives

    (A,B,C)=(15360,0,-65),
    H_m=15360*3^m-65.

The checker constructs the source and next clone graph and independently counts46015 and138175 Hamilton paths. This recovers the earlier physical mixed-cancellation example by an inverse arithmetic calculation, rather than discovering it through a finite profile census.

### Missing-prime obstruction

Keep the bank(K,W), but take P=K1,R=Tr3. Now kappa=7/2. Every bank ratio has zero7-adic valuation, while kappa has valuation one. Thus no allowed U,V, at any order, can cancel the binary mode. A finite search cutoff is unnecessary.

### An integrality obstruction despite full prime support

Return to P=K1,R=Tr4 and kappa=5/2, but use the bank(Tr3,W). Its matrix is

    M=[[3,5],[0,-1]], t=[-1,1].

Its rational span contains t, but its integer lattice does not. Every generated valuation vector satisfies

    v_2+5v_5=0 modulo3,

whereas the target gives -1+5=4, not divisible by3. Equivalently z_2=-1 would require3z_1=4. This is a genuine Smith divisibility obstruction, not a missing-prime-support test.

The zero cases are physical too: P=R=Tr2 cancels for every U,V; P=R=K1 never cancels; P=Tr2,R=K1 has exactly one zero e and never cancels.

## 6. Scope and reproduction

The algorithm decides pure ternary cancellation for fixed physical P,R and a finite bank of allowed ordered-join factors for U,V. It does not decide whether the resulting A,C match a separately prescribed pair, nor whether any output is prime or7p. In the positive example every count is a composite multiple of5. The unit-offset and prime-offset rigidity theorems remain separate restrictions.

Run `python3 verify_ratio_lattice.py` and `python3 -O verify_ratio_lattice.py`. The checker verifies all physical bank profiles, Smith certificates, positive graph recipes, two different all-order impossibility certificates, the nonempty common-padding rule, and the zero-sign branches. It stress-tests the general Smith routine on exact random rectangular matrices and known integer targets, checking E M F=D, unimodular determinants, the divisibility chain, and every returned witness or obstruction. No floating-point linear algebra is used.

The inherited physical coefficient formula is `UNIT_OFFSET_RIGIDITY.md`, frozen SHA256 `3c116086bfd697eaeaf5fcbf46308f1e31cdf653eb2a651053f01bce46ce374b`. Its transitive quotient source and full-bin linearity dependencies are traced there. The new result here is the exact multiplicative-response/lattice equivalence and its constructive all-order certificates for this module language.

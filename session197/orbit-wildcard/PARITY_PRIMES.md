# Which primes can be sensitive to module-order parity?

Separate refinement of the frozen `PROOF.md`, Session197. No change to that snapshot. Its OCF resolution, exact parity-context gcd theorem and inherited maximal-run substitution are the inputs.

For a fixed rooted context let a_1,...,a_D be its nonnegative run coefficients, with a_1,a_D>0. Let G_odd,G_even be its exact universal divisors over all nonempty input tournaments of the chosen parity, across all orders. Let G_all be the divisor with both parities allowed, and put g=gcd(a_1,...,a_D).

## 1. Complete unconditioned divisor

The transitive clone sequence is h_j=sum_k a_k k!S(j,k). The D by D matrix (S(j,k)), j,k=1,...,D, is integer triangular with diagonal1. Both it and its inverse have integer entries. Therefore

    G_all = gcd_(k=1)^D (k! a_k)
          = gcd(G_odd,G_even).                          (1)

The first equality uses the same physical transitive inputs and positive resolution as the main theorem. Equivalently, it follows by applying its unrestricted-parity version. This is a dual statement to the S193 input-profile content theorem: the fixed object here is the context, not the input.

## 2. Parity anomalies have only small prime factors

**Theorem.** For every prime ell>2D-1,

    v_ell(G_odd)=v_ell(G_even)=v_ell(G_all)=v_ell(g).     (2)

In particular any prime factor of G_odd/G_all or G_even/G_all is at most2D-1.

**Proof.** Write h_j=sum_(b=1)^D beta_b b^j. The integer triangular change from a to beta has diagonal1, so gcd(beta_1,...,beta_D)=g.

For e=1 or2 the first D parity terms are V_e beta, where

    (V_e)_(t,b)=b^e(b^2)^t, t=0,...,D-1.

Its determinant, up to the harmless positive ordering sign, is

    Delta_e=(product_(b=1)^D b^e)
            product_(1<=i<j<=D)(j-i)(j+i).             (3)

Every positive factor in (3) is at most2D-1. Thus V_e is invertible over integers localized at ell. Multiplication by V_e preserves the minimum ell-adic valuation of the coordinates. The finite parity-gcd theorem then proves the first two equalities in (2). Since ell>D, none of1!,...,D! has an ell factor, and (1) gives the other equality. QED.

An entirely integral quantitative form, with no localization terminology, is

    g divides G_e and G_e divides g*Delta_e.             (4)

The first assertion follows from the integer matrix V_e. Multiplying its output by its integer adjugate shows that G_e divides Delta_e*beta_b for every b. Taking the gcd over b proves the second. In particular G_e/G_all divides Delta_e, since g divides G_all. Formula(2) is also immediate by taking prime valuations in (4) for primes not dividing Delta_e.

The bound is not merely a large-number screen. For the triangle context D=2, G_odd=3 and G_even=1, so the endpoint prime2D-1=3 really occurs. For the carousel-five context D=3, G_odd=15 and G_even=1, so the endpoint prime5 occurs. No sharpness assertion for every D is made.

## 3. Geometric and computational form

The inherited S196 clone matching theorem identifies D=1+nu_r, where nu_r is the maximum matching size between root out-neighbors and root in-neighbors along directed triangles through the root. Thus a parity-specific prime obstruction must satisfy

    ell <= 2nu_r+1.                                    (5)

A large prime cannot mysteriously appear only because the input module has odd rather than even order: beyond this geometric cutoff, it divides every output exactly to the valuation already shared by all run coefficients.

The full universal divisors themselves are witnessed by finitely many physical graphs. Their largest necessary positive input length is2D, so on a rooted quotient of order q the required substituted graphs have at most

    q-1+2D = q+2nu_r+1 <= 2q

vertices. This is a finite certificate bound for universal divisibility, not a bound on witness order for an arbitrary target value and not a coverage theorem.

`verify_parity_primes.py` checks the integral identities on a complete small formal coefficient universe and the named physical contexts. Formal coefficient vectors are used only as algebra controls; they are not asserted to be realizable tournament contexts.

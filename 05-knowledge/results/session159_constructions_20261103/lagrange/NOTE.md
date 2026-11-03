# Branched continued fraction and explicit strong-tournament gadgets

Status: all displayed infinite identities are proved below. Finite calculations independently check small cases. This builds on the path-cover transfer in the [companion path-cover operator note](../profile_operator/NOTE.md); no novelty or exhaustive spectrum claim is made.

## Definitions and operator

For a finite tournament X, let p_k(X) count unordered covers by k nonempty directed paths, and set P_X(z)=sum_k p_k(X)z^k. Let G(X) be the tournament on X and two new vertices a,b with X -> a -> b -> X. For nonempty X, G(X) is strong. The positive path-cover transfer proved in the [companion path-cover operator note](../profile_operator/NOTE.md) is

  p_k(G(X)) = k(k+1)p_{k+1}(X) + k(k+2)p_k(X)
              + (2k-1)p_{k-1}(X) + p_{k-2}(X).                  (1)

All out-of-range terms vanish. Equivalently, with primes denoting d/dz,

  P_G = (z+z^2)P_X'' + (3z+2z^2)P_X' + (z+z^2)P_X.       (2)

The Hamilton-path count is h(X)=p_1(X). Starting from X_0=T_1 and iterating X_{m+1}=G(X_m), (1) gives the linear production rules on levels j>=1:

  j -> j-1 with weight j(j-1),   j -> j with weight j(j+2),
  j -> j+1 with weight 2j+1,    j -> j+2 with weight 1.

The weight of a walk is the product of its steps' weights. In particular h(X_m) equals the total weight of length-m walks from level 1 back to level 1. The lower boundary is automatic because the 1 -> 0 weight is zero. This yields a concrete lattice-walk certificate for each count, not merely an asymptotic analogy.

## Exact branched continued fraction

Let E_j(t) be the generating series for excursions from j to j staying at or above j, using those same step weights. Then, as formal power series,

  E_j = 1/[1 - A_j t - B_j t^2 E_{j+1}
               - C_j t^3 E_{j+1} E_{j+2}],               (3)
  A_j=j(j+2),  B_j=j(j+1)(2j+1),  C_j=j(j+1)^2(j+2).

Proof: An excursion is a sequence of primitive returns to j. A primitive either is a flat step of weight A_j t; starts with a one-level up-step, has an arbitrary excursion at j+1, and returns by a down-step, contributing B_j t^2 E_{j+1}; or starts with a two-level up-step, has an excursion at j+2 until its first descent to j+1, then an excursion at j+1 until its first descent to j, contributing C_j t^3 E_{j+2}E_{j+1}. These possibilities are disjoint and exhaustive. Summing sequences gives (3). The infinite system is unambiguous coefficientwise: a return that visits level j requires at least j-1 up-height and j-1 down-steps, so only finitely many levels affect a fixed coefficient. The desired series is E_1; its first terms are

  1 + 3t + 15t^2 + 123t^3 + 1479t^4 + 23955t^5
    + 496479t^6 + 12725547t^7 + ... .

At the first level,

  E_1 = 1/[1 - 3t - 6t^2 E_2 - 12t^3 E_2 E_3].

Every coefficient in (3) is nonnegative. Thus it can also carry extra nonnegative weights, without cancellation; it is a branched continued fraction, not an asserted ordinary J-fraction or a Pfaffian formula.

## Lagrange inversion for the walk skeletons

Forget the height-dependent step weights and assign indeterminates a,b,c to primitive flat, one-rise and two-rise blocks, respectively. Their total lengths are 1,2,3. A primitive block contains respectively 0,1,2 recursively nested excursions. Let E be their generating series. The same decomposition gives

  E=1/[1-a t-b t^2 E-c t^3 E^2].

More finely, introduce an independent variable s counting primitive blocks and put W=E-1. Then

  W=s(1+W)[a t + b t^2(1+W) + c t^3(1+W)^2].

The formal Lagrange inversion formula gives, with r=i+j+k>0 and m=i+2j+3k,

  [s^r a^i b^j c^k t^m] E
    = (1/r) multinomial(r;i,j,k) binomial(m,r-1).            (4)

Proof: [s^r]W=(1/r)[u^{r-1}]((1+u)[a t+b t^2(1+u)+c t^3(1+u)^2])^r. Choosing i,j,k factors leaves (1+u)^{r+j+2k}=(1+u)^m, yielding (4). The empty skeleton has coefficient 1. These are counts of unweighted excursions/skeletons with the specified numbers of flat, one-rise and two-rise primitives; their length is m. Formula (4) is not itself the weighted tournament count because the actual weights depend on level.

## Parameterized strong gadgets

Start instead with the n-vertex transitive tournament T_n, n>=1, and define Y_{n,m}=G^m(T_n). It has n+2m vertices and is strong for m>=1. A path in T_n is uniquely determined by its vertex set, so p_j(T_n)=S(n,j), a Stirling number of the second kind. Consequently (1) gives the exact general construction

  h(Y_{n,m}) = sum_{j=1}^{m+1} (M^m)_{1j} S(n,j),          (5)

where M is the positive transfer in (1); entries beyond j=m+1 cannot reach level 1 in m steps. In particular, for all n>=1,

  h(Y_{n,1}) = 1+2^n,
  h(Y_{n,2}) = 2*3^n + 5*2^n - 1,
  h(Y_{n,3}) = 6*4^n + 28*3^n + 13*2^n - 11.           (6)

For an explicit derivation, rows 1 of M, M^2 and M^3 are respectively

  (3,2), (15,22,12), (123,266,312,144).

Substituting S(n,1)=1, S(n,2)=(2^n-2)/2, S(n,3)=(3^n-3*2^n+3)/6, and S(n,4)=(4^n-4*3^n+6*2^n-4)/24 proves (6), including n=1. For example, n=2,m=2 gives an explicit 6-vertex strong tournament with 37 Hamilton paths; n=2,m=3 gives an 8-vertex one with 389. These are constructive witnesses; no nonexistence or optimality claim is made.

## A finite exponential formula at every nesting depth

For every m>=0 there are integer coefficients c_{m,q}, independent of n, such that

  h(Y_{n,m}) = sum_{q=1}^{m+1} c_{m,q} q^n,              (7)
  c_{0,1}=1,
  c_{m+1,q} = (q-1)c_{m,q-1} + q^2 c_{m,q} - (q+1)c_{m,q+1}.

All out-of-range coefficients are zero; apply the recurrence only for q>=1. The proof uses the ordered r-slot cover counts F_r(T_n)=r^n and the exact slot transfer F_r(G(X))=rF_{r+1}(X)+r^2F_r(X)-rF_{r-1}(X). In particular, the coefficient of 0^n can be discarded throughout, because n>=1 and that term remains zero under later iteration. Reading off the coefficient of q^n gives (7). This supplies a quadratic-in-m coefficient construction for every n, and is independent of expanding Stirling numbers individually.

The top coefficient obeys c_{m,m+1}=m!, since each recurrence step has just one way to reach the new top index. It follows that for fixed m, h(Y_{n,m}) ~ m!(m+1)^n as n tends to infinity. Moreover the sequence in n is annihilated by the shift polynomial product_{q=1}^{m+1}(E-q), yielding an exact constant-coefficient recurrence for this explicit gadget family. Signs in (7) may be negative; the actual path counts are nevertheless positive.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`. The scripts compare (3) coefficientwise through t^15 with independent forward walk iteration; compare the first seven coefficients against subset/endpoint Hamilton-path DP on the literal tournament; compare (4) with direct unweighted step enumeration through length 12; and check all three parameterized formulas (6) by literal tournament DP for 1<=n<=4. Explicit conditions remain active under Python -O. These are controls, not substitutes for the general proofs.

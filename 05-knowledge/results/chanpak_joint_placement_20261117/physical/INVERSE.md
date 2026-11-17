# Exact joint inversion also works with a changing strong flank

The [fixed-core inverse](../inverse/PROOF.md) extends to the
[variable-height Farey flank](PROOF.md). For each pair (m,G) there
is at most one physical member with given source count N and
signed shift Delta. Searching all heights and scales takes
polynomial bit time: the extra scale does not create a hard
inverse search within this grammar.

Use h,P,A,C as in the physical proof, let S=u+w,d=u-w, and put
q=2^m, l=2q-1. The constant of the physical lower count is
K_m=96q+71+142l, by collecting its explicit formula. Thus

    H = P(S²-d²)/4 + (A+C)S/2 + (A-C)d/2 + K_m + h(F(u)+F(w)).

Let a=4PG+h/(32G), b=3h/(32G), and
c=b(1-10G)-Delta-16PG²-6h+4hG+4G(A-C). Direct substitution of
S+1=G(32y+16z+10) and d=16Gz+2r+1 into the shift law gives

    a d = b S + c - 8h rev(r) + h(2r+1)/(32G).

For even m>=2, P/h>=320 (use 3^m/2^m>=9/4 and 2^m>=4
in the explicit strong-flank formulas), hence the difference
of d from p(S)=(bS+c)/a is less than 1/100. Also A,C are
nonnegative and |A-C|<=P by their displayed formulas. Let
Phi be the quadratic above, replacing d by p(S) and omitting
the F term. Every physical solution satisfies

    |p(S)| <= S+1/100,    |Phi(S)-N| <= PS/30.

For the second bound, charge at most .006PS for the quadratic
change, .005PS for the linear change, and 3hS <= .009375PS
for F (S>=1). On the interval for S cut out by the first
inequality, p has slope b/a<1/10000, so Phi'(S)>=.49PS.
If two values of S satisfy the band inequality, integration
gives .245P(S2-S1)(S1+S2) <= P(S1+S2)/30, hence they
differ by less than one. Integer bisection with Phi±PS/30
finds the only candidate. Recover d by rounding, then z,r,
then x,y as in the fixed-core proof. Verify the exact N, Delta
and the graph conditions. Approximate elimination does not
replace the exact check.

In the positive-address branch all terms of the explicit
physical read before collecting are nonnegative, so
H>=Puw, and u,w>=G. Therefore m=O(log N), G<=sqrt(N/P),
and S<=N. For signed lower *source* N, enumerate O(log² N)
pairs of m,G and do O(log N) rational comparisons at each;
require also any advertised shift/scale gates if desired.
An upper or destination target is handled by subtracting 2
or the prescribed shift. This is exact inversion for the
growing grammar, not density or existence at every target.

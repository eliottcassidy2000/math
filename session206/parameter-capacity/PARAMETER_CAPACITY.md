# Exact parameter capacity of an O/Z tail source

This is an all-order necessary condition for finite-bank coverage, not a covering theorem and not a general finite-bank impossibility theorem.

Let a physically attained unit-start source have state s=(h,R,e,d), with d>0 and g=R-d>0. The established tail formula is

H_s(m)=h+Rm-dQ(m)-e(m mod 2), m>=0,

where Q(0)=0, Q(2m)=Q(m), and Q(2m+1)=Q(m)+m+(m mod 2). Define the parameter count, with multiplicity,

A_s(X)=#{m>=0:H_s(m)<=X}.

Then, as X tends to infinity,

A_s(X) = (X/d) log(R/(R-d)) + O_s(sqrt(X)).

In particular the coefficient is strictly smaller than the earlier reciprocal-width bound 1/(R-d). For any finite physical bank B, eventual coverage of all odd integers requires

sum_{s in B} log(R_s/(R_s-d_s))/d_s >= 1/2.

This is only a necessary condition. Distinct m can yield the same scalar, and different seeds can overlap. For the rooted edge alone the parameter coefficient is log 2 > 1/2, so this result does not even prove that the rooted-edge tail image omits infinitely many odds.

## Digit identity

Write m=sum_{j>=0} b_j 2^j, b_j in {0,1}, with finitely many nonzero digits. Unrolling the Q recurrence gives exactly

Q(m)=sum_{j>=0} b_j floor(m/2^(j+1)) + sum_{j>=0} b_j b_{j+1}.

Put u(m)=sum_{j>=0} b_j/2^(j+1). If ell(m) is the binary length, then

|Q(m)-m u(m)| <= ell(m).

Indeed the error equals the sum of -b_j frac(m/2^(j+1))+b_j b_{j+1}; each summand is in [-1,1] and there are at most ell(m) terms. This identity is valid at m=0 with ell(0)=0.

For L>=1 let u_L(m) retain only j<L. Then 0<=u(m)-u_L(m)<2^(-L). Consequently, uniformly for 0<=m<=C X,

H_s(m)=m(R-d u_L(m))+O_{s,C}(log(2+X)+X/2^L).

Crucially, u_L depends only on m modulo 2^L; as that residue runs over its complete set, its values are exactly j/2^L, j=0,...,2^L-1, in bit-reversed order.

## Residue-class count

Set K=2^L and F_L(m)=m(R-d u_L(m)). Every slope is at least g. For Y>=0, counting separately in each residue class modulo K gives

#{m>=0:F_L(m)<=Y}
 = (Y/K) sum_{j=0}^{K-1} 1/(R-dj/K) + O(K).

The constant in O(K) is absolute: each residue-class count differs from its interval-length/K estimate by at most a constant.

For large X choose K to be a power of two between sqrt(X)/2 and sqrt(X). Both H_s(m)<=X and F_L(m)<=X+O_s(sqrt(X)) restrict m to m<=C_s X: the former follows from the established inverse bound H_s(m)>=h+g m, and the latter from F_L(m)>=g m. On this whole range the digit error is O_s(sqrt(X)). Sandwiching the H-count between F-counts at X plus or minus that error therefore gives

A_s(X) = (X/K) sum_{j=0}^{K-1} 1/(R-dj/K) + O_s(sqrt(X)).

The summand is a continuously differentiable function on [0,1], with denominator bounded below by g>0. Its left Riemann sum differs from its integral by O_s(1/K). Hence

A_s(X)=X integral_0^1 du/(R-du)+O_s(sqrt(X))
      =(X/d) log(R/g)+O_s(sqrt(X)).

All constants depend on the fixed source, not X. Summing over a fixed finite bank proves the asserted capacity condition, since a cofinite odd set has X/2+O(1) members and the union of scalar images is bounded above by the sum of parameter counts.

## Scope and limits

Physical attainability remains necessary for using a seed as a tournament construction. The arithmetic asymptotic itself only needs fixed real parameters d>0 and R>d, with a harmless fixed h,e correction; the inverse bound used above is available for physical states (and a weaker bound suffices arithmetically). No image-density asymptotic is asserted. No finite computation is used in the proof. This does not settle whether some finite bank covers all sufficiently large odd scalars.

# Exact rooted slot reciprocity: the source of the factor 3

The ordinary weighted path-end congruence is the first specialization of a
polynomial reflection about the half-integer slot number (-1/2).
The resulting modulus applies to fixed rooted signatures, not to every
possible tournament repair.

## Squarefree path identity

Work with commuting variables (x_v) and (x_v^2=0). Let (P) be 1 plus
the sum of the monomials of all nonempty directed paths. Let (R_v^+)
and (R_v^-) be the sums restricted respectively to paths starting and
ending at (v); sums count each distinct path, not only its vertex set.
Then

~~~text
R_v^+(x)=-P(x)R_v^-(-x).
~~~

For a proof, induct on the vertex set and use
(R_v^-(-x)=-x_v(1+Σ_{u→v}R_u^-(-x))).
Multiplying by (x_v) discards terms involving (v), so the induction
applies to each (u) in the deletion (W-v). The right-hand side above
reduces to (x_v(P-Σ_{u→v}R_u^+)), which is
(x_v(1+Σ_{v→u}R_u^+)=R_v^+), because the tournament chooses
exactly one direction for each pair.

Since (Σ_v R_v^+=Σ_v R_v^-=P-1), summing this same identity
also proves (P(x)P(-x)=1). This is the signed two-slot cancellation
in the [emptying proof](EMPTYING_EXCLUSIONS.md), now proved by rooted
path splitting alone. It holds in every induced subtournament.

## Exact polynomial reflection and lifted endpoint congruences

For a nonnegative integer (k), the coefficient of the squarefree
vertex-set monomial of (A) in (P^k) is the number (Q_k(A)) of
ordered (k)-slot path covers allowing empty slots. It extends to all
integer (k) by nilpotence, and the [odd-cycle packing expansion](https://arxiv.org/abs/2307.05569)
gives the integral polynomial
(Q_k(A)=Σ_D2^c k^{|A|-s+c}), where the packing
D contains (c) pairwise vertex-disjoint odd directed cycles
on (s) vertices. To see the counted weight, assign each selected cycle
and each unoccupied vertex to a path slot. Thus (Q_k(-1)) or any
negative evaluation means polynomial evaluation, not a literal
negative number of slots.

Define (S_v(k)=Σ_{L starts v}Q_k(W-V(L))), and (T_v(k))
by ending at (v), as in the question when (k=1). Multiply the rooted
identity by (P^k), then replace (x) by (-x), using
(P(-x)^{k+1}=P(x)^{-k-1}). At order (m) this yields

~~~text
S_v(k)=(-1)^{m+1}T_v(-k-1).
~~~

The displayed equality is polynomial, since it holds for all
nonnegative integers (k). Therefore (T_v(k)-(-1)^{m+1}S_v(k)) is
divisible by (2k+1) as an integral polynomial. In particular:

- At (k=1), T_v ≡ (−1)^(m+1) S_v (mod 3).
- For every odd (M), taking (k=(M-1)/2) gives this endpoint relation
  modulo (M), but for **(k)-slot weights**, not Hamilton counts if
  k≠1.
- For Hamilton weights the first stronger defect is
  ~~~text
  [T_v(1) − (−1)^(m+1) S_v(1)]/3 ≡ T'_v(1) (mod 3).
  ~~~
  The derivative marks one assigned odd-cycle block or uncovered
  singleton in the path complement. It is not in general zero.
  Already the source of TT2 has (S=2,T=1), a defect of 3.

The symmetric combination (T_v(k)+(-1)^{m+1}S_v(k)) and the defect
divided by (2k+1) are both fixed by k ↦ −k−1; thus both are
polynomials in (k(k+1)) with integer coefficients. Indeed any invariant
integral polynomial has even degree and may subtract a multiple of
([k(k+1)]^{d/2}) to lower degree inductively. This distinguishes
what the reflection forces (one factor of 3 at one slot) from what
requires further information.

## A literal same-count caution and one design rule

For the two-wing family in the [emptying note](EMPTYING_EXCLUSIONS.md),
put (X=2^{a-1}) and (Y=2^{b-1}). The common root has
(S_r=2XY+2Y-1) and, by converse, (T_r=2XY+2X-1).
The common scalar (H(W)) does not record this directional bias.
After attaching the outlet, swapping the wing sizes changes the
count by (2(Y-X)). Its 3-adic valuation is zero if (a-b) is odd;
if (a-b) is nonzero even it is
(1+v_3((a-b)/2)), by applying the geometric-sum valuation to (4^d-1).
In particular the (H=15) deleted cores of wing sizes ((2,3),(3,2))
yield (37) and (33), a difference of four. This is not a contradiction
to the endpoint identity: the identity retains the complementary
orientation and order parity for the **same** root, whereas these two
repairs exchange which response the outlet reads. A step-of-three
rule therefore needs a matching-signature hypothesis; it is not a
universal step rule for physical repairs.

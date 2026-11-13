# Unequal modules do not neutralize the fixed H25 bridge

**PROVED IN NOTES:** every one- or two-port nontrivial inflation fails;
the exact three-port response criterion below. **FINITE-EXACT:** coefficient
and literal graph certificates. This is a theorem about the specified bridge,
not about all exact-count repairs or all three-slot substitutions.

## The inherited pair and the sharper question

Use the earlier pair

    A=(44,17,2,6,13,30,63),
    B=(110,16,2,6,13,30,62).

Both have H=25. A is a singleton source before a strong six-core; B is
strong. The two changed arcs are {0,1} and {0,6}. The earlier
[H25-BRIDGE](../../strategy_fibres_absorption_20261113/fibres/H25-BRIDGE.md)
proved that no single vertex admits a nontrivial arbitrary-module inflation
preserving the equality, and tested identical transitive modules at two ports.
The missing question was whether **different** modules could balance the
positive and negative response changes.

## Exact two-module response and physical constraints

Let X,Y be nonempty tournaments substituted at vertices p<q. Put
x_i=K_i(X), y_j=K_j(Y), for ordered nonempty path covers. The induced
two-vertex quotient is acyclic, so between successive exterior vertices a
run word can visit each module at most once. There are finitely many legal
run words, and

    H(B[X,Y])-H(A[X,Y]) = sum D_ij x_i y_j.            (1)

The complete matrices D for all 21 port pairs are in the receipt. They are
computed by a finite exterior-word recurrence and independently compared
with literal subset DP for 378 substituted graphs.

Every physical nonempty module satisfies:

* H=K₁ is a positive odd integer.
* A singleton has K_i=0 for i≥2.
* For a nonsingleton, K₂≥2H and K₂≡2 (mod4).
* K_i is divisible by i!, since ordering its nonempty disjoint paths is a
  free permutation action. In particular K₃ is even and K_i≡0 (mod4) for i≥4.

For the K₂ lower bound, cutting Hamiltonian paths gives
K₂≥(|M|−1)H when |M|≥3; order2 has K₂=2H directly. For the congruence,
P₂=Σ_U H(U)H(Uᶜ) is divisible by4 at every order at least two. Pair U with
Uᶜ; there are 2^(n−1) pairs, each with odd weight. Therefore
K₂=P₂−2H≡2 modulo4.

## The all-module two-port exclusion

**Theorem.** For each pair p<q and every nonempty X,Y,

    H(A[p←X,q←Y]) = H(B[p←X,q←Y])

holds if and only if X and Y are both singletons.

**Proof.** The finite response matrices divide into four complete categories.

1. At (0,1),(0,2),(0,3),(0,5),(1,2),(1,3),(1,5),(2,3),(2,5),(3,5),
   every nonzero coefficient in (1) is negative. Both D₁₂ and D₂₁ are
   negative. Thus any nonsingleton module gives a strictly negative change.
2. At (4,6), all nonzero coefficients are positive, with D₁₂,D₂₁>0.
3. At (2,4),(2,6),(3,4),(3,6),(5,6), replacing y₂ by 2y₁+L, L≥0,
   makes every coefficient nonnegative and the x₁y₁ coefficient positive.
   Hence the change is positive whenever Y is nonsingleton. If Y is a
   singleton, the remaining term is −2x₂. For (4,5), the same argument uses
   x₂=2x₁+L; if X is a singleton, the remaining term is −2y₂.
4. At (0,4),(0,6),(1,4),(1,6), the complete matrix gives a difference
   congruent to2 modulo4 whenever X is nonsingleton. Every other term
   vanishes modulo4 using the divisibilities above. If X is a singleton,
   the remaining expression is either 22y₂+12y₃ or 10y₂, so it is positive
   for nonsingleton Y.

These 10+1+6+4 cases exhaust all pairs. The singleton/singleton case is the
original H25 equality. ∎

For example, at the genuinely mixed pair (0,6),

    ΔH=10x₁y₂−9x₂y₁+10x₂y₂−6x₃y₁.

The relaxed positive response cone has points on this kernel: for example
(x₁,x₂,x₃;y₁,y₂)=(1,2,12;1,3). But the physical congruence K₂≡2 modulo4
excludes that point and, by the proof, every actual nontrivial module pair.
This is the useful Diophantine connection: an algebraic response equation
and positivity do not characterize the realizable integer response set.
The missing coordinate is a small-adic congruence, not a larger size cutoff.

## A three-port escape produces a genuine ternary child

The exclusion does not extend to every operation on three ports. Put a
singleton at0, TT_b at1, and an arbitrary nonempty M at6, b≥1. Let
h=H(M), p=P₂(M). Exact run enumeration gives

    source = (14·3^b−17)h,
    target = 2(4·2^b−3)p+(3·2^b−1)h.                (2)

A remains an ordinal source M followed by the inflated six-core; B is
strong. The scalar difference from the ternary target is

    target−3source
      = 2(4·2^b−3)p+(3·2^b+50−42·3^b)h.           (3)

For b=2, take M to be the old strong six-core

    (44,17,2,6,13,30),      (H,P₂)=(25,304).

Equation (3) becomes 2(13p−158h)=4. Thus the actual two thirteen-vertex
graphs have counts

    2725 -> 8179 = 3·2725+4,

and 8179 is prime. The literal rows and independent subset counts are in the
receipt. The transformation reverses eight arcs: two from the inflated
{0,1} pair and six from {0,6}. It is not count-neutral and does not preserve
the old graph as an induced subgraph. It is an exact target-specific core
replacement beyond the two-port obstruction.

Formula (3) is a necessary and sufficient response criterion for this entire
three-port chart. Its existence premise is still physical: an arithmetic
pair (h,p) must be realized by M. For fixed b and fixed digit, p/h is bounded,
so the path-cut inequality bounds |M|. Consequently the one positive module
does not give an unbounded fixed-b family. Varying b or changing the core is
necessary for further growth. This boundary motivated the distinct induced
family in [TERNARY.md](TERNARY.md), which succeeds for every m with unit
module count rather than attempting to renew this fixed mixed response.

Reproduce both the coefficient theorem and the literal three-port success
with [check.py](check.py). The theorem uses physical factorial divisibility
and a complete finite response table; it does not extrapolate from small
module samples.

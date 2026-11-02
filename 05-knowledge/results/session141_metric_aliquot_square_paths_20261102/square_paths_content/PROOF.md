# Binary content from reversal, free orientations and joint multiplicity

**Status: PROVED elementary refinement; exact n25/n26 application.** This
note adds no new census. It uses the already verified joint-union profiles.
The content/reversal statement below was not found in the three explicitly
inspected inherited dependencies of the square-path packet; no global
repository or literature-priority claim is made.

## 1. Integer translation gives the exact content

Fix an undirected guard graph `G` on `n>=2` labelled vertices, and let `P`
be its complete directed Hamilton-path bank. More generally, the theorem
holds for any finite bank of distinct directed spanning paths closed under
reversal. Set `N=binom(n,2)` and

    Z(z) = sum_T z^(H_P(T)) = sum_(j>=0) M_j (z-1)^j,
    M_j = sum_(compatible S subset P, |S|=j) 2^(N-|U(S)|).

In particular `M_0=2^N`. These are **binomial moments**
`M_j=sum_T binom(H_P(T),j)`, or falling-factorial moments divided by `j!`.
The convention matters already at `j=2`.

**Theorem.** Let `content(Z)` be the positive gcd of all its ordinary
integer coefficients. Then

    content(Z) = gcd(M_0,M_1,...) = 2^e,
    e = min_(M_j != 0) v_2(M_j).                            (1)

**Proof.** Each coefficient of `Z(z)` is an integer combination of the
`M_j`, by expanding `(z-1)^j`. Conversely, substitute `z=w+1`; every
`M_j=[w^j]Z(w+1)` is an integer combination of the coefficients of `Z`.
Thus each gcd divides the other. The common gcd divides `M_0=2^N`, so
it is a power of two, and its exponent is the stated minimum. This is
integer translation, not a division modulo a possibly nonunit modulus.
QED.

## 2. The guaranteed binary factor

Assume `P` is nonempty and put

    Umax = max_(nonempty compatible S subset P) |U(S)|.

**Theorem.** Every coefficient of `Z` is divisible by

    2^(N-Umax+1).                                         (2)

**Proof.** Reversing every member maps a compatible subset to another
compatible subset of the same size and the same union-arc cardinality.
This involution has no fixed nonempty compatible subset: a fixed subset
would contain both some path `p` and its reversal, prescribing opposite
directions on its `n-1>=1` path edges. Thus, for every `j>=1` and union
size `u`, the number `c_(j,u)` of such subsets is even. Their paired
contribution `c_(j,u)2^(N-u)` is divisible by `2^(N-Umax+1)`. This holds
for all nonconstant binomial moments. Also `Umax>=1`, so the exponent in
(2) is at most `N` and divides `M_0=2^N`. Apply (1). QED.

The bank-empty case is separate: `Z(z)=2^N` and its content is `2^N`;
`Umax` is undefined. The `n>=2` and reversal-closure hypotheses are real:
at `n=1`, the self-reverse empty-arc path gives `Z=z`; a bank containing
only one direction of the two-vertex path gives `Z=1+z`. Neither has an
extra factor of two.

The bound (2) need not be sharp. Profiles `c_(j,u)`, including their
multiplicities, determine the exact moments, and (1) then determines the
exact additional binary factors. Pairwise compatibility alone does not
supply these weights.

## 3. Why the factors at 25 and 26 differ

At order 25, the nonempty profiles are

    c_(1,24)=20, c_(2,26)=4.

Thus `N=300`, `Umax=26`, and reversal guarantees only `2^275`. The exact
binomial moments have valuations

    v2(M0)=300,
    M1=20*2^276, v2(M1)=278,
    M2=4*2^274,  v2(M2)=276.

Formula (1) proves `content(Z25)=2^276`, one binary factor beyond the
universal reversal bound. The primitive polynomial is

    Z25(z) = 2^276 (16777197 + 18z + z^2).

The extra factor is visible in the FOUR compatible pairs, not merely in
the existence of an involution. The primitive constant and quadratic
coefficients are odd.

At order 26 the profiles are

    c_(1,25)=24,
    c_(2,27)=8, c_(2,29)=2, c_(2,30)=2,
    c_(3,29)=2.

Here `N=325`, `Umax=30`, so (2) guarantees `2^296`. In guard-orientation
coordinates the moments are `2^34,12288,1120,64`; multiplying by the
`2^291` free outside orientations gives full moment valuations

    325, 303, 296, 297.

Therefore this time the reversal bound is SHARP, and

    Z26(z) = 2^296 (536870561 + 320z + 29z^2 + 2z^3).

In particular the quadratic coefficient 29 and the constant coefficient
are odd. The change is explained by free-orientation dimensions, reversal
orbits and exact joint multiplicities. It is not inferred from an
unrelated appearance of 1,2,4,20 or a numerical difference such as 10.

## 4. A leaf gives a concrete reversal half

If the guard has a leaf `ell`, every spanning Hamilton path uses it as an
endpoint. Partition the bank into paths starting at `ell` and paths
ending there. A compatible family cannot mix the two halves because they
prescribe opposite directions on the unique leaf edge. Reversal bijects
the halves, preserving subset size and union cardinality.

Consequently, for `j>=1`, each `c_(j,u)` is twice the corresponding count from the
leaf-starting half. This realizes the free involution with a distinguished
source rather than an abstract orbit choice. A unique leaf is sufficient
but unnecessary; any specified leaf gives this partition. At n25 and n26,
18 is a leaf with neighbour 7, and the displayed source-preserving block
carriers begin at 18. This explains their reversed partner families without
postulating a new arithmetic symmetry.

## 5. Two distinct overlap statistics, and a finite probability reversal

At n26, the maximum compatible-family size is three, but `Umax=30` is
attained by compatible PAIRS; the triples use only 29 arcs. Maximum
response and maximum prescribed-edge support are therefore different
statistics. The latter controls the reversal bound. Indeed the 30-arc
pairs contribute `2*2^295=2^296` to `M2`; its other contributions are
divisible by `2^297`, leaving `M2=35*2^296` and forcing exact content.

For uniformly chosen full labelled tournaments, the already established
coefficients give

| order | maximum guarded response | expected response | probability of positive response |
|---:|---:|---:|---:|
|25|2|`5/2^22`|`19/2^24 = 608/2^29`|
|26|3|`3/2^22`|`351/2^29`|

Thus the new three-path overlaps coexist with both a smaller mean and a
smaller probability of any allowed path. These are exact consequences of
the two finite distributions, not an asymptotic or monotonicity law.

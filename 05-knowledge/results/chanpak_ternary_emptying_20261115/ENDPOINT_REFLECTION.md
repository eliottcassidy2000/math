# An exact reflection of the rooted endpoint polynomial

The user's mod-3 endpoint relation is the fixed point of an
integer-polynomial identity, not a universal restriction on repairs.

## Definitions and the induced-subtournament recurrence

For a tournament U let Q_U(a) be the number of ordered covers
by a path slots, allowing empty slots, evaluated first at a
nonnegative integer a. Expand the Hamilton path count of each
slot using the odd-cycle packing formula of
[Grinberg and Stanley](https://arxiv.org/abs/2307.05569).
An odd-cycle packing with c nontrivial cycles and k blocks
(including unused singleton vertices) contributes 2^c a^k:
each block independently chooses one slot. Thus Q is a polynomial
in Z[a], Q_empty=1, and Q_U(−a)=(-1)^|U| Q_U(a), since every
block is odd.

For vertex v in T, let S_v(a),T_v(a) sum Q_(T minus P)(a) over
all nonempty paths P starting or ending at v respectively.
When a=1 these are exactly the weighted path start and end.
At a=0 they give ordinary Hamilton path endpoint counts.
Delete v and put U=T−v. Then
```
S_v(a) = Q_U(a) + sum_(v→w) S_w(U;a),
T_v(a) = Q_U(a) + sum_(w→v) T_w(U;a).
```
The singleton path supplies Q_U(a); every longer path has
a unique next or previous vertex. An extra slot gives
`sum_w S_w(U;b)=Q_U(b+1)−Q_U(b)` by distinguishing its
nonempty path. These are identities on **every** induced
subtournament.

Let m=|T| and sigma=(-1)^(m+1). Induction with b=−a−1
gives the exact relation:
```
T_v(a) = sigma S_v(−a−1).
```
Indeed for each w in U, the induction sign is −sigma, so
T_v(a)/sigma = Q_U(−a)−sum_(w→v) S_w(U;b).
Now Q_U(−a)=Q_U(b)+sum_w S_w(U;b), which proves the formula.
The singleton base has both endpoint weights one.

There is also an exact signed version. If a rooted path is
weighted by epsilon^(|P|−1), epsilon in {+1,−1}, throughout
both definitions, then T_v^epsilon(a)=sigma S_v^epsilon(−a−epsilon).
In the deletion recurrence each prepended vertex contributes
epsilon, and the sum of the signed (epsilon=−1) start weights
is Q_U(b)−Q_U(b−1): let the distinguished block color have
weight −1 per vertex, which is −1 per odd block. Equivalently,
for either sign, epsilon times the sum is Q_U(b+epsilon)−Q_U(b).
The same induction proves the signed identity. Its a=1,
epsilon=−1 case reads the *ordinary reverse endpoint*, as in
the inherited fan identity, not the weighted endpoint.

## Optimal modulus and the first obstruction after reduction

For every nonnegative integer a,
```
T_v(a) ≡ sigma S_v(a) mod (2a+1).
```
This modulus is **optimal at every fixed order m>=2** for
unsigned endpoints over all tournaments and vertices. In a
transitive tournament with v at position i (from the source),
S_v(a)=a^(i−1)(a+1)^(m−i) and T_v(a) is the expression with
i−1 and m−i swapped. For m even, take i=m/2; for m odd m>=3,
take i=(m−1)/2. The resulting difference is, up to sign,
(2a+1) times a power of a(a+1); the latter is coprime to
2a+1. To exclude an extra common factor outside 2a+1, use the
first vertex: its difference is coprime to a(a+1).
For a=0 the source itself gives difference of magnitude one.
The m=1 identity is exact, a separate boundary.

Use Hasse derivatives S^[j], defined by
S_v(a+b)=sum_j S_v^[j](a)b^j. Then the exact discrepancy at
a=1 is
```
[T_v(1)−sigma S_v(1)]/(-3 sigma)
       = S_v^[1](1) −3 S_v^[2](1) +9 S_v^[3](1) − ...
```
The first lift to mod 9 is therefore governed by S_v^[1](1)
mod 3. It counts a rooted path, a packing of its complement,
and one distinguished block (possibly a singleton), weighted
by the packing. Higher lifts distinguish several complementary
blocks. This is a literal emptying-side statistic: the counted
blocks live **outside** the rooted path. The residue relation
discarded that information and cannot be lifted without it.
For two complementary slots the analogous exact reflection
modulus is 5, a different readout, not a claim about all H=5
tournaments or about a physical five-step repair.

## The window step it really constrains

When comparing two same-order rooted cores at the same
unsigned S value, the unsigned T difference is divisible
by 3 at a=1. If the signed end totals also agree, the unsigned
difference is even (unsigned and signed totals have the same
parity), hence the difference is divisible by 6. If an ordinal
reader sees the unsigned end linearly and all its other relevant
coordinates are held fixed, its repairs are therefore multiples
of three **odd-count half-steps**. Changing the reader or the
data held fixed invalidates that argument.

For a simple physical demonstration, take the strong order-7
core with outgoing rows (76,9,98,84,71,91,2), and root 0.
It has H=81, equal ordinary start and end 9, and weighted
S=131,T=137. An apex beating just the root gives H+S−9=203;
an apex losing just to the root gives H+T−9=209. Both strong,
their difference is 6. The two interfaces need not agree:
this example proves a scalar difference, not matched continuation.
The sharper matched regime is treated separately.

In contrast, [primitive two-tape repairs](../chanpak_viability_shape_20261115/PRIMITIVE_REPAIR_LANGUAGE.md)
give actual +2 shifts. The cyclic transport spectrum has an
exact ternary count of its own, derived in the companion note.
These mechanisms are different and consistent.

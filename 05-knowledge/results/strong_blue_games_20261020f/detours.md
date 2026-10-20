# Two-block repairs and arc-disjoint squares for strong blue edges

Status: **PROVED** local identities and sufficient blackness conditions;
**FINITE-EXACT** small presentation controls. These give black detours for
specified blue lines, not all-order black-core connectivity. The new two-block
repair closes the explicit two-chord family in section 3 at every order.

## 1. Inheritance, conventions, and the missing coordinate

Write `T^P` for reversal of exactly the arcs of a directed Hamiltonian path
`P` in `T`. The repository's complement operation is `K_P(T)=(T^P)^op`;
these give the same move after merging converse. A path is **blue** precisely
when its positional reversal `r_P` is an anti-automorphism of its tournament.
A black path need not move to a different merged node.

The inherited all-order [outer-component theorem](../spectrum_basins_20261020e/outer.md)
bypasses every mixed blue line incident to a nonstrong tournament. Together
with the [pure-blue hubs](../spectrum_basins_20261020e/pure_hubs.md), this
reduces the remaining connectivity problem exactly to strong-to-strong mixed
blue lines. Descent alone does not join its terminal SCC basins.

The closest earlier square result is section 5 of
[the paired-black bridge packet](../paired_black_bridges_order_nine_20261020/BRIDGES_THROUGH_NINE_AND_LOCAL_CYCLE_LAWS.md).
It warns that an all-black square may fold in the quotient. Here the desired
consequence is a path between the two endpoints of a *blue* edge, so repeated
merged vertices do not invalidate a black detour. No claim about a specific
paired-black edge being nonbridging is inferred from the square.

The least-used sidecar is the actual positional anti-involution. Changing a
Hamiltonian presentation nonautomorphically does **not** imply that it becomes
black in an arbitrary strong tournament. For example, in the circulant
`C5` with forward steps `{1,2}`, both `01234` and `02341` are blue, although
the latter permutation is not an automorphism. The separated outer SCC in the
earlier proof was essential to its stronger inference. This hostile is replayed
in `detours_square.py`.

## 2. A two-block endpoint repair

**Theorem (PROVED).** Let `P=(L,a,R)` be blue in an n-vertex tournament `T`,
where `L` and `R` are nonempty directed path blocks. Put `x=first(L)` and
`z=last(R)`. Assume:

1. `Q=(L,R,a)` is directed. In addition to the original path arcs, this asks
   for `last(L)->first(R)` and `z->a`.
2. In `T[L]` after reversal of its displayed path, there is a Hamiltonian
   path `H_L` with the same first and last vertices as `L`. The analogous
   endpoint-preserving path `H_R` exists in `T[R]` after reversal of `R`.
3. `d_T(x)+d_T(a)` belongs to neither `{n-1,n}`.

Then the following three **black R moves** end exactly at `T^P`:

```
Q = (L,R,a)
U = (H_L,a,H_R)
V = (a,reverse(H_R),reverse(H_L)).
```

The equivalent K words are `Q,reverse(U),V`; they end exactly at `K_P(T)`.

**Proof.** After reversing `Q`, the induced tournaments on the two blocks
have precisely the path reversals stipulated in assumption 2. The two
connecting arcs in `U`, namely `last(L)->a` and `a->first(R)`, are original
P arcs absent from Q and are unchanged. Hence U is directed. Reversing U
makes both reversed H blocks directed. The other two V joins are
`a->z` and `first(R)->last(L)`; these are reversals of the two Q joins and
were not used by U. Thus V is directed.

Each internal block edge used by H occurs in both U and V and cancels.
Internal L and R edges occur only in Q and remain reversed. The two joins
of Q are canceled by V, leaving precisely the two joins of P from U.
Consequently `E(Q) xor E(U) xor E(V)=E(P)`, proving the exact endpoint.

A path reversal changes only the scores of its endpoints: its start loses
one and its end gains one. Since P is blue, `d_T(x)+d_T(z)=n-1`. The three
endpoint score sums, in their respective input tournaments, are

```
d_T(x)+d_T(a),   n-2,   d_T(x)+d_T(a)-1.
```

Every blue path has endpoint score sum n-1. Assumption 3 therefore makes
all three moves black. QED.

**Stronger predicate version.** Assumption 3 may be replaced by: Q is
actually black, and `d_T(x)+d_T(a) != n`. The middle move is always black;
the last is then certified by scores. If the last sum is balanced, a direct
positional-reversal check can still establish blackness, but it is not
implied by the theorem.

This construction preserves the full labelled target and requires no
isomorphism search, no strongness hypothesis, and no realization of an
abstract wiring. The lost information in a score-only application is exactly
the two endpoint-preserving Hamiltonian paths, which must be supplied as
sidecars.

### Transitive blocks supply the missing paths

Let `J_m` be the transitive tournament on `0,...,m-1` with all consecutive
arcs reversed. For every `m>=4`, `J_m` has a Hamiltonian path `h_m` from
0 to m-1. Bases and a recursion are

```
h4 = (0,2,1,3)
h5 = (0,3,2,1,4)
h6 = (0,2,1,4,3,5)
h_(m+3) = (h_m,m+1,m,m+2).
```

The three new joins are respectively a forward distance-two arc, a reversed
consecutive arc, and a forward distance-two arc. Thus the recursion is
directed and preserves the required endpoints. This proves the assertion
for all sizes; the computation merely checks its implementation. The boundary
`m=2,3` has no such endpoint-preserving path. A fixed two-swap word works
at m=4 but is **not** a valid substitute for this recursion at larger sizes.

## 3. An infinite hostile family now has black detours

Let `T_k`, `k>=2`, be `TT_(2k+1)` with exactly the two arcs `(0,k)` and
`(k,2k)` reversed. The standard path `P=(0,...,2k)` is blue: reflection
`i -> 2k-i` exchanges the two exceptional pairs and reverses every other
arc. The source is strong, because it contains the Hamiltonian cycle
`(0,1,...,k-1,k+1,...,2k,k,0)`. Its blue partner is strong as well: P is
open in T, so complementation closes P. The certificates below themselves
show both endpoints are black-active.

**Corollary (PROVED).** The blue line given by P has a black detour of
length at most three in the converse-merged metagraph, for every `k>=2`.

For `k>=4`, apply section 2 with `L=(0,...,k-1)`, `a=k`, and
`R=(k+1,...,2k)`. Both blocks are transitive of size k and both extra Q
joins have the required orientation. Here

```
d(0)=2k-1, d(k)=k, d(2k)=1,
three endpoint sums = 3k-1, 2k-1, 3k-2,
n-1 = 2k.
```

All three sums differ from 2k for k>=4. The paths are explicit using h_k.
`detours_family.py` replays the exact labelled identity and blackness for
k=4,...,100, i.e. orders 9,...,201.

The two smaller boundaries require the isomorphism/converse quotient:

- k=2: from T use the black R path `(1,3,4,2,0)`; from `T^P` use
  `(0,4,3,2,1)`. Their resulting tournaments are isomorphic. This is a
  two-black detour.
- k=3: from T use black R paths `(0,1,2,4,5,6,3)` and
  `(0,2,1,3,4,6,5)`; from `T^P` use `(3,2,1,0,6,5,4)`.
  The source result maps to the converse of the target result by the
  vertex map `(4,5,3,6,0,2,1)`. This is a three-black detour.

The small certificate file retains the exact map at k=2 as well. These are
finite checks at two specified orders, not an induction inferred from a
census. The all-order part is the two-block theorem and h recurrence.

## 4. The exact colour law for an arc-disjoint square

Let P,Q be arc-disjoint directed Hamiltonian paths in T, and assume P is
blue. Reversals commute, and each path remains directed after reversing the
other. Then

```
P is blue in T^Q  iff  r_P(E(Q))=E(Q)  iff  r_P=r_Q.
```

The first equivalence follows by toggling the Q edge set in the
anti-automorphism equation for r_P. The only automorphisms of an undirected
path are its identity and reversal. Since r_P is not the identity for n>=2,
preservation of the Q path graph is exactly `r_P=r_Q`.

Two useful consequences are therefore **PROVED**:

- If Q is black both in T and in `T^P`, then `Q,P,reverse(Q)` are three
  black R moves ending exactly at `T^P`. Indeed, Q black in T implies
  `r_Q != r_P`, making the middle P move black. The third has the same
  colour as Q in `T^P`. K words are `Q,reverse(P),reverse(Q)`.
- If Q is also blue in T and `r_P != r_Q`, the two partners `T^P,T^Q`
  are joined by the other **two black** sides of the square. If their
  reversals coincide instead, the square is all blue.

A convenient sufficient certificate for the first consequence uses endpoint
scores. Write `P=(a,...,z)`, `Q=(u,...,v)`, and `D=d_T(u)+d_T(v)`. It is
enough that both

```
D != n-1,
D - [u=a] - [v=a] + [u=z] + [v=z] != n-1.
```

No extra endpoint-unpaired condition is necessary: Q black already implies
different positional reversals. Actual blackness checks are stronger than
these sufficient score inequalities.

**Rigid-partner refinement (PROVED).** If `Aut(T^P)` is trivial, any
P-arc-disjoint Q that is black in T also stays black in `T^P`; the first
square consequence then applies. Indeed, r_P remains an anti-automorphism
of `T^P`. If Q became blue there, its distinct reversal r_Q would be a
second anti-automorphism, and their product a nonidentity automorphism.
This contradiction removes the second blackness test under the stated
rigidity hypothesis. It does not guarantee an arc-disjoint black Q exists:
the entire T_k family has just one arc-disjoint Q, always blue, as proved in
[the family counting note](double_return_counts.md).

## 5. Finite scope and reproducibility

`detours_probe.py` exhausts standard-blue presentations of orders 4..7,
orients to the open side, retains strong mixed endpoint pairs, and removes
ordinary-isomorphism loops. It tries a fixed-boundary common hub, the exact
square, and literal suffix/extraction templates. The following counts overlap
across mechanisms:

| n | Nonloop presentations | Square covered | All four tests uncovered |
|---:|---:|---:|---:|
|5|3|0|3|
|6|25|11|10|
|7|239|209|8|

The new two-block construction is not included in that earlier small probe;
its transitive-block application starts at order 9. The root's independent
rotation-hub incidence test, combined with these tests, leaves only two
presentations at order 5, seven at order 6, and one at order 7. The sole
order-seven residual is T_3 above. The uniform family corollary repairs that
mechanism failure rather than inferring disconnection from it.

At order five, the strong blue presentation `(18,4,9,19,6)` has no
P-arc-disjoint Hamiltonian path at all. The two other open nonloop examples
have a unique arc-disjoint Q, blue on both relevant sides with the same
positional reversal. Thus neither availability nor blackness of an
arc-disjoint replacement can be silently assumed.

`detours_square.py` independently checks the exact colour law for every
standard-blue presentation of orders 3..7 and every arc-disjoint Hamiltonian
path, including pure sources and loops. It includes both successful detour
types and the nonautomorphic-blue-tail hostile. All validation uses explicit
checks, not Python assertions.

Reproduction, from the repository root:

```
python 05-knowledge/results/strong_blue_games_20261020f/detours_probe.py
python 05-knowledge/results/strong_blue_games_20261020f/detours_square.py
python -O 05-knowledge/results/strong_blue_games_20261020f/detours_square.py
python 05-knowledge/results/strong_blue_games_20261020f/detours_family.py
python -O 05-knowledge/results/strong_blue_games_20261020f/detours_family.py
```

**OPEN:** guaranteeing an appropriate two-block repair, square, or hub for
every strong-to-strong mixed blue line. The T_k corollary closes an explicit
infinite residual family; it does not settle that universal existence claim.

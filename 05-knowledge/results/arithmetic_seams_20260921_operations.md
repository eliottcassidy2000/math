# Shifted operations, diagonal forests, and the missing remainder coordinate

**Date:** 2026-09-21. **Status:** PROVED elementary statements; FINITE-EXACT
controls. The user's f and g do not yet have unambiguous definitions.
The operations below are explicitly proposed precise models, not a claim
to have uniquely decoded that notation. No Collatz convergence result or
universal explanation of exceptional prime patterns is asserted.

## Inheritance and concept board

The equal-parent deletion mechanism is already proved in
[THM-2422, operation-fibres-summand-closure-and-twin-center-ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md),
[THM-2433, operation-fibre-deletion-incidence-and-startup-scar](../../01-canon/theorems/THM-2433-operation-fibre-deletion-incidence-and-startup-scar.md),
and the [2026-09-17 summand synthesis](arithmetic_braids_20260917_summand.md).
We use it rather than relabeling it as new. The hostile is an operation
shadow that loses its companion or ambient relation. The corrected near
miss is treating graph complement, arrow reversal, and inverse operation
as the same construction. The underused coordinate here is the remainder
missing from a divisibility incidence address.

The [hyperoperation-grid hypothesis HYP-3087](../hypotheses/HYP-3087-lrc14-hyperoperation-grid-address.md)
is an operation-address synthesis, not a proved cross-domain theorem.
The present work retains its useful instruction to specify a coordinate
map and its lost data, but proves the following claims independently.

| Board object | Operation / retained coordinate | Cheap hostile |
|---|---|---|
| Shifted product | Transport by h(a)=a+1 | Identity is 0; additive zero requires -1 |
| Additive shadow | Distinct parent pairs and ascending ambient arcs | Allow equal parents and the missing forest vanishes |
| Multiplicative shadow | Divisibility ambient arcs, unit loops removed | Pair(2,3) is not in either direction of divisibility |
| Plane address | Parent pair versus arc pair | Multiplicative incidence has zero planar density |
| Residue scaling | Both multiplier and modulus | Multiplication by m annihilates every residue modulo m |
| Sequence diagonal | Affine conjugacy with doubling/squaring | Prime support is not preserved by adding one |

Anchor: repair the operation/graph identities. Niche: distinguish genuine
plane coordinates from an embedding of a graph. Wildcard: make the
micro-remainder/macro-quotient idea exact through Euclidean division.

## 1. A coherent arithmetic loop obtained by shifting the origin

On `D={-1,0,1,2,...}` define

```text
a boxplus b = a+b+1,
a boxtimes b = ab+a+b.
```

The bijection `h:D->N_0`, `h(a)=a+1`, gives

```text
h(a boxplus b)=h(a)+h(b),
h(a boxtimes b)=h(a)*h(b).                                (O1)
```

**PROVED.** These operations are associative, commutative, and distributive;
their additive zero is -1 and multiplicative identity is 0. Their entire
semiring structure is ordinary nonnegative arithmetic transported through
h. Restricting to nonnegative a,b keeps both operations closed but removes
the additive identity. A new symbol is not a new arithmetic invariant.

There is nevertheless a useful combinatorial realization. For nonnegative
integers a,b, pointed sets with a,b nonbasepoint elements have a+1,b+1
elements. The product of two pointed
sets therefore has `ab+a+b` nonbasepoint elements, naturally partitioned
into an a-axis, a b-axis, and an ab-interior. This gives the exact product
underlying the [tournament edge decomposition](arithmetic_seams_20260921_tournaments.md).
The distinguished root is a real coordinate, not a negligible endpoint.

Both diagonal maps now occur in one explicit structure:

```text
a boxplus a = 2a+1,
a boxtimes a = a^2+2a,
h(2a+1)=2h(a),              h(a^2+2a)=h(a)^2.             (O2)
```

Starting from 0, repeated boxplus-diagonal gives `M_j=2^j-1`.
Starting from 1, repeated boxtimes-diagonal gives `2^(2^r)-1`.
Consequently the Fermat numbers satisfy

```text
F_r=M_(2^r)+2.
```

This relates the two sequences by a precise coordinate change and index
subsequence. It does **not** transport primality or primitive-divisor
behavior through the added constant. Those require the separate
[prime-factor analysis](arithmetic_seams_20260921_primes.md).

Repeated boxtimes also has a well-defined exponentiation:
`a^[r]=(a+1)^r-1`, with zeroth power 0. For example
`a^[r+s]=a^[r] boxtimes a^[s]`. This is one coherent interpretation of
an arithmetic hierarchy containing addition inside multiplication and
exponentiation. It does not satisfy the user's literal self-scaling
equation unless the domain or equality relation is changed.

## 2. Complements of the operation graphs: the inherited precise statement

Use positive integer vertices and suppress loops. Inside the ascending
ambient relation `0<x<z`, the distinct-summand graph is

```text
G_plus={(x,z):0<x<z, z!=2x}.
```

Its missing arcs are exactly `x->2x`. The omitted companion would be
`z-x=x`, so this is precisely deletion of equal-parent witnesses.

For multiplication the ambient relation changes to proper divisibility:

```text
D={(x,z):0<x<z, x divides z},
G_times={(x,z) in D:z!=x^2}.
```

The complement **inside D** is `x->x^2`, x>=2. It is not the complement
inside all ascending pairs: for example (2,3) is absent for lack of any
integer cofactor, rather than because it lies on the square diagonal.
The raw unit pair {1,z} also supplies a loop at z; suppressing loops is
necessary for this stated simple-graph convention.

The removed arcs form forests: odd roots for doubling, and nonsquare roots
at least 2 together with isolated vertex 1 for squaring. Multiplicative ancestry follows by writing
`n=product p^e_p` and extracting the largest common power of two from
the exponents. These are exact instances of swap-fixed parent deletion,
already proved in the cited work.

## 3. Why these graphs do not embed as planar tilings

**PROVED.** Each underlying undirected operation graph contains a copy
of K_m for every positive m.

For addition take the m vertices `m+1,...,2m`. No two have ratio two,
so every pair is an edge. For multiplication fix a prime p and take
`p^(3^0),...,p^(3^(m-1))`. Each smaller element divides every larger
one, and no larger exponent is twice a smaller exponent. Thus no edge
is removed by the square diagonal. Both graphs consequently contain
every finite simple graph as a subgraph after edge deletion, and hence
as a minor. Unrestricted graph-minor membership has no discriminatory
power between these full graphs or between their number-theoretic labels.

Small explicit nonplanarity witnesses are already available:

```text
G_plus:  K_(3,3) between {1,2,4} and {3,5,6};
G_times: K_(3,3) between {1,2,3} and {6,12,18}.
```

All nine cross edges exist in each case. K_(3,3) is nonplanar by the
Euler bound `e<=2v-4` for simple bipartite planar graphs: 9>8.
Therefore neither full graph admits a crossing-free embedding in the
plane. This does not preclude drawing its *incidence matrix* in the
plane; a two-dimensional array and a planar graph embedding are distinct
objects. Nor does ordinary undirected minor theory preserve directed
Hamiltonian paths or the numerical labels of operation fibres.

The removed doubling and squaring forests are planar, but passing to
them discards every unequal-parent incidence. Their planarity says
nothing by itself about the original graph or a selected Collatz orbit.

## 4. There is also an exact density obstruction to a finite-copy lattice tiling

Interpret an arc (x,z) as a lattice position. In `[1,N]^2`, the counts are

```text
|G_plus| = binomial(N,2)-floor(N/2),
|G_times| = sum_(d=1)^N floor(N/d)-N-floor(sqrt(N))+1.      (O3)
```

The first graph has asymptotic density one half in the positive square.
The second has `O(N log N)` positions and thus density zero. Its count
comes from all divisibility pairs, minus N diagonal pairs, minus the
`floor(sqrt(N))-1` removed square arcs.

**PROVED.** Finitely many fixed invertible affine-linear images of
G_times, its reversed arcs, and its removed square forest cannot cover
the integer lattice plane. The preimage of a box of side O(N) under
each fixed inverse affine map lies in a box of side O(N). Each image
therefore contributes only O(N log N) positions; a finite union still
has o(N^2) positions. This theorem concerns lattice incidence sets and
fixed affine copies, not arbitrary nonlinear transformations or an
infinite family of row-dependent copies.

If “inverse graph” instead means the complement inside *all* lattice
positions, coverage is true by definition, but that complement is not
just the square forest. Naming the ambient universe is indispensable.

## 5. The strongest tiling repair retains the parent or the remainder

Every positive parent pair (a,b) is a valid input for both operations.
The square grid of parents therefore has two total evaluation maps:

```text
(a,b) -> (a,a+b),
(a,b) -> (a,ab).                                         (O4)
```

The first is an integer shear with determinant one and inverse
`(a,z)->(a,z-a)`, on its positive wedge. The second has inverse
`(a,z)->(a,z/a)` only on divisibility positions. Its Jacobian is a,
and its row-a output is spaced in multiples of a. This is exactly where
the claimed common plane tiling loses information.

Restoring the remainder gives a total address for every a>=1,z>=0:

```text
z=a*b+r,             b>=0, 0<=r<a.                       (O5)
```

Existence and uniqueness are Euclidean division. For a fixed row a,
the a possible remainders translate its multiplication sites to cover
all nonnegative targets. The number of copies depends on a, consistently
with the finite-copy obstruction. This is a rigorous version of the
user's bounded microscopic remainder and unbounded quotient coordinates.
For a fixed modulus m=AB, write `x=m*N+q`, `0<=q<m`. The quotient N
is unbounded; it is not another finite microscopic coordinate.

## 6. Repairing the self-scaling equation gives modular dynamics

Over integers or rationals, `q=m*q` with m!=1 forces q=0. Over
`Z/MZ`, with integer M>=1, however, the exact classification is

```text
q=m*q  iff  q=j*(M/g), 0<=j<g, where g=gcd(m-1,M).        (O6)
```

To prove it, divide `M | (m-1)q` by g; the remaining multiplier and
modulus are coprime. There are exactly g solutions. For m>=2 and r>=1:

- with modulus M=m, only q=0 survives;
- with modulus M=m-1, every residue survives;
- with modulus M=m^r-1, multiplication by m has period dividing r
  on every residue, because m^r=1 modulo M.

For base ten and M=10^r-1, the last operation rotates r-digit words,
allowing leading zeroes. The all-nine word and the all-zero word
represent the same residue, an explicit quotient ambiguity that must
not be discarded. This is a genuine cyclic arithmetic action behind
digit patterns. The [append-digit prime analysis](arithmetic_seams_20260921_primes.md)
uses modular orders for a different affine recurrence; it retains the
starting phase needed to decide which index is divisible by a prime.

## 7. What can and cannot be transported

The coordinate shift h transports algebraic equalities, not primality.
The parent-to-arc map transports an operation incidence, but multiplication
requires a remainder to reach the full target grid. The graph shadow
forgets the operation's numerical state when passed to an unlabelled minor.
The residue quotient transports divisibility tests, not ordinary height.
Each of these maps is useful precisely because its loss is now explicit.

**OPEN research direction:** build a labelled, directed carrier retaining
both the Euclidean remainder and a selected dynamics' height change, then
ask for a graph operation that preserves its descent predicate. The full
operation graphs are minor-universal, while their diagonal forests are
too small to retain this predicate. A workable reduction would have to
specify a restricted class and a lawful contraction operation. No such
global Collatz carrier is constructed here.

## 8. The append-digit gcd law also applies within each inverse Collatz fibre

The incoming concurrent trunk note at origin/main commit78d2b2c7c
suggested retaining an explicit common-image test. The inverse-fibre
mechanism itself is inherited from the
[2026-09-17 Collatz note](arithmetic_braids_20260917_collatz.md);
the following elementary calculation needs no claims from the incoming
note's other investigations.

For any positive odd x_0 and j>=0 put

```text
x_j=4^j*x_0+(4^j-1)/3.
3x_j+1=4^j*(3x_0+1).                                  (O7)
```

All x_j are positive odd and have the same image under accelerated
Collatz; their halving exponents differ by2j. Moreover, for i>=0,t>=1,

```text
gcd(x_i,x_(i+t))=gcd(x_i,(4^t-1)/3).                    (O8)
```

Indeed `x_(i+t)=4^t*x_i+(4^t-1)/3`; subtract its first term inside the
gcd. Thus consecutive members of this inverse family are coprime. This
is the same repunit-return mechanism used for append-digit families,
now on an actual Collatz carrier. The members are **siblings** having a
common target, not consecutive nodes along a forward trajectory. Hence
(O8) gives no successive-orbit factor novelty or descent theorem.

## Reproduction

Run `python 04-computation/experiments/arithmetic_seams_20260921_operations.py`
and the same command with `python -O`. The [script](../../04-computation/experiments/arithmetic_seams_20260921_operations.py)
uses exact integer arithmetic. Its [JSON](../../04-computation/experiments/arithmetic_seams_20260921_operations.json)
declares finite universes: semiring variables -1..12; graph cutoffs 1..160;
clique controls of size2..8; multipliers2..30 and moduli1..60; and
Euclidean addresses for A,B<=20,x<200; and15600 inverse-sibling gcd
controls from100 positive odd starts up to199. The all-size statements above
are proved algebraically, not inferred from these checks.

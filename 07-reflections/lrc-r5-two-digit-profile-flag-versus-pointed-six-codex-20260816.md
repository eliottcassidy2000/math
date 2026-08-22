# The two-digit period-three pattern is a profile flag, not a `C_3` action

**Status: DERIVED FINITE-EXACT PROFILE-FLAG CONSEQUENCE OF THE INDEPENDENTLY
AUDITED TWO-DIGIT OBJECT AT `abf64fa44`; CANDIDATE AND CLEAN-ROOM PROFILES,
BOUNDARIES, AND OWNER-VISIBLE MASK AGREE ENTRYWISE; THE FOUR-WAY POINTED LIFT
REMAINS UNDER AUDIT; NO PHYSICAL CURRENT AND NO LRC(14) CONCLUSION.**

The conditional endpoint ranks

```text
(4,3,3,4,3,3,4,3,3,4,3,3,4)
```

do not come from a three-cycle representation on `F_13`.  Before endpoint
integration, the exact nested-window profiles form a reflection-stable flag

```text
L_3 inside M_6 and B_10 inside G_12,
M intersect B = S_6,
```

where `S_r` is the span of the thirteen `r_1` profiles at fixed `r_0=r`.
Every interior nonmultiple of three gives the same `L`.  Each of
`S_3,S_6,S_9` is `L` plus one line, and those three quotient lines are
independent.  Endpoint integration preserves each individual interior slice
but maps the three extra lines onto one line.  Its two-dimensional kernel on
`M/L` is one trivial plus one sign direction.

This gives the precise verdict on rank six.  The six-dimensional endpoint
relation image is genuinely the already identified pointed quotient `W`, but
six is not a ceiling on the two-digit ancestry carrier.  The owner-visible
weighted profile carrier has rank twelve, and its interior six-space has state-block
ranks `(1,3,1,3)`, not the pointed grading `(1,2,1,2)`.  Thus six is a genuine
**projected response channel**, not a global or pre-integration carrier bound.

The incoming nested source/current result at `2d52215a3` makes this separation
decisive on a second joint ancestry object.  Its double-character profile
space has rank seventeen, while endpoint transport has rank four and kernel
thirteen; nevertheless its relation rows equal the pointed spaces statewise
with dimensions `(1,2,1,2)` and globally with rank six.  Rich ancestry
amplitudes and the pointed relation carrier are therefore two different tensor
axes, not competing estimates of one dimension.

## Typed objects and map

Work over the same split field `k=F_p` as the exact endpoint computations.
Let `C` be the 33 open source cells cut out by the 34 pinned profile
boundaries.  The independently reconstructed chamber/state predicate marks
16 of them owner-visible; the remaining 17 coordinates are set to zero.  Let

```text
P_6=((0,0),(1,0),(1,6),(3,6),(3,12),(2,12))
```

be the six realized `(Boolean state,source root)` incidences.  For each
two-digit address `(r_0,r_1)`, define

```text
p_(r_0,r_1) in A=Fun(C x P_6,k)
```

to be its exact weighted nested-window source profile after this visibility
projection, zero outside the matching state block.  The candidate and
clean-room profile arrays, boundaries, and masked cell tuples agree entrywise;
their common masked-cell hash is
`56fe6e86a6d49ae4d7bdb0cdde5a5b70b0b2d8364c06dc7c22c7d8fafa3df4f5`.
Put

```text
S_r = span{p_(r,r_1): r_1 in F_13}.
```

The endpoint pipeline is linear in this left profile once the owner base,
right profile, literal guards, and endpoint factors are fixed.  It therefore
defines a linear response map

```text
E : G -> Fun(V_4(state) x F_13(relation),k),
G = sum_r S_r.
```

This `E` is the combined-address flattening.  The relation flattening is its
transpose-shaped companion: it spans response rows over relation while
state and address index the columns.  Their ranks `4` and `6` refer to
different typed flattenings and must not be conflated.

The six `P_6` entries are path incidences.  The digits `(r_0,r_1)` are radix
tree branches.  An ordered root pair `(u,q)` is a cut arc and
`s=u-q` is its colour.  The reflection below is a path/radix involution.
None of these sets is equipped here with a complete antisymmetric pairwise
relation, so no tournament is present.

## Finite profile-flag theorem

Use the canonical integer representatives `0,...,12` and define

```text
U_1={1,4,7,10},  U_2={2,5,8,11},
A_3={3,6,9},     D={0,12}.
```

The exact rank certificates imply the following finite theorem on the common
candidate/clean-room owner-visible profile table.

1. For every `r` in `U_1 union U_2`, the spaces `S_r` are one common space
   `L` of dimension three.

2. Each of `S_3,S_6,S_9` has dimension four and contains `L`.  Every pair has
   union rank five, while

   ```text
   M=S_3+S_6+S_9
   ```

   has dimension six.  Hence the three lines `S_r/L` are independent and
   `dim(M/L)=3`.

3. The boundary space `B=S_0+S_12` has dimension ten; each boundary slice has
   dimension six.  One has `L subset B`, `dim(M+B)=12`, and

   ```text
   M intersect B = S_6.
   ```

   Indeed `B+S_6` has rank ten, whereas `B+S_3` and `B+S_9` have rank eleven;
   also `dim(M intersect B)=6+10-12=4=dim S_6`.  Each individual boundary
   slice meets `L` in dimension two, even though their sum contains all of
   `L`.

4. The full weighted profile space `G=M+B` has dimension twelve.  Replacing
   every nonzero amplitude by one reduces the full rank to three.  Amplitudes,
   not support geometry alone, carry nine of the twelve directions.

The source ranks by fixed `r_0` are therefore

```text
(6,3,3,4,3,3,4,3,3,4,3,3,6).                 (1)
```

Equation (1) already separates the two boundary fours in the endpoint record
from the three interior fours: the former begin as six-dimensional slices
and lose two dimensions, while the latter begin and remain four-dimensional.

## Exact reflection representation

The typed reflection is

```text
(cell,state,u,r_0,r_1)
 -> (reversed cell,state XOR 2,12-u,12-r_0,12-r_1).
```

On the six pointed incidences it is the three-pair reversal

```text
(0,1,2,3,4,5) -> (5,4,3,2,1,0).
```

Since the field has odd characteristic, every `C_2` module splits into plus
and minus eigenspaces.  Exact ranks give

```text
             dimension   plus   minus
L                 3         2      1
M                 6         3      3
B                10         5      5
G                12         6      6
S_6               4         2      2
S_3+S_9           5         3      2.
```

Consequently

```text
M/L = 1 trivial + 2 sign.
```

The paired lines `S_3/L` and `S_9/L` contribute one trivial and one sign
direction.  The reflection-fixed digit `r_0=6` contributes the remaining
**sign** line `S_6/L`.  This is the exact representation-theoretic mechanism
behind the three interior rank jumps.

There is also an exact hostile to an untyped identification with the pointed
permutation module.  Projection of the profile spaces to the four Boolean
state blocks gives

```text
                 state 0  state 1  state 2  state 3
L                   1        2        1        2
M                   1        3        1        3
B                   3        3        3        3
G                   3        4        3        4.
```

The pointed carrier itself has state grading `(1,2,1,2)`.  Although `M` and
the six-pointed permutation module both have abstract character `(3,3)`, `M`
has a third state-local direction in each doubleton state.  No lossless
state-graded identification of these two six-spaces exists.  Character
agreement is therefore necessary but not a typed factorization.

## Endpoint projection explains the visible `4,3,3`

Now add the independently audited endpoint facts from `abf64fa44`.  The conditional
dimensions `dim E(S_r)` are

```text
(4,3,3,4,3,3,4,3,3,4,3,3,4),                 (2)
```

and `dim E(G)=4`.  Comparing (1) and (2) proves:

```text
E is injective on S_r for 1 <= r <= 11;
dim ker(E|S_0)=dim ker(E|S_12)=2;

dim E(L),E(M),E(B),E(G) = (3,4,4,4);
dim ker(E|L),ker(E|M),ker(E|B),ker(E|G) = (0,2,6,8).
```

In particular, the three independent lines of `M/L` map onto the single line
`E(M)/E(L)`.  Reflection equivariance and injectivity on `S_6` show that the
surviving line is the sign line `S_6/L`.  Hence

```text
ker(E|M) = 1 trivial + 1 sign,
E(M)     = 2 trivial + 2 sign.                     (3)
```

Thus the endpoint map kills one regular two-point `C_2` packet among the
three interior extras.  Every multiple of three still displays rank four
because its own extra line survives, but the three extras are identified in
the global endpoint image.  This is a projection theorem, not a guessed
recurrence.

The boundary entries in (2) reach the same number four by a different
mechanism: each loses two of its six source-profile directions.  The apparent
fivefold repetition of `4` therefore splices an interior flag mechanism to a
boundary collapse.

## The incoming `17 -> 4` result gives the exact two-rank factorization

Let `C_17` be the selected basis of the complete double-nontrivial
`b_source x r_owner` profile sector from `2d52215a3`, and let `E_sr` be its
common endpoint operator.  The exact ranks state

```text
dim C_17=17,
dim im(E_sr)=4,
dim ker(E_sr)=13.                                  (4)
```

For each Boolean state `i`, let `P_i` denote the pointed response rows over
relation.  The exact union-rank certificates are

```text
(number of pointed rows, joint rank, pointed rank, union rank)
  state 0: (1,1,1,1)
  state 1: (2,2,2,2)
  state 2: (1,1,1,1)
  state 3: (2,2,2,2),
  global:  (6,6,6,6).                              (5)
```

Equations (4)--(5) give a typed factorization

```text
C_17 --quotient by K_13--> A_sr,4
     --injective endpoint map--> direct_sum_i P_i,
```

whose image has channel dimension four but whose relation-mode span is the
entire pointed six.  In tensor language these are different flattenings:

```text
profile/channel rank = 4,
relation-carrier rank = 6.
```

A four-dimensional subspace of `V_4 x Fun(F_13,k)` can have relation-mode rank
six; there is no rank contradiction.  Calling the kernel in (4) a
“six-carrier kernel” would be a type error: it is the kernel of the
17-dimensional source-profile amplitude map.

This new result explains how to read the two-current conditionals.  Their
fixed-`r_0` ranks `3` or `4` measure how each radix cylinder samples a
four-dimensional **amplitude quotient**.  They do not measure the dimension of
the relation carrier.  The source profile flag says:

```text
nonmultiple interior cylinders -> one common amplitude hyperplane (rank 3),
multiple interior cylinders    -> add one line and reach rank 4,
boundary cylinders             -> start at rank 6 but project to rank 4.
```

Meanwhile the relation flattening remains the same pointed six.  The nested
source/current package supplies an exact statewise realization of this
`large profile -> four amplitudes inside six relation rows` architecture, with
`17 -> 4` replacing the two-current package's `12 -> 4`.

What is not yet proved is that the two rank-four output subspaces are the same
subspace of `Fun(V_4 x F_13,k)`.  Equality of their dimensions is insufficient.
Once both endpoint bases are available, the cheapest discriminator is one
RREF stack:

```text
rank(A_two-current + A_source-current).
```

Rank four proves a common amplitude quotient; rank five through eight proves
distinct four-charts inside the shared pointed relation carrier.  This basis
comparison needs no new endpoint sweep or spectral census.

## Why this is not a `C_3` or translation representation

Multiplication or addition by three in `F_13` cannot produce a three-cycle.
Translation by three is one 13-cycle:

```text
0,3,6,9,12,2,5,8,11,1,4,7,10.
```

Along this single orbit the profile and endpoint ranks are respectively

```text
(6,4,4,4,6,3,3,3,3,3,3,3,3),
(4,4,4,4,4,3,3,3,3,3,3,3,3).
```

An equivariant group action would preserve rank around its orbit, so both
records are exact hostile counterexamples.  “Divisible by three” refers only
to the ordered representatives `0,...,12`; it is not a subgroup, quotient,
or character of `F_13`.  Carries and the two end strata are load-bearing.

The honest interpretation is therefore a radix-cylinder flag on an ordered
interval, not a period-three recurrence.  A recurrence still requires a
typed operator at a third digit or a lawful clock transport.

This is not a no-go theorem for general autonomous recurrences.  As repaired
by the independent audit, the fixed-`r_0` rank test by itself excludes only a
scalar state/relation-independent factor `K(r_0,r_1)T_1`; hidden-state,
matrix-valued, state-dependent, nonlinear, and temporal transitions remain
untested.  The flag theorem additionally refutes the specific `C_3`/translation
explanation, not those broader mechanisms.

## Where the genuine rank six lives

Let `R_2` be the span of the two-digit endpoint rows over the relation
coordinate.  Its `r_1` marginal is the audited one-digit relation space `W`.
Marginalization makes `W` a subspace of `R_2`, while the exact ranks are

```text
dim W = dim R_2 = 6.
```

Therefore `R_2=W`.  The prior pointed audit identifies

```text
W = mu(P_6),
```

where `mu` sums source-root difference and is injective on the six-pointed
response space.  The second current digit adds no endpoint relation direction
after this marginalization.  That equality makes six a genuine typed quotient
channel rather than a numerical coincidence.

The nested source/current result strengthens this from a global equality to a
statewise one for a different lawful fibre product: its joint rows equal
`P_i` in every state and span all six pointed rows globally.  It is therefore a
positive exact instance of pointed closure **after** endpoint transport.
Because that probe also marginalized `s=u-q`, it does not settle closure in
the expanded `(difference,relation)` target.

It does **not** make six a pre-integration ceiling.  The present owner-visible
profile space has dimension twelve; its state grading is incompatible with `P_6`; and
the independently audited unrestricted root-difference response has relation
rank thirteen raw and twelve centred.  The strongest sound statement is:

> rank six is the pointed endpoint quotient after root-pair and profile
> marginalization; it is not a global ancestry-carrier ceiling.

Whether the two-digit response already factors through the pointed bundle
while source-root difference is retained remains open because `e9131a159`
aggregates away the right root before address expansion.

## Cheapest decisive lifted test

Do not repeat the full joint Fourier census.  Reuse the compressed endpoint
event pass and change only the information discarded by its aggregation.
The current code stores

```text
key=(cell_index,selected_u_mask),
scalar=(sum_q v_q)*jump.
```

This destroys `q`, so `s=u-q` cannot be reconstructed.  Instead store

```text
key=(cell_index,selected_mask),
scalar=jump,
```

or an equivalent `q`-resolved vector.  During address expansion, emit

```text
address_value[u,r_0,r_1] * v_value[q] * scalar
```

for every selected `u,q`, binned by `s=u-q`.  Work either in relation space
or in the invertibly equivalent `tau`-frequency space.

For each state `i`, precompute an exact RREF annihilator of the pinned pointed
space

```text
P_i subset Fun(F_13(s) x F_13(t),k),
dim(P_i)=(1,2,1,2).
```

Stream the `169` address rows for that state through the annihilator and stop
at the first nonzero residual.  No four-way Fourier census, projective
marginal atlas, or storage of Gibbs's full joint tensor is needed.

Mandatory controls are:

```text
sum_s          -> the pinned two-digit tensor at e9131a159;
sum_r1         -> the one-digit current-branch x root-difference tensor;
sum_r0,sum_r1  -> the audited Boolean-square x root-difference tensor;
same-root, literal guards, and reflection -> their pinned values.
```

The outcome is decisive:

- if all residuals vanish, every depth-two row lies statewise in `P_i`; the
  one-digit marginal already spans `P_i`, so this constructs the desired
  pointed-bundle factorization and proves a six-dimensional ceiling for this
  exact lifted response;
- the first nonzero residual gives a minimal counterexample
  `(state,r_0,r_1,s,t,pivot)` and proves that the observed relation rank six
  arose only after lossy root-difference marginalization.

This is the cheapest exact discriminator between a genuine depth-two pointed
carrier and a projection artifact.

## Reproduction and exact ledger

The profile-only probe deliberately does not execute the endpoint event
sweep, character bank, relation inversion, or spectral census.  It imports
both source reconstructions and first requires exact equality after the
independently derived owner-visible mask:

```text
python -B 04-computation/lrc_r5_two_digit_period3_pointed_profile_probe_20260816.py
python -B -O 04-computation/lrc_r5_two_digit_period3_pointed_profile_probe_20260816.py
```

Its weighted/support profile hashes are
`f7d1ac1bba79a25b232dd4f0539b9a236f047b5224e49050a44c70f4d2544b68`
and `6fd60baa0f82c5f234164b2869cb369468fcbcf128e995603c49e5c07ed7ea68`;
the semantic hash is
`0d5798935986c13b1aa70fc1db7abac162994dabaf7dd3ab5f323ffeb1d1d63e`.

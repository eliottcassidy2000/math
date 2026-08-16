# The common rank-six quotient is the six-pointed carrier

**Status: FINITE-EXACT REPRESENTATION SYNTHESIS ON THE CANONICAL `r=5` OWNER
BASE; THE FOUR-WAY PHYSICAL CANDIDATE STILL AWAITS INDEPENDENT AUDIT;
POINTED-BUNDLE CLOSURE REMAINS A HYPOTHESIS; LRC(14) remains OPEN.**  The
source-time high digit `b_source=floor(n/13)` and the current-leg last digit
`r_owner=a mod 13` do not merely give two tables of rank six.  Their exact
endpoint-relation row spaces are equal.  The common six-space carries the same
reflection in both presentations and decomposes as

```text
W = 3 trivial + 3 sign.
```

The incoming current-branch x root-difference table supplies the missing
carrier.  Statewise, its branch rows span exactly the corresponding rows of the
six pointed pairs, and difference marginalization is injective on that pointed
space.  Thus the common endpoint six-space is the image of the six-pointed
carrier, whose path reversal consists of three transpositions.  Its
four-dimensional Boolean-square parent is `2 trivial + 2 sign`, so the branch
refinement adds exactly one even and one odd response direction.

This does not identify the two radix digits.  Equality is obtained only after
quotienting two different 52-dimensional radix/state carriers by their
different 46-dimensional response kernels.  In particular, the source-time
digit has not yet been factored through the pointed carrier before root
difference is marginalized.

The broad statement that every one-owner endpoint response has rank at most
six is false.  The already audited root-difference and pointed-tail tables on
the same base have raw relation rank thirteen and centred relation rank twelve.
Six is therefore a typed channel dimension, not a global endpoint ceiling.

## Inheritance pass

The closest finite mechanism is the exact current-sheet reflection from the
reflection frozen at `origin/main` commit `085a3b7ab`:

```text
1-X_(u,a)(y) = X_(12-u,13^5-1-a)(1-y).
```

On the retained last digit and owner-visible state it induces

```text
(r_owner,state) -> (12-r_owner,state XOR 2).
```

The concurrent independently audited source sheet writes `n=13b+r`.  Source
reflection sends `n` to `168-n`, hence

```text
(b_source,r) -> (12-b_source,12-r).
```

After the low digit is marginalized, it has the same formal reversal on
`(b_source,state)`.  The canonical hostile is MISTAKE-417: full Fourier support
can be a rank-one delta-cell lift.  The present test therefore compares exact
row spaces, kernels, and involution graphs rather than support counts.

The decisive incoming mechanism is the four-way current-branch x
root-difference candidate at `origin/main` commit `80d7fed8e`, merged into the
rolling frontier at `23b6e7981`.  Its pointed parent is independently audited;
the four-way integrand itself is not.  The least-used next sidecar is the
source-time branch jointly retained with root difference.

## Typed response maps

Let `k=F_p` for

```text
p = 755373809845391722745761,
H = Fun(F_13(relation),k).
```

The two domain carriers are

```text
D_source = k[F_13(b_source) x V_4(state)],
D_owner  = k[F_13(r_owner)  x V_4(state)].
```

They have the same dimension `52` but different temporal meanings.  The exact
endpoint integrations and relation inversion define

```text
Phi_source : D_source -> H,
Phi_owner  : D_owner  -> H.
```

The image of a basis vector is the corresponding 13-entry endpoint-response
row.  Summing the radix label in either domain gives the same audited
four-state Boolean-square map with image `W_0`.

This is the relevant use of “representation.”  The carrier points in `H` are
endpoint relation residues `(1,0,t)`.  They are not radix branches, path
vertices, source roots, cut arcs, or tournament contestants.

## Exact common-image result

The exact rank ledger is

```text
dim im(Phi_source)                         = 6,
dim im(Phi_owner)                          = 6,
dim(im(Phi_source)+im(Phi_owner))          = 6,
dim W_0                                    = 4.
```

Therefore

```text
im(Phi_source) = im(Phi_owner) =: W,
dim(W/W_0) = 2.
```

The equality is not inferred from dimensions.  Independent reduced-row-echelon
bases of the two `52 x 13` tables have the same SHA-256

```text
6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4.
```

The source table digest is the one reproduced by the concurrent clean-room
audit.  The current table is rebuilt from the disjoint `13^4`-fold then
`r_owner`-window construction.  Both exact parents agree entrywise before any
rank comparison.

## The incoming four-way table identifies the carrier

Let

```text
K = Fun(F_13(source difference) x F_13(relation),k)
```

and let `P_6` be the six realized pointed cut-state incidences

```text
((0,0),(1,0),(1,6),(3,6),(3,12),(2,12)).
```

The incoming joint tensor retains `(state,r_owner,u-q,t)`.  For each fixed
state `i`, compare in `K`

```text
M_i = span of its thirteen r_owner rows,
P_i = span of the pointed rows (i,u).
```

The exact certificates

```text
state 0: (number of points,rank M_i,rank P_i,rank(M_i+P_i)) = (1,1,1,1),
state 1:                                                        (2,2,2,2),
state 2:                                                        (1,1,1,1),
state 3:                                                        (2,2,2,2)
```

prove `M_i=P_i` statewise.  Hence the current-branch image in `K` is exactly
the six-dimensional pointed space `P=P_0+P_1+P_2+P_3`.

Now let

```text
mu : K -> H,       mu f(t)=sum_s f(s,t),
```

be difference marginalization.  A separate derived audit obtains

```text
rank P                         = 6,
rank mu(P)                     = 6,
rank im(Phi_owner)             = 6,
rank(mu(P)+im(Phi_owner))      = 6.
```

Thus `mu|P` is injective and

```text
mu(P) = im(Phi_owner) = W.
```

The canonical RREF digest of `mu(P)` is exactly the previously pinned common
endpoint digest

```text
6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4.
```

Since `im(Phi_source)=im(Phi_owner)`, the source-time table also reaches
`mu(P)`.  This last equality occurs after difference marginalization; it does
not yet construct a source-time map into `P` inside `K`.

## The reflection descends, but not as a residue permutation

Let `J_source` and `J_owner` reverse the path state and the appropriate radix
digit:

```text
J(d,state) = (12-d,state XOR 2).
```

For a map `Phi`, `J` descends to the image exactly when `ker(Phi)` is
`J`-stable.  The audit tests this without choosing a basis: horizontally join
the response matrix with its `J`-permuted copy.  For both maps the joined rank
remains six.  Thus there are induced involutions

```text
j_source, j_owner : W -> W.
```

Their graph spaces also agree: stacking the two six-dimensional graphs still
has rank six, with common canonical graph digest

```text
67850c2586ac9a4f242ada7173eba2fe1b3c5e6d1ae20a5e2b6d7ad0dc7d3f88.
```

Hence `j_source=j_owner=:j`.  Since `p` is odd, the `C_2` representation is
semisimple, and the exact plus/minus ranks give

```text
W^+ = 3,       W^- = 3,
W_0^+ = 2,     W_0^- = 2,
(W/W_0)^+ = 1, (W/W_0)^- = 1.
```

The pointed path reverses by

```text
(0,1,2,3,4,5) -> (5,4,3,2,1,0).
```

Its permutation module is therefore three two-point permutation modules,
namely `3 trivial + 3 sign`.  The derived audit compares involution graphs and
gets joint graph rank six with the same canonical graph digest above.  Hence
the abstract quotient involution `j` is exactly pointed-path reversal transported
through the injective marginal `mu|P`; it is no longer an unexplained kernel
symmetry.

Likewise the four Boolean parent states form two reversal pairs, giving
`2 trivial + 2 sign`.  The quotient `W/W_0` is the one remaining pointed
reversal pair, `1 trivial + 1 sign`.  This is the representation-theoretic
content behind the repeated rank jump `4 -> 6`.

There is an important hostile qualification.  Under every affine action

```text
t -> c*t+d,       c != 0,
```

on the thirteen endpoint residues, the stabilizer of `W` is only the identity.
No affine reflection preserves the raw or pure six-space, and no entrywise
identity is obtained by combining path/radix reversal with an affine target
map and a sign.  Thus `W` is not “the six odd Fourier modes,” and `j` is not a
cosmetic permutation of endpoint labels.  It is an abstract quotient action
created by the response kernels.

## Why the pure three-way rank is also six

The two pure branch/state/relation ANOVA tables also have one common row space
`W_pure`:

```text
dim W_pure(source) = dim W_pure(owner) = 6,
dim(W_pure(source)+W_pure(owner))      = 6.
```

Their common canonical row-space digest is

```text
ae86ad2a3fd03bea95c823b2454b78f244581aa048cb2da63a03a75f484cc596.
```

Let `Q` subtract the mean on the relation axis.  The constant vector is not in
`W`, so `Q|W` is injective.  Exact ranks show

```text
W_pure = Q(W),
dim W_pure = 6,
dim(W intersect W_pure) = 5.
```

Thus ANOVA does not reveal another independent six-space.  It transports the
raw quotient isomorphically into the relation augmentation space.  The raw and
pure spaces together have rank seven, exactly as this rank-one gauge change
predicts.

## Path, cut, tournament, and tree ledger

These four structures remain distinct:

| structure | typed objects | relation/operation | status in this audit |
|---|---|---|---|
| owner path | six pointed `(state,u)` incidences over four visible cut states | reversal of the six-term pointed path | the actual rank-six carrier |
| source cut arcs | ordered root pairs crossing a state cut; colour `u-q` | directed partial cut relation, with within-side pairs missing | retained in `K`, then marginalized injectively on `P`; rank-13 outside `P` |
| tournament | vertices with exactly one oriented edge for every pair | complete antisymmetric orientation | absent; neither radix labels nor endpoint residues carry this relation |
| radix tree | inverse-word children labelled by base-13 digits | parent/child extension and digit reversal | `b_source` and `r_owner` are children at different stages of different trees |

The path involution is not a tournament orientation.  The six pointed
incidences are not six tournament vertices.  A cut arc is not a path edge, and
a radix child label is not a pairwise relation.  Equal `F_13` cardinality does
not identify the source and current trees.

## Hostile boundary: the broad ceiling is refuted

Restoring source-root difference gives exact flattening ranks

```text
raw:   (state,difference,relation) = (4,12,13),
ANOVA:                              (3,12,12).
```

Marking the absolute tail gives

```text
raw:   (pointed,difference,relation) = (6,12,13),
ANOVA:                                (5,12,12).
```

These are on the same owner base and use the same endpoint relation carrier.
Consequently the statement

```text
all one-owner endpoint-response images have dimension at most six
```

is **REFUTED**.  The missing hypothesis is not reflection but the ancestry
quotient: root difference and absolute tail must already have been
marginalized.

## Two-digit update: six is a projected carrier, not the ancestry ceiling

The object at `e9131a159`, independently hostile-audited at `abf64fa44`,
retains two lawful current digits
`a=r_0+13r_1+169c`.  Its relation flattening has rank six and its `r_1`
marginal is the one-digit table with relation image `W`.  Marginal row spaces
are contained in the child row space, so equal dimensions prove that the
two-digit relation image is again exactly `W`.

The owner-visible profile carrier before endpoint integration is much larger.
The derived exact sidecar, now gated entrywise against the independent source
reconstruction,
`lrc-r5-two-digit-profile-flag-versus-pointed-six-codex-20260816.md` finds a
weighted profile space `G` of dimension twelve.  Its interior-multiple space
`M` has dimension six and the same abstract reflection character `(3,3)` as
the pointed carrier, but its state-block ranks are `(1,3,1,3)`, versus the
pointed grading `(1,2,1,2)`.  Thus there is no lossless state-graded
identification: the shared character is not a typed factorization.

Endpoint integration maps the three independent quotient lines at
`r_0=3,6,9` onto one odd line and kills one trivial plus one sign direction.
The combined-address rank four is therefore an explicit projection artifact.
The relation rank six remains a genuine pointed **quotient** because it equals
`W=mu(P)`, but it is not a pre-integration or global carrier ceiling.

The nested source/current candidate at `2d52215a3` supplies the complementary
positive datum.  Its complete double-character profile sector has dimension
seventeen, but the endpoint channel map has rank four and kernel thirteen.
At the same time, its relation rows equal the pointed rows statewise with
dimensions `(1,2,1,2)` and globally with record `(6,6,6,6)`.  Hence it factors
exactly as

```text
profile_17 -> amplitude quotient_4 -> direct_sum_state P_i,
```

while spanning all six relation directions.  The numbers four and six are
Tucker ranks of different tensor axes: four counts independent channel
amplitudes in `state x relation`, while six counts the relation-mode carrier.
The thirteen-dimensional kernel belongs to the profile-to-amplitude map, not
to the pointed relation map.

This also locates the two-current `3/4` conditional pattern.  It describes
sections of a four-dimensional amplitude quotient: nonmultiple cylinders see
one common rank-three hyperplane, interior multiples add one line, and the two
boundary profile six-spaces collapse to rank four.  It does not describe the
six-dimensional relation carrier, which remains `W`.

## Precise surviving conjecture

The earlier one-digit “endpoint ceiling” formulation is **SUPERSEDED** by a
typed factorization statement.

**Pointed-bundle closure conjecture (HYPOTHESIS).**  Fix the canonical `r=5`
owner base, literal guard order, endpoint factors, and relation inversion.  For
an inverse-tree digit or finite inverse word `w`, retain source-root difference
before integration and form the response rows over `(difference,relation)`.
For each Boolean state `i`, the `w`-resolved row space is contained in the
pointed-tail space

```text
P_i = span{ response(i,u,-,-) : (i,u) is one of the six pointed states }.
```

Equivalently, all such ancestry words act through the fixed six-dimensional
pointed bundle `P`.  The current last digit proves equality for one word level.
After difference marginalization this implies containment in `W=mu(P)`, but
that rank-six endpoint statement is a consequence, not the primitive claim.

The conjecture does not yet cover the source-time high digit: its endpoint
marginal equals `W`, but its row space over `(difference,relation)` has not been
constructed.  The two-current-digit endpoint marginal now also equals `W`,
but that computation summed the right root before address expansion, so it
does not test the conjecture in `(difference,relation)`.  The nested
source/current fibre product does prove statewise pointed equality after
endpoint transport, but it too sums root difference.  Other clocks, rows,
primes, physical currents, and LRC(14) remain outside its scope.

## Cheapest decisive computation

The cheapest immediate mechanism test now reuses the compressed two-digit
event pass without repeating its Fourier census.  Its current aggregation key
`(cell_index,selected_u_mask)` multiplies by `sum_q v_q`, irreversibly losing
the right root.  Retain `selected_mask` and the unsummed jump instead; during
address expansion emit each selected pair `(u,q)` into `s=u-q`.

For every state, stream the `169` address rows over `(s,t)` through an exact
annihilator of the pinned pointed space `P_i`, whose dimensions are
`(1,2,1,2)`.  Stop at the first nonzero residual.  The mandatory marginals are

```text
sum_s          -> the pinned two-digit tensor;
sum_r1         -> the one-digit current-branch x difference tensor;
sum_r0,sum_r1  -> the audited square x difference tensor.
```

All residuals zero construct the depth-two pointed-bundle factorization; one
nonzero residual gives the minimal hostile
`(state,r_0,r_1,s,t,pivot)` and proves that root-difference marginalization
created the apparent ceiling.  This test needs neither the full joint tensor
ledger nor the projective Fourier marginals already computed by Gibbs.

The orthogonal next test remains the lawful source-time table
`T_source(state,b_source,u-q,t)`, with the same statewise annihilator.  It asks
for the missing source factorization rather than duplicating the two-digit
current probe.

## Connection contract

| field | exact answer |
|---|---|
| source | six pointed `(state,u)` incidences, plus two distinct radix/state presentations |
| expanded target | `K=Fun(F_13(u-q) x F_13(relation),F_p)` |
| endpoint target | `H=Fun(F_13(relation),F_p)` |
| map | pointed response into `K`, then the injective-on-`P` difference marginal `mu:P->W`; source/current endpoint maps land in the same `W` |
| preserved | owner base, Boolean state, literal guard order, pointed tail in `P`, cut colour in `K`, endpoint factors, path reversal |
| destroyed by `mu` globally | root difference outside `P`; the source/current presentations additionally lose 46 kernel directions and temporal stage data |
| positive gate | statewise current rows equal pointed rows; nested source/current gives profile `17 -> 4` with statewise pointed ranks `(1,2,1,2)` and global `(6,6,6,6)`; equal RREF and involution-graph hashes |
| hostile | affine target stabilizer is trivial; unrestricted root-difference/pointed tables have relation rank `13/12` |
| boundary | both joint positive results marginalize root difference; two-current and source-time factorizations in `K`, physical current, row exclusion, and LRC(14) remain open |

## Reproduction

```text
python -B 04-computation/lrc_r5_endpoint_response_rank_six_representation_audit_20260816.py
python -B -O 04-computation/lrc_r5_endpoint_response_rank_six_representation_audit_20260816.py
python -B 04-computation/lrc_r5_pointed_carrier_rank_six_representation_audit_20260816.py
python -B -O 04-computation/lrc_r5_pointed_carrier_rank_six_representation_audit_20260816.py
python -B 04-computation/lrc_r5_two_digit_period3_pointed_profile_probe_20260816.py
python -B -O 04-computation/lrc_r5_two_digit_period3_pointed_profile_probe_20260816.py
```

The common-image, pointed-carrier, and profile-flag semantic digests are
`7baeb128a4c4d5998342611fdcf821d002ffb4622692dd72c03e6f11c8d9825a`
`5e491136809fc164bbbfc7aeb9c272b7aff05020992985a39f291bf31297903e`,
and `0d5798935986c13b1aa70fc1db7abac162994dabaf7dd3ab5f323ffeb1d1d63e`.

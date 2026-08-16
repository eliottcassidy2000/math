# The two radix sheets have one six-dimensional endpoint quotient

**Status: FINITE-EXACT CROSS-AUDIT ON THE CANONICAL `r=5` OWNER BASE;
SCOPED ONE-DIGIT CEILING REMAINS A HYPOTHESIS; LRC(14) remains OPEN.**  The
source-time high digit `b_source=floor(n/13)` and the current-leg last digit
`r_owner=a mod 13` do not merely give two tables of rank six.  Their exact
endpoint-relation row spaces are equal.  The common six-space carries the same
abstract reflection in both presentations and decomposes as

```text
W = 3 trivial + 3 sign.
```

Its four-dimensional Boolean-square parent is `2 trivial + 2 sign`, so the
branch refinement adds exactly one even and one odd response direction.  This
does not identify the two radix digits: equality is obtained only after
quotienting two different 52-dimensional source carriers by their different
46-dimensional response kernels.

The broad statement that every one-owner endpoint response has rank at most
six is false.  The already audited root-difference and pointed-tail tables on
the same base have raw relation rank thirteen and centred relation rank twelve.
The only surviving ceiling is therefore explicitly scoped to one retained
radix digit after cut-arc colour and absolute tail have been marginalized.

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

The least-used nearby sidecars are the source-root difference and pointed
absolute tail.  They are decisive hostile boundaries because they already
break any unscoped rank-six ceiling.

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

This is the representation-theoretic content behind the repeated rank jump
`4 -> 6`: each radix presentation reaches the same additional trivial/sign
pair.

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
| owner path | four visible cut states, inherited from the five-state physical support path | reversal `state -> state XOR 2` | retained and used in `J` |
| source cut arcs | ordered root pairs crossing a state cut; colour `u-q` | directed partial cut relation, with within-side pairs missing | marginalized; exact rank-13 hostile when restored |
| tournament | vertices with exactly one oriented edge for every pair | complete antisymmetric orientation | absent; neither radix labels nor endpoint residues carry this relation |
| radix tree | inverse-word children labelled by base-13 digits | parent/child extension and digit reversal | `b_source` and `r_owner` are children at different stages of different trees |

The path involution is not a tournament orientation.  A cut arc is not a path
edge.  A radix child label is not a pairwise relation.  Equal `F_13`
cardinality does not identify the source and current trees.

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

## Precise surviving conjecture

**Scoped one-digit endpoint-image conjecture (HYPOTHESIS).**  Fix the canonical
`r=5` owner base, literal guard order, endpoint factors, and relation inversion
used above.  Retain exactly one base-13 inverse-tree digit from one source or
current stage, and marginalize source-root difference, pointed absolute tail,
all other inverse digits, deep/horizon labels, and exact relation addresses.
Then its endpoint response map has image contained in the canonical
six-dimensional `C_2` module `W`.  If the digit is nonflat, its image equals
`W`; its branch marginal is `W_0`; and the induced reflection quotient is

```text
W       = 3 trivial + 3 sign,
W/W_0   = 1 trivial + 1 sign.
```

This quantifies only the fixed finite owner construction.  It is not a theorem
about other clocks, rows, primes, physical currents, or LRC(14).

## Cheapest decisive computation

The decisive two-stage hostile is the lawful joint table

```text
T_joint(b_source,r_owner,state,t),
```

formed before either digit is marginalized.  It must recover the source table
after summing `r_owner` and the owner table after summing `b_source`.  A formal
product of the two final marginals is inadmissible.

Use the seven-dimensional annihilator of the pinned RREF basis of `W` as a
streaming test.  Contract every joint response row against those seven
functionals during relation inversion:

```text
first nonzero contraction  => finite counterexample to the ceiling;
all contractions zero      => im(T_joint) is contained in W;
rank six plus both exact marginals => the joint refinement still reaches only W.
```

This is cheaper than first analyzing all Fourier modes and is logically
decisive.  If the joint rank exceeds six, the first surviving annihilator and
its `(b_source,r_owner,state)` row are the minimal hostile witness.  If it
remains six, the result still does not identify the two digits; it says only
that the endpoint response forgets their difference at this quotient.

## Connection contract

| field | exact answer |
|---|---|
| source | two distinct `52`-dimensional radix/state carriers, one before source transfer and one after the depth-`13^5` current fold |
| target | `H=Fun(F_13(relation),F_p)` |
| map | actual endpoint integration followed by the canonical relation inversion |
| preserved | owner base, Boolean state marginal, literal guard order, endpoint factors, relation residue, abstract path/radix reflection quotient |
| destroyed | 46 kernel directions in each carrier, low/older digits, root difference, pointed tail, deep/horizon labels, chronology, exact address |
| positive gate | equal RREF row spaces, equal involution graphs, `3+3` split, common `4`-dimensional parent |
| hostile | affine target stabilizer is trivial; root-difference and pointed tables have relation rank `13/12` |
| boundary | no source-to-current digit map, joint ancestry table, physical current, row exclusion, or LRC(14) conclusion |

## Reproduction

```text
python -B 04-computation/lrc_r5_endpoint_response_rank_six_representation_audit_20260816.py
python -B -O 04-computation/lrc_r5_endpoint_response_rank_six_representation_audit_20260816.py
```

The semantic digest is
`7baeb128a4c4d5998342611fdcf821d002ffb4622692dd72c03e6f11c8d9825a`.

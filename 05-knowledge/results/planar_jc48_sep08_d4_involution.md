# A complete fixed-sheet Euler test for the finite D4 reflection model

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED for the finite theorem; geometric consumer CONDITIONAL.**
The finite theorem concerns the literal even signed permutation group of
order 192. Its geometric use requires an actual monodromy map into this
group, the six declared access pairs, and equality of every retained
sheet set with the complete corresponding fixed set. These hypotheses
are not inferred from the three-cusp numerical scout.

## 1. Precise finite theorem

Let `W` be the group of even signed permutations of four coordinates,
acting on the eight letters `+1,...,+4,-1,...,-4`. Write a reflection
`s_ij` for simultaneous swaps `i<->j`, `-i<->-j`; write `s_i,-j`
for swaps `i<->-j`, `-i<->j`. Fix

    a=s12, b=s23, c=s34, d=s3,-4, e=bcb^-1.

These four reflections generate all 192 elements. Leaves `a,c,d`
commute pairwise and each braids with central generator `b`.
The standard Weyl-group interpretation of this explicit signed model
is **CITED** from Noam Elkies's author-hosted
[Lie Groups and Lie Algebras notes, Lecture 21](https://people.math.harvard.edu/~elkies/M222.23/index.html).
The proof below uses the literal finite group, not an unverified claim
that a particular geometric monodromy factors through it.

For any transitive finite `W`-set `Omega` put `d0=|Omega|` and
`a0=|Fix(b)|>0`. For a pair `p,q` put `f(p,q)=|Fix(p) intersect Fix(q)|`.
Retain the six original access pairs, in this order:

    (b,c), (b,d), (a,e),
    (e^-1ae,b), (a,d), (c,d).                           (1)

Define

    E(Omega)=3d0-8a0+sum_(six pairs) f(p,q).             (2)

**Finite theorem.** If `E(Omega)=1`, then `d0=1` or `d0=4`.
For `d0=4`, necessarily `a0=2`; the first three pair counts are all
one, and the last three counts are a permutation of `(2,0,0)`.
There are six labelled stabilizers in the fixed normalization below
giving those degree-four cases. No conjugacy-class count is asserted.

## 2. Why only 26 subgroups need to be checked

The element `z=-I` is a central involution in `W`. In a transitive
action a central element either fixes every point or fixes no point:
if it fixes one point, commuting with the group transports that fixedness
to the whole orbit. If `z` acts freely, all fixed sets and simultaneous
fixed sets in (2) are invariant under this free involution. Their sizes
and `d0` are even, so `E(Omega)` is even. Thus `E=1` forces `z`
to act trivially.

Because `a0>0`, choose a point fixed by `b`. Its stabilizer `H`
then contains both `z` and `b`. The action is the literal left-coset
action `W/H`, with the original generators and access pairs retained.
No change of meridian labels or forgotten conjugating path is involved
in choosing this base point.

It remains to enumerate all overgroups of `K=<z,b>`. Start from `K`
and repeatedly adjoin every element outside each discovered subgroup,
closing under multiplication. This procedure is exhaustive: any larger
subgroup can be reached by adjoining its elements one at a time. One
representative of each coset of the form `gH` suffices during construction,
since `<H,gh>=<H,g>`. The final verification then checks **every**
element extension of **every** discovered subgroup again, independently
of that construction shortcut. Exactly 26 overgroups result.

For each such subgroup the source explicitly constructs all left cosets,
the actions of all generators needed in (1), and the six full fixed-set
intersections. Every subgroup is rechecked from its displayed generating
set, and every coset has the required size. Exactly seven labelled
subgroups have `E=1`: the whole group, and six index-four subgroups.
The listed degree and fibre counts follow by exact enumeration.

This is a finite proof on a fully specified group, not a search up to
an arbitrary degree cutoff. The central-involution argument and the
positive fixed-point hypothesis explain both inherited filters.

## 3. The actual geometric consumer and its boundaries

For a whole irreducible support normalized by `A1`, with exactly three
ordinary cusps and three ordinary nodes, the inherited specialization
ledger is

    1=-2a0+sum_cusps n_i+sum_nodes omega_j.             (3)

The actual retained/deleted-sheet passage is supplied by
[the two-cusp passport, Sections 2--3](planar_jc48_sep08_two_cusp_passport.md),
with the number of singular points retained in this ledger. Assume an
actual monodromy map from `W`, with all six pairs in (1) as the actual
cusp and node pairs, and assume that every actual retained set is the
complete fixed set of its meridian. Then the cusp counts are the first
three fixed-set intersections. At a node the two deleted sets are the
complements of those fixed sets, so

    omega_j=d0-2a0+f(p_j,q_j).

Substitution in (3) is exactly (2). Hence a nontrivial such cover has
degree four, which is excluded for polynomial Keller maps by the
already cited degree-four geometric theorem in the inherited passport.
This is a **CONDITIONAL geometric consumer** of a proved finite result.
It is not a new proof of the classical degree-four theorem.

The full-fixed-set hypothesis can be supplied by saturation of a genuine
local cusp injection, as in the inherited passport. Saturation does not
force meridians to be involutions: longer even cycles are an explicit
hostile. Conversely involutive monodromy alone does not say which
inertia-fixed letters are actually retained. Neither shortcut is used.
The finite calculation is therefore not a general three-cusp exclusion.

The six scout words' full abstract presentation was subsequently
identified with Artin D4. That group identification preserves the words
but does not by itself pay their geometric realization, finite quotient,
or retained subsets. The present finite theorem deliberately consumes
the original six pairs, including their conjugates, so it need not
identify arbitrary retained subsets after a presentation simplification.

## 4. Controls, provenance, and reproduction

The closest successful mechanism is the two-cusp finite passport's
retained-set equality in the saturated case. The named hostile is the
nonabelian `S4` quotient from
[the six-word group note](planar_jc48_sep08_three_cusp_group.md): its
natural degree-four action has `E=1`, so parity alone cannot exclude
every action. The classical degree-four consumer is essential.

Another exact hostile is the natural signed eight-letter action. Its
generic fixed-set size is four, its cusp counts are `(2,2,2)`, and
its node fixed-set counts are `(0,0,4)`. Its Euler expression is two,
not one. Here the central involution acts freely, exactly as predicted
by the parity reduction. This transitive representation is not a Keller
cover with the required Euler characteristic.

The live concepts are the full fixed-sheet sets, a central involution,
point stabilizers, the six access pairs, and the actual Euler ledger.
The map from a transitive action to `W/H` retains all these fixed-set
counts. Passing from actual sheets to fixed sets would destroy deletion
data without the explicit equality sidecar. The cheapest decisive test
was the eight-letter hostile followed by the complete 26-overgroup
enumeration; it found the surviving degree-four consumer rather than
mistaking a finite representation for a counterexample.

Reproduce from the worktree root:

```sh
python3 -B 04-computation/planar_jc48_sep08_d4_involution.py
python3 -B -O 04-computation/planar_jc48_sep08_d4_involution.py
```

The [source](../../04-computation/planar_jc48_sep08_d4_involution.py)
uses only exact integer permutations and always-active checks. Both
modes pass **5,631 gates** and reproduce
[the frozen output](planar_jc48_sep08_d4_involution.out). The [independent audit](planar_jc48_sep08_d4_involution_audit.md)
accepts the reduction, exhaustive universe, original pairs and both
replays. Its separate signed-kernel/lift classifier reproduces all26
subgroups and their fixed-count characters without the primary BFS
or literal coset-action computation.

Frozen source SHA-256: `752ef3690dd6e608a98c9bdcb9632f71911588f3e1e4a9523f4fb71278975c2b`
(5158 bytes). Frozen output and both producer replays SHA-256:
`7e7b10f6d3cd4e0ed116fa74c44d001e09815077a1f92a74184a510562ad2384`
(439 bytes).

Independent audit:12530 bytes, SHA-256
`ce45b6d2eda2d01e96047b0025ddcf9e0795c7c3590ad15348e5546b2d77fa87`.

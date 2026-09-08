# Every H4 meridian of type (3)(2) has the same image

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** September 8, 2026.
The minimal power-transfer guess is REFUTED. The repaired theorem is
uniform in ambient degree, including arbitrarily many fixed letters.

## 1. Precise theorem, scope and inheritance

Let a,b,c,d be permutations of a finite set, each with exactly one
three-cycle and one disjoint transposition, and all other letters fixed.
Suppose

    aba=bab,       bcb=cbc,       cdcdc=dcdcd,
    [a,c]=[a,d]=[b,d]=1.                              (1)

Then

    a=b=c=d.                                         (2)

Thus their generated action is never transitive: the common permutation
has distinct nontrivial orbits of lengths three and two. No full-retention
assumption, Euler identity, involution quotient or mapping-degree cutoff
is used in this abstract conclusion.

The [actual mixed-cusp braid supplier](planar_jc48_sep08_mixed_cusp_braid.md)
pays a four-positive-meridian H4 quotient with chain `(z,y,w,x)` and labels
3,3,5. Relabel it `(a,b,c,d)` in (1). Its connected actual complement cover
must have transitive monodromy. Therefore a hypothetical whole-support
Keller map for that curve cannot have positive meridians of type (3)(2).
The same exclusion transports through the proved connected intrinsic good
boundary family of
[the (5,3,3) cusp geometry](planar_jc48_sep08_three_cusp_boundary.md).
This does not assert that the complement is equal to the Artin group or
that any abstract representation realizes a Keller source.

The closest proved mechanism is
[the H4 single-cycle theorem](planar_jc48_sep08_h4_single_cycle.md):
its all-ambient conclusion is used only after the transposition factors
have been proved identical. The other recently paid type is the
[H4 involution exclusion](planar_jc48_sep08_h4_involution.md), but it is
not a dependency of (2). The
[odd-retention gate](planar_jc48_sep08_odd_retention.md) supplies a separate
proved mapping-degree floor; no step here uses that floor.

Root proposed cubing the mixed generators to recover a Coxeter action.
The first test refuted that map at the fifth edge, on the smallest possible
support. The successful repair retains the cross-intersections between
the two cycle lengths instead. The concept board is canonical factors,
odd local relations, cross-length support cells, centralizer orbit checks,
and the surrounding three commuting relations. The corrected near miss
is assuming that taking a generator power preserves a braid relation.
The least-used sidecar is the direction of a cross-length intersection:
one vanishing direction at the fifth edge forces both to vanish.

Targeted recovery found the general odd-run mechanism and the single-cycle
H4 supplier, but no existing (3)(2) equality theorem. This records a local
statement comparison, not an external priority claim. No additional
literature input or exceptional-group classification is required.

## 2. The smallest hostile to cubing a braid relation

On five letters take

    sigma=(1 2 3)(4 5),       tau=(1 2 4)(3 5).         (3)

The alternating words of length five agree. However

    sigma^3=(4 5),           tau^3=(3 5).

The product of these two transpositions has order three, so their
alternating length-five words differ. Thus the proposed power map does
not preserve the fifth relation, even with the exact same cycle type and
no extra fixed letters. Five is the minimal ambient degree for this type.

The first failed implication is a local one, before any transitivity or
actual retained-page predicate: a generator power is not a homomorphism
on noncommutative words. The ordinary length-three relation does survive
cubing for this particular cycle type, as proved below. This strongest
local survivor must not be upgraded to the fifth edge without the other
H4 relations. The pair (3) is not an H4 tuple or a Keller realization.

## 3. A complete pair theorem with a finite-to-unbounded bridge

For a permutation s of the declared type define its canonical factors

    P_s=s^4,    Q_s=s^3,    s=P_s Q_s=Q_s P_s.         (4)

Here P_s is its single three-cycle and Q_s its unique transposition.
They have disjoint supports. For a pair sigma,tau record the four cells

    (i,j,k,l) = (
       |supp(P_sigma) intersect supp(P_tau)|,
       |supp(P_sigma) intersect supp(Q_tau)|,
       |supp(Q_sigma) intersect supp(P_tau)|,
       |supp(Q_sigma) intersect supp(Q_tau)| ).        (5)

The complete accepted cell tables are as follows. Counts refer to fixed
`sigma=(1 2 3)(4 5)` on ten labelled letters; all other sigma letters
are fixed.

| Relation | Cell (i,j,k,l) | Number of partners |
|---|---|---:|
| braid3 |(3,0,0,2)|1|
| braid3 |(3,0,0,1)|10|
| braid3 |(2,0,0,2)|15|
| braid3 |(2,0,0,1)|120|
| braid5 |(3,0,0,2)|1|
| braid5 |(2,1,1,1)|6|
| braid5 |(2,1,1,0)|30|
| braid5 |(1,0,0,2)|60|

In particular there are 146 ordinary partners, all with ordinary-braided
transposition factors. There are 97 fifth-edge partners, only 61 of which
retain braid5 on their transposition factors. The hostile (3) belongs to
the six-element `(2,1,1,1)` class; the other thirty failures are retained,
not discarded from the table.

### Why the ten-letter computation is a uniform proof supplier

Each pair has at most ten moved letters in its union. Restrict to that
union, discard any jointly fixed labels, and add fresh jointly fixed
labels if needed to reach ten letters. All braid identities, canonical
cycle factors and cells (5) are unchanged. A simultaneous relabelling
puts sigma in the displayed form. There are exactly

    2 * binom(10,3) * binom(7,2) = 5040                (6)

possible tau of type (3)(2). The source enumerates every one: choose its
three-element support, one of its two orientations, and a two-element
support in the complement. It applies no support, parity, transitivity,
retention or preliminary relation filter to this pair universe. The two
literal alternating words are then compared exactly.

This is a finite-reduced theorem about all pairs, not a conjecture based
on a scan through ambient degree ten. An arbitrary four-generator tuple
may move more than ten labels; only each individual pair is embedded in
the complete universe (6). The proof in §4 applies the resulting pair
lemmas separately and retains their actual supports in the original
ambient set. It never embeds a general full tuple into ten letters.

The entire raw universe has a second exact check. The centralizer of the
fixed sigma is `C3 x C2 x S5`: the different nontrivial cycle lengths are
preserved and the five fixed labels may be freely permuted. Its standard
cycle generators and four adjacent transpositions on the fixed labels
partition all 5040 partners into **64 complete orbits**, without a braid
filter. Independent orbit closure, predicate invariance and weighted
cell counts recover the same table. Each literal word is also evaluated
by a separate rightmost-first action on individual letters. The source
records semantic hashes of both complete raw and complete orbit data.

For direct inspection the positive counts have simple size checks. In
the ordinary table they are respectively `1`, `2*5`, `3*5`, and
`3*5*2*4`; the fifth table has `1`, `3*2`, `3*2*5`, and
`3*binom(5,2)*2`. The orientation predicates are decided by the exact
word computation, not by these counts alone.

### Two consequences needed by the H4 reduction

First, every ordinary pair has both cross-length cells zero. Its two
three-cycle supports and its two transposition supports therefore form
disjoint invariant unions. Restricting its ordinary braid relation to
these two unions gives

    P_sigma P_tau P_sigma = P_tau P_sigma P_tau,
    Q_sigma Q_tau Q_sigma = Q_tau Q_sigma Q_tau.       (7)

Second, for a fifth-edge pair, if even one cross cell, say
`supp(P_sigma) intersect supp(Q_tau)`, is empty, the two mixed rows in
the table are impossible. The remaining rows have the other cross cell
empty and l=2. Since a transposition is uniquely determined by its support,

    one empty cross cell at braid5 => Q_sigma=Q_tau. (8)

Statement (8) does not assert the power map on arbitrary fifth-edge
pairs. Its hypothesis is exactly the information supplied globally below.

## 4. The uniform H4 reduction

We first record an elementary single-three-cycle fact. If p and r commute,
and q ordinarily braids with both, then p=r. Indeed commuting single
three-cycles have either disjoint or equal supports. If disjoint, the
ordinary half-support bound makes q meet each in at least two letters,
which its three-element support cannot do. If their supports agree,
either r=p or r=p^-1. In the latter case q cannot have that same support,
since commuting plus the odd relation forces equality. If q has a
four-letter union with p, conjugate by a power of p and relabel the
outside letter to put

    p=(1 2 3),  q=(1 4 2).

This is the only orientation with the shared pair {1,2} satisfying braid3.
The inverse p^-1 fails braid3 with that q: at letter 2 the two alternating
words have images 4 and 3. This exhausts the remaining case and excludes
r=p^-1. The tiny two-orientation calculation is checked independently
in the exact source.

Apply (7) to the ordinary pairs (a,b) and (b,c). The three-cycle factors
P_a and P_c commute because a and c do. The elementary fact therefore gives

    P_a=P_c.                                         (9)

The relation [a,d]=1 implies [P_a,Q_d]=1. A three-cycle commuting with a
transposition has disjoint moved support: that two-element support is
invariant under the three-cycle and cannot contain a nontrivial orbit of
length three. Hence, by (9),

    supp(P_c) intersect supp(Q_d) is empty.

Apply the genuinely conditional fifth-pair consequence (8) to (c,d):

    Q_c=Q_d.                                         (10)

Now [b,d]=1 gives [Q_b,Q_d]=1, and the ordinary pair (b,c) gives braid3
between Q_b and Q_c. By (10) these two transpositions both commute and
ordinarily braid. In any group the commuting odd braid relation forces
equality, so Q_b=Q_c. Likewise [a,c]=1 and the ordinary pair (a,b) give
Q_a=Q_b. Thus all four generators have a common canonical transposition Q.

This already proves intransitivity: every generator preserves the same
two-letter support of Q, and each also moves its disjoint three-letter
support. To obtain the stronger equality (2), note that Q commutes with
every P_i and each generator is P_i Q. In an alternating word of odd
length L, the common Q contributes the same final Q on both sides. It
cancels, leaving the same length-L relation on the P_i. The commuting
relations similarly reduce to those among P_i. Therefore the four
single-three-cycle factors satisfy the full H4 presentation (1).

The already PROVED
[H4 single-cycle theorem](planar_jc48_sep08_h4_single_cycle.md) now gives
`P_a=P_b=P_c=P_d`. Combining with the common Q proves (2). The dependency
is invoked after its hypotheses are paid, not through an assumed power
map. No actual retained sheet was substituted for a fixed sheet, and no
Euler or Coxeter-group assertion entered the argument.

## 5. Consequence, equality control and what remains open

The conclusion is sharp as an abstract equality theorem: taking all four
generators equal to any `(123)(45)`, with arbitrary fixed padding,
satisfies every H4 relation. Their cyclic image has distinct orbits of
length three and two and is never transitive. The exact source retains
these controls at ambient sizes 5,6,11 and25.

An additional bounded check on ten letters fixes b=sigma, keeps every
one of its 146 ordinary partners for a and c, and every one of its 62
commuting type-(3)(2) partners for d. It then applies the two remaining
commutators and the fifth relation. Exactly one tuple remains: all four
generators equal. This is a supplemental hostile/control of the complete
proof, not its finite-to-unbounded justification.

For the actual marked H4 boundary family, (3)(2) monodromy is therefore
excluded at every mapping degree, beyond the earlier purely numerical
degree floor. The surviving cycle types are not inferred from a proper
power of a braid relation: for example (4)(2) and (3)(3) are not settled
by this argument. No general mixed-cycle closure or Jacobian conjecture
conclusion is claimed.

The source-to-target map first extracts canonical factors by powers and
retains all four cross-support cells of every local pair. It preserves
cycle type and the exact original relation. It destroys the relation on
separate factors unless cross-support separation is first proved. The
minimal hostile supplies that failure, and the surrounding H4 graph
repairs it through one empty cross cell. The next test should therefore
retain cycle-length incidence at the remaining mixed types; simply taking
powers is not a valid supplier.

## 6. Reproduction and independent audit status

Source:
[planar_jc48_sep08_h4_mixed32.py](../../04-computation/planar_jc48_sep08_h4_mixed32.py).
Output:
[planar_jc48_sep08_h4_mixed32.out](planar_jc48_sep08_h4_mixed32.out).
Run from the actual worktree:

```bash
python3 04-computation/planar_jc48_sep08_h4_mixed32.py
python3 -O 04-computation/planar_jc48_sep08_h4_mixed32.py
```

The complete universe and each filter are specified in §§3 and5. Every
check raises explicitly and remains active under optimized Python. The
source performs no numerical approximation, ambient extrapolation or
unrecorded orbit thinning. Normal/optimized replay and exact pins are
recorded below. The [independent complete audit](planar_jc48_sep08_h4_mixed32_audit.md)
accepts the ambient-unbounded support reduction, every pair table, common
transposition and single-cycle stripping. A separate literal-action census
reconstructs all5040 partners and both tables. All earlier frozen artifacts
are unchanged.

Frozen normal and optimized runs are byte-identical: **51,345 always-active
gates**, 3,607 output bytes.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source |8661|`9f83b3345e0d39da9f7abf1bede93b45deecd3234501d5331ebf28a4b640f0f6`|
| Output |3607|`399694228f01eb1409c1476d2a42e25eda8eae32811b3367d9a65fb7f5da44ef`|

Complete raw-pair data SHA256:
`9d33769e2dd066673ed4c828be3d280694ec702f1d364581f00111ded1559031`.
Complete centralizer-orbit data SHA256:
`554bcdaaf35b135645f4c5e4b215e3fca68e56f76b7c1b70c99305972589fadd`.
These hashes pin reproducible finite universes; the source retains and
checks every element before either semantic hash is formed.

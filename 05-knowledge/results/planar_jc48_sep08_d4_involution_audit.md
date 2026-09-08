# Independent audit of the finite D4 fixed-sheet Euler classifier

**Status: independent analytic/source audit PASS; normal, optimized,
and frozen output agree.** The finite result is accepted on its
literal group and complete fixed sets. Its whole-Keller consumer
remains conditional on the declared actual quotient, access pairs,
and retention equalities. This audit does not certify the numerical
three-cusp paths or promote their geometric presentation.

Audited:

- [Primary proof](planar_jc48_sep08_d4_involution.md).
- [Frozen source](../../04-computation/planar_jc48_sep08_d4_involution.py).
- [Frozen output](planar_jc48_sep08_d4_involution.out).
- The actual subset and Euler suppliers in
  [the two-cusp passport, Sections 2--4](planar_jc48_sep08_two_cusp_passport.md).
- The six common-access prefix pairs in
  [the three-cusp group proof, Section 2](planar_jc48_sep08_three_cusp_group.md).

## 1. Literal group, pair labels, and the parity reduction

The source universe is exactly `4!*2^3=192` signed coordinate
permutations with an even number of sign changes. Composition is
rightmost first on the eight signed letters; the stored eight-entry
permutations and their multiplication table preserve that convention.
The four declared reflections are distinct involutions and generate
the whole universe. The three leaves commute pairwise and braid
with the central reflection `b`.

I also read the cited primary author source,
[Elkies, Math 222, Lecture 21](https://people.math.harvard.edu/~elkies/M222.23/index.html),
specifically its root-system examples and Weyl-group paragraph.
It identifies the `D_n` reflections with ordinary coordinate swaps
and simultaneous signed swaps, and gives the even signed permutation
group of order `2^(n-1)n!`. This supports the name and interpretation;
the finite proof independently defines and generates the group.

All the declared meridians are conjugate in this literal group.
Coordinate permutations conjugate the ordinary swap reflections;
an even sign change using a spare coordinate conjugates a positive
swap to the negative swap. The additional symbols `e=bcb^-1`
and `e^-1ae` are specified conjugates. Their individual fixed-set
sizes therefore equal `a0=|Fix(b)|` in any action.

The six simultaneous fixed sets consumed are precisely

    (b,c), (b,d), (a,e),
    (e^-1ae,b), (a,d), (c,d).

The last pair comes from the original prefix pair `(e,bdb^-1)`
by simultaneous conjugation with `b^-1`. Such a common conjugation
preserves the intersection cardinality. No separate relabelling of
the two meridians or replacement of an arbitrary retained subset
is involved. The other five pairs agree literally with their
declared prefixes.

The central element `z=-I` is an involution different from the
identity. On a transitive action, if it fixes one point it fixes
all points; otherwise it acts freely. In the latter case every
full fixed set and every simultaneous fixed set is `z`-invariant,
since `z` commutes with each group element. Their cardinalities
and the total degree are even. Hence

    3d0-8a0+sum six intersections

is even. Euler value one forces `z` to be trivial on the whole
action. This argument would not apply to arbitrary retained
subsets of fixed sets; the primary correctly retains that
distinction.

The assumption `a0>0` supplies a point fixed by the original `b`.
Its stabilizer contains both `b` and `z`. Choosing this point
does not change a generator or an access path. Thus all actions
under consideration are among `W/H` with `<z,b> <= H`.
No action-faithfulness assumption is imposed.

## 2. Primary exhaustive enumeration is valid

The producer begins with `<z,b>` and closes successively under
adjoining every group element, storing the resulting subgroup and
an explicit generating set. Closure under positive multiplication
is enough in a finite group: inverse powers are eventually positive
powers. Every subgroup containing the base is obtainable in this
way by adjoining its elements, so the universe is not an arbitrary
degree cutoff.

During construction the code tries only one element from each set
`gH`. This is sound because `<H,gh>=<H,g>` for `h in H`.
In standard terminology these are left cosets, although the
original source comment called them right cosets; their literal
formula and construction are unambiguous. The parent was notified
to use “cosets of the form gH” in the primary prose. This is a
terminology clarification, not a change to the finite engine.

After discovering the 26 subgroups, the producer checks all
`26*192` single-element extensions without the coset shortcut.
It also reconstructs each subgroup from its recorded generators.
Therefore the resulting finite family is demonstrably upward
closed under every possible extension and contains the base;
this supplies a second completeness check.

The action is correctly the left action on cosets `gH`, implemented
by the assignment of `qg` to its coset. Representatives do not
alter fixedness. All cosets have size `|H|`, exhaust `W`, and have
the central involution trivial. The source counts full fixed-set
intersections directly, retaining every one of the six labels.

## 3. Independent semidirect classification and character counts

I independently rebuilt the quotient

    W/<-I> = V semidirect S4,
    V={even vectors in F2^4}/<1111>,       |V|=4.

This uses no primary multiplication table or eight-letter subgroup
search. For an explicit implementation, represent `V` by the bit
masks `0,3,5,6`, reducing a mask `v` to `min(v,v xor 15)`.
Permutations act by permuting its four bits. Multiplication is

    (v,p)(w,q)=(v+p(w),p q),

with addition followed by that reduction. The positive reflections
have zero first coordinate; the negative `d` has mask three and
the same projected swap as `c`. This gives the 96-element quotient
and the exact six images used by the primary.

For an arbitrary subgroup `Hbar` containing the prescribed lift
`b=(0,(23))`, put `P=projection(Hbar)` and `K=Hbar intersect V`.
Then `P` is a subgroup of `S4` containing `(23)`, and `K` is a
`P`-invariant subspace. There are exactly six possible projected
subgroups, of orders

    2, 4, 6, 6, 8, 24.

These were enumerated in the 24-element permutation group alone.
For each, I tested the zero subspace, the three lines, and all of
`V` for invariance. Fix a generating set of `P` starting with
`(23)`. For every other generator, choose each representative
of `V/K` as its lift; include the fixed lift of `b` and every
element of `K`, and close in the semidirect product. Retain a
candidate exactly when its sign kernel is the prescribed `K`.
This is exhaustive: any subgroup with that projection and kernel
contains such a lift of each generator, and changing a lift by
`K` leaves the generated subgroup unchanged.

There are **32** such invariant-kernel/lift attempts, giving exactly
**26** distinct accepted subgroups. The `(order P, order K)` counts
are

| `(order P, order K)` | Count |
|---|---:|
| `(2,1)`, `(2,2)`, `(2,4)` | 1 each |
| `(4,1)`, `(4,2)` | 2 each |
| `(4,4)` | 1 |
| `(6,1)` | 4 |
| `(6,4)` | 2 |
| `(8,1)` | 4 |
| `(8,2)` | 2 |
| `(8,4)` | 1 |
| `(24,1)` | 4 |
| `(24,4)` | 1 |

For this independently constructed list I counted fixed points
using the character formula, rather than a literal coset action:

    |Fix(g1) intersect ... intersect Fix(gk)|
       = |{x in Wbar: x^-1 gi x in Hbar for every i}|/|Hbar|.

The denominator divides each numerator exactly. All 26 resulting
profiles, including their multiplicities, agree with the primary.
The complete independent profile bank is below; `C` is the ordered
triple of cusp counts and `N` the ordered triple of node counts.

| Multiplicity | `d0` | `a0` | `C` | `N` | Euler |
|---:|---:|---:|---|---|---:|
| 1 | 1 | 1 | `(1,1,1)` | `(1,1,1)` | 1 |
| 1 | 3 | 1 | `(0,0,0)` | `(1,1,1)` | 4 |
| 2 | 4 | 2 | `(1,1,1)` | `(0,0,2)` | 1 |
| 2 | 4 | 2 | `(1,1,1)` | `(0,2,0)` | 1 |
| 2 | 4 | 2 | `(1,1,1)` | `(2,0,0)` | 1 |
| 3 | 6 | 2 | `(0,0,0)` | `(2,2,2)` | 8 |
| 2 | 12 | 2 | `(0,0,0)` | `(0,0,2)` | 22 |
| 2 | 12 | 2 | `(0,0,0)` | `(0,2,0)` | 22 |
| 2 | 12 | 2 | `(0,0,0)` | `(2,0,0)` | 22 |
| 1 | 12 | 4 | `(0,0,0)` | `(4,4,4)` | 16 |
| 4 | 16 | 4 | `(1,1,1)` | `(0,0,0)` | 19 |
| 1 | 24 | 4 | `(0,0,0)` | `(0,0,4)` | 44 |
| 1 | 24 | 4 | `(0,0,0)` | `(0,4,0)` | 44 |
| 1 | 24 | 4 | `(0,0,0)` | `(4,0,0)` | 44 |
| 1 | 48 | 4 | `(0,0,0)` | `(0,0,0)` | 112 |

Only after completing this separate construction and character
calculation did I compare its profile multiset to the primary's
stored rows; they match exactly. In particular the seven
Euler-one labelled stabilizers are the whole group and the six
index-four cases. This is not a count of conjugacy classes of
subgroups or isomorphism classes of transitive actions.

The natural eight-letter action independently gives `a0=4`, cusp
counts `(2,2,2)`, node counts `(0,0,4)`, and Euler value two.
The central element is free there. The natural degree-four quotient
is a valid positive control for the finite statement, not a
polynomial Keller cover.

## 4. Independent geometric ledger and exact conditional scope

Let the whole irreducible support have normalization `A1`, exactly
three ordinary cusps and three ordinary nodes, with no other
singularities. A node lowers Euler characteristic by one relative
to the normalization; an ordinary cusp does not. Thus

    chi(S)=-2,
    chi(S smooth)=-8,
    chi(C^2 minus S)=3.

The generic off-curve degree is `d0`, and the actual smooth-curve
fibre has `a0` points. The actual specialization theorem gives
the cusp fibre as the intersection of its two re-accessed retained
sets. At a node its fibre is `2a0-d0+omega`, where `omega` is
the intersection size of the two deleted sets. These are actual
sets supplied by the inherited local normalization/purity
argument; permutation fixedness alone is insufficient.

Under the explicit hypothesis that each retained set is its full
meridian fixed set, the node formula becomes its full fixed-set
intersection count, since `omega=d0-2a0+f`. Therefore direct
Euler integration gives exactly

    1=3d0-8a0+sum six f.

Equivalently, substituting the three node formulas first gives
`1=-2a0+sum cusp n+sum node omega`, as in the primary.
The coefficient `-2a0` is correct for three cusps; it was not
copied from the two-cusp coefficient `-a0`.

A genuine nonautomorphic Keller cover has `d0>1` by the inherited
birational-etale argument. The finite result then forces mapping
degree four. The primary's already audited classical degree-four
Keller exclusion is the final consumer. Nothing here substitutes
the source normalization's degree for the map's geometric degree.

This conclusion is conditional on an actual action of the literal
finite group on the sheets, carrying the six declared pairs and
the retained-set equalities. An abstract Artin `D4` presentation
does not force involutivity. Saturation of the ordinary-cusp
injection forces full fixed sets and even nontrivial cycles but
also does not force transpositions; the inherited four-cycle
hostile remains valid. Conversely a finite involutive action
does not identify retained subsets unless that extra sidecar is
paid. The primary keeps all these conditions explicit, so it
does not overstate a general three-cusp or JC(2) exclusion.

## 5. Frozen replay pins and disposition

Independent executions from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_d4_involution.py
python3 -B -O 04-computation/planar_jc48_sep08_d4_involution.py
```

Both pass **5,631 always-active gates** and reproduce all **439
bytes** of the frozen output. The source is **5,158 bytes**, SHA256
`752ef3690dd6e608a98c9bdcb9632f71911588f3e1e4a9523f4fb71278975c2b`.
The output and both replays have SHA256
`7e7b10f6d3cd4e0ed116fa74c44d001e09815077a1f92a74184a510562ad2384`.

The separate semidirect/kernel-lift character computation was
performed in a temporary standalone audit script without
importing the primary's mathematical implementation. Its only
use of the primary program was the final, after-the-fact profile
comparison. The complete algorithm, universe, count checks and
resulting bank are recorded above for reconstruction.

No mathematical or source correction remains. The minor `gH`
coset terminology was reported to the owner; the corrected primary
now says “coset of the form gH” and has been accepted.
Status promotion and any downstream use remain with the parent;
this audit does not edit the frozen producer or the geometric
path status.

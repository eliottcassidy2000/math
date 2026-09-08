# Independent audit of the D4 partial-retention obstruction

**Status: independent analytic/source audit PASS; all 12,509 normal,
optimized and frozen gates agree.** The finite result and its geometric
consumer in [the primary note](planar_jc48_sep08_d4_retention.md) are
accepted. No full-fixed-set or central-involution hypothesis has been
silently retained. The geometric conclusion excludes involutive
monodromy under the declared actual marking; it does not exclude all
Artin D4 monodromy or settle the whole Keller problem.

## 1. Actual retained-set bounds

Let the generic degree be `D`, the actual number of retained sheets be
`k`, and the common meridian fixed count be `a0`. These are distinct
quantities, with `1<=k<=a0` and `k<D`. The common value of the fixed
count follows from conjugacy of positive meridians of the whole
irreducible curve.

At an actual cusp, write

    A subset Fix(sigma), B=(sigma*tau)A subset Fix(tau),
    |A|=|B|=k, n=|A intersect B|,
    C=Fix(sigma) intersect Fix(tau), |C|=f.

Every element of `A intersect C` is fixed by `sigma*tau`, so belongs
to `B` as well. Since the complement of `C` inside `Fix(sigma)` has
size `a0-f`,

    n>=k-a0+f.

Together with nonnegativity and the proved actual cusp injection
`2n>=3k-D`, this gives exactly the primary's cusp lower bound. No
involutivity, saturation or equality of actual retained sets with full
fixed sets is used in this step.

At an actual node, let the two full fixed sets be `F1,F2`, each of
size `a0` and with intersection size `f`. For arbitrary retained
subsets `A subset F1`, `B subset F2` of size `k`, at most
`2(a0-k)` points of `F1 intersect F2` can be lost from `A intersect B`.
Thus

    |A intersect B|>=max(0,f-2(a0-k)),
    omega=D-2k+|A intersect B|
         >=max(0,D-2k,D-2a0+f).

Here `omega` is the exact intersection cardinality of the two deleted
sets, so its nonnegativity is independent of either algebraic lower
bound. Choosing retained subsets independently at distinct points
relaxes the necessary condition and cannot create a false exclusion.

The actual Euler identity is also correctly typed for arbitrary
`N>=3`. The whole irreducible curve has normalization `A1`, exactly
three ordinary cusps and `N` ordinary nodes. Its Euler characteristics
are `chi(S)=1-N` and `chi(S_smooth)=-2-2N`; the complement has
Euler characteristic `N`. A node's actual fibre count is
`2k-D+omega`, and a cusp's is `n`. Therefore Euler integration gives

    1=D*N+k*(-2-2N)+sum_(three cusps) n
      +sum_(N nodes)(2k-D+omega)
     =-2k+sum_(three cusps) n+sum_(N nodes) omega.

Keeping only the three declared distinct node contributions and their
lower bounds gives `1>=B`. Each extra node adds a nonnegative term.
The proof neither uses a three-node equality in an N-node application
nor discards unrecorded affine singularities: the hypothesis expressly
excludes any other affine singularities.

## 2. The marked-pair sidecar is sufficient

The marked finite pairs are

    (b,c), (b,d), (a,e), (e^-1ae,b), (a,d), (c,d),
    e=bcb^-1.

For the literal curve, the full
[three-cusp braid audit](planar_jc48_sep08_three_cusp_braid_audit.md)
checks actual common-access pairs through six certified stems and
uniform complex-disk two-root clusters. This is stronger information
than the global word relations alone. Two stems differ from the formal
prefix pairs by an inverse Hurwitz move entirely inside the colliding
pair. That move preserves the local subgroup, hence its joint fixed
set and `f`. On actual sheets its transport

    (A,B) -> (B,tau^-1 A),      B subset Fix(tau),

preserves the retained intersection cardinality; it also preserves
`D-2k+|A intersect B|`. Thus the bounds derived on the actual local
pair apply with exactly the formal pair's finite count `f`.

This argument does not reassert a reaccess equation after an arbitrary
abstract change of generators. The cusp bound is proved first at the
actual geometric pair, and only the valid count invariants are then
transported. The global marked D4 quotient, the local access supplier
and this count transport are distinct load-bearing inputs.

## 3. Full finite universe and primary implementation

A transitive W-set with a point fixed by the literal reflection `b`
is `W/H` with `b in H`. Choosing that point changes neither `b` nor
the other marked group elements. The full group is the 192 even signed
permutations on four coordinates. The primary source constructs its
complete multiplication table, the four specified reflections and all
six pair elements, checks generation and the marked relations, and
starts its subgroup search at `<b>`.

For each subgroup it adjoins every outside element up to the harmless
shortcut of one representative from each set `gH`. This shortcut is
valid because `<H,gh>=<H,g>` for `h in H`. After construction, it
separately checks all 192 one-element extensions of every subgroup,
without the shortcut. Starting from `<b>`, these checks pay the
completeness argument: every larger subgroup can be reached by adding
its elements successively. Exactly 53 subgroups occur.

There is no filter requiring `-I in H`, no parity filter, no index
cutoff and no faithfulness assumption. The previous 26-subgroup
central-quotient computation is not used as this universe. The source
constructs the left action on all cosets of the form `gH` and checks
both their partition and the common meridian fixed count. It considers
every integer `1<=k<=min(a0,D-1)`, giving exactly 211 trials. The
trivial degree-one subgroup contributes no trial, as required by the
finite statement's `D>1` hypothesis. Its data are retained in the
53-subgroup classification.

The ceiling expression `(3*k-D+1)//2` is correct for positive or
negative integer numerator. All other operations in the bounds are
integer arithmetic and exact finite-set counts. Only the six labelled
degree-four stabilizers have `B<=1`, and all have

    D=4, a0=k=2, B=1,
    cusp f=(1,1,1), node f a permutation of (2,0,0).

Full retention is a conclusion on these survivors, not an entry
assumption. Every other trial has the strict integer obstruction
`B>=2`.

## 4. Independent semidirect classification of every subgroup

I independently reconstructed the entire universe as

    W=V semidirect S4,
    V={v in F2^4 : sum(v)=0}, |V|=8,
    (v,p)(w,q)=(v+p(w),pq).

The sign vector is indexed by target coordinate. In this convention
`a,b,c` have zero sign vector and the usual adjacent transpositions;
`d` has sign vector with bits three and four set and permutation
`(34)`. This is the full sign group, without quotienting by its
central all-one vector.

For an arbitrary subgroup `H` containing `b`, set `P=projection(H)`
and `K=H intersect V`. There are six projected subgroups of `S4`
containing `(23)`, of orders `2,4,6,6,8,24`. The vector space `V`
has sixteen subspaces, enumerated by closure under xor. Only the
`P`-invariant ones can be `K`.

For every such `(P,K)`, I fixed a generating set of `P` beginning
with `(23)`, fixed its lift to be the literal `b`, and tried every
lift modulo `K` of the other generators. Together with `K`, those
lifts generate all possible `H`: any subgroup with that kernel and
projection contains such lifts, and any element differs from a word
in them by an element of `K`. I rejected precisely the trials whose
generated kernel enlarged `K`. There were 89 lift trials and exactly
53 distinct accepted subgroups. This algorithm does not perform the
primary's 192-element overgroup search.

I computed their fixed counts using the character formula

    |Fix_(W/H)(q1,...,qr)|
      = |{g in W : g^-1 qi g in H for every i}| / |H|,

rather than constructing the primary's coset permutations. Divisibility
by `|H|` was checked. The complete independent profile multiset was:

| Multiplicity | D | a0 | Three cusp counts | Three node counts |
|---|---:|---:|---|---|
| 1 | 1 | 1 | 1,1,1 | 1,1,1 |
| 1 | 3 | 1 | 0,0,0 | 1,1,1 |
| 2 each | 4 | 2 | 1,1,1 | each permutation of 2,0,0 |
| 3 | 6 | 2 | 0,0,0 | 2,2,2 |
| 2 each | 8 | 4 | 2,2,2 | each permutation of 4,0,0 |
| 2 each | 12 | 2 | 0,0,0 | each permutation of 2,0,0 |
| 1 | 12 | 4 | 0,0,0 | 4,4,4 |
| 4 | 16 | 4 | 1,1,1 | 0,0,0 |
| 1 | 24 | 2 | 0,0,0 | 0,0,0 |
| 3 each | 24 | 4 | 0,0,0 | each permutation of 4,0,0 |
| 3 | 24 | 6 | 0,0,0 | 4,4,4 |
| 4 | 32 | 8 | 2,2,2 | 0,0,0 |
| 4 | 48 | 4 | 0,0,0 | 0,0,0 |
| 1 each | 48 | 8 | 0,0,0 | each permutation of 8,0,0 |
| 1 | 96 | 8 | 0,0,0 | 0,0,0 |

After this independent construction was complete, I compared all 53
literal subgroup sets through the explicit signed-coordinate
identification. Every set matched the primary. All fixed-count
profiles and all 211 tuples of retained size, individual local bounds
and total bound matched as well. Thus the comparison is stronger than
agreement only on the survivor count or degree ceiling.

## 5. Coxeter quotient and actual geometric consequence

I read the primary course source
[Crain and Clement, Lecture 2: Coxeter Groups](https://www.math.ucdavis.edu/~anne/WQ2009/MAT280-Lecture2.pdf),
January 7, 2009: Definition 2 on PDF page 1 and the signed-permutation
and type-D examples on page 2. It supplies the classical presentation
identification of the involutive Artin quotient with even signed
permutations. Reversing coordinate indices sends its `s3,s2,s1,s00`
to the exact `a,b,c,d` used here. Therefore involutive marked sheet
images give a homomorphism from literal `W` to the sheet group;
merely checking relations in a chosen matrix group would not establish
that direction.

Under the primary's whole-curve geometric hypotheses, the cover action
is transitive and `1<=k<D`. All positive meridians are conjugate. If
one were the identity, all generating meridians would be identities,
contradicting transitivity in degree greater than one. If their common
cycle type contained only one- and two-cycles, their images would be
involutions and hence would give a transitive W-set with positive
reflection fixed count.

The actual bounds then force degree four. The classical geometric
degree-four Keller exclusion, already cited and source-audited in
[the three-cusp passport](planar_jc48_sep08_three_cusp_passport.md),
is essential to remove this survivor. The finite action itself exists
and has Euler bound one; the computation does not refute it. The
conclusion is exactly that a putative whole-curve Keller realization
under these access hypotheses must have a meridian cycle of length
at least three.

The result does not assume or prove that arbitrary Artin D4 images
are involutions. It does not replace longer even cycles by
transpositions using cusp saturation. It does not prove existence of
any longer-cycle realization, nor exclude all degrees. An individual
component of a larger nonproperness set, other singularity types or a
parameter family without paid marking transport is outside the theorem.

## 6. Controls and frozen reproduction

The natural signed eight-letter action has common fixed count four
and joint profile `(2,2,2;0,0,4)`. For actual retained sizes
`1,2,3,4`, the independent bounds are respectively `16,8,5,2`.
The primary enumerates every retained subset at each cusp under the
actual reaccess `B=(sigma*tau)A`, and every pair of retained subsets
at each node. Its proper-retention control with `k=1,n=0` is real;
it prevents silently identifying actual and fixed sheets. The six
degree-four survivors are the necessary sharp positive controls.

I read the complete source and independently ran

    python3 -B 04-computation/planar_jc48_sep08_d4_retention.py
    python3 -B -O 04-computation/planar_jc48_sep08_d4_retention.py

Both completed with **12,509 always-active gates** and the exact same
575 bytes as the frozen output. Frozen pins:

- [Source](../../04-computation/planar_jc48_sep08_d4_retention.py),
  8,306 bytes, SHA256
  `94928ce8d0a1099f604fd57234932a8394053879c5e135c704dfde58b32575a3`.
- [Output](planar_jc48_sep08_d4_retention.out), 575 bytes, SHA256
  `bb250893c2f7afbc9f07a8294eae18d23b0db0f79ae6ae6bb01edfb973626d67`.
- Semantic digest
  `8d1dd3ceafa209c6450a09de2b991dbe3638652de9d76993b9624dcbfddbe660`.

No correction to the mathematics, source or frozen output was needed.
The candidate is accepted for owner-controlled promotion with the
explicit actual-access and involutive-monodromy scope above.

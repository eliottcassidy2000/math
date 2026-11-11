# A genuine pentad inside the path compiler: its five indecomposable modules

**PROVED IN THIS NOTE, elementary linear algebra over any field.** The
application to the inherited path compiler is exact. The five modules are not
identified with the five selected tournament templates. Chan–Pak remains OPEN.

## 1. The actual algebra, and where its five objects come from

Session252's [algebra proof](../../session252_fivefold_bridges_20261111/algebra/PROOF.md)
identifies the polynomial path-interface algebra as

    A = { [[d,o,j],[0,e,i],[0,0,d]] }.

It is the endomorphism algebra of the dual-number module
`k[epsilon]/epsilon² ⊕ k`. The terminology *Auslander algebra* for this
endomorphism construction follows the definition in the primary textbook
chapter [Finite representation type](https://www.cambridge.org/core/books/abs/representation-theory-of-artin-algebras/finite-representation-type/5678CAD8D8CCDBDDFD7005D1197F05E3).
The classification below is proved directly; no classification theorem is
imported, and no literature-novelty claim is made.

Let U be the outer two-coordinate space, and V the middle coordinate. The
algebra has two idempotents and two arrows

    a: V -> U,       b: U -> V,       b a = 0.

On its natural three-dimensional module, a sends the middle coordinate to
the output coordinate, and b sends the input coordinate to the middle.
The composite `a b` sends input to output; it is nonzero but square-zero.
The basis is precisely the two idempotents, a, b, and ab.

**Theorem.** Every finite-dimensional A-module is a direct sum of exactly
five possible indecomposable types:

| Name | dim U | dim V | Arrow data |
|---|---:|---:|---|
| S_U | 1 | 0 | both zero |
| S_V | 0 | 1 | both zero |
| M_a | 1 | 1 | a=1, b=0 |
| M_b | 1 | 1 | a=0, b=1 |
| P | 2 | 1 | a(v)=u_0; b(u_1)=v; b(u_0)=0 |

Moreover, the multiplicities are uniquely determined by five elementary
integers. Write

    u=dim U, v=dim V, r=rank a, s=rank b, t=rank(ab).

Then

    mult(P)   = t,
    mult(M_a) = r-t,
    mult(M_b) = s-t,
    mult(S_U) = u-r-s,
    mult(S_V) = v-r-s+t.                              (1)

This holds over every field, including characteristic two.

### Proof

First split the nonzero composites. Let X complement
`im b intersect ker a` in `im b`. Its dimension is t. Choose a basis x_i of X
and lifts y_i in U with b(y_i)=x_i. Then the triples

    U_i=span{a(x_i),y_i},       V_i=span{x_i}

have the displayed P form. They are independent because applying b detects
the y_i coefficients, after which independence of a(X) detects the rest.

These copies are direct summands, not merely submodules. Choose a complement
Y of X in V containing ker a and with a(Y) complementary to a(X) in im a.
Then `Y intersect im b = im b intersect ker a`. Split ker b into a(X), a(Y),
and a further complement, and choose lifts of `im b intersect Y` in U.
Those latter spaces form an invariant complementary representation. On it
both composites ab and ba are zero.

With both composites zero, im a lies in ker b and im b lies in ker a.
Choose bases of the two images and separate preimages. Each nonzero a-arrow
then splits as one M_a, each nonzero b-arrow as one M_b, and unused basis
vectors are S_U or S_V. This proves existence. Removing t copies of P leaves
r-t and s-t nonzero arrows, yielding (1). The remaining simple multiplicities
are the residual dimensions. Their nonnegativity also follows directly from
`im a subset ker b` and
`dim(im b intersect ker a)=s-t`.

The four modules of dimension at most two have only scalar endomorphisms.
For P, endomorphisms are

    U: [[lambda,mu],[0,lambda]],       V: [lambda].

Their only idempotents are zero and identity, over any field. Thus all five
types are indecomposable. Formula (1) gives uniqueness without appealing to
a separate decomposition theorem. QED.

## 2. This pentad is visible in the compiler's natural flag

There is a canonical invariant flag

    0 < F1=<output> < F2=<output,middle> < F3=<output,middle,input>.

Its six nonempty contiguous subquotients have only five distinct types,
because F1 and F3/F2 are both S_U. Explicitly:

| Subquotient | Module | Information retained under every ordinal gluing |
|---|---|---|
| F1, or F3/F2 | S_U | d=H(C) |
| F2/F1 | S_V | e=H(C)+K_C(z) |
| F2 | M_a | d,e,O_C(z) |
| F3/F1 | M_b | d,e,I_C(z) |
| F3 | P | d,e,O_C(z),I_C(z),J_C(z) |

Thus these five objects classify actual partial observers of the physical
path compiler. Every one is a subquotient of its natural three-state
representation. No arbitrary five-set has been imposed.

The two two-dimensional modules have the same dimension and retain opposite
endpoint extensions. Duality exchanges M_a and M_b and fixes the other three
types. This differs from the Platonic-solids duality pattern (one fixed solid
and two pairs). It also does not send the two H9 core templates to these two
modules: **every** coloured core acts on the same natural module P. The core
is an algebra element, whereas the five types classify representations of the
whole algebra. This type distinction is exactly what makes the connection
usable rather than a chosen bijection.

Keeping both two-dimensional observers recovers d,e,O,I but still loses J.
The missing space is the square-zero central ideal k E13. Multiplication is

    J(C▷D)=d_C J_D+O_C I_D+J_C d_D.                     (2)

The product O_C I_D is the path passing from the first block's outgoing side
to the second block's incoming side. Thus J is a compositional interaction,
not a redundant endpoint statistic. Formula (2) is inherited from the exact
path decomposition; the module pentad explains its observer hierarchy.

## 3. Pentagon versus pentagram: an exact central defect

A particularly concrete example lies just outside the five chosen cores.
Let R be the regular cyclic tournament on five vertices, with row masks

    (6,12,24,17,3),         i -> i+1,i+2 modulo5.

The ten two-element apex cuts split into two C5-orbits of five:
`{i,i+1}` are the pentagon edges; `{i,i+2}` are its diagonals, forming the
pentagram. Compare O={0,1} with O'={0,2}.

The **known** unit switch H=41 to43, its boundary-balanced interpretation,
and its repeated-block progressions are already in
[Session192, cut sumsets §5](../../session192_surprising_approaches_20261105/cut_sumsets/PROOF.md).
The arbitrary-context protected-closure version and three-value next layer
are already in
[Session229, dynamic §2–3](../../session229_combined_constructions_20261109/dynamic/PROOF.md).
Those are recovered mechanisms, not new claims here.

The full polynomial calculation adds a sharp observation to the current
five-dimensional-algebra picture. For either cut,

    K(z)=25z+20z²+20z³+25z⁴+15z⁵,
    O(z)=2K(z)/5,       I(z)=3K(z)/5.

The two J polynomials are

    J_edge(z)=7z²+7z³+6z⁴+4z⁵,
    J_diag(z)=5z²+5z³+9z⁴+5z⁵.

Consequently

    T_diag(z)-T_edge(z)
      = z²(z-1)(z²+4z+2) E13.                         (3)

This is an identity of polynomial interfaces, not just equality at one
fugacity. At z=1 even the full positive path matrix misses the switch. At
z=-1 the defect is +2, producing the actual adjacent odd counts. Both
proper two-dimensional module observers miss it at every z.

The pentagon supplies golden geometry, while this exact counting defect has
the different quadratic factor with roots -2±sqrt(2). That is a concrete
reason not to expect a single distinguished quadratic irrational to control
all observations on a pentagonal object. No Pell recurrence or graph
spectral statement is inferred merely from those roots.

For arbitrary common ordinal contexts L,R, equation (3) becomes

    T_(L▷diag▷R)-T_(L▷edge▷R)
      = H(L)H(R) z²(z-1)(z²+4z+2) E13.               (4)

This follows from `T E13=E13 T=d E13`. In particular arbitrary transitive
flanks preserve the +2 Hamiltonian difference. The exact two-flank family is

    F_edge(A,B)=41+15(g(A)+g(B)-2)+120AB
                  +63A+42B+9(A mod2)+6(B mod2),
    F_diag(A,B)=F_edge(A,B)+2.                        (5)

All graphs in (5) are actual apex attachments with the indicated cuts.
The standalone check evaluates the polynomial coefficients by literal
permutation paths and counts small actual children independently.

## 4. What this can and cannot do for the missing-values proof

The positive contribution is a finite, complete classification of linear
subquotient observers and an exact place to put the interaction sidecar.
A proposed composition proof can now say which of the five observers it
uses. A proof using only S_U sees H; a proof using M_a and M_b sees both
endpoint responses; neither can recover the central switch (3). P restores
it with one additional polynomial coordinate.

This gives a precise *sufficient implementation rule*: a target-preserving
rewrite proved in the full P interface transports through arbitrary ordinal
contexts, and equations (2) and (4) calculate its exact scalar defect.
Checking the scalar or the two proper extensions alone is insufficient.
This rule is an inherited algebra law made explicit through the pentad, not
an all-target rewrite or a new coverage theorem.

A cheap positive-construction probe asked for whole-response ternary lifts

    F_D(A,B)=3F_C(A,B)+delta,       delta in {-4,-2,2,4},

using unchanged ordinary binary addresses. All 62,422 core/cut records from
the 532 tournament classes of orders1..7 were compared on their exact
coefficient signatures. There were 518 scaled matches, all with delta=0.
This is **FINITE-EXACT**, only for that literal fixed-address architecture.
It says to change the address/carry state when seeking ternary digits, not
that the true branches or Chan–Pak fail. The modules classify linear
observers over a field; this is not a five-type classification of integral
modules or modules over k[z]. The natural flag and its exact maps remain
valid over the integer polynomial ring, but field decomposition does not
classify integer orbits, target occupancy, or a terminating inverse.

The next useful question is therefore specific: construct a rewrite which
changes the flank addresses and carries the central extension defect, with
a proper integer height decreasing. The full module P is a faithful small
place to express that rewrite. The five indecomposables themselves supply no
such decrease.

## 5. Reproduction and source recovery

Run from the repository root:

    python3 -B 05-knowledge/results/five_core_web_20261111/construction/check_pentad.py
    python3 -B 05-knowledge/results/five_core_web_20261111/construction/search_scaled.py

`check_pentad.py` independently enumerates F2 representations with dimensions
0..2 on each vertex, computes change-of-basis orbits and nontrivial commuting
idempotents, and finds exactly the five displayed indecomposables in that
finite universe. The proof above covers every dimension and field. It also
checks (3) from permutation paths and (5) from literal child DP for A,B=0..4.
`search_scaled.py` derives all classes through order7 by vertex addition,
canonicalizing only within score classes, and compares coefficients. It
imports the inherited compiler only for its subset count helper; class and
cut enumeration are local. No graph-order census is misread as an all-order
spectrum exclusion.

The exploratory `search_central.py`/`central_search.json` recovers the known
central pairs and triples; its discoveries are explicitly credited above.
The orthogonal historical sources consulted include the May31 reflection
[If fifteen were prime](../../../../07-reflections/if-fifteen-were-prime-structure-shapes-truth-s17.md),
whose useful move is to distinguish additive, multiplicative and observer
structure; and Session41's
[additive-support/golden/Gaussian comparison](../../session41_owned_attachments_fibonacci_20261023/additive/RESULTS.md),
whose actual quadratic-ring sidecar prevents an unwarranted transport from
shared four-state cardinality. Neither supplies a tournament existence theorem.

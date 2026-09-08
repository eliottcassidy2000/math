# Independent audit: Artin D4 and the conditional degree-eight floor

**Status: VERIFIED / independent analytic and exact-source audit PASS.**
Root referee, September 8, 2026. The audited primary is
[three_cusp_passport](planar_jc48_sep08_three_cusp_passport.md), produced
by `three_ray_geometry`. Root read its complete proof and source,
independently reconstructed the group and degree arguments, and replayed
both modes against the frozen output. The motivating geometric words
remain **HEURISTIC** in this bundle. Acceptance does not certify paths.

## 1. Group equivalence is reversible

The six prefix/core pairs agree with the independently audited
[earlier group note](planar_jc48_sep08_three_cusp_group.md). Writing
`e=bcb^-1`, the relation `bcb=cbc` gives both
`e b e^-1=c` and `c e c^-1=b`. Conjugating the original node0
commutator by `e` turns it into `[a,c]=1`. With that relation,
conjugating the original cusp-two braid relation by `c` turns it into
`aba=bab`. The other four relations are unchanged, apart from the
simultaneous conjugation of node2 giving `[c,d]=1`.

Every step reverses, so the presentation is exactly Artin D4 on the
same four marked letters: pairwise commuting leaves `a,c,d`, each
braiding with `b`. Neither a Coxeter square relation nor a finite-group
assumption has been added. The independent source also verifies that
the differences between the two key conjugacy expressions and their
targets are exactly the defining braid relator and its inverse in the
free group. The source does not use finite representations to prove
this arbitrary-group equivalence.

## 2. The retained-sheet consumer has the right hypotheses

The geometric statement assumes whole irreducible support, normalization
`A1`, precisely three ordinary cusps and `N>=3` ordinary nodes, a
connected degree-`D` covering with four positive meridian generators,
and actual common-access interpretations of all six core pairs at the
declared three cusps and three distinct nodes. Extra nodes and relations
are allowed. These hypotheses are not supplied by a numerical word or
by an arbitrary isomorphism of abstract presentations.

Let `k` be the actual generic retained count. The inherited local
ordinary-cusp theorem applies with no one-cusp upper bound on special
counts: `2n_i>=3k-D`. If one inequality is an equality, it identifies
the full fixed set with the retained set and forces every nontrivial
cycle to have even length. Global conjugacy of positive meridians then
identifies every retained set with the corresponding full fixed set;
commuting node monodromy preserves each support, so every deleted
support intersection is a union of even cycles and has even size.
This uses the actual equality theorem, not parity of an unlabelled
permutation statistic.

The source's Euler ledger was independently reconstructed. Removing
three cusp preimages and `2N` node preimages gives smooth-stratum Euler
characteristic `-2-2N`. Restoring target singular points gives the
curve characteristic `1-N`, so its plane complement has characteristic
`N`. Integrating actual fibre counts yields exactly

    1=-2k+n1+n2+n3+W,
    omega_p>=max(0,D-2k),  W=sum omega_p.

The two-cusp scalar table is not imported. Positivity `k>=1` and
strict inequality `k<D` have their actual quasi-finite and nonproperness
suppliers. No deleted-sheet number is treated as the number of boundary
components.

## 3. Support and low-degree implications

For transposition images, fix `b=(12)`. A leaf braiding with `b`
is either `b` or shares one endpoint with it. If one leaf is `b`,
commutation of leaves forces all to equal `b`. Otherwise distinct
commuting leaves are disjoint, so at most two distinct leaves can occur,
using at most four labels together with `b`. This proves the support
bound in every ambient degree. The finite transposition bank checks
representative sizes but is not the source of the quantifier.

At degree five, `k=1` contradicts `W>=9` and `W<=3`.
For `k=2`, cusp inequalities give `sum n_i>=3`, hence `W<=2`,
whereas nodes give `W>=3`. Thus `k>=3`; every nonidentity meridian
image with that many fixed labels is a transposition. The uniform
support bound excludes transitivity.

At degree six, node bounds exclude `k=1,2`. Nontrivial images fixing
at least three letters are transpositions or single three-cycles.
Transpositions are excluded. Three-cycles force `k=3`, so actual
retained sets are full fixed sets. Cusp inequalities give `n_i>=2`
and `W<=1`. Commuting single three-cycles have identical or disjoint
supports, giving node overlaps three or zero. Hence every node overlap
is zero. The original node2 pair, conjugated simultaneously by `b^-1`,
makes supports of `c,d` disjoint. Each of the original cusp pairs
`(b,c)`, `(b,d)` requires at least two labels of `supp(b)` in the
corresponding support. Four required labels cannot fit into three.

At degree seven, `k<=2` is immediately incompatible with the node
bound. For `k=3`, each cusp has at least one retained point and
`3<=W<=4`. A cusp count one is a saturated injection; every node
overlap must then be even and at least one, forcing `W>=6`.
If no cusp count is one, `sum n_i>=6` and `W<=1`. Thus `k>=4`.
Again only transpositions or single three-cycles are possible;
the latter force `k=4` and full fixed sets. Then `n_i>=3`, `W=0`,
and the same two disjoint node supports require four labels in the
three-label support of `b`.

All these bounds use `N>=3`, not `N=3` as an equality. Additional
nodes preserve the argument. The known mapping-degree two, three and
four exclusions are credited to their inherited primary references;
they are not coordinate-degree statements. Thus the conditional
consumer gives `D>=8` with its precise access hypotheses.

## 4. Hostiles and the boundary of presentation changes

The named `S4` tuple satisfies the complete six words. With full fixed
sets it has cusp counts `(1,1,1)`, node deleted overlaps `(0,2,0)`,
and Euler value one. It supplies actual local re-access and an abstract
common sheet passport, but not a polynomial Keller realization. Its
proper `S3` cusp-subgroup images prohibit inferring cusp-local global
generation. The cited degree-four theorem is therefore essential.

The natural eight-letter even signed permutation action is transitive,
has four retained fixed letters per meridian, and gives cusp counts
`(2,2,2)` with node overlaps `(0,0,4)`. Its Euler value is two.
This is one failed higher-degree control; it does not exclude all
eight-sheet passports or impose involutivity on general images.

Root specifically checked that the low-degree support arguments consume
the original cusp pairs `(b,c)`, `(b,d)` and the original node2 pair
with a simultaneous conjugation of both sets. They never replace an
arbitrary retained subset by the fixed set merely because the abstract
group presentation has been simplified. The future actual-braid consumer
must still identify the displayed core pairs with the actual local
pairs. A full loop word alone supplies global fixedness, and does not
automatically supply that extra retained-sheet identification.

## 5. Independent replays

Root read all source code and checked the declared finite universe:
free-word identities; all transposition leaf triples for ambient sizes
two through seven with `b` fixed by conjugacy; scalar rows for
`D=5,6,7`, `N=3,4,5`; all named three-cycle support controls; and
the full `S4` and signed eight-letter actions. The scalar enumeration
is explicitly a relaxation and its surviving rows are not asserted to
have geometric realizations. Always-active gates survive optimization.

Both root runs passed **789 gates**, matching all **1712 output bytes**.

| Artifact | SHA-256 |
| --- | --- |
| [Source](../../04-computation/planar_jc48_sep08_three_cusp_passport.py), 7333 bytes | `a77300a0dc7b9f376c97dff000b4b3ae0441429da6dd379430fe4c5778837049` |
| [Frozen output](planar_jc48_sep08_three_cusp_passport.out), normal and optimized | `b8050ea0c25e8ffbfcc3ab86e36a7ae99a9ef125ce4ec36b4cacda730e887ef3` |

```sh
python3 -B 04-computation/planar_jc48_sep08_three_cusp_passport.py
python3 -B -O 04-computation/planar_jc48_sep08_three_cusp_passport.py
```

**PASS** for the exact Artin-D4 identification and the conditional
degree-eight floor. The primary may be promoted with the geometric
word status still marked **HEURISTIC**. No source or output was changed
during this independent audit.

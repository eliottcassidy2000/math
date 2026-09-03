---
id: THM-4381
title: "Exceptional quartic seminormalization and conductor-fibre classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, in the exceptional
  restriction-curve scope stated below. For S=K[B,C,E] inside its
  normalization K[x], let
  c=x^2(x^2-1)^2 h_172 be the THM-4034 conductor and put
  r=x(x^2-1)h_172.  The reduced conductor-fibre algebra has exact dimension
  87.  Over an algebraic closure the exceptional curve consequently has
  exactly one ordinary triple point, with normalization fibre {-1,0,1}, and
  exactly 86 ordinary nodes.  Its seminormalization is S+K r, its
  seminormal defect has length one, and its conductor in K[x] is rK[x].
  Thus dim K[x]/S=89=86+3 while dim K[x]/S^sn=88=86+2.  This classifies one
  restricted exceptional curve; it proves no chart entry, Keller map, or
  consequence for JC(2) or DC(2).
source: root + quartic_niche / JC2 continuation session, 2026-09-03
audit: >
  PASS for the exact characteristic-zero certificate, with an independent
  proof/type audit and alternate residual-minor rank check. The primary verifier
  imports the independently frozen THM-4034 conductor, uses the THM-3703
  monic Apery basis, and computes the reduced conductor-fibre rank 87 over
  the quartic field itself.  It verifies r is not in S, r^2=ch_172, the exact
  retained-fibre gcd, and the three distinct tangent rows.  Rational-fibre
  censuses at two split reductions modulo 137 are hostile controls only.
  Independent normal and optimized replays byte-match the frozen transcript;
  the script has no Python assert nodes or float literals. No full second
  conductor reconstruction is claimed.
depends_on:
  - THM-4034-exceptional-quartic-global-conductor-degree-178
  - THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction
related:
  - THM-3703-russell-cylinder-exceptional-quartic-sagbi-module
  - THM-3737-russell-cylinder-exceptional-quartic-jacobian-image-hyperplane
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-4063-finite-graph-period-connection-and-ramification-no-go
script: 04-computation/jc2_exceptional_quartic_seminormalization_conductor_fibre_thm4381.py
output: 05-knowledge/results/jc2_exceptional_quartic_seminormalization_conductor_fibre_thm4381.out
script_sha256: 76ff209a959e23015a84fe857a10787917a700cea80454c3621459227ab5e0e4
output_sha256: 495f3267c72b5e0e8a244966267513dd981f986eabe615d0f28b2d5215e73b0f
hash_basis: raw LF bytes
---

# THM-4381 -- the exceptional conductor fibre is one triple plus 86 nodes

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, for the single
exceptional restriction curve below.** The degree-178 conductor of THM-4034 does more than bound a finite
defect.  Its radical sees every geometric branch above the singular locus,
and the restriction algebra identifies those 175 branches in exactly 87
fibres.  This forces a complete fibre classification and leaves precisely one
seminormal class missing from the curve.

All rings have characteristic zero.  Geometric singularity statements are
made after scalar extension to an algebraic closure of the quartic field.

## 1. Conductor, radical, and the exact rank-87 certificate

Retain the exceptional field and restriction algebra of THM-3703/4034,

```text
K=Q[alpha]/(F_6(alpha)),
S=K[B,C,E] subset N=K[x],
normalization(S)=N.                                      (1)
```

Put

```text
L=x(x^2-1).
```

THM-4034 proves

```text
c=(S:N)=L^2 h N,          deg h=172,                    (2)
```

where `h=h_172` is monic, squarefree, and coprime to `L`.  It also proves

```text
dim_K(N/S)=89,       dim_K(S/cN)=89.                    (3)
```

Define the monic radical conductor polynomial

```text
r=Lh,                  deg r=175.                       (4)
```

Then `r` is squarefree, `c=Lr`, and

```text
r^2=ch.                                                  (5)
```

Consider the reduced restriction map

```text
pi:S/cN -> N/rN.                                        (6)
```

The companion computes its rank exactly over `K`, not from a rank failure in
a finite field:

```text
boxed: dim_K im(pi)=87.                                 (7)
```

Here is the small exact certificate behind `(7)`.  THM-3703's monic Apery
presentation supplies 89 canonical basis elements of `S/cN`, one in every
leading-semigroup degree below 178.  Exactly 86 have degree below 175.  Their
images modulo the monic degree-175 polynomial `r` remain independent because
they are monic with distinct degrees.  The remaining three have degrees

```text
175,176,177.                                             (8)
```

Reduce each first modulo `r`, then by the exact THM-3703 normal form.  The
resulting three vectors on the 89 gap coordinates have exact `K`-rank one.
Every reduced representative has degree below 175, and the monic filtered
basis identifies `S` in that range with the span of the 86 low elements.
Thus these gap normal forms compute precisely the quotient by the low image,
and the rank is `86+1=87`.  The companion freezes the pivot and serialization
hash of those three exact gap vectors.  A reduction at `(p,alpha)=(137,44)`
also has top rank one, but this is only a hostile control.

The same exact normal form gives

```text
r notin S.                                               (9)
```

This matters: replacing the conductor by its radical inside `S` would erase
the theorem's one surviving class.

The kernel of `(6)` is exactly the nilradical of `S/cN`.  Indeed, `(5)` makes
every element represented by a multiple of `r` square-zero modulo `cN`.
Conversely a nilpotent maps to zero in the reduced algebra `N/rN`.  Therefore
`(3)` and `(7)` give

```text
dim_K Nil(S/cN)=2,
dim_K (S/cN)_red=87.                                    (10)
```

## 2. Eighty-seven reduced fibres force `1+86`

After extending scalars to an algebraic closure `k` of `K`, squarefreeness of
`r` gives

```text
N/rN tensor_K k ~= k^175.                              (11)
```

A finite reduced subalgebra of `k^175` is the algebra of functions constant
on the blocks of a unique partition of those 175 points.  By `(7)`, the image
of `(6)` is `k^87`; hence the 175 conductor-support points form exactly 87
normalization fibres over the reduced singular locus.

The distinguished fibre is exactly

```text
{-1,0,1} -> (B,C,E)=(0,0,-3).                           (12)
```

It is not merely a known three-point subset.  The exact polynomial check is

```text
monic gcd(B,C,E+3)=x(x^2-1)=L.                          (13)
```

There is also a direct compiler proof.  From

```text
q=Q_alpha(x),
B=(D-1)(D+2)^2,       C=xD(D+2),       E=q(D+3),        (14)
```

the equation `C=0` gives `x=0`, `D=0`, or `D=-2`.  The middle case has
`B=-4`; the first gives `x=0`; and `D=-2`, together with `E=-3`, gives
`q=-3` and hence `x^2=1`.  This proves `(13)` without elimination.

It remains to classify the other 172 roots of `h`.  Every one is a simple
conductor branch by THM-4034.  Such a branch cannot be a singleton fibre.
Indeed, over `k`, a one-branch completed local subring `A subset k[[t]]` with
conductor exponent one contains both `k` and `t k[[t]]`; it is therefore all
of `k[[t]]`, contradicting that the point belongs to the conductor support.

Nor can a normalization point outside the conductor support share one of
these singular fibres.  In the product of completed normalization branches,
a branch on which the conductor is the unit ideal contributes its branch
idempotent to the conductor and hence to the local ring.  If another branch
lay over the same point, that would be a nontrivial idempotent in a local
ring.  Thus every branch in a multiple fibre occurs among the 175 points of
`(11)`, and the preceding non-singleton lemma applies to the entire fibre.

After removing the one block of size three in `(12)`, equations `(7)` and
`(11)` leave

```text
172 points in 86 blocks,             every block size >=2. (15)
```

Consequently every block has size exactly two.  This proves at once both
that there are 86 double fibres and that no larger residual fibre exists:

```text
boxed: conductor fibres = one triple plus 86 pairs.     (16)
```

This counting step is why both the exact rank and the non-singleton lemma are
load-bearing.  A finite census of rational points would prove neither.

## 3. The local singularities are ordinary

Let one of the 86 pair fibres have normalization

```text
k[[t_1]] direct_sum k[[t_2]].                           (17)
```

The conductor exponent is one on both branches, so the completed local ring
contains `t_1k[[t_1]] direct_sum t_2k[[t_2]]`.  Its residue image is diagonal
because the ring is local.  Hence it is exactly

```text
{(f_1,f_2):f_1(0)=f_2(0)},                              (18)
```

the ordinary-node ring.  Thus all 86 pair fibres are ordinary nodes; no
separate tangent genericity assumption is being smuggled into the count.

At the retained fibre, the conductor exponent is two on each of the three
branches.  Exact differentiation gives

```text
B'= (0,0,0),
C'= (3,3,3),
E'= (-9,4,9)                 at x=(-1,0,1).             (19)
```

The two latter rows are independent, and their three branch columns

```text
(3,-9), (3,4), (3,9)                                  (20)
```

are pairwise nonproportional.  All elements of `S` have one common retained
value, and their derivative rows lie in the plane spanned by `(19)`; the two
generators `C,E+3` realize that whole plane.  At every exponent-one node,
`cN` already equals `rN` locally, so that node contributes zero to the kernel
of `(6)`, hence zero to `Nil(S/cN)`.  Thus the full two-dimensional
nilradical in `(10)` is supported at the retained fibre.  The completed local
ring consequently consists exactly of the common value, this two-dimensional
first-jet plane, and all branch terms of order at least two.  Pairwise
nonproportionality of the three branch columns in `(20)` says that this plane
cuts out three distinct tangent lines.  Finally, the ambient relation

```text
C^2E=B(B+4)                                             (20a)
```

has derivative `4` with respect to `B` at `(B,C,E)=(0,0,-3)`.  Formal
implicit elimination therefore expresses `B` uniquely as a series in
`C,E+3`; the completed curve is genuinely planar in those two coordinates.
Its three smooth normalization branches have the three distinct tangents in
`(20)`.  The singularity is therefore an ordinary plane triple point.

This distinction is essential.  An ordinary **plane** triple point is not
seminormal: its graph value-equalizer permits all three first derivatives,
whereas the plane embedding retains only the two-dimensional subspace in
`(19)`.  THM-4067 records this as the one-class defect
`binom(3-1,2)=1`.  Calling the triple “ordinary” does not make graph incidence
period-complete.

## 4. The explicit seminormalization

Nodes in `(18)` are already seminormal.  At the triple, seminormalization
replaces the two-dimensional derivative plane by the full three-branch
first-jet space while retaining equality of values.  It follows locally and
globally that

```text
length_K(S^sn/S)=1.                                    (21)
```

The polynomial `r` in `(4)` vanishes at every one of the 175 reduced
conductor points.  At each node it is already in the local conductor, while
at the retained triple it has a simple zero on each branch.  Hence `r` lies
in the graph value-equalizer, and therefore in `S^sn`.  Equation `(9)` shows
that its class is nonzero.  It spans the one-dimensional quotient `(21)`:

```text
boxed: S^sn=S+K r.                                      (22)
```

Equivalently, `(22)` is the smallest subring of `N` obtained by restoring
the single derivative direction destroyed by the planar triple embedding.
Equation `(5)` is a useful closure check: the square of the new class already
lies in the old conductor.

The conductor of the seminormalization has exponent one at every branch of
every node and at every branch of the seminormal triple.  Its monic generator
is therefore exactly the radical:

```text
boxed: (S^sn:N)=rN=x(x^2-1)h_172 K[x].                 (23)
```

The length ledgers are

```text
dim_K N/S      =89=86*1+3,
dim_K N/S^sn   =88=86*1+2,
dim_K S^sn/S   = 1.                                    (24)
```

Thus the original ordinary triple has delta three, while its seminormal
value-equalizer has normalization delta two.  The `88=86+2` equality is the
seminormalized normalization defect, not a revised value for the original
curve's delta invariant.

## 5. Graph-equalizer consequence and boundary

For the complete normalization-branch graph of this curve, THM-4067's
function equalizer is precisely `S^sn`.  Equations `(21)--(22)` therefore
identify its hidden graph-gluing term exactly:

```text
S_Gamma/S ~= K·( r mod S ),                            (25)
```

The right side is a one-dimensional `K`-vector space.  Under exact edgewise
primitives, the fixed period kernel has this one extra class.  The mixed
cokernel is **not** automatically one-dimensional: THM-4067's exact sequence
still requires intersecting this class with the first-output transgression.
No connection quotient is silently identified with `(25)`.

The exact connection contract is:

| field | contract |
|---|---|
| source | the exceptional restriction algebra `S=K[B,C,E]` and its normalization `K[x]` |
| target | the reduced conductor fibre and the graph value-equalizer |
| map | reduction `S/cN -> N/rN`, followed by the partition of geometric conductor points |
| preserved predicate | equality of target values on normalization fibres; conductor multiplicity; local first-jet span |
| destroyed information | stable-order lifts, convergence, global two-variable maps, behaviour outside this exceptional curve, injectivity, and properness |
| restored sidecar | the single class of `r=Lh_172`, which restores the missing retained derivative direction |
| cheapest decisive tests | exact rank `87`, `r notin S`, gcd `(13)`, tangent rows `(19)`, and the non-singleton exponent-one lemma |

## 6. Exact replay and hostile controls

The primary verifier imports the frozen THM-4034 exact conductor certificate
and performs the following new checks:

1. reconstruct the 89 canonical classes of `S/cN`;
2. prove the 86 low classes independent modulo `r` and compute the exact
   `K`-rank-one top gap residual;
3. verify `r notin S`, `r^2=ch`, and replay THM-4034's good-reduction
   squarefree/coprime certificate for `r`;
4. verify the exact retained gcd and tangent rows; and
5. at the split prime 137, compare the embeddings `alpha=44` and `alpha=92`.

At `alpha=44`, the only rational multiple fibre is the retained triple.  At
`alpha=92`, a rational node `{44,134}` also splits, and all tangents in both
multiple fibres are nonproportional.  These changing rational counts are a
hostile reminder that rational-point enumeration is not the geometric fibre
classification.  The proof of `(16)` is the characteristic-zero rank and
length argument above.

Reproduce with

```bash
python3 -B 04-computation/jc2_exceptional_quartic_seminormalization_conductor_fibre_thm4381.py
python3 -O -B 04-computation/jc2_exceptional_quartic_seminormalization_conductor_fibre_thm4381.py
```

Both commands byte-match the stored LF transcript.  The raw-byte hashes are
pinned in the frontmatter.

This theorem classifies the singular fibres and seminormalization of one
exceptional rational curve inside the Russell-cylinder construction.  It
does not construct a global polynomial pair, enter an arbitrary planar map
into this chart, control the stable lift past THM-4046's obstruction, prove a
Keller counterexample, or imply anything about `JC(2)` or `DC(2)`, which
remain **OPEN**.  **QED.**

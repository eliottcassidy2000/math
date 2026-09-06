# Independent audit: nodal nonproperness and actual-page overlap

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.** September 6, 2026.
This audits [the primary note](planar_jc48_sep06_infinity.md), its standalone
exact source, and its frozen output. The assertion concerns the **whole**
nonproperness support of an actual polynomial Keller map `C² -> C²`.
There is no reduction of a general Keller map to nodal or irreducible support,
and `JC(2)` remains open. The finite controls do not prove the geometric
theorem.

## 1. Exact inherited scope and the source correction

I read the complete Jelonek version-3 primary paper, including Theorem 1.1
and its proof and the local normal form in Theorem 3.1. The smooth and
cuspidal endpoints used here fall within its no-self-intersection theorem.
Its connected Euler conclusion belongs to its **smooth-support** theorem;
it is not a theorem about arbitrary connected singular support.
See [Jelonek, arXiv:2011.03472v3](https://arxiv.org/pdf/2011.03472v3).
The [official current record](https://arxiv.org/abs/2011.03472) is withdrawn
and explicitly distinguishes the incorrect versions 4–5 from versions 1–3.
No stronger conclusion from the withdrawn abstract is used.

I also read the actual downloaded public-manuscript PDFs, rather than their
web summaries. Paper III's Lemmas 5.1 and 6.2 are prior actual-page and exact
node-specialization statements. There is a still more direct antecedent in
[Paper II](https://supermindai.github.io/Jacobian-Conjecture/publication/02-one-place-curves-and-page-density.pdf),
Section 5, Corollary 5.2: its global ledger is
`d-1=sum delta_i+Phi-K`, and its following paragraph explicitly sets `Phi=0`
at ordinary normal crossings. Hence the current nodal ledger is a direct
specialization. With one component it forces `a=1`, and Paper II's
Corollary 4.2 already excludes that sector. This is an inherited implication
of that repository manuscript, not an external peer-review or priority claim.

The current proof independently derives the same nodal exclusion and closes
its last case by `d=2` and purity. This avoids Paper II's separate curve
singularity bound and topological-degree-at-most-five dependency. The exact
antecedent was reported to the author and is now explicitly included in the
primary note's header and Section 1. I checked the repaired text. No
mathematical correction was needed.

The local canon dependencies were checked in
[THM-3578, sheet debt](../../01-canon/theorems/THM-3578-zariski-main-boundary-rank-and-sheet-debt.md)
and [THM-4122, intrinsic asymptotic width](../../01-canon/theorems/THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction.md).
Only the constructible Euler identity and the `A¹` normalizations are used;
the class-group lower bound for a different target supplies no extra degree
obstruction here.

## 2. Finite envelope and local sheets

Write the finite normal envelope as `A² --j--> Xbar --pi--> A²`, and let
`B=Xbar\j(A²)`. The proof that `B` has no isolated component is valid:
remove its divisorial components to obtain a normal open `W`; both source
coordinates extend over the remaining codimension-two complement. They give
`r:W->A²`. The maps `j r` and `id_W` agree on the dense source, so separatedness
forces agreement everywhere. Thus `W=j(A²)`. In particular an isolated
omission cannot occur in a local degree-one normalization sheet whose axes
are both present. This uses the affine source and normal envelope, not merely
an arbitrary open surface.

The generic actual cardinality on a component satisfies `1<=a_i<d`.
For the lower bound, a nonconstant component equation pulls back to a
nonconstant polynomial; each source curve in its zero locus has dense image
in the target component by quasi-finiteness. At its generic DVR the finite
normalization has total degree `d`. The deleted contribution is
`delta_i=sum_D e_D f_D=d-a_i`, a sum of positive integers. Here `f_D` is the
residue degree. This is weighted generic degree, not a special-fibre boundary
point count.

At a smooth point of the whole support, the local complement has cyclic
fundamental group. A connected finite normal extension has form
`(u,v)->(u^e,v)`. An actual affine point is etale, so an actual sheet must
have `e=1`. Its axis can be wholly deleted or retained. Pure boundary prevents
deleting only the special point of a retained axis. This proves actual
cardinality constancy on each component minus the singularities of the whole
support. It does not require the global inertia fixed set to equal the actual
sheet set.

At an ordinary node, choose one common nearby labelled fibre `Omega`.
The two local meridians commute. Membership in the actual subset `A` on the
first punctured branch implies fixation by the first meridian; membership in
`B` implies fixation by the second. A label in `A intersect B` is therefore
a singleton orbit under the entire local complement group. Its finite normal
extension is a degree-one bidisk. Both axes are retained there; every possible
deleted divisor would map onto one of those axes, so none remains. Pure
boundary retains the origin. Conversely, an actual point over the node has a
local inverse and supplies a label on both branches. Thus

```text
#F^(-1)(node) = |A intersect B|.
```

This is an equality, including at a self-node. The transport of raw labels
must be shared on the two branches. Equal cardinalities on a self-node do not
make the two subsets equal. The two-copy deleted-axis hostile correctly
separates actual subsets from inertia fixed sets; the punctured-bidisk hostile
correctly identifies the pure-boundary hypothesis.

## 3. Global degree and the purity endpoint

If `b_i` counts the branches above all nodes on the normalization of `C_i`,
then `chi_c(C_i\Z)=1-b_i`. A self-node contributes two to its own `b_i`.
Constructible Euler integration and the local intersection equality give

```text
d-1 = sum_i delta_i(1-b_i) + sum_p delta(p)
    = sum_i delta_i - sum_p omega_p,
omega_p = |(Omega\A_p) intersect (Omega\B_p)| >= 0.
```

The signs and each branch multiplicity cancel exactly. Disconnected supports
and crossings between distinct components require no alteration of this
formula. It supplies a necessary condition there, not a general exclusion.

For a single component, `delta<=d-1` forces `delta=d-1`, `a=1`, and every
`omega_p=0`. If there is at least one node, the two missing subsets there have
size `d-1`; their intersection has size at least `d-2`. Hence `d=2` for a
nonautomorphic map. Its missing generic length is exactly one. There is
therefore exactly one deleted prime over the component, with `e=f=1`.
This assertion is about the **deleted** prime; it does not infer arbitrary
residue degrees from an unweighted count. Actual primes are already etale.
Thus no codimension-one ramification remains anywhere in the finite envelope.

The application of [Stacks Project, Lemma 58.21.4, tag 0BMB](https://stacks.math.columbia.edu/tag/0BMB)
has all its hypotheses: the source is normal, the target is regular, the map
is finite, corresponding local dimensions agree, and it is unramified at the
height-one specializations. It follows that the envelope is etale. The
connected finite analytic cover of simply connected `C²` has degree one,
contradicting `d=2`.

With zero nodes the argument using `omega_p>=d-2` cannot be invoked. The
separate cited smooth endpoint is necessary and is correctly retained.

## 4. Intrinsic `(2,3)` consequence

For intrinsic normalization degrees `(2,3)`, a translation of the parameter,
nonzero target scalings, and affine target shears produce
`U=t²`, `V=t³+c t`. These are admissible coordinate operations. In particular,
the equation `V²=U(U+c)²` is a relation, not an asserted target coordinate
change.

For `c!=0` there is exactly one affine singular point, `(-c,0)`; its two
parameters satisfy `t²=-c`. Their tangent determinant is `-8ct!=0`, so it is
an ordinary node. For `c=0` the cusp has bijective normalization and falls
under the separately cited endpoint. Consequently a **sole** intrinsic
`(2,3)` component is excluded in every mapping degree. In the notation of
THM-4122 this removes `rho=1` when the reduced pair is `(2,3)` and there is
one component. Cover-inflated parametrization degrees, an arbitrary small
source width, and the occurrence of that reduced pair have not been inferred.

## 5. Source and frozen replay

I read the complete standalone Fraction source and independently ran normal
and optimized Python; both outputs equal the frozen output byte for byte.
There are **4,223 always-active gates**. The exact universe is 1,136 ordered
pairs of nonempty proper actual subsets for `d=2,...,5`, 777 debt assignments
on seven specified incidence graphs for `d=2,...,7`, finite one-component
arithmetic for `d=2,...,17`, and the named cubic/node controls. The symbolic
topology, purity, and all-degree quantifiers come from the analytic proof.

```sh
python3 -B 04-computation/planar_jc48_sep06_infinity.py
python3 -B -O 04-computation/planar_jc48_sep06_infinity.py
```

```text
source SHA256 a1a90e2827316bf1da3fe4c99570c093101dda892b387c27e1ef193bc9e82e21
output SHA256 d725ea2e517b82788b6a5faf9a8ca10079263003027012264fde98dac6a1b646
raw-labelled pair digest 90e1f46065dcec63638d3897f3b3b43b71c651d9af1ebbc75d7cc50aaee2607f
```

The sole surviving finite formal passport `(d,delta,sum omega)=(2,1,0)` is
correctly preserved as a hostile to arithmetic-only exclusion. The proof's
purity step, not the finite enumeration, eliminates it. Primary mathematical
statements and source are accepted; the direct public-manuscript overlap was
the only requested scope clarification and has been incorporated. This audit
is frozen.

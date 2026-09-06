# Nodal nonproperness: an overlap debt and a one-component exclusion

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED; RECOVERED COROLLARY
with an independent shorter proof. Bounded controls are FINITE-EXACT.
JC(2) remains OPEN.**
September 6, 2026. This note
concerns actual polynomial Keller maps over `C`, in every degree. It imposes
a condition on their **whole** nonproperness support; it does not prove that
an arbitrary map has that support.

## 1. Inheritance, source checks, and target

The inherited mechanisms are:

- [THM-3578, Zariski-main boundary rank and sheet debt](../../01-canon/theorems/THM-3578-zariski-main-boundary-rank-and-sheet-debt.md),
  especially its exact constructible Euler identity, not its class-group
  bound (which is vacuous for target `A2`).
- [THM-4122, asymptotic width and resonant shear contraction](../../01-canon/theorems/THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction.md):
  normalization of every nonproperness component is `A1`; its **intrinsic**
  pole pair is distinct from a covering-inflated parametrization degree.
- [THM-4124, integral-degree-ratio all-vertex shear](../../01-canon/theorems/THM-4124-planar-keller-integral-degree-ratio-all-vertex-shear.md)
  and [THM-4126, extremal-target vertical nonproperness](../../01-canon/theorems/THM-4126-jc23-extremal-target-vertical-nonproperness.md)
  retain their stated chart and actual-automorphism restrictions.

The closest prior page lemmas are in the complete primary manuscript
[Paper III, *Disjoint components and Suzuki pencils*](https://supermindai.github.io/Jacobian-Conjecture/publication/03-disjoint-components-and-suzuki-pencils.pdf),
Lemmas 5.1 and 6.2. Its source directory is
[`manuscripts/03-disjoint-components-suzuki-pencils`](https://github.com/SuperMindAI/Jacobian-Conjecture/tree/main/manuscripts/03-disjoint-components-suzuki-pencils).
The published-on-repository manuscript states proved smooth-support page
constancy and exact node specialization, and an all-degree exclusion for
its specified two-nodal-cubic support. These are antecedents, not external
peer-review claims. We read the actual PDF proofs; the local arguments
below are given again so this note need not import a public summary.

[Paper II, *One-place curves and page density*](https://supermindai.github.io/Jacobian-Conjecture/publication/02-one-place-curves-and-page-density.pdf),
Section 4, gives the more general one-component density bound
`(d-3a)R <= a+1`, with arbitrary singularities, and excludes its unique
actual-page sector using the low-degree theorem. Its source is
[`manuscripts/02-one-place-curves-page-density`](https://github.com/SuperMindAI/Jacobian-Conjecture/tree/main/manuscripts/02-one-place-curves-page-density).
More directly, its Section 5, Corollary 5.2 already proves
`d-1=sum_i(d-a_i)+Phi-K`; the immediately following paragraph gives
`Phi=0` at ordinary normal crossings. Therefore (5) below is the exact
wholly-nodal specialization of that **inherited** ledger. With one component
it forces `a=1`, and Paper II's Corollary 4.2 already excludes that sector.
The present note recovers this consequence explicitly and gives a short
independent finish through `d=2` and branch purity, avoiding the general
branch budget and topological-degree-at-most-five theorem. Section 4 connects
the recovered exclusion to intrinsic target poles. It neither extends the
nodal calculation to arbitrary singularities nor claims a new external
theorem or replaces Paper II's general singularity budget. This direct
antecedent was identified by the independent audit before promotion.

For the smooth or cuspidal endpoint we use **CITED** Jelonek, *A note on the
Jacobian Conjecture*, Colloq. Math. 170 (2022), 85–90,
[DOI 10.4064/cm8671-12-2021](https://doi.org/10.4064/cm8671-12-2021),
Theorem 1.1 as read in the correct
[author version 3](https://arxiv.org/html/2011.03472v3).
It excludes a whole nonproperness curve without self-intersections. Its
Theorem 3.1 also supplies the usual local smooth-discriminant normal form.
**Version guard:** the [current arXiv record](https://arxiv.org/abs/2011.03472)
is withdrawn and says versions 4 and 5 are incorrect while versions 1–3
are correct. We do not use the stronger disconnectedness assertion in the
withdrawn landing-page abstract.

The control board is: normalized target poles; actual versus deleted
unramified pages; node branch labels; Euler debt; purity; polynomial target
normalization. The principal hostile is a deleted inertia-fixed page. The
corrected near miss is to add generic debts without the node overlap. The
least-used sidecar is the actual subset of the local proper-locus labels.

## 2. The finite envelope and the actual-page local lemmas

Let

```text
F=(P,Q): A2_C -> A2_C,       Jac(P,Q) in C*,
d=[C(x,y):C(P,Q)],
S = the nonproperness set of F.
```

Suppose `F` is nonautomorphic, so `d>1` and `S` is a nonempty curve.
Let `C_1,...,C_s` be its irreducible components. Zariski Main factors `F` as

```text
A2 --j(open)--> Xbar --pi(finite degree d)--> A2,
B=Xbar\j(A2),                                      (1)
```

where `Xbar` is the normal finite normalization of the target in the source
function field. The finite map is etale away from `S`; every branch divisor
lies over `S`, because `pi` agrees with `F` away from the deleted boundary.

**Pure deleted boundary.** `B` is the union of its prime divisors, with no
additional isolated deleted points. Indeed, remove those divisors from
`Xbar` to get a normal open surface `W`. The complement of `j(A2)` in `W`
has codimension at least two. Normal Hartogs extends the source coordinate
functions `x,y` to `W`, defining `r:W->A2`. The maps `j r` and `id_W` agree
on the dense open `j(A2)`, and hence everywhere since `W` is separated.
It follows that `W=j(A2)`. Finiteness preserves dimension, so boundary
divisors map onto curves, and the finite-envelope properness criterion
identifies their union of images with `S`.

Let `a_i` denote the generic number of **actual affine points** over `C_i`
and put `delta_i=d-a_i`. Then

```text
1 <= a_i <= d-1,       1 <= delta_i <= d-1.             (2)
```

For the first lower bound choose an irreducible target equation `h_i` of
`C_i`. Dominance makes `h_i(P,Q)` nonconstant. Its zero locus in the affine
source contains a curve, and quasi-finiteness makes the image of that curve
dense in `C_i`. Positivity of `delta_i` follows from nonproperness. More
precisely, at the generic point of `C_i` the deleted boundary contribution
is a sum of **lengths** `e_D f_D`, where `e_D` is ramification index and
`f_D` residue degree; actual primes have `e_D=1`. Their sum is `delta_i`.
These are not counts of geometric boundary points at a special target value.

We give the local page proof in the analytic topology. Choose a small
bidisk at a point of `S`, with local axes as stated below, and label a
proper-locus fibre by a set `Omega` of `d` letters. All comparisons use
transport into this one local fibre; no preferred global labeling is needed.
Connected components of a finite cover of the axis complement correspond
to monodromy orbits. Their finite normal extension is unique. In particular
a singleton orbit extends as a degree-one copy of the bidisk.

**Smooth support.** At a smooth point of the whole reduced `S`, its local
equation is one axis and the complement has fundamental group `Z`. Each
connected finite-cover component is cyclic, with normal form
`(u,v)->(u^e,v)`. An actual page along the axis has `e=1`. In such a copy,
the axis is either deleted or retained, and if retained its special point
cannot be deleted alone by purity of `B`. Thus the actual fibre cardinality
is locally constant on the smooth locus of the whole support. In particular
it equals `a_i` everywhere on `C_i` outside the singularities of `S`.

**An ordinary node.** Suppose the whole reduced `S` has local equation
`uv=0`. The two meridians give commuting permutations `sigma,tau` on
`Omega`, since the complement has fundamental group `Z^2`. Let `A_p,B_p`
be the actual-page subsets over the two punctured branches. An actual page
in `A_p` is fixed by `sigma`; a page in `B_p` is fixed by `tau`. If a letter
is in both, it is fixed by the group generated by both meridians. Its
normalization component is therefore a degree-one copy of the bidisk.
Both axes are retained in that copy. Any deleted divisor there would have
to map into one of the axes, so there is none. Pure deleted boundary now
forces the origin to be retained. Conversely, an actual point over the
node is etale and gives a local inverse; its label belongs to both actual
branch subsets. Distinct such points give distinct labels. Therefore

```text
#F^(-1)(p)=|A_p intersect B_p|.                         (3)
```

This proves the equality both for a crossing of two components and for a
self-node. The sizes on a self-node's branches are the same generic `a_i`,
but the two subsets need not be equal. Inertia-fixed letters alone would
not prove (3): a degree-one normalization page can be deleted along one
axis. Actual-page membership is an indispensable extra datum.

Write the missing subsets as `M_p=Omega\A_p`, `N_p=Omega\B_p`. Set

```text
omega_p=|M_p intersect N_p| >= 0.
```

If the two local branches belong to `C_i,C_j` (allowing `i=j`), (3) gives

```text
delta(p)=d-#F^(-1)(p)
        =|M_p union N_p|=delta_i+delta_j-omega_p,
omega_p >= max(0,delta_i+delta_j-d).                    (4)
```

## 3. The nodal overlap identity and the one-component obstruction

Assume every singular point of the **whole reduced** `S` is an ordinary
node. Let `Z` be its finite node set, and let `b_i` count the preimages of
all nodes in the normalization `A1` of `C_i`. A self-node contributes two
to its component; a crossing contributes one to each component. Thus
`chi_c(C_i\Z)=1-b_i`. The Euler debt from THM-3578 is

```text
d-1 = sum_i delta_i(1-b_i) + sum_(p in Z) delta(p).
```

Substitute (4). Every branch debt subtracted in `delta_i b_i` appears once
in the node sums and cancels. We obtain the exact necessary condition

```text
d-1 = sum_i delta_i - sum_(p in Z) omega_p.             (5)
```

This is the wholly-nodal instance of Paper II's Corollary 5.2, independently
rederived above. It includes disconnected supports, self-nodes, intersections between
components, and the case `Z` empty. It does not assert that arbitrary sets
satisfying (5) arise from a connected finite envelope or a polynomial map.

**Irreducible nodal exclusion.** If `S` is irreducible and has at least one
node, (2) and (5) force

```text
delta_1=d-1,       a_1=1,       omega_p=0 for all p.    (6)
```

At any one node, both missing subsets have size `d-1`. Equation (4) then
gives `0=omega_p >= d-2`, hence `d<=2`. Since `d>1`, only `d=2` remains.
Its generic missing **length** over the unique component is `delta_1=1`.
Each boundary prime contributes the positive integer `e_D f_D`, so the
total one forces precisely one such prime with `e_D=f_D=1`. All actual
primes are etale. Consequently `pi` is unramified at every codimension-one
point over `S`, as well as away from `S`.

Purity of the branch locus for a finite dominant generically separable map
from a normal surface to the regular surface `A2` now makes `pi` etale
everywhere; the exact local hypotheses are those of
[Stacks, Lemma 58.21.4, tag 0BMB](https://stacks.math.columbia.edu/tag/0BMB).
Normal source, regular target, quasi-finiteness, equal local dimensions and
unramified height-one specializations all hold here. Its analytification
is a finite covering of simply connected
`C^2`. The envelope is connected, so the covering has degree one, contrary
to `d=2`. This proves the exclusion.

If `S` is irreducible with **zero nodes and no other singularities**, its
normalization `A1` is the smooth curve itself. This is the separately
**CITED** smooth case of Jelonek's version-3 Theorem 1.1; the preceding
node argument is not applied without a node. Combining the two cases:

> A hypothetical nonautomorphic planar Keller map cannot have an
> irreducible whole nonproperness curve with only ordinary nodal
> singularities (including the smooth case).

For a reducible nodal support, (5) remains a necessary condition and does
not supply a general exclusion. For a cusp, tangency, triple point, or
singular local branch the two commuting meridians used in (3) are not the
local model; no corresponding union rule has been inferred here.

## 4. Intrinsic pole pair (2,3) cannot be the sole component

Suppose the whole `S` is one component whose normalization has polynomial
coordinates `U(t),V(t)` of exact degrees `(2,3)` (or the swapped pair).
These are **intrinsic normalization degrees**, not degrees after a finite
cover of the parameter line. Affine reparametrization and invertible
triangular affine target changes put them into

```text
U=t^2,            V=t^3+c t,
V^2=U(U+c)^2.                                           (7)
```

Indeed, complete the square in the quadratic `U`, scale its leading term,
then subtract the quadratic and constant terms of `V` and scale its leading
cubic coefficient. All target changes used here have nonzero constant
Jacobian; no change `V^2-U^3` is used as a target coordinate.

If `c!=0`, the only affine singular point of (7) is `(-c,0)`. Its two
normalization preimages satisfy `t^2=-c`, and their tangent vectors are
`(2t,-2c)` and `(-2t,-2c)`, with nonzero determinant. It is an ordinary
node, excluded by Section 3. If `c=0`, (7) is the cusp; its normalization
is bijective and the curve has no self-intersection. It is excluded by
the separate cited Jelonek theorem. Thus the sole-component `(2,3)`
sector is excluded in every mapping degree.

In the intrinsic notation `(rho d_0,rho e_0)` of THM-4122, this implies:
if a nonautomorphic Keller map has **one** nonproperness component and
reduced coordinate-degree ratio `(d_0,e_0)=(2,3)`, then `rho>=2`.
No assertion is made that a degree-minimal representative has that reduced
ratio, that a width-six chart is attained, or that raw parametrization
multiplicity is invariant. THM-4126's one-nodal-fibre branches are already
closed by later canon; the present statement has no elliptic-pencil or
maximum-weight-eight hypothesis.

There is no new width bound here: THM-4122 already gives width at least six
unconditionally in its reduced `(2,3)` sector, through its equation (16)
and THM-3544. The explicit connection retained here is the conditional
intrinsic `rho>=2` statement. A next geometric frontier must retain
Paper II's nonzero conductor-specialization term `Phi` at singularities
beyond ordinary nodes; repeating a `Phi=0` support does not supply that step.

## 5. Hostiles, source/target map, and reproducibility

The source object is the actual finite envelope with its deleted primes.
The first map retains actual subsets in the local proper-locus fibre and
forgets internal ramification data. It preserves (3) at ordinary nodes.
The second map takes cardinalities and Euler weights, producing (5), and
forgets monodromy, connectedness and realization by polynomials. The missing
ramification data must be restored for the final degree-two purity step.

The named controls are:

1. `Omega={0,1}`, both meridians trivial, `A={0}`, `B={1}`. Then the node
   has no actual page, even though both inertia fixed sets have size two.
   It is realized **locally** by two trivial copies, deleting one opposite
   axis from each. It is not a connected global Keller counterexample.
2. The same formal passport has `d=2`, `delta=1`, `omega=0` and satisfies
   the one-component Euler identity for any positive number of self-nodes.
   Thus Euler and subset arithmetic alone do not exclude degree two.
3. A punctured bidisk would retain both generic axes but omit their origin.
   It defeats (3) without the affine-source/pure-boundary hypothesis. It is
   a local open-surface hostile, not a polynomial map from `A2`.
4. The Danielewski examples in THM-3578 have a different target Euler
   characteristic. Replacing `d chi_c(Y)-chi_c(X)` by `d-1` there is invalid.

The source enumerates all nonempty proper actual subsets on `d=2,...,5`
letters, preserving their raw labels; tests a declared list of finite
node-incidence graphs (including loops); and checks the cubic
normalization and its named node/cusp endpoints by exact polynomial
arithmetic. These are finite controls of the formulas, not substitutes for
the topology, purity or all-degree proof. Every check remains active under
`python3 -O`.

```sh
python3 04-computation/planar_jc48_sep06_infinity.py
python3 -O 04-computation/planar_jc48_sep06_infinity.py
```

Normal and optimized replay are byte-identical: **4,223 always-active
gates**, all 1,136 ordered actual-subset pairs on `d=2,...,5`, and 777 debt
assignments on the seven declared incidence graphs. The sole surviving
formal one-component passport is `(d,delta,sum omega)=(2,1,0)`, which the
analytic purity step excludes. Frozen files:

```text
source SHA256 a1a90e2827316bf1da3fe4c99570c093101dda892b387c27e1ef193bc9e82e21
output SHA256 d725ea2e517b82788b6a5faf9a8ca10079263003027012264fde98dac6a1b646
```

The output is [planar_jc48_sep06_infinity.out](planar_jc48_sep06_infinity.out).
The complete independent
[analytic/source/primary-scope audit](planar_jc48_sep06_infinity_audit.md)
passes the local page argument, Euler cancellation, weighted degree-one
boundary contribution, purity, separate zero-node/cusp citation, intrinsic
normalization and normal/optimized/frozen replay. It identified the direct
Paper II antecedent, which is now credited in Section 1. Source/output were
not changed during audit. No canonical ID is reserved by this note. Root
owns integration.

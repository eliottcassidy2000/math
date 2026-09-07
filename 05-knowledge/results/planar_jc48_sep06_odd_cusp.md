# Odd cusp passports have a finite divisor spectrum

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED. Exact controls
are FINITE-EXACT. Primary local topology and low mapping degree are CITED.
The necessary passports do not assert realization. JC(2) remains OPEN.**
September 6, 2026.

**Current target status:** the concrete five-node curve proposed in §5
is now [excluded by its abelian affine complement](planar_jc48_sep06_global_curve.md).
The [four-node equality successor](planar_jc48_sep06_boundary_plumbing.md)
is also excluded by its marked infinity relations. The general theorem
below is unchanged. The [three-node successor](planar_jc48_sep06_next_infinity.md)
passes that necessary boundary test; its global realization remains open.

## 1. Exact target and inheritance

Let a nonautomorphic polynomial Keller map of the complex affine plane
have **whole** nonproperness set an irreducible curve whose normalization
is `A1`, with exactly one analytic cusp `v^2=u^m`, where `m>=3` is odd,
and `N>=1` ordinary nodes, with no other singularities. Write `d` for its
mapping degree, `a` for the generic actual affine fibre count along the
smooth support, `delta=d-a`, and `n` for its actual cusp fibre count.

Then, before the cited degree-three/four exclusions, a necessary spectrum is

```text
q>=3 an odd divisor of m,       n in {0,1},
d=q+n,       a=(q-1)/2+n,       delta=(q+1)/2.             (1)
```

If `N>=2`, necessarily `n=1`, so

```text
d=q+1,             a=delta=(q+1)/2.                       (2)
```

The cited low-degree results remove `q=3` in both cases. Thus every
geometric survivor has `q>=5`. Formula (1) is a necessary condition,
not a geometric realization theorem. In particular the ordinary cusp
`m=3` is excluded, while this argument does not exclude all higher cusp
cases. The displayed degree spectrum uses only the cited degree-two
through degree-four filter; it is not an assertion that every listed
degree is unresolved by other literature.

The immediate proved supplier is
[the ordinary cusp passport, §§2–3 and §5](planar_jc48_sep06_cusp_passport.md):
its actual-subset re-access lemma, node Euler identity, pure boundary,
and correctly typed low-mapping-degree citations. It inherits
[the nodal infinity note](planar_jc48_sep06_infinity.md) and
[THM-3578, Zariski-main boundary rank and sheet debt](../../01-canon/theorems/THM-3578-zariski-main-boundary-rank-and-sheet-debt.md).
The ordinary cusp result is already independently audited. The extension
here changes the marked meridian word and proves a stronger sparse
permutation lemma that does not need an Artin-relation hypothesis.

The concept board is: **actual retained subsets; meridian access paths;
odd product cycles; the Euler node budget; and mapping-degree exclusions**.
The inherited hostile is the nontransitive ordinary four-sheet passport.
The corrected near miss is replacing actual retained subsets by all
inertia-fixed labels. The underused sidecar is the simultaneous access
path and transported subset: it strengthens a conjugacy-class statement
to the equation `B=gA`. A targeted tracked-canon/results search for odd
cusp/divisor-spectrum statements found no earlier version of this exact
passport theorem; this is a local duplicate audit, not a literature
priority claim.

## 2. Marked meridians and the actual Euler budget

Use left permutation actions and rightmost-first composition. The primary
marked-meridian input is Ichiro Shimada,
[*Lectures on Zariski van-Kampen theorem*](https://www.math.sci.hiroshima-u.ac.jp/shimada/LectureNotes/LNZV.pdf),
§6, pp. 24–25. In its notation, take `p=m=2r+1` and `q=2`. The two
fibre meridians generate the local group and satisfy

```text
m_group=ell_1 ell_0,
ell_(j+2)=m_group ell_j m_group^-1,      ell_(j+m)=ell_j.
```

Consequently `ell_(2r+1)=m_group^r ell_1 m_group^-r=ell_0`.
With `sigma=ell_1`, `tau=ell_0`, this supplies the actual local loop

```text
g=(sigma tau)^r,                g sigma g^-1=tau.         (3)
```

These are marked normal meridians generating the full local group, not
merely unspecified generators of an abstract torus-knot group. The
author PDF's actual passage was read locally; the web-tool timeout was
not treated as evidence. The same presentation is commonly written
`(sigma tau)^r sigma=tau(sigma tau)^r`.

Apply the ordinary supplier's re-access construction unchanged. If `c`
accesses a smooth transversal of this irreducible cusp and `A_eta` is
the actual retained endpoint subset, let `A=T_c^-1 A_eta`. Re-access
that same smooth point by `c'=c g^-1`. Then

```text
A subset Fix(sigma),       B=gA subset Fix(tau),
|A|=|B|=a,                n=|A intersect B|.              (4)
```

The last equality retains geometry: a common label is fixed by both
generators, hence gives a singleton orbit of the entire local group.
Its normal finite extension is a copy of the local target. Its generic
cusp point is retained, and pure deleted boundary prevents removing
only the cusp origin. Conversely an actual cusp point supplies an etale
local inverse and therefore such a common label. No transitivity or
preferred walk along the cusp is assumed.

Normalization `A1` and a unibranch cusp give exactly the same compactly
supported Euler counts as in the ordinary case:
`chi(S_smooth)=-2N` and `chi(A2\S)=N`. For each node `p`, the exact
actual-subset formula has an integer overlap `omega_p` satisfying

```text
#F^-1(p)=d-2delta+omega_p,
omega_p>=max(0,d-2a),         n+sum_p omega_p=1.           (5)
```

Thus `n<=1`. Put `h=|Omega\(A union B)|=d-2a+n`. For any node,
`omega_p>=d-2a`, so `(5)` gives `h<=1`. If `N>=2`, the integer `d-2a`
cannot be positive: otherwise the sum in (5) is at least two. Hence

```text
N>=1 => h<=1;           N>=2 => d<=2a and h<=n.           (6)
```

## 3. The sparse re-access lemma, without an Artin assumption

Let `sigma,tau` be arbitrary permutations of a finite `Omega`. Suppose
`A subset Fix(sigma)`, `B subset Fix(tau)`, `|A|=|B|=a`,
`1<=a<d=|Omega|`, and `B=(sigma tau)^r A` for an integer `r>=0`.
Write

```text
I=A intersect B,   A0=A\B,   B0=B\A,
n=|I|<=1,         k=|A0|=|B0|,         h=|Omega\(A union B)|<=1.
```

Then `h=1`. If `k=0`, the only possibility is the formal identity
passport `d=2,a=n=1`. If `k>=1`, after relabelling,

```text
Omega=I disjoint_union {o,a_1,...,a_k,b_1,...,b_k},
sigma=(o,b_1,...,b_k),       tau=(o,a_1,...,a_k),
w=sigma tau=(o,a_1,...,a_k,b_1,...,b_k),
q=2k+1,                     r=k mod q.                  (7)
```

Both permutations fix `I` pointwise. Conversely these patterns give
all the required equations. In particular `(7)` automatically implies
the odd Artin relation `w^r sigma=tau w^r`.

**Proof.** The common set `I` is fixed pointwise by both generators.
If `h=0`, `sigma` acts only on `B0` and `tau` only on `A0`; both preserve
`A`. Thus `w^r A=A=B`, and `h=0` gives `Omega=A`, contrary to `a<d`.

Now `h=1`; denote the outside label by `o`. The permutation `sigma`
acts only on `B0 union {o}`. Any of its cycles contained entirely in
`B0` is also fixed pointwise by `tau`, hence is invariant as a set under
`w` and `w^-1`. Since that cycle lies in `B=w^r A`, it would lie in
`A`, a contradiction. This argument excludes one-cycles as well as
longer cycles away from `o`. Thus `sigma` is one cycle through every
point of `B0 union {o}`. The identical argument for a `tau`-cycle in
`A0`, using its invariance and `w^r A=B`, gives the other cycle.

Choose the cycle labels in (7). Direct composition gives the displayed
`q`-cycle `w`. Along it, `A0` occupies positions `1,...,k` and `B0`
occupies positions `k+1,...,2k`. Translating the first block by `k`
gives the second. A nonempty proper consecutive block on a cycle has
exactly one entry boundary; any translation preserving the block must
fix that boundary and is zero modulo `q`. Therefore `w^r A0=B0`
holds exactly when `r=k mod q`.

For the converse, at `r=k` the conjugate of the `sigma`-cycle is
`(a_k,o,a_1,...,a_(k-1))`, which is the same cycle as `tau`; adding
multiples of `q` does not change the action. All equations follow.
Finally `k=0` gives `A=B=I`; `1<=a=n<=1` and `h=1` force the two-label
identity passport. This completes the abstract classification.

## 4. Geometric spectrum, equality, and limits

Apply (7) to (3)–(6). For `k>=1`, the congruence in (7) is equivalent to
`2r+1=0 mod (2k+1)`. Hence `q|m`, and counting gives exactly (1).
Since `h=1`, the strengthened inequality `h<=n` for `N>=2` forces `n=1`
and yields (2). All node overlaps then vanish. For `N=1,n=0`, the
unique overlap is `omega=1`, and the actual fibre at that node is zero;
this is allowed by the necessary ledger. Thus no unproved positivity
of every special actual fibre has been inserted.

The formal `k=0` passport has `delta=1`. It is excluded geometrically:
the generic missing length one forces a single deleted boundary divisor
of ramification and residue degrees one. Retained divisors are already
unramified because the original map is etale. The finite normal envelope
is therefore unramified in codimension one over the smooth target;
purity makes it finite etale everywhere, hence degree one over `C^2`.
This contradicts `d=2`. This is precisely the inherited `delta=1`
purity argument, and is not a consequence of permutation algebra.

For `q=3`, the only degrees in (1) are three and four. The primary
mapping-degree exclusions of Orevkov and Domrina are stated, linked,
and verified in [the ordinary supplier, §1](planar_jc48_sep06_cusp_passport.md).
They eliminate these two cases; the degree-two theorem is an alternative
way to exclude the formal identity case. These results concern generic
preimage degree, not the coordinate degrees of a polynomial map.

| Cusp exponent | Divisors before low-degree exclusions | Degrees left by the cited 2–4 filter, `N=1` | Degrees left by the cited 2–4 filter, `N>=2` |
|---|---|---|---|
| `m=3` | `q=3` | none | none |
| `m=5` | `q=5` | `5,6` | `6` |
| `m=9` | `q=3,9` | `9,10` | `10` |
| `m=15` | `q=3,5,15` | `5,6,15,16` | `6,16` |

Proper divisors are essential. In particular the local `(2,9)` group
has both the `q=3` and `q=9` passports before the cited geometric filter.
The abstract `q=5,m=5,n=1` control has

```text
Omega={0,1,2,3,4,5},  sigma=(0,3,4),  tau=(0,1,2),
A={1,2,5},           B={3,4,5},      g=(sigma tau)^2,
d=6, a=delta=3,      I={5},         O={0}.                (8)
```

It attains `h=1` and every `N>=2` node budget with all `omega=0`.
It proves that this local permutation-plus-Euler argument cannot alone
exclude all higher cusps. It does not provide a Keller realization,
a full complement representation, or global affine source geometry.
For `n=1` the local action is a `q`-label orbit plus a singleton and
must not be forced to be transitive.

## 5. Connection contract, hostile controls, and exact reproduction

The source is the actual finite envelope of a Keller map with the stated
**whole** support. The target is the labelled permutation passport
`(Omega,sigma,tau,A,B)` together with integer node overlaps. The map is
local path transport plus integration of actual fibre counts. It
preserves actual retention, access orientation, fixed sheets, the Euler
identity and mapping degree. It loses global complement monodromy,
affine-source realization and coefficient/degree information. The
necessary sidecars are pure boundary and both actual transported subsets.
Source pole degrees, intrinsic target parametrization degrees, cusp
exponent `m`, and mapping degree `d` are distinct quantities.

The finite source contains five cheap failures of stronger implications:

1. Identity permutations on three labels satisfy the Artin relation and
   fix distinct singleton choices `A,B`, but fail `B=gA`. Inertia data
   alone cannot supply the actual-subset bridge.
2. Add an extra unretained fixed sheet to the `q=3,n=1` passport. Local
   re-access and Artin still hold, but `h=2` violates the node Euler
   budget. Dropping that budget destroys the degree classification.
3. Add two common fixed sheets instead. The local equations hold with
   `n=2`, but Euler rejects them.
4. In (8), replacing `sigma tau` by `tau sigma` while leaving the
   transported subset formula unchanged fails. A convention change must
   transport the subset and conjugator together.
5. The two-label identity passport passes every local equation for every
   `m`. Its exclusion needs the geometric `delta=1` purity sidecar.

The strongest survivor after these failures is (1), with the exact
cycles (7), and the `q=5` positive passport (8). The next meaningful
test is an additional global obstruction in mapping degree six for a
whole support having one `(2,5)` cusp and several nodes. Another possible
extension would let `h>1`: the same invariant-cycle proof still requires
every nonretained meridian cycle to meet the outside set, but it no
longer forces a single cycle. Neither extension is proved here.

A concrete next **OPEN** whole-support realization target is the curve
parametrized by `(t^4+t^2,t^6+t^5+t^2)`, supplied by the concurrent root
and geometry lane. Its detailed curve classification is a separate
control, not an input proved by this permutation producer. No Keller
map with this support is being proposed or constructed here.

The standalone [producer](../../04-computation/planar_jc48_sep06_odd_cusp.py)
uses the standard library and explicit failure gates, including under
`python -O`. Its [frozen output](planar_jc48_sep06_odd_cusp.out) records:

- All `400` independent supported permutation pairs for `k=1,2`,
  `r=0,...,9`, with no Artin prefilter; exactly `11` satisfy re-access,
  and all satisfy the predicted full-cycle/divisibility classification.
- `882` canonical rows with `q=3,5,...,19`, `m=3,5,...,99`, `n=0,1`;
  `116` divisibility survivors. Direct tuple composition is checked
  against independent cyclic-position translation.
- Nine named positive passports, the five hostile mechanisms above, and
  the pre/post-low-degree spectra in the table. These are finite controls
  for the analytic proof, not a census-based theorem.

Reproduce from the repository root:

```sh
python3 04-computation/planar_jc48_sep06_odd_cusp.py
python3 -O 04-computation/planar_jc48_sep06_odd_cusp.py
```

Both modes pass **4,072 gates** with byte-identical output. Pins:

```text
source SHA-256:
91f487bac1217da20f235ec0746fd7d587b8d5fa3e4ce8f304d19472df790879
output SHA-256:
18368ecf7ee29d02c11d5daca31d1879582ac5d5bd30c1aecc19e6191aac3cbc
semantic SHA-256:
20873175031fe168cff89a75e43a8fc7f2c1068a64c8a888e3246846702526dd
```

The [independent audit](planar_jc48_sep06_odd_cusp_audit.md) passes the marked-meridian
input, actual re-access and intersection, the stronger Artin-free cycle
lemma, proper divisors and converse, the `q=1` purity boundary, both node
regimes, and the stated low-degree-filter scope. It also rereads the
source and reproduces normal/optimized/frozen output byte for byte.
Audit SHA-256:
`8a16a7d433f56fb007478223341ece4299a649a645a7b17341c7c48bdd3ae7b0`.
All three owned files are frozen; no geometric realization or general
Jacobian-conjecture conclusion is claimed.

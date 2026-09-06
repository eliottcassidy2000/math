# A sparse cusp passport excludes one cusp with ordinary nodes

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED.
The degree-two/three/four exclusions are CITED; the permutation controls
are FINITE-EXACT. JC(2) remains OPEN.** September 6, 2026.

## 1. Target, inheritance, and primary inputs

The conclusion is the following precise support exclusion.

> A nonautomorphic polynomial Keller map `F:A2_C -> A2_C` cannot have
> **whole** nonproperness set an irreducible curve whose normalization is
> `A1`, and whose singularities consist of exactly one ordinary cusp and
> `N>=1` ordinary nodes, with no other singularities.

The map has arbitrary polynomial degree. Throughout, `d` means its mapping
degree `[C(x,y):C(F_1,F_2)]`, or equivalently the number of preimages of a
generic target point. The argument first forces `d` into `{2,3,4}` and then
uses the correctly typed low-mapping-degree theorems. It does not show that
every Keller map enters this support class.

The closest proved mechanism is the exact ordinary-node actual-subset
intersection and pure deleted boundary in
[the nodal infinity note](planar_jc48_sep06_infinity.md), inheriting
[THM-3578, Zariski-main boundary rank and sheet debt](../../01-canon/theorems/THM-3578-zariski-main-boundary-rank-and-sheet-debt.md).
The immediate supplier is §5 of
[the exact curve probe](planar_jc48_sep06_curve_probe.md), which retains the
actual fibre at its single nonnodal point. The general conductor/overlap
ledger in Paper II, *One-place curves and page density*, §5, is an
antecedent as documented in the infinity note. This is an explicit further
specialization with local cusp monodromy, not an external priority claim.
A targeted current-canon/results search found the three cusp parameters
still recorded as OPEN in the frozen curve probe; that note is unchanged.

The control board is: actual retained subsets; re-access of a smooth point;
the cusp braid group; the Euler overlap budget; low mapping degree; and
intrinsic target parametrizations. The named hostile is the four-sheet
`S3` block plus singleton in §5. The corrected near miss is to regard all
inertia-fixed labels as actual, or to identify subsets using an unproved
preferred peripheral path. The least-used sidecar is the access path
together with the actual retained subset, not merely a conjugacy class.

We use these verified primary inputs:

- **CITED local meridian presentation.** Ichiro Shimada,
  [*Lectures on Zariski van-Kampen theorem*](https://www.math.sci.hiroshima-u.ac.jp/shimada/LectureNotes/LNZV.pdf),
  §6, pp. 24–25, gives the local complement presentation for `x^p-y^q` in
  fibre meridians `ell_j`. At `p=3,q=2`, set `m=ell_1 ell_0` and use
  `ell_(j+2)=m ell_j m^-1`, `ell_(j+3)=ell_j`. The relation at `j=1`
  gives `ell_1 ell_0 ell_1=ell_0 ell_1 ell_0`; the two meridians generate
  the local group. This is the meridian-marked cusp presentation
  `B3=<sigma,tau | sigma tau sigma=tau sigma tau>`, not only an abstract
  isomorphism with a torus-knot group. Its Theorem 6.0.4 and Corollary 6.0.5
  give the equivalent torus-knot presentation. The actual author PDF was
  read; a web-tool fetch failure was not treated as source evidence.
- **CITED degrees two and three.** S. Yu. Orevkov,
  [*On three-sheeted polynomial mappings of C2*](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  §1, Theorem 1.1. The definition immediately preceding the theorem is the
  number of generic preimages. Under a nonzero constant Jacobian, mapping
  degrees two and three are excluded. We read the author's primary PDF.
- **CITED degree four.** A. V. Domrina,
  [*On four-sheeted polynomial mappings of C2. II. The general case*](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=273&what=fullteng),
  *Izvestiya: Mathematics* **64:1** (2000), 1–33,
  [DOI 10.1070/im2000v064n01ABEH000273](https://doi.org/10.1070/im2000v064n01ABEH000273).
  Page 1 defines topological degree by generic preimages and states the
  exclusion for degree four. This is the general-case theorem, not Part I's
  additional dicritical restriction and not a theorem about polynomial
  coordinate degree. The publisher PDF was read directly.

No degree-five theorem or withdrawn connectedness claim is used.

## 2. The global ledger retains the actual cusp fibre

Suppose the stated support occurs. Let `a` be the generic actual affine
count over its smooth locus, and let `delta=d-a`. The finite-envelope
facts proved in the infinity supplier give

```text
d>1,                  1<=a<=d-1,
A2 --j(open)--> Xbar --pi(finite degree d)--> A2,
B=Xbar\j(A2) is pure divisorial.                         (1)
```

Here `Xbar` is normal. The finite map is etale off the whole support `S`.
The deleted generic contribution is a sum of positive integer lengths
`e_D f_D`; it is not a count of special boundary points. Briefly, purity
of the deleted set follows by removing its divisors, extending the two
affine source coordinates over the remaining codimension-two set by
normal Hartogs, and obtaining an inverse to the open inclusion there.

The actual count on the whole smooth stratum is `a`, as proved in the
supplier. This constancy includes retained versus deleted unramified
pages; it is stronger than inertia-cycle constancy. At each ordinary node
`p`, the exact local intersection of the two retained subsets yields

```text
#F^(-1)(p)=d-2delta+omega_p,
omega_p >= max(0,d-2a),       omega_p integer.             (2)
```

Write `n=#F^(-1)(z)` at the unique cusp. Its normalization has one point
over the cusp and two over each node. Thus

```text
chi_c(S_smooth)=-2N,      chi_c(A2\S)=N.
```

Integrating actual fibre counts over these strata gives

```text
1=dN-2aN+sum_nodes(d-2delta+omega_p)+n
 =n+sum_nodes omega_p.                                    (3)
```

In particular `n` is zero or one. This step inserts no cusp-intersection
formula. The cusp's actual affine count remains an independent coordinate
until the local argument below proves exactly what it counts.

## 3. Actual cusp subsets and re-access

Take a sufficiently small local neighborhood `D` of the cusp containing
no other singularity. Its complement of the cusp has a `d`-sheet cover
from (1). Let `Omega` be a labelled fibre and let its local monodromy be
the permutation action of

```text
G=<sigma,tau | sigma tau sigma=tau sigma tau>,
g=sigma tau,            g sigma g^-1=tau.                 (4)
```

The action is not assumed transitive. The cited presentation allows
`sigma` to be a based normal meridian at a smooth point of this same local
cusp. Its `a` actual retained pages define
`A subset Fix(sigma)`, `|A|=a`.

Define the second subset by **re-accessing that same smooth point**.
This construction is important: we do not equate an arbitrary meridian
conjugator with a preferred geometric walk along the cusp. Here is one
fully specified convention. Write path composition with the rightmost
path traversed first. Let `c` run from the basepoint to the nearby point
of a normal transversal, and let `m` be its small normal meridian. If
`T_c` denotes path transport and `A_eta` the actual subset at its endpoint,
then `A=T_c^-1 A_eta` and `sigma=T_(c^-1 m c)`. Now use

```text
c'=c g^-1,
T_(c')^-1 A_eta = g A =: A',
T_((c')^-1 m c') = g sigma g^-1 = tau.                    (5)
```

Thus `A'` is an actual transported subset of size `a` and lies in
`Fix(tau)`. Formula (5) simultaneously transports the subset and meridian;
it does not assume a peripheral-centralizer invariant. With another
standard path convention the action can be defined by inverse lifting,
giving the same simultaneous conjugation and subset transport.

We claim the exact equality

```text
                    n = |A intersect A'|.                (6)
```

If a label lies in the intersection, it is fixed by both generators in
(4). Its orbit under the full local group is therefore a singleton. The
finite normal extension of this degree-one complement cover is a copy of
`D`: a finite birational map onto the normal target is an isomorphism.
One can equivalently use uniqueness of normalization of the local cover.
On that copy, the generic point of the irreducible cusp is retained,
because the label belongs to `A`. Any deleted divisor on the copy would
have to lie over the cusp, so would be the whole cusp divisor. That is
impossible since its generic page is retained. The pure deleted boundary
in (1) then excludes deleting only the origin. The singleton sheet
therefore gives an actual affine cusp point.

Conversely, an actual affine point over the cusp has an etale local
inverse. After shrinking `D`, it supplies a singleton local-cover sheet
whose generic cusp page is retained. It is fixed by both generators and
belongs to both transported subsets. Distinct points give distinct
labels. This proves (6). Notice that being inertia-fixed alone would
not establish retention: the proof needs membership in the actual `A`.

## 4. A complete sparse braid lemma

**Abstract lemma.** Let permutations `sigma,tau` on a finite set `Omega`
of size `d` obey (4). Let
`A subset Fix(sigma)`, `1<=a=|A|<d`, and `A'=gA`. Put

```text
I=A intersect A',   U=A\A',   V=A'\A,
O=Omega\(A union A'),       n=|I|,       h=|O|.
```

If `n<=1` and `h<=1`, precisely the following numerical possibilities
survive, up to relabelling of the displayed permutations and subset:

| `d,a,n` | `sigma,tau` | Actual subset pattern |
|---|---|---|
| `2,1,1` | both identities | one common retained label, one outside label |
| `3,1,0` | two transpositions generating `S3` | `A` the fixed label of `sigma` |
| `4,2,1` | those transpositions on three labels and a fixed fourth | `A` contains that fourth label and the fixed label in the three-label block |

All three patterns occur abstractly. No local transitivity is required.

**Proof.** By (4), `A' subset Fix(tau)`. Every label in `I` is fixed
pointwise by both generators. Since `g` fixes `I` and takes `A` to `A'`,
it takes `U` bijectively to `V`. Also

```text
                 h=d-2a+n.                              (7)
```

If `h=0`, `sigma` acts only on `V` and `tau` only on `U`. They commute,
and their braid relation then forces `sigma=tau`. Since their supports
are disjoint, both are identities. Thus `A=A'` and `Omega=A`, contrary
to `a<d`.

Hence `h=1`; write `O={o}`. The permutation `sigma` fixes `I union U`
and acts only on `V union {o}`; `tau` fixes `I union V` and acts only on
`U union {o}`. If `x in U` and `tau x in U`, applying the braid relation
to `x` gives `tau x=tau^2 x`, so `tau x=x`. The label is then fixed by
both generators and by `g`. Since it belongs to `A`, it belongs to `gA`,
contradicting `x in U`. Consequently `tau U subset {o}` and `|U|<=1`.
The same holds for `V`, or follows from `gU=V`.

If `U` is empty, `a=n=1` and `d=2`; both permutations fix both labels.
If `U` has one element, `a=n+1` and `d=n+3`. Each generator swaps the
outside label with its unique nonfixed retained label, fixing `I`.
As `n` is zero or one, these are exactly the degree-three and degree-four
patterns in the table. This proves the lemma in every finite degree.

Apply it to the actual cusp data. By (3) and (6), `n<=1`. Choose any one
of the `N>=1` nodes. Equations (2) and (3) give

```text
d-2a <= omega_p <= 1-n,
h=d-2a+n <=1.                                            (8)
```

Thus `d` is two, three, or four. Orevkov's Theorem 1.1 excludes the first
two degrees, and Domrina's general-case theorem excludes the third.
This proves the support exclusion in §1. If `N>=2`, the integer ledger
already implies `d<=2a`: otherwise each node contributes at least one
to (3). This removes the three-sheet pattern before the literature step.
The two-sheet pattern can alternatively be excluded by `delta=1` and
weighted-boundary branch purity as in the infinity supplier. None of
these alternatives removes the genuine abstract four-sheet pattern.

## 5. Hostiles, the explicit family, and the connection's scope

The sharp abstract four-sheet control is

```text
Omega={1,2,3,4},  sigma=(12),  tau=(23),  g=(123),
A={3,4},        A'={1,4},   I={4},       O={2}.            (9)
```

It satisfies the braid relation, actual-subset cardinalities, and the
numerical cusp ledger with `n=1` and zero node overlaps. Its local action
has a three-sheet orbit and a singleton. Removing the common fixed label
gives the degree-three control with `a=1,n=0`, which is compatible with
the integer ledger for exactly one node. Hence neither the scalar Euler
identity nor the sparse braid lemma by itself closes the problem.
The cited low-mapping-degree theorems are essential to the stated finish.
These are abstract local passports, not assertions of globally realizable
polynomial Keller maps or actual global node subsets.

The frozen curve probe established the complete singularity inventory
for the intrinsic birational parametrization

```text
U(t)=t^4+t,        V_lambda(t)=t^6+t^2+lambda t.
```

Away from its three exceptional polynomial loci, it has six nodes. Its
four triple-point parameters and six tacnode parameters were excluded
by the single-exception ledger. Its remaining parameters are exactly

```text
C(lambda)=128lambda^3-288lambda-283=0,
```

each with one ordinary cusp and five ordinary nodes. The result here
excludes all three as whole nonproperness supports. Combined with that
already audited classification, **every parameter in this explicit
one-parameter family is excluded as a sole support**. This note supersedes
only the three-parameter OPEN boundary of that frozen family result.
It excludes neither all intrinsic `(4,6)` parametrizations nor arbitrary
Keller targets. The cubic is irreducible over `Q` since its reduction
modulo seven has no root, but no Galois transport is needed for this
uniform conclusion.

The connection can be recorded without losing its data:

| Coordinate | Content |
|---|---|
| Source | whole-support Euler ledger and the finite normal envelope |
| Target | a finite permutation action of the local cusp group |
| Map | transport actual smooth-cusp pages to a labelled complement fibre, then re-access the same point by `g^-1` |
| Preserved predicate | actual retention, generator fixation, subset cardinality, and the exact cusp count (6) |
| Lost information | global monodromy, source polynomial coordinates, all actual node subset labels |
| Needed sidecar | pure boundary, the whole-support node ledger, and the correctly typed low-degree theorems |
| Cheapest hostile | the four-sheet pattern (9), which survives every abstract local condition |

The local mechanism can be reused only when the singularity's marked
meridian generators and simultaneous actual-subset transport are known.
It does not license replacing an arbitrary singularity by a cusp group,
discarding boundary multiplicities, or assuming a generic target chart.

## 6. Exact controls and audit state

The source is
[planar_jc48_sep06_cusp_passport.py](../../04-computation/planar_jc48_sep06_cusp_passport.py).
It enumerates all ordered permutation pairs on `2<=d<=5`, retains exactly
the braid pairs, and visits every nonempty proper subset of the first
generator's fixed labels. The filters `n<=1` and `h<=1` are the declared
abstract lemma hypotheses. All labelled survivors are printed, including
the two identity, six three-sheet, and twenty-four four-sheet controls.
The expected labelled counts are also independently checked from their
structural descriptions. No local transitivity or quotient representative
is imposed. All gates remain active under `python -O`.

An independent integer-ledger control uses `2<=d<=40`, `1<=a<d`,
`1<=N<=8`, and `n in {0,1}`. The only remaining overlap mass is zero or
one, so `N max(0,d-2a)<=1-n` exactly decides whether integer allocations
exist. It checks (8), and `d<=2a` for `N>=2`. These bounded controls do not
prove the all-degree topological or algebraic-geometric statements.

```sh
python3 04-computation/planar_jc48_sep06_cusp_passport.py
python3 -O 04-computation/planar_jc48_sep06_cusp_passport.py
```

Normal and optimized outputs are byte-identical: **16,815 always-active
gates**, with **32 raw labelled sparse passports** and **6,419 eligible
integer-ledger rows**. The latter are a separate relaxation: they need not
have nonnegative `d-2a+n`, and are not claimed to be actual permutation
passports. Source and output are frozen.
The scope of the proof remains §1, and the scope of its parameter
application remains the displayed family.

```text
source SHA256 8a7d40c0f001a2063f976b949d49033087adc3ae1616072c1c644172918e3d6f
output SHA256 6c330f724ac39a6960e22d13bee1d54efcd3cd13ff891acb470de41eccd04d95
```

The saved output is
[planar_jc48_sep06_cusp_passport.out](planar_jc48_sep06_cusp_passport.out).

The complete independent
[analytic and source audit](planar_jc48_sep06_cusp_passport_audit.md)
passes the same-point re-access convention, actual intersection equality,
normal-extension and pure-boundary step, every-degree sparse lemma,
`N>=1` ledger, primary theorem degree types, and exact family application.
Its independent normal/optimized/frozen replay matches all 16,815 gates
and both source/output hashes. No mathematical correction was required.

```text
independent audit SHA256 308094dd7e3df9dff412b024e7aec2df042d0d08983bd33782f4aec482e5ff2f
```

# THM-2892 structural-compression audit

**Status:** scratch audit of a concurrent `RESERVED / UNPROVED PROOF
CANDIDATE`; no THM-2892 result is claimed here.

## 1. Inheritance and exact mechanism

The closest proved mechanism is
`THM-2893-complement-cap-finite-core-flag-lemma.md`, specialized to

```text
(k,s,ell)=(4,2,3).
```

The canonical near miss is an `H`-only `K4` census: THM-2893 forces only
three of the four hypothetical covering speeds into the finite high core.
The fourth speed may lie outside `H`, so the theorem-bearing object is a
literal residual behind a heavy triangle, not a `K4` inside `H`.

The candidate's closure mechanism is:

1. THM-2888 gives a global pair cap `B<5h/7` on each honest first-apex
   carrier, after a forced-edge/avoiding-edge dichotomy on 26 carriers.
2. With `theta=h-B`, the symmetric heavy graph has edge `xy` precisely when
   `U_C({x,y})>=theta`.  Its finite high core is
   `H={x:c_C(x)>=theta/2}`.
3. Any four-cover induces a heavy `K4`; every heavy edge meets `H`; hence
   at least three vertices of the `K4` lie in `H` and form a heavy triangle.
4. The concurrent exact atlas enumerates 372209 such triangles.  Every
   literal residual is nonempty and no single further danger comb contains
   it.  Thus no heavy triangle extends to a four-cover.
5. All 10202 post-THM-2888 carriers become terminal.  Returning them to the
   weighted root gate gives a positive strict margin on all 2508 active
   roots.  THM-885 supplies the disjoint 495 low roots.

This explains *why* the finite verdict closes the rung: the cap forces a
three-vertex flag, and the fourth-slot question is exactly singleton
containment of the literal triangle residual.

### Hypergraph-link formulation

The proof can be stated more compactly at the consequence-object level.
For each body `E`, define its exact-cover hypergraph

```text
X_E={Q subset {15,16,...}: |Q|=5 and G_(E union Q) is empty}.
```

THM-2885 says every edge of `X_E` meets its global top-ten apex bank.
Every already terminal apex has empty link in `X_E`.  For a remaining apex
`a`, an edge in `link_(X_E)(a)` is a four-cover of the literal carrier
`G_E minus D_a`.  THM-2893 maps such an edge to a heavy triangle in the
finite core.  The singleton-residual certificate proves that no heavy
triangle extends to an edge of this link.  Hence every remaining apex link
is empty, and therefore `X_E` itself is empty.

The 26 normalized carriers simply split one apex link into:

```text
four-sets containing the forced pair  -> closed by THM-2888,
four-sets avoiding the forced pair    -> governed by the deleted pair cap.
```

This link language exposes the recursion without changing any quantifier:
the computation empties all possible links of the exact-cover hypergraph,
not merely a scalar proxy.

## 2. A smaller analytic singleton certificate

The discrepancy estimate is not needed to seal the final singleton tail.
There is a strictly smaller geometric certificate.

### Longest-component tooth lemma

Let

```text
R = union_j I_j,       I_j=(a_j,b_j),       m=|R|>0,
ell=max_j(b_j-a_j).
```

For the LRC danger comb

```text
D_w = union_(n in Z) [n/w-1/(14w), n/w+1/(14w)],
```

put

```text
W_R=floor(1/(7 ell)).
```

If `R` is contained in `D_w` up to endpoints, every connected interval
`I_j` lies in one tooth, because distinct teeth are separated by a positive
gap.  A tooth has length `1/(7w)`, so necessarily

```text
w <= W_R.                                                (A)
```

For a finite candidate `w<=W_R`, containment of `I_j` in one tooth is
equivalent to the existence of an integer `n` with

```text
w*b_j-1/14 <= n <= w*a_j+1/14.                           (B)
```

Thus `I_j` is an exact noncontainment witness precisely when

```text
ceil(w*b_j-1/14) > floor(w*a_j+1/14).                    (C)
```

Consequently a triangle residual is singleton-closed if and only if, for
every allowed

```text
15 <= w <= W_R
```

outside the selected labels, some component satisfies `(C)`.  Speeds above
`W_R` are excluded analytically by `(A)`.

This certificate uses one failed interval/tooth incidence rather than a
full scalar coverage.  It also has the universal comparison

```text
ell >= m/r
W_R <= floor(r/(7m)),
```

whereas the candidate's discrepancy horizon is

```text
floor((99/70)r/(6m)) = floor(33r/(140m)).
```

The geometric coefficient is `1/7=20/140`, only `20/33` of the discrepancy
coefficient before floors.

### Exact controls already obtained

The independent 50-carrier pilot checked 212 triangle residuals:

```text
empty residuals                 0
singleton covers                0
maximum geometric tail first   53
maximum discrepancy tail first 165
```

It checked the interval predicate against both literal subtraction and
scalar full coverage at every positive candidate and at scan boundaries.
An independent primary/alternate-path comparison also reproduced all 279
high vertices, 434 heavy edges, 212 triangles, and 212 singleton closures
under byte-identical normal and optimized runs.  In this sample all 212
triangle residuals are distinct even after keying by carrier; literal
residual deduplication gives no compression there.

On the candidate's nominal minimum-margin row

```text
E=(1,5,6,7,8,9,11,13), a=19, L=(37,121,130),
m=1099999/24444420, r=42,
```

the exact longest interval is

```text
ell=911/220220,
W_R=34,
```

instead of discrepancy horizon 219.  The exact maximum over the geometric
head is still speed 25 with

```text
coverage=7111/459800,
m-coverage=137171599/4644439800.
```

Hence the printed `1/285184900` discrepancy-tail margin is a weakness of
that scalar certificate, not evidence of an almost-covering fourth comb.

### Theorem-ready formulation

The candidate may replace its singleton paragraph by:

> For every enumerated heavy triangle `L`, the literal residual `R_L` is a
> nonempty finite union of rational intervals.  Let `ell_L` be the longest
> interval and `W_L=floor(1/(7 ell_L))`.  Exact integer arithmetic verifies
> `(C)` for at least one component of `R_L` for every allowed
> `15<=w<=W_L`.  For `w>W_L`, the longest component is wider than a danger
> tooth, so `R_L` cannot be contained in `D_w`.  Thus no heavy triangle
> extends to a four-cover.

This removes the THM-735/THM-2885 discrepancy estimate from the *second*
triangle-residual tail.  The discrepancy estimate remains load-bearing in
the earlier construction of the finite high core `H`.

A full 372209-row replay of `(C)` would make this an independent replacement
certificate.  Until that replay is frozen, the current scalar atlas remains
the proved candidate evidence and this lemma is a proof compression, not a
substitute transcript.

## 3. Why a finite phase atlas remains load-bearing

The longest-component lemma makes the tail analytic, but it does not decide
the finite head.  The statistics

```text
(mass, number of components, component lengths)
```

cannot do so: translating an interval of length below `1/(7w)` from the
center of a tooth to the center of a gap preserves all those statistics but
changes whether it lies in `D_w`.  The missing coordinate is endpoint phase
modulo `1/w`, exactly the data in `(B)` and `(C)`.

Likewise, vertex coverages do not determine heavy edges:

```text
U_C({x,y})=c_C(x)+c_C(y)-|C intersect D_x intersect D_y|
```

retains an overlap term which scalar ranking forgets.  After a triangle is
chosen, the residual endpoint phases retain still more information than the
weighted heavy graph.

Therefore no currently proved statistic-level lemma makes the 372209
triangle rows redundant.  What can be compressed is the certificate for
each row:

```text
weighted heavy graph
  -> heavy flag
  -> literal rational-interval residual
  -> finite list of failed interval/tooth incidences.
```

An analytic removal of the atlas would require a new phase theorem
classifying the endpoint residues of all heavy-triangle residuals.  Neither
THM-2888 nor THM-2893 supplies it.

## 4. Equality and scope audit

The load-bearing boundaries are:

| Gate | Correct boundary |
|---|---|
| Pair cap | `B<5h/7`; equality does not make the high core finite. |
| Heavy edge | `U_C({x,y})>=theta`; equality is included. |
| High vertex | `c_C(x)>=theta/2`; equality is included. |
| Triangle residual | `m>0`; `m=0` is an extending three-cover, not a closure. |
| Geometric horizon | scan `w<=floor(1/(7ell))`; equality may align with one tooth. |
| Finite singleton gate | full containment is equality `c_R(w)=m` and is an actual extension; closure needs strict noncontainment. |
| Weighted root gate | the complement margin must be `>0`; equality remains open. |

The candidate correctly reports zero empty triangle residuals.  This is
load-bearing in light of the current THM-2893 wording, which permits
recursion only on nonempty literal residuals.

In the current discrepancy proof, choosing

```text
W=floor((99/70)r/(6m))
```

already makes the numerical tail cap at `W+1` strictly below `m`.  A
non-strict source estimate would still suffice at that singleton step.
The source estimate's strictness is load-bearing earlier when the finite
`H` cutoff lands at equality.  Candidate wording should therefore say that
the **positive tail-cap margin** in its equation (13) is essential, not that
the first strict discrepancy sign is itself essential there.

Before promotion, `depends_on: []` must be repaired.  Minimal dependencies
depend on the stated conclusion:

- for the exact statement `E in C({1,...,14},8)` and five external speeds:
  THM-2888 and THM-2893.  THM-2888's proved statement already imports the
  THM-885 low chamber; because the candidate prose invokes THM-885 directly,
  either list THM-885 as a direct dependency or rewrite that sentence as an
  inherited THM-2888 conclusion;
- for the larger composed “at least eight speeds in `{1,...,14}`” sector:
  also THM-741 for the nine-or-more-low chamber.

THM-2885 is inherited through THM-2888 and may remain `related` unless the
proof cites its gate directly.  The 26 forced-edge normalizations are
proved THM-2888 branches, not a fresh unconditional pair cap.

The honest non-consequence remains:

```text
v_8>=15  (at most seven normalized speeds in {1,...,14})
```

is outside the theorem, so THM-2892 would not prove unrestricted LRC(14).

## 5. Tournament validity gate

There is no intrinsic tournament here.

```text
vertices:            allowed external speeds
pairwise observable: U_C({x,y})
relation:            undirected threshold U_C({x,y})>=theta
ties:                exact threshold equality is a heavy edge
preserved target:    a four-cover induces a heavy K4
lost by orientation: overlap mass, literal residual, and support symmetry
needed sidecar:      the full residual interval set
```

Orienting by `c_C(x)>c_C(y)` produces only a preorder (with genuine ties)
and forgets the overlap term deciding heaviness.  A tie-break would be
cosmetic.  THM-2894 proves the stronger structural reason: literal deletion
factors through the commutative idempotent finite-support semilattice, so it
contains no intrinsic order, and every union-multiplicative group image is
trivial.

The faithful combinatorial object is therefore the **residual-decorated
flag complex** of the heavy graph:

```text
simplex/flag L  |->  C minus union_(w in L) D_w.
```

Inclusion of flags reverses inclusion of residuals.  This order-reversing
semilattice map preserves the covering target exactly; a tournament quotient
would add information not present in the carrier and discard information
the singleton decision needs.

## 6. Recommended disposition

1. Promote the finite-exact THM-2892 candidate only after its locked normal,
   optimized, stored-output, hash, and independent scope audits pass.
2. Repair dependencies and the singleton strictness wording above.
3. Add the longest-component tooth lemma either to the theorem as an
   analytic explanation or as a small independent companion.
4. If replacing the scalar triangle tail, freeze a full exact replay of
   condition `(C)` on all 372209 residuals.  The 50-carrier pilot and the
   nominal sharp row are positive controls, not the theorem universe.
5. Do not add a tournament.  Retain the undirected weighted heavy graph,
   its flag complex, and the literal residual sidecar.

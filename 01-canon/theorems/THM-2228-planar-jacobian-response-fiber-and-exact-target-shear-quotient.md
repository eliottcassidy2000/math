---
id: THM-2228
title: "Planar Jacobian response fibers are exactly target-shear orbits"
status: >
  PROVED + CITED. Over any characteristic-zero field, if P has one
  constant-Jacobian mate, then ker(Q -> Jac(P,Q)) is exactly k[P].
  Consequently every nonempty Jacobian-response fiber is one affine
  target-shear orbit, the response induces a linear isomorphism from the
  complement space modulo k[P] onto its image, and all constant-response
  complements are explicitly alpha B+k[P]. The kernel theorem is the
  Cheng--McKay--Wang theorem over characteristic-zero fields as proved in
  Moskowicz's Theorem 2.3; the fiber, quotient, scalar-response splitting,
  and normal-form consequences are proved here. This is an exact gauge
  theorem, not a proof of planar JC.
source: codex-2026-07-24-planar-jacobian-response-fiber
depends_on:
  - "Moskowicz, The two-dimensional Centralizer Conjecture, arXiv:1802.04685v2, Theorem 2.3"
related:
  - THM-2084-cubic-fiber-low-complement-gauss-manin-gate
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
  - THM-2181-exact-square-prefix-compression-and-monic-depressed-quartic-closure
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
---

# THM-2228 -- the Jacobian response is the exact target-shear quotient

Let `k` be a field of characteristic zero and put

```text
V=k[x,y],
D_P:V -> V,                    D_P(Q)=Jac(P,Q).       (1)
```

The point of this record is not a new centralizer theorem. It is the exact
response interpretation of that theorem: once `P` has a Jacobian mate, the
usual reduction of its complement modulo polynomial target shears loses
precisely the information invisible to `D_P`, and no more.

## 1. Exact response-fiber theorem

Assume there are `B in V` and `kappa in k*` such that

```text
D_P(B)=kappa.                                        (2)
```

Then:

1. the zero-response fiber is

   ```text
   ker D_P=k[P];                                     (3)
   ```

2. for every `R in im D_P` and every one choice `Q_R` with
   `D_P(Q_R)=R`,

   ```text
   D_P^(-1)(R)=Q_R+k[P];                             (4)
   ```

3. regarding `k[P]` as a `k[P]`-submodule of `V`, the induced response

   ```text
   Dbar_P: V/k[P] -> im D_P,
   [Q] |-> D_P(Q)                                    (5)
   ```

   is a `k[P]`-module isomorphism;

4. every constant-response fiber is explicit:

   ```text
   {Q:D_P(Q)=lambda}=(lambda/kappa)B+k[P]
                                            for every lambda in k. (6)
   ```

5. after normalizing `R_0=B/kappa`,

   ```text
   {Q:D_P(Q) in k}=k[P] direct_sum k R_0,             (6a)
   D_P(H(P)+aR_0)=a.
   ```

In particular, for a fixed nonzero `lambda`, two complements `Q_1,Q_2`
of `P` have the same constant Jacobian if and only if

```text
Q_2=Q_1+H(P)                         for some H in k[T]. (7)
```

Thus the polynomial target shears fixing the first output,

```text
(P,Q) |-> (P,Q+H(P)),                               (8)
```

act simply transitively on every response fiber. The polynomial `H` is unique:
(2) makes `P` nonconstant, so substitution `k[T] -> k[P]` is injective.
Equivalently, `D_P(Q)` is a complete invariant of the target-shear orbit of
the complement once `P` has a mate.

## 2. Proof

The Cheng--McKay--Wang centralizer theorem over characteristic-zero fields
states that, under (2),

```text
D_P(w)=0  implies  w in k[P].                        (9)
```

Moskowicz proves this field form directly. The reverse inclusion is elementary:
the chain rule gives

```text
D_P(H(P))=H'(P)D_P(P)=0.                             (10)
```

Equations (9) and (10) prove (3).

If `D_P(Q)=D_P(Q_R)=R`, linearity of the Jacobian in its second entry gives

```text
D_P(Q-Q_R)=0.
```

By (3), `Q-Q_R in k[P]`, proving one inclusion in (4); (10) proves the
other. The map in (5) is therefore well-defined and injective, and it is
surjective by the definition of `im D_P`. It is `k[P]`-linear because

```text
D_P(F(P)Q)=F(P)D_P(Q).                               (5a)
```

Finally,
`D_P((lambda/kappa)B)=lambda`, so (4) gives (6). This also proves (7).
For (6a), every scalar-response `Q` belongs to the fiber in (6), and the sum
is direct because applying `D_P` to `aR_0 in k[P]` forces `a=0`.
QED.

## 3. What complement minimization now means

Many planar-Jacobian arguments in this repository fix `P` and define a
reduced complexity such as

```text
min_(H in k[T]) deg_y(Q-H(P)).                        (11)
```

Under (2), (4) says that (11) minimizes over the **entire** collection of
complements with the same Jacobian response. It is not merely a convenient
family of legal normalizations. Conversely, there is no second,
non-target-shear complement with the same response that a reduced-degree
argument has silently omitted.

Moreover, if `D_P(Q)=D_P(B)=kappa`, then (7) gives

```text
k[P,Q]=k[P,B].                                       (11a)
```

Thus replacing a mate by a shear-reduced mate changes neither the generated
subalgebra nor whether the pair is a polynomial automorphism. For a different
nonzero scalar response, rescaling the second output and then shearing gives
the analogous statement.

This applies directly to the cubic and quartic source-fiber programs:
their Faber and pole descents choose representatives in the exact fiber (4).
A theorem uniform over the shear-reduced representatives is therefore
uniform over every mate of the fixed `P` with that Jacobian. What remains
open is whether such a `P` must be a coordinate; response-fiber uniqueness
does not imply polynomial invertibility.

The quotient in (5) is a quotient of `k[P]`-modules (and hence of
`k`-vector spaces), not a quotient ring: `k[P]` is generally not an ideal of
`k[x,y]`. This type distinction is load-bearing when the response quotient
is compared with other context compressions in the repository.

## 4. Equality and failure boundaries

Two small examples isolate both hypotheses.

### Positive control

For `P=x`, `B=y`, and `kappa=1`,

```text
D_x(Q)=Q_y.
```

Hence

```text
{Q:D_x(Q)=lambda}=lambda y+k[x],
```

exactly as (6) predicts.

### A mate is necessary

Over characteristic zero, take `P=x^2`. Then

```text
D_(x^2)(Q)=2x Q_y,
ker D_(x^2)=k[x]  strictly contains  k[x^2].          (12)
```

There is no constant-Jacobian mate because every response in (12) is
divisible by `x`. Thus the existence of the unit response in (2), not merely
the formal derivation `D_P`, collapses the centralizer to `k[P]`.

### Characteristic zero is necessary

Over a field of characteristic `p>0`, take `P=x` and `B=y`. Although
`D_x(B)=1`, one has

```text
D_x(y^p)=0,                  y^p notin k[x].           (13)
```

So (3)--(7) fail in positive characteristic.

### The field scope is deliberate

For a characteristic-zero integral domain `D`, Moskowicz's proved domain
statement places a centralizing polynomial only in `Frac(D)[P]`; membership
in `D[P]` is her two-dimensional Centralizer Conjecture over `D`. Thus the
field hypothesis here cannot be silently weakened to an arbitrary domain.

## 5. Response-quotient perspective

The general response pattern is:

```text
object  --response-->  observable,
fiber of response = information discarded by the observable.       (14)
```

In the coin extractor of THM-2225 one must retain the image of a cyclic
checksum response, not only its variance. In the tournament transport of
THM-2221 one must retain the cut-metric sidecar. Here the favorable planar
Jacobian phenomenon is stronger: (3) computes the response kernel exactly,
so (5) is already the minimal lossless quotient for this observable.

This suggests a precise audit rule for future planar normal forms: whenever
a proposed reduction of the complement is coarser than target-shear
equivalence, it must name the extra information it destroys, because the
Jacobian response itself does not justify that additional quotient.

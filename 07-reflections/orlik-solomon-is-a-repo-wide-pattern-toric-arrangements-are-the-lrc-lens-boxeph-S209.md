# The LRC arrangement lens is phase-height, labelled, and oriented

*Scope-repaired audit of boxeph-2026-07-21-S209 / HYP-8830.*

> **Correction:** MISTAKE-224 retracts the identification of the Fourier
> relation lattice with a toric layer poset, of `G_delta` with an ordinary
> toric complement, and of its Fourier series with an arithmetic-Mobius sum.
> The braid and Shi controls survive. The useful repair is a different finite
> toric arrangement in `(time,height)`.

## Four objects that must stay distinct

For a primitive speed vector `v=(v_1,...,v_n)`, define

```text
phi_v : T -> T^n,       t |-> (v_1 t,...,v_n t),
H_v   = image(phi_v),
L_v   = {k in Z^n : k.v=0},
F_delta = {x in T : ||x|| >= delta}.
```

Then `L_v` is the character annihilator of the orbit subtorus `H_v`, and

```text
G_delta(v) = phi_v^{-1}(F_delta^n).
```

Character orthogonality on `H_v` gives THM-1820's exact Fourier identity

```text
|G_delta(v)| = sum_{k in L_v} product_j ghat_delta(k_j).
```

This is an infinite, sinc-weighted, `delta`-dependent series. Its summation
index is not the layer poset of a finite toric arrangement. A toric layer is a
connected component of an intersection of finitely many character hypertori.

The distinction is not cosmetic. For `S={1}`, `|G_delta|=1-2delta`, whereas
the ordinary complement of the character kernel `{1}` in the circle has Haar
measure `1`. Thus thickening a resonance point changes the measured object.

The **full embedded** lattice `L_v subset Z^n` is also stronger than S209 said:
its rational orthogonal complement recovers the line spanned by `v`, hence the
primitive vector up to scale and coordinate relabelling. What loses the LRC
predicate is an abstract lattice isomorphism type, a coefficient truncation,
a small-relation count, or a Betti/Mobius summary—not necessarily `L_v`
itself.

## The exact constructive replacement

Put coordinates `(t,delta)` on the cylinder `T x [0,1/2]` and use the finite
character list

```text
X_S = {(v,+1),(v,-1) : v in S} subset Z^2.
```

Its walls are

```text
v t + sigma delta in Z,   sigma in {+1,-1}.
```

The exact LRC feasible region is the closure of the selected inequality cells

```text
E_S = {(t,delta) : 0 <= delta <= 1/2 and delta <= ||v t|| for every v in S}.
```

It is **not** the complement of all the walls. The height functional is part
of the object: `M(S)=max{delta:(t,delta) in E_S}`. Each constraint wall also
carries its runner owner, sign, and the cellwise selected side in the strip.
The bottom `delta=0` is ambient boundary, and a periodized torus wall has no
global left/right side.

After periodizing the height coordinate, two independent full-torus
characters `(v,sigma)` and `(w,tau)` have intersection index

```text
|det((v,sigma),(w,tau))| = |v tau-w sigma|.
```

Consequently:

- same-sign walls contribute `|v-w|`;
- opposite-sign walls contribute `v+w`;
- the two walls of one runner contribute `2v`.

This full-torus index is arithmetic, not the number of intersections retained
by `0<=delta<=1/2`; dependent characters have a non-proper intersection.

At a strict top of the lower envelope, an increasing and a decreasing active
branch meet, or one runner is at its cusp. These are precisely the opposite-
slope cases, so the denominators divide `v+w`, with the cusp as `v=v`. This is
the geometric content of the proved pair-sum denominator lemma in
[THM-1002](../01-canon/theorems/THM-1002-pair-sum-denominator-bound-and-the-bounded-gap-case.md).
The arrangement lens therefore explains an existing exact theorem without
replacing it by topology.

Deleting a speed removes its **pair** of signed walls. That makes paired
deletion/restriction a natural experiment for the AP-core supplier. But a
faithful state must retain:

```text
wall owner + sign + selected side + attained height + deletion identity.
```

Ordinary Orlik--Solomon cohomology or an arithmetic Tutte polynomial is not
currently known to recover that selected state. THM-2047 now proves this
phase-height dictionary, the top/pair-sum mechanism, and exact paired
deletion; its localization route toward an AP core is open, not a proved LRC
step.

## What survives from the original computation

Two classical control families remain exact and useful:

- the braid complement has Poincare polynomial
  `product_(k=1)^(n-1)(1+k t)`, with unsigned Stirling Betti numbers; its real
  chambers are the `n!` total orders;
- the Shi arrangement has characteristic polynomial
  `q(q-n)^(n-1)` and `(n+1)^(n-1)` real regions.

Neither statement supplies an LRC map. In particular, braid cohomology belongs
to one configuration space, not to each tournament, and the Shi walls
`x_i-x_j in {0,1}` are not the signed safety walls `v t +/- delta in Z`.

The script also exactly finds that `(1,2,3)` uniquely maximizes a `B=2`
relation count among 72 hard-coded primitive triples. That is a finite scout,
not a theorem about Betti or Mobius mass. The natural one-circle hypertorus
control goes the other way: the union for `(1,2,3)` has four points, while
`(2,3,4)` has six.

## Live question

Can paired deletion/restriction in the labelled phase-height complex, with a
chamber-height sidecar, force the AP maximum-deletion core required by
THM-1017? A successful argument must show how the top cell and its owner data
survive localization. A failure should identify the first pair of speed sets
with the same proposed quotient but different top-cell feasibility.

Primary background: [De Concini--Procesi on finite toric arrangements](https://arxiv.org/abs/math/0505351),
[Moci on the multiplicity Tutte polynomial](https://arxiv.org/abs/0911.4823),
and [Stanley's arrangement survey](https://math.mit.edu/~rstan/papers/nas.pdf).

Artifacts: corrected HYP-8830; MISTAKE-224;
`04-computation/orlik_solomon_across_the_repo_boxeph_S209.py` and its frozen
output.

---
id: THM-965
title: T4 strip-line decomposition — exact lattice fibers and H-free atom product, open uniform tail estimate
status: PARTIAL EXACT REDUCTION + FINITE FLOATING-POINT PROBES — the k-fibers and H-free product of two-pole masses are exact, but the claimed T4(H) = O(L^3/H) estimate is not proved; the probes cap the coefficient and k ranges, have no tail bound, and their three advertised dissociated examples all have short exact relations
source: kind-pasteur-2026-07-17-S128 cont.41/44; rigor correction after upstream audit
depends_on:
  - THM-946 (corrected two-pole atom and the strip/slab problem statement)
  - THM-952 (punctured near-pole congruence mechanism at support three)
related:
  - THM-935/948 (the E_s frame and exact packet audits)
  - THM-968 (the still-open T5 plane scheme)
scripts: 04-computation/t4_strip_referee_kps_S128c41.py -> 05-knowledge/results/t4_strip_referee_kps_S128c41.out; 04-computation/t4_composed_atom_referee_kps_S128c44.py -> 05-knowledge/results/t4_composed_atom_referee_kps_S128c44.out
---

# THM-965 — the T4 strip-line scheme

## Exact structural reduction

For fixed outer coordinates `(u,t)`, put

```text
k = v3*u + v4*t.
```

For each fixed integer `k`, the solutions of this equation are either empty or
one affine coset of the kernel of `(u,t) -> v3*u+v4*t`.  Writing
`g34=gcd(v3,v4)`, the primitive kernel direction is

```text
(v4/g34, -v3/g34).
```

Consequently the strip `|k| <= K0` is a finite family of congruence-coset
lines.  This fiber statement is elementary and exact.  It is useful because it
identifies the arithmetic carrier that any proof must price.

## Exact H-free composed-atom reduction (cont.44)

There is a further exact separation which survives the audit.  Write the
four-variable relation as

```text
v1*h1 + v2*h2 = -k,
v3*u  + v4*t  =  k.
```

For fixed `k` the solution fibre is the Cartesian product of the two affine
two-pole fibres.  If `S12(-k)` and `S34(k)` denote their punctured absolute
Fourier masses, then dropping the height cutoff gives the rigorous envelope

```text
T4(H) <= sum_k S12(-k) * S34(k).                    (1)
```

Applying THM-946's proved two-pole atom to both factors yields a summable
H-free majorant: each atom is `O(log(2+|k|)/|k|)` away from its pole, so their
product is `O(log(2+|k|)^2/k^2)`.  This is genuine structural progress over a
naive planar estimate and explains why the correct carrier is the shared
relation value `k`, with two affine pole fibres attached.

The cont.44 script checks a finite `BOX=200` truncation against a correspondingly
truncated composed bound.  Its displayed bounds are between roughly `465` and
`2063` for true masses below `2e-4`, so this is a coarse structural sanity check,
not evidence for a sharp constant or for `1/H` decay.

## What is not proved

The desired estimate

```text
T4(H) <= C4 * L^3 / H
```

does not follow merely from the line decomposition or from (1).  The cutoff
`max(|h1|,|h2|,|u|,|t|)>H` does not factor into the product of the two H-free
atoms.  Under `H`-dissociation the two fibres cannot both contain a point in
the height box, but converting that exclusion into a uniform tail price still
requires a shifted two-pole tail lemma, including near-coordinate-zero collars,
gcd and zero-coordinate cases.  Invoking an `H`-dissociation floor is legitimate
only after the speed tuple is actually shown to have no nonzero relation of the
required height.  That hypothesis was not checked in either finite probe.

Thus the strip/complement recursion is an **open proof scheme**, not a
coarse-constant theorem.  The T5 slab discussion formerly attached here is
likewise only a scheme; it is recorded separately as THM-968.

## Exact scope of the finite probe

`t4_strip_referee_kps_S128c41.py` computes truncated floating-point sums only:

- every coefficient is restricted to `[-250,250]` (`BOX=250`), so there is no
  bound for the omitted infinite tail;
- `K0` is replaced by `min(K0,4000)`, changing the proposed strip for the first
  and third test tuples;
- the Fourier weights are ordinary floating-point values, with no interval or
  rounding certificate;
- at `H=160`, the reported mass is only the shell `160 < max|h_i| <= 250`.

The three tuples called “dissociated” in the script are not uniformly
dissociated over the tested horizons.  They satisfy the exact relations

```text
 5*307  - 3*425  + 2*541  - 2*671  = 0,
 3*800  + 5*944  - 2*1413 - 2*2147 = 0,
13*541  - 7*1087 + 2*1943 -   3310 = 0.
```

The first two already violate height-10 dissociation; the third violates the
height-40 and height-160 tests.  Therefore the reported bounded truncated
envelopes and strip percentages are finite exploratory observations only.
They neither establish the dissociation-floor premise nor control the tail.

## Remaining proof obligations

- prove a uniform shifted-line estimate with explicit dependence on gcd and
  successive-minimum data;
- exploit the exact `k`-fibre product (1) to prove that when the two factors
  cannot both meet the height box, one factor pays a uniform `1/H` tail;
- separate genuinely `H`-dissociated tuples from the structured-relation
  branch before using the `1/H` floor;
- bound the complement and the coefficient tail beyond `BOX`;
- assemble the uncapped `k`-sum with constants uniform in the speeds.

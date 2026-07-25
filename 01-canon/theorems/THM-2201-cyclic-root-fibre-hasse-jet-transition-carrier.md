---
id: THM-2201
title: "Cyclic root-fibre Hasse jet is a faithful triangular transition carrier"
status: >
  PROVED + VERIFIED-EXACT. Over F_13, the group algebra of the thirteen root sheets is the
  truncated local algebra F_13[epsilon]/(epsilon^13). The binomial/Hasse jet
  is an invertible Pascal transform of sheet incidence, deck translation is
  triangular unipotent, and labelled endpoint deletion/insertion is an exact
  affine jet update. Augmentation alone detects THM-2198's first-image
  occupancy because it is at most ten; the full thirteen-sheet fibre is
  instead epsilon^12 and is invisible to augmentation. Five labelled unit
  mask polynomials retain owner counts without modular wrap, but the rooted
  guard block and event/winding word are still needed to choose the updates.
  This supplies a faithful finite transition representation, not the missing
  transition theorem and not a proof of LRC(14).
source: codex-2026-07-24-root-fibre-hasse-jet
depends_on:
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2200-convex-semigroup-and-finite-place-support-hole-trichotomy
  - THM-840
  - THM-853
script: 04-computation/lrc14_root_fibre_hasse_jet_thm2201.py
output: 05-knowledge/results/lrc14_root_fibre_hasse_jet_thm2201.out
script_sha256: 0988095205ac3488946d1e1314ef941ce88d31d25006c57027714d69b508bd22
output_sha256: a5c1facfb2eed6aa4023d03e86903999b5adc6ceda78a185b7e96841511b4613
hash_basis: working-tree bytes (LF)
---

# THM-2201 -- the cyclic root-fibre Hasse jet

THM-2198's first image remembers only whether a root fibre is occupied.
The first place that this can fail modulo thirteen is not mysterious: it is
the socle of the cyclic group algebra.

Let `C_13=<g>` label the thirteen rooted sheets and put

```text
R=F_13[C_13]=F_13[g]/(g^13-1).                       (1)
```

In characteristic thirteen,

```text
g^13-1=(g-1)^13.
```

Writing `epsilon=g-1` therefore gives an algebra isomorphism

```text
R isomorphic F_13[epsilon]/(epsilon^13).              (2)
```

## 1. The full Hasse jet reconstructs sheet incidence

For a sheet-weight vector `f=(f_0,...,f_12) in F_13^13`, encode

```text
Phi(f)=sum_(k=0)^12 f_k g^k in R.                     (3)
```

Since `g=1+epsilon`,

```text
Phi(f)=sum_(j=0)^12 J_j(f) epsilon^j,
J_j(f)=sum_(k=0)^12 f_k binom(k,j).                   (4)
```

Here `binom(k,j)=0` for `j>k`. The Pascal matrix

```text
(binom(k,j))_(0<=k,j<=12)                             (5)
```

is triangular with diagonal entries one. It is invertible over `F_13`.
Consequently

> **Hasse-jet faithfulness.** The thirteen coordinates
> `(J_0(f),...,J_12(f))` reconstruct the entire rooted sheet vector `f`.

This is a change of basis, not a quotient.

## 2. Translation and endpoint moves are triangular

Translation of the rooted sheet cycle by `a in F_13` multiplies (3) by
`g^a=(1+epsilon)^a`. Hence

```text
J_j(g^a Phi(f))
 =sum_(r=0)^j binom(a,j-r)J_r(f).                    (6)
```

The update is triangular unipotent: its diagonal coefficient is
`binom(a,0)=1`.

Let `delta_u` be the unit sheet vector at `u`. Deleting an endpoint at `u`
and inserting one at `v` sends

```text
f -> f-delta_u+delta_v.                               (7)
```

Equation (4) gives the exact affine jet update

```text
J_j -> J_j-binom(u,j)+binom(v,j),       0<=j<=12.    (8)
```

Thus every named endpoint transport has a finite exact representation in
the Hasse coordinates.

## 3. The top augmentation coordinate is the full-fibre flag

The augmentation is

```text
aug(Phi(f))=Phi(f)|_(g=1)=J_0(f)=sum_k f_k.           (9)
```

For a subset `Z subset F_13`, take `f=1_Z`. If `|Z|<=10`, then

```text
J_0(1_Z)=0 in F_13             iff             Z=emptyset. (10)
```

This is precisely THM-2198's first-image situation: its guard restriction
gives root-fibre occupancy at most ten, so reduction modulo thirteen cannot
erase a nonempty fibre.

For the full fibre, however,

```text
N=1+g+...+g^12
  =(g^13-1)/(g-1)
  =epsilon^12 in R.                                  (11)
```

Therefore

```text
J_0(N)=...=J_11(N)=0,             J_12(N)=1.         (12)
```

The occupancy-thirteen obstruction is exactly the top
augmentation-filtration coordinate. Augmentation forgets it; the full jet
does not. This explains algebraically why THM-2198's mod-thirteen support
test is exact on the first image and can fail later.

## 4. The labelled owner carrier

Encode the five unit masks by a labelled tuple

```text
(Phi(f_1),...,Phi(f_5)) in R^5.                       (13)
```

At each sheet the owner count `sum_i (f_i)_k` lies in `{0,...,5}`, so
reduction modulo thirteen does not wrap it. From (13), or equivalently its
five Hasse jets, one can recover:

1. every labelled mask--sheet incidence;
2. every sheet's owner count;
3. the deficiency set inside a separately supplied rooted guard-safe block;
4. the result of every **named** translation, deletion, or insertion via
   (6) and (8).

The label tuple is necessary. THM-2197's hostile pair

```text
({x},{x,y})             versus             ({x},{y}) (14)
```

has the same covered union, hence the same aggregate deficiency polynomial.
Deleting the first labelled mask gives different unions. By THM-840's
operation-kernel criterion, no aggregate-union or deficiency quotient can be
deletion-Markov.

Nor does the algebra decide which endpoint move occurs next. The complete
transition address remains

```text
five labelled Hasse jets
+ rooted guard block
+ active bits of the two shallower deep combs
+ cyclic event order and signed coefficient winding.        (15)
```

Given the last three sidecars, (6)--(8) update the incidence component
exactly. Without them, (13) is a faithful snapshot but not an autonomous
dynamical state.

## 5. Exact computational consequence and boundary

At any fixed finite depth, once the named event maps on (15) are derived and
the winding is replaced by its sufficient finite residue/event address, they
act on a finite state space. The coarsest future-coverage congruence can then
be computed by Moore refinement, as in THM-853. The Hasse basis makes deck
motion triangular and singles out the sole first augmentation failure
`epsilon^12`; this is a more structured carrier than a raw list of 512
deficiency subsets.

The theorem does not prove that the private deepest piece in THM-2198 keeps
occupancy below thirteen after an owner expires, derive the event word from
the coefficients, or eliminate any deeper valuation profile. Its cheapest
decisive use is now precise: compute the first exhausted-owner transfer in
the labelled jets and test whether the `epsilon^12` full-fibre coordinate can
occur. If it cannot, THM-2198's no-cancellation support pump iterates. QED.

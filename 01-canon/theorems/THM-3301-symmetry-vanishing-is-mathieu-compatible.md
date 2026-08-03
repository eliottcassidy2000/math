---
id: THM-3301
title: "Symmetry vanishing is Mathieu-compatible"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.
  If a compact group `G` preserves the moment functional and `P` is an
  eigenvector for a character `chi` of infinite order, then `L(P^m)=0` for
  EVERY `m>=1` -- and yet for every fixed `Q` one also has `L(Q P^m)=0` for all
  but finitely many `m`, the exceptional set being exactly the `m` with
  `chi^(-m)` in the isotypic support of `Q`.  So the Mathieu/vanishing
  conclusion HOLDS.  **No counterexample to a vanishing conjecture of
  GMC/GVC/FC/Image type can be produced by a symmetry argument alone.**  A
  refutation must instead use a mechanism that is uniform in `m` on `P` while
  failing to be uniform on `QP^m`; THM-3290's degree-versus-flatness collision
  is such a mechanism, and its `P` is verifiably NOT an eigenvector -- it
  carries six distinct torus weights.  For a finite group of exponent `e` the
  argument degrades further: the surviving set is the progression `e Z`.
audit: >
  The exact companion verifies on the `rho=t^2+xy` Gaussian functional that a
  pure-weight eigenvector kills every moment (24 rows), that for six test
  polynomials `Q` and three weights the nonvanishing exponents are EXACTLY the
  predicted weight-cancelling ones -- at most one per pair, so the Mathieu
  conclusion holds -- with a hostile control showing a weight-zero element need
  not have vanishing moments at all; that THM-3290's `P` has torus weight set
  `{-2,0,2,4,6,8}` and so is not an eigenvector, while violating the Mathieu
  conclusion at `m=1,2,3,4`; and that a finite order `e` leaves exactly the
  multiples of `e` alive.  Normal and `-O` replay are byte-identical.
  Independent immutable audit is pending.
source: death-star-fc-archimedes-port-2026-08-03
depends_on: []
related:
  - THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3204-parabolic-continuant-single-gate-and-jacobi-smith-obstruction
script: 04-computation/symmetry_vanishing_mathieu_thm3301.py
output: 05-knowledge/results/symmetry_vanishing_mathieu_thm3301.out
script_sha256: 266c019a7d1f2b7bd115dd33d829d524ac368038e945a743c3f60cbe5ea0b66a
output_sha256: e8bcce93f2287be41732eb84b04a05057d373035b2b7d829040c348d91226f05
hash_basis: LF-normalized bytes
---

# THM-3301 -- symmetry vanishing is Mathieu-compatible

**PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.**

THM-3300 showed that the Archimedes mechanism cannot be transported to the
Factorial lane because the factorial test functions are exactly the
torus-invariants.  The natural next move is to run the mechanism in reverse:
use a *nonzero* weight deliberately, since a nonzero weight kills every moment
for free.  This theorem shows that move can never work -- for any conjecture in
this family -- and that is why the flatness collision of THM-3290 was not a
convenience but a necessity.

## 1. Setting

Let `L` be a linear functional on a polynomial algebra, invariant under a
compact group `G` acting by algebra automorphisms: `L(f o g) = L(f)`.  The
Gaussian moment functional with `G = SO(n)`, and the factorial functional with
`G = U(1)^n` (THM-3300), are the cases of record.  Call a **Mathieu/vanishing
conjecture** the assertion

```text
L(P^m)=0 for all m>=1   =>   for every Q, L(Q P^m)=0 for all m >> 0.   (1)
```

GMC, GVC and the Image Conjecture have this shape; FC replaces the conclusion
by `P=0`, and the theorem below applies to its `(1)`-shaped weakening.

## 2. The theorem

**Theorem 1.**  Let `chi` be a character of `G` and let `P` satisfy
`P o g = chi(g) P` for all `g` in `G`.

(a) If `chi^m != 1` then `L(P^m)=0`.  In particular, if `chi` has infinite
order then `L(P^m)=0` for **every** `m>=1`.

(b) For every fixed polynomial `Q`, decompose `Q = sum_psi Q_psi` into its
finitely many `G`-isotypic components.  Then

```text
L(Q P^m) = 0   unless   psi * chi^m = 1  for some psi occurring in Q.   (2)
```

Since only finitely many `psi` occur and `chi` has infinite order, the set of
exceptional `m` is finite.  **Hence `L(Q P^m)=0` for all `m >> 0`: conclusion
`(1)` holds.**

*Proof.*  `P^m` is a `chi^m`-eigenvector, so `L(P^m) = chi(g)^m L(P^m)` for
every `g`; choosing `g` with `chi(g)^m != 1` gives (a).  For (b), `Q_psi P^m`
is a `psi chi^m`-eigenvector, so the same argument kills it unless
`psi chi^m = 1`.  `Q` is a polynomial, so it lies in a finite-dimensional
`G`-stable space and has finitely many isotypic components; each contributes at
most one exceptional `m` because `chi` has infinite order.  QED

**Corollary (no symmetry refutation).**  A counterexample to `(1)` can never be
a `chi`-eigenvector for a character of infinite order of any compact group
preserving `L`.  Symmetry produces total vanishing of `L(P^m)` and total
vanishing of `L(QP^m)` *together*; a counterexample must separate them.

## 3. The finite-order degradation, and where each lane sits

If `chi` has finite order `e`, (a) only gives `L(P^m)=0` for `e` not dividing
`m`, and the surviving obligation is the progression `e Z`.  This is the exact
statement behind THM-3300 section 3: the simplex has affine automorphism group
`S_n`, which is finite, so the cyclic route on `Delta_(n-1)` never even reaches
hypothesis (a) for all `m`.

The three lanes therefore line up as follows.

```text
lane            group / character        vanishing reach     refutes (1)?
--------------------------------------------------------------------------
GMC by symmetry  U(1) subset SO(n)       all m               NO  (Thm 1b)
HFC by symmetry  S_n on Delta_(n-1)      m not in eZ         NO  (not even 1a)
THM-3290         no group; order clash   all m               YES
```

The companion confirms the middle column exactly: for `P=x^w` and each of six
test `Q`, the nonvanishing exponents are precisely the weight-cancelling ones,
at most one per pair.

## 4. THM-3290 is not an eigenvector, and had to not be

With `rho=t^2+xy` and the circle `x -> e^(i theta) x`, `y -> e^(-i theta) y`,
the torus weight of `x^a y^b t^c` is `a-b`.  Then

```text
weights(A) = {0,2},   weights(C) = {-1,1,3},
weights(P) = {-2,0,2,4,6,8}.                                           (3)
```

So THM-3290's `P` carries six distinct weights and is as far from an
eigenvector as its degree allows, while `L(P^m)=0` and `L(x^2 P^m) != 0` for
every tested `m`.  By Theorem 1 this had to be so: the counterexample cannot be
weight-homogeneous.  Its mechanism -- a prefactor of degree `k-nu` annihilated
by a derivative of strictly higher order `k-delta`, the `t`-antiderivative being
flat -- is uniform in `m` on `P` precisely because the *gap* `nu-delta` is
uniform, and it fails on `QP^m` precisely because `Q` shifts `delta`.  That
asymmetry is unavailable to a character, whose action on `Q P^m` is the same
multiplication as on `P^m`.

## 5. Reading

The useful form is a search rule.  **When hunting a counterexample to a
vanishing conjecture, first compute the isotypic decomposition of the
candidate.**  If it is concentrated in one isotypic component, Theorem 1 says
the candidate satisfies the conjecture and the search is wasted; the candidate
must mix components, and the mechanism must be an order or degree collision
rather than a character.  Conversely, symmetry remains the cheapest way to
produce elements with *all* moments zero -- it just cannot produce
counterexamples.

This also explains a pattern visible across the repository's refutation
lane: THM-3204's parabolic/non-parabolic dichotomy is the same phenomenon in
multiplicative dress.  A semisimple (diagonalizable) transfer has an eigenvalue
ratio of finite multiplicative order, giving periodic vanishing; the parabolic
(unipotent) transfer has no such ratio and gives the single congruence.
Torsion in the acting character is what turns total annihilation into periodic
gaps, in both settings.

## 6. Scope

Proved: Theorem 1 and its corollary, for any compact `G` preserving `L`.
Verified-exact: `(3)` and the isotypic exceptional-set prediction on the stated
bank.

**Not proved:** any instance of `(1)`, any refutation of it, `FC(n)`, `HFC(n)`,
`GMC(n)` for any `n`, or `JC`.  Theorem 1 constrains *mechanisms*, not truth
values.  It does not say counterexamples are rare or hard to find, only that
they are never eigenvectors.  Non-compact groups, and functionals with no
invariance group at all, are outside its scope.

Run

```text
python 04-computation/symmetry_vanishing_mathieu_thm3301.py
python -O 04-computation/symmetry_vanishing_mathieu_thm3301.py
```

and compare LF-normalized bytes with the declared output.  Exact rational
arithmetic only; no floating point, random sampling, imported executable, or
assertion-sensitive test.

**QED.**

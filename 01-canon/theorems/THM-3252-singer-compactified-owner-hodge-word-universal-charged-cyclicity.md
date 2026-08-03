---
id: THM-3252
title: "Singer-compactified owner Hodge word universal charged cyclicity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Adjoin to THM-3246's exact 168-owner second-corrector word the unique
  completion value making its THM-3234 Singer-plane sum zero.  In every one
  of the 8,064 cyclic Singer relabellings, the resulting 13 by 13
  rational coefficient matrix is nonsingular.  It is therefore a cyclic
  pointed current in every nonzero central block of THM-3250.  A rational
  central contrast combines the twelve blocks into one signed packet whose
  orbit spans the full 2,028-dimensional charged subspace.  This is a linear
  frame for a signed asymptotic corrector, not a positive physical current.
source: root/2026-08-03
audit: >
  The assertion-independent exact companion pins THM-3234, promoted
  THM-3246 and promoted THM-3250; independently reconstructs the complete
  rational owner word from the Bernoulli endpoint formulas; builds the
  deterministic F_169 Singer orbit; and checks all 8,064 gauges by exact
  integer determinants modulo 1,000,000,007.  Both the zero and zero-sum
  completions are nonsingular in every gauge; exact Bareiss determinants at
  the identity gauge and the rational central-contrast Fourier transform
  provide separate controls.  Independent implementations reconstructed the
  word from a common-breakpoint interval integral and an explicit Singer
  matrix, checked all 16,128 determinants with separate elimination code and
  additional primes, matched both exact identity determinants, and verified
  the central contrast and twelve rank-169 blocks directly over a second
  finite field.  Fresh normal and optimized runs byte-match the stored
  transcript.  The audits also repaired stale dependency pins and the
  coefficient-labelling, characteristic, external-gauge, digest-order, and
  physical-typing boundaries below.
depends_on:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word
  - THM-3250-charged-heisenberg-blowup-address-intertwiner-and-pointed-multiplicity-gate
related:
  - THM-3224-complete-lrc-orbit-bernoulli-gcd-carry-and-owner-hodge-splitting
  - THM-3240-exact-address-heisenberg-clutch-on-carrier-imbalance
script: 04-computation/lrc_owner_hodge_charged_cyclicity_thm3252.py
output: 05-knowledge/results/lrc_owner_hodge_charged_cyclicity_thm3252.out
script_sha256: 4471da35fa0fd63dc5c920ef7b695936be722ccbfbb12aeb8cbd55e2607d15c9
output_sha256: 6576214230219b9759646d50f88636ec7a35eab459b19ea722813302744e9d99
hash_basis: LF-normalized bytes
---

# THM-3252 -- Singer-compactified owner Hodge word universal charged cyclicity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3246 promotes an exact signed word on the 168 cyclic owners.  THM-3234
identifies those owners, after one-point completion, with the punctured
affine plane over `F_13`.  THM-3250 then asks a concrete question of any
pointed charged current on that plane: is its `13` by `13` multiplicity
matrix invertible?

For the owner word the answer is uniformly yes.  It stays yes under every
cyclic Singer relabelling, so the determinant gate is not the missing
obstruction.
The remaining obstruction is physical: the word and the central contrast
used below are signed.

## 1. The zero-sum owner compactification

Let

```text
q=(q_0,...,q_167)                                      (1)
```

be THM-3246's all-dilation second-corrector word.  Thus

```text
sum_j q_j=1/24696,
#{q_j>0}/#{q_j<0}/#{q_j=0}=156/12/0.                  (2)
```

Adjoin the Hodge completion

```text
q_*=-sum_j q_j=-1/24696.                               (3)
```

The completed 169-point word has total sum zero.  This is a normalization
of the abstract completion point in THM-3234; it is not a claim that the
point is a physical asymptotic head.

Use the deterministic model

```text
F_169=F_13[u]/(u^2-2),          alpha=1+2u,             (4)
```

where `alpha` has order 168.  For

```text
a in (Z/168Z)^*,       b in Z/168Z,                    (5)
```

define a rational function on `F_13^2` by

```text
A^(a,b)(0,0)=q_*,
A^(a,b)(alpha^(b+aj))=q_j              (0<=j<168).     (6)
```

Regard `(6)` as a `13` by `13` matrix, with the first affine coordinate as
row and the second as column.  This uses a chosen labelling of those
coordinates by THM-3250's coefficient indices `(s,t)`.  It does not identify
the THM-3234 affine `H_13`-set with THM-3250's address carrier as an
`H_13`-set.

## 2. Universal Singer-gauge determinant

For every pair `(a,b)` in `(5)`,

```text
det A^(a,b) != 0.                                      (7)
```

Here `(a,b)` externally relabels the coefficient plane by a primitive Singer
ordering and phase.  It is not a Singer action on `G_delta`, and no
Singer-equivariance of the THM-3250 source or its nonvertical target stratum
is asserted.

This is an exact finite theorem.  Clear the common denominator by

```text
L=32006016000,
B^(a,b)=L A^(a,b).                                     (8)
```

The 168 punctured entries of `B` sum to `1296000`, while its completion
entry is `-1296000`; the complete integer word is primitive.  The companion
performs Gaussian elimination modulo the prime

```text
ell=1000000007                                         (9)
```

on all

```text
phi(168)*168=48*168=8064                               (10)
```

matrices.  Every determinant is nonzero modulo `ell`.  Hence every integer
determinant is nonzero, proving `(7)` over the rationals.  The residues form
168 distinct values, each repeated 48 times.  Their ordered digest is

```text
c4a504a86966c5dc7c0d375e962ae2419fe99271c79b8ae7ddc2f25a5d82d730. (11)
```

The digest order takes increasing integer unit representatives `a` in the
outer loop, `b=0,...,167` in the inner loop, and matrix rows and columns in
the order `0,...,12`.  The same convention is used in `(13)`.

At the identity gauge the exact fraction-free determinant of the cleared
matrix is

```text
-41129257721095723275557032050508281832453443425542278710444911947. (12)
```

As a robustness boundary, replacing `(3)` by the zero completion also makes
all 8,064 matrices nonsingular.  Its ordered modular digest is

```text
e8391c1b3f0a4584f6df15953aa89ed49a8b7dedaf8707ccbfa3689b8fa229e2. (13)
```

Thus zero-sum normalization is natural but is not responsible for
cyclicity.

## 3. Cyclic pointed currents in all charged blocks

Let `K` be a characteristic-zero splitting field containing a primitive
thirteenth root `zeta`, with the conventions of THM-3250.  Fix a nonzero
central character `kappa in F_13^*` and use that theorem's basis

```text
E_(s,t)^kappa
 =sum_(delta in F_13) zeta^(-kappa delta)[s,t,delta].   (14)
```

For any gauge `(a,b)`, put

```text
v_kappa^(a,b)=sum_(s,t) A^(a,b)(s,t) E_(s,t)^kappa.    (15)
```

THM-3250 proves that the orbit span of `(15)` has dimension

```text
13 rank(A^(a,b)).                                      (16)
```

Equations `(7)` and `(16)` therefore give

```text
span_K(H_13.v_kappa^(a,b))=K[G_delta]_kappa,
dim=169                                                 (17)
```

for every nonzero `kappa` and every cyclic Singer relabelling.  The explicit
charged intertwiner of THM-3250 transports `(15)` to a cyclic vector in the
regular nonvertical blowup-flag block as well.

This conclusion is stronger than full support or nonzero energy.  It clears
the exact multiplicity-frame determinant for this signed owner word.

## 4. One rational packet spanning all charged sectors

The twelve currents in `(15)` are the central Fourier blocks of one rational
signed packet.  Define

```text
c(delta)=12,   delta=0;
c(delta)=-1,   delta!=0,                               (18)

W^(a,b)=sum_(s,t,delta) A^(a,b)(s,t)c(delta)[s,t,delta]. (19)
```

Its raw central Fourier transform is

```text
sum_delta c(delta)=0                    at kappa=0,
sum_delta c(delta)zeta^(kappa delta)=13 at kappa!=0.   (20)
```

After scalar extension from `Q` to `K` and the inverse-Fourier factor `1/13`,
the exact Fourier identity is

```text
W^(a,b)=sum_(kappa in F_13^*) v_kappa^(a,b).           (20a)
```

Central idempotents separate the inequivalent blocks, so `(17)` yields

```text
dim span_K(H_13.W^(a,b))=12*169=2028.                  (21)
```

Thus `(19)` saturates the entire charged subspace and sidesteps the neutral
sector in which THM-3250 isolates an exact 156-dimensional rank defect.  It
does so by centering and signing, not by constructing a full positive carrier.

## 5. Scope and the surviving physical gate

The uniformity in `(7)` means that changing the primitive Singer generator
or phase does not destroy cyclicity.  It does not make the matrix itself
canonical.  THM-3234 shows that Singer and Heisenberg together generate the
full affine group and that the completion point is mobile; `(6)` still
requires a chosen affine basis and cyclic identification.

The determinant and orbit-span conclusions are characteristic-zero
statements over `Q` and its splitting extension `K`.  The one modular prime in
the companion certifies that the integer determinants are nonzero; it does not
assert nonsingularity after reduction in every positive characteristic.

More importantly, the `q_j` are signed second asymptotic correctors, not
positive cell masses or endpoint currents.  The completion `(3)` is negative,
and the central contrast `(18)` has both signs.  Neither `(15)` nor `(19)` is
constructed or typed here as a Boolean packet, positive measure, Markov
kernel, lawful LRC current, or Singer-equivariant physical clutch.  No
owner-to-endpoint map, coefficient-current covariance, row exclusion,
arbitrary-radial NC2 theorem, or `LRC(14)` decrement follows.

The exact gain is a gate removal: after passing to the linear charged model,
the THM-3246 owner Hodge word already supplies a full multiplicity frame in
every cyclic Singer relabelling.  Positivity and physical typing, not cyclic
rank, are the remaining debts for this bridge.

## 6. Exact companion

Run

```text
python3 04-computation/lrc_owner_hodge_charged_cyclicity_thm3252.py
python3 -O 04-computation/lrc_owner_hodge_charged_cyclicity_thm3252.py
```

and compare LF-normalized bytes with the declared output.  The companion
uses exact rational, integer, finite-field and modular arithmetic only.  It
contains no randomness, floating point, discovery cache or optimization-
sensitive assertions.

QED.

---
id: THM-2022
title: "Frobenius amplification of the lowest balanced Wick face proves NC2 and GMC(2)"
status: >
  PROVED. For every finite exact support in C[Z,W], a complex torus point of
  the Gaussian moment nullcone descends to an algebraic torus point. The
  lowest balanced face supplies a nonzero Laurent constant term Q by the
  project-internal effective compound-root theorem THM-2095, whose proof is
  Galois-free beyond the small-root logarithmic identity. At a suitable good
  prime p, divide the moment of order p*m0 by the common factorial (p*A0)!.
  Kummer and Lucas identify the surviving residue layer with the p-fold
  dilation of the face channels, and Frobenius makes its residue exactly
  Q^p. Thus no support whose charge convex hull contains zero can be null;
  the nullcone consists exactly of the two strict one-sided charge loci.
  The earlier exposed-two-vertex gap>1 candidate reserved under this ID is
  subsumed and no longer needed.
source: codex-2026-07-21-NC2-followup
supersedes_reservation: "exposed two-vertex factorial face with gap greater than one"
depends_on:
  - THM-2095-effective-compound-root-bound-for-one-variable-constant-terms
  - THM-1540-gmc2-reduced-to-the-nullcone-structure-theorem
related:
  - THM-2067-galois-orbit-product-closes-one-variable-dvdk
  - THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2
  - THM-1645-gmc2-angular-layer-is-dvdk-the-gap-is-purely-radial
  - THM-2019-gmc2-affine-height-supports
  - THM-2020-gmc2-finite-place-channel-separation
  - THM-2033-the-nc2-wall-is-the-confluent-transitivity-vandermonde
  - THM-2040-the-de-factorialization-principle
  - THM-2041-frobenius-stability-of-exact-period-projectors
  - HYP-8800-lrc14-face-carry-frobenius-transfer
  - HYP-8765-gmc2-radial-channel-return-tower
scripts:
  - 04-computation/gmc2_frobenius_lowest_face_codex_20260721.py
  - 04-computation/gamma_radial_frobenius_face_codex_20260721.py
outputs:
  - 05-knowledge/results/gmc2_frobenius_lowest_face_codex_20260721.out
  - 05-knowledge/results/gamma_radial_frobenius_face_codex_20260721.out
formalization:
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2Reduction.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2WickChannels.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2NormalizedMoment.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2LowestFacePackage.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2FaceHeightFloor.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2FaceSeed.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2FaceSeedDescent.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2IntegralFaceSeedDescent.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2ChannelDilation.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2FrobeniusResidue.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2ResidueAssembly.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2NormalizedResidue.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2SupportFaceBridge.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2NC2.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2NC2Capstone.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2HeightWitness.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2Formalization.lean
formalization_root_imported: true
formalization_status: >
  PARTIAL at one explicit analytic interface. The concrete
  normalized three-case residue, exact support-face seed transport, and the
  finite-field zero/nonzero contradiction are kernel-checked. GMC2NC2 derives
  NC2 and GMC(2) from DvdK1. GMC2HeightWitness proves and root-imports
  `heightWitnessSupplier_holds`, removing that former interface. The Lean
  proposition `DvdK1` is now the sole explicit premise; THM-2095 proves its
  mathematical content internally and effectively (THM-2067 is an alternate
  historical route), but the
  small-root/compound proof has not been formalized.
---

# THM-2022 -- Frobenius amplification of the lowest balanced Wick face

Let `Z` be a circular complex Gaussian and put `W=Zbar`, so

```text
E[Z^A W^B] = A! if A=B, and 0 otherwise.
```

For a nonzero polynomial, combine equal monomials and delete zero
coefficients:

```text
P = sum_(i=1)^k c_i Z^(a_i) W^(b_i),       c_i != 0,
q_i = a_i-b_i.
```

The theorem proves the two-dimensional nullcone conjecture NC2:

> **Theorem.** If `E[P^m]=0` for every `m>=1`, then either every `q_i>0`
> or every `q_i<0`. Conversely, either strict one-sided condition makes all
> positive moments vanish. Consequently GMC(2) holds.

The proof does not isolate first-return circuits. It preserves the whole
lowest balanced face, including every balanced-channel collision on that
face, and amplifies
its nonzero constant term by Frobenius.

## 1. Exact Wick channels

For `m>=1`, let

```text
R_m = {r in N^k : |r|=m and sum_i q_i r_i=0},
A(r) = sum_i a_i r_i = sum_i b_i r_i       (r in R_m).
```

Direct multinomial expansion and Wick balance give

```text
M_m(c):=E[P^m]
 = sum_(r in R_m) binom(m;r) A(r)! c^r,                 (1)
binom(m;r)=m!/prod_i r_i!.
```

The crucial two-dimensional feature is that a balanced channel has one
nonnegative radial exponent `A(r)` and hence one monotone Wick factorial
`A(r)!`.

## 2. Reduction to algebraic coefficients on a fixed exact support

Fix the monomial support and regard the moments in (1) as polynomials over
`Q` in the variables `c_i`. Put

```text
I=(M_1,M_2,...) subset Q[c_1,...,c_k].
```

If the complex null locus of `I` meets the coefficient torus, then the
localized algebra

```text
R=Q[c_1,...,c_k,(c_1...c_k)^(-1)]/I
```

is nonzero. (Noetherianity first replaces the infinite moment list by a
finite generating sublist.) A maximal-ideal residue field of `R` is a field
finitely generated as a `Q`-algebra, hence a finite extension of `Q` by
Zariski's lemma. Therefore a complex torus null point would imply a torus
null point

```text
c=(c_1,...,c_k) in (K^*)^k                              (2)
```

for some number field `K`. It is enough to exclude (2).

There is a shorter equivalent descent, now used by the Lean development.
Starting directly from a complex torus point `z`, take the subalgebra

```text
A=Q[z_1,...,z_k,z_1^(-1),...,z_k^(-1)] subset C.
```

It is nonzero, finite type over `Q`, all `z_i` are units, and every rational
moment relation remains zero in `A`. A maximal quotient of `A` is finite over
`Q`; the quotient map preserves the relations and cannot kill a unit. This
avoids both the infinite moment ideal and its Noetherian finite-sublist step.
It produces exactly the point (2).

### 2.1. Shorter integral specialization: let the quotient choose the prime

For the formal proof there is a still shorter route.  Once the lowest-face
constant term `Q(c)` is known to be nonzero, form the finite-type integral
subalgebra

```text
A_Z=Z[c_1,...,c_k,c_1^(-1),...,c_k^(-1),Q(c)^(-1)] subset C.
```

Every maximal quotient of `A_Z` is a finite field of some prime
characteristic `p`: finite type over the Jacobson ring `Z` makes the field
finite as a `Z`-module, and integrality rules out characteristic zero.  The
coordinates and `Q` remain nonzero because they are units.  More strongly,
the quotient map preserves **every** integral polynomial relation that
vanishes at `c`.

The order of quantifiers matters.  The normalized moment polynomial used
below depends on the residue characteristic `p`, so it should be selected
only after the quotient reveals `p`.  Universal preservation of zero
relations permits exactly this.  Since the non-dilated multinomial-isolation
lemma is valid for every prime (it does not require `p>m0`), an arbitrary
maximal quotient works.  Thus the number-field and good-prime construction
remains a valid independent route, but it is no longer required on the
formalization's critical path.

## 3. The lowest balanced face retains a nonzero toral address

Assume the support is not strictly one-sided. Equivalently,

```text
0 in conv{q_1,...,q_k}.
```

Consider the rational linear program

```text
delta = min {sum_i a_i x_i : x_i>=0, sum_i x_i=1,
                                  sum_i q_i x_i=0}.          (3)
```

Linear-programming duality supplies `lambda in Q` such that

```text
a_i-lambda q_i >= delta                                    (4)
```

for every `i`. Let `F` be the equality set in (4). An optimal
point in (3) is supported on `F`, so

```text
0 in conv{q_i:i in F}.                                     (5)
```

Define the face Laurent polynomial

```text
f_F(u)=sum_(i in F) c_i u^(q_i).                            (6)
```

There is no hidden projection collision in (6). If two points on `F` had
the same charge, (4) at equality would give the same `a`; then
`b=a-q` would also agree, contradicting distinctness of the exact monomial
support.

By (5) and the effective compound-root theorem
(`THM-2095-effective-compound-root-bound-for-one-variable-constant-terms.md`),
the constant terms of all positive powers of `f_F` cannot vanish. Choose
`m0>=1` such that

```text
Q:=CT_u(f_F^m0) != 0.                                      (7)
```

If the exact negative and positive charge widths of `f_F` are `M,N`, the same
theorem gives the effective choice

```text
m0 <= binom(M+N,a),
a=min(M,N).                                                (7a)
```

This bound is generally non-sharp and does not make the later good-prime
choice coefficient-uniform.  Beyond its small-root logarithmic identity,
THM-2095 is Galois-free.  For bare existence one may instead cite THM-2067's
Galois orbit-product proof or the stronger published DvdK theorem, THM-1630;
neither is needed in this effective route.

Every multiplicity vector contributing to (7) is balanced, has length
`m0`, and is supported on `F`. For all of them, (4) gives the same integer

```text
A0=A(s)=delta*m0 in N.                                     (8)
```

This also covers endpoint cases. If the optimal face meets charge zero only
at a neutral monomial, then its nonzero coefficient already gives a nonzero
constant term (and one may take `m0=1`).

## 4. Kummer isolates the dilated face layer

Choose a rational prime `p` satisfying

```text
p>m0,                                                       (9)
```

outside the finite set at which one of the algebraic numbers `c_i` or `Q`
fails to be a local unit. (One may also exclude ramified rational primes,
although ramification is not used.) Let `O_K` be the ring of integers of `K`
and fix a prime ideal `pfrak` of `O_K` above `p`. Put

```text
M=p*m0.
```

For every `r in R_M`, (3)-(4) imply

```text
A(r)>=delta*M=p*A0,                                        (10)
```

with equality exactly when `r` is supported on `F`.

Kummer's carry formula for a multinomial coefficient says

```text
v_p(binom(M;r)) = total base-p carry count (with multiplicity). (11)
```

Because `p>m0`, the base-`p` digits of `M` are `(0,m0)`. Thus
there are no carries in (11) if and only if every units digit of every
`r_i` is zero, equivalently

```text
r=p*s,       |s|=m0.                                       (12)
```

Balance of `r` is equivalent to balance of `s`.

Work termwise in the localization `O_(K,pfrak)`. By (10), every channel
factorial in (1) is divisible as a rational integer by `(p*A0)!`, so this
common integer factor can be cancelled *before* reduction. We do not invert
its residue—it is generally zero modulo `pfrak`. The resulting normalized
sum is `pfrak`-integral, and its terms fall into exactly three cases after
reduction.

1. If `r` is not divisible componentwise by `p`, (11) makes the multinomial
   coefficient divisible by `p`, so the normalized term vanishes.
2. If `r=p*s` but `s` is not supported on `F`, then `A(s)>A0`.
   Both are integers, so `A(s)>=A0+1`; the factorial quotient
   `(p*A(s))!/(p*A0)!` contains the factor `p*(A0+1)` and vanishes.
3. If `r=p*s` and `s` is a balanced face channel, then
   `A(r)=p*A0`, so the factorial quotient is exactly one.

Therefore the complete residue layer of `M_M/(p*A0)!` is precisely

```text
{p*s : |s|=m0, q dot s=0, supp(s) subset F}.                (13)
```

It may contain many channels. No unique-channel assertion is used.

## 5. Lucas and Frobenius prevent cancellation inside the layer

For a channel `p*s` in (13), the multinomial Lucas congruence gives

```text
binom(p*m0;p*s_1,...,p*s_k)=binom(m0;s_1,...,s_k) mod p.     (14)
```

Writing bars for residues in the characteristic-`p` residue field, (7),
(13), and (14) yield the exact normalized identity

```text
M_(p*m0)(c)/(p*A0)!
 = sum_s binom(m0;s) c^(p*s)                    mod pfrak
 = (sum_s binom(m0;s) c^s)^p                    mod pfrak
 = Qbar^p                                       mod pfrak. (15)
```

The middle equality is Frobenius: coefficients from `F_p` are fixed by the
`p`-th power map. By the choice of `p`, `Qbar` is nonzero. Thus (15) is
nonzero, and so

```text
M_(p*m0)(c) != 0,
```

contradicting the assumed algebraic torus null point. Section 2 then excludes
every complex torus null point on every support whose charge convex hull
contains zero.

It follows that a null polynomial has all charges strictly positive or all
strictly negative. Conversely, charges add under multiplication, so a strict
one-sided polynomial has no charge-zero monomial in any positive power and
all its moments vanish. This proves NC2.

The implication recorded in
`THM-1540-the-two-dimensional-nullcone-conjecture.md` now gives GMC(2): if,
say, every charge of `P` is at
least one, then for a fixed polynomial `Q_0`, every charge of `Q_0 P^m` is
nonzero once `m` exceeds the finite negative-charge range of `Q_0`; hence
`E[Q_0P^m]=0` eventually. The negative-charge case is identical.

## 6. Why the transfer works, and why it is two-dimensional

The mechanism was found by importing two seemingly unrelated repository
ideas: retain a side channel under Frobenius twist, and attach the base-`p`
carry state to a tropical wall. Here the lower balanced face is the exact
wall, `Q` is the retained side channel, and Kummer says every non-dilated
allocation crosses a carry wall. Frobenius preserves the *whole* face sum
`Q`, so colliding primitive circuits never need to be separated.

This does not prove an analogous statement in higher Gaussian dimension.
For one complex variable, balance leaves the single scalar Wick factor
`A(r)!`, monotone in the lower-face height. With several complex variables
the Wick factor is a product of coordinate factorials; lowering one scalar
functional need not increase every coordinate valuation. That is precisely
where the above minimum-layer argument stops.

For Tournament Analysis, take balanced channels as vertices and compare
their divisibility after the normalization `(p*A0)!`, with lexicographic
order only as a tie
path. The selected quotient preserves the minimum residue layer but
forgets its residue sum. Formula (15), not transitivity or arbitrary
tie-breaking, restores that missing coordinate. The challenged assumption
is therefore explicit: noncancellation does not require a dominant channel;
an entire tied face can survive as one Frobenius power.

## 7. Scope and supersession

The proof uses arbitrary finite monomial support and arbitrary complex
coefficients. It includes arbitrary radial polynomials, arbitrary many
charges, arbitrary many primitive return circuits, neutral terms, and all
degree-resonance bands. It therefore closes the residual left by THM-2017,
THM-2018, THM-2019, HYP-8765, and HYP-8766, without invalidating their finer
effective and asymptotic statements.

The former THM-2022 candidate required an exposed two-vertex face with
factorial gap greater than one. That archimedean estimate is unnecessary:
the lowest balanced face always exists, and the carry/Frobenius argument
handles a many-vertex face with no metric gap beyond strict face separation.

## 8. Exact extension to rational Gamma radial laws

The proof uses a prime-block property of the scalar radial weight, not the
identity `A!=Gamma(A+1)` itself. This gives a genuine new family.

Let `alpha=h/k>0` be rational in lowest terms. Let `T` have the Gamma law of
shape `alpha` and unit scale, let `U` be uniform on the unit circle and
independent of `T`, and put `Z=sqrt(T)U`, `W=Zbar`. Then

```text
E_alpha[Z^A W^B] = 0                 if A!=B,
E_alpha[Z^A W^A] = (alpha)_A         if A=B,             (16)
```

where `(alpha)_A=alpha(alpha+1)...(alpha+A-1)`.

> **Gamma-radial corollary.** For every positive rational `alpha`, the
> moment nullcone of `E_alpha` on `C[Z,W]` is exactly the union of the two
> strict one-sided charge loci. Consequently its polynomial Mathieu
> implication holds: if `E_alpha[P^m]=0` for every positive `m`, then
> `E_alpha[Q_0P^m]=0` for every fixed `Q_0` and all sufficiently large `m`.

The algebraic descent and constant-term face seed in Sections 2--3 are unchanged,
because all moment polynomials have rational coefficients. For the finite
place, choose `p>m0` with `p` not dividing `k` and outside the same finite bad
set. Normalize the moment of order `pm0` by `(alpha)_(pA0)`. If `n>=pA0`,

```text
(alpha)_n/(alpha)_(pA0)
 = product_(j=pA0)^(n-1) (h+kj)/k                         (17)
```

is `p`-integral. More importantly, for every integer `A'>A0`, the interval
`pA0 <= j < pA'` consists of `A'-A0` complete blocks of length `p`. Since
`k` is a unit modulo `p`, each block contains exactly one solution of
`h+kj=0 mod p`. Hence

```text
(alpha)_(pA')/(alpha)_(pA0) = 0 mod p.                   (18)
```

Kummer kills the non-`p`-divisible channels, (18) kills the dilated off-face
channels, and face channels have quotient one. Lucas and Frobenius therefore
give the same residue `Q^p` as (15). This proves the corollary.

More generally, the same proof works for a scalar weight `w(A)` valued in a
fixed number field (or a fixed rational model), provided the moment
polynomials admit the same algebraic descent and face seed, `w(pA0)` is
nonzero, and infinitely many finite places are simultaneously good for the
coefficients, seed, and ratios. At each such place require
`w(n)/w(pA0)` to be integral on the chosen side of the face and
`w(pA')/w(pA0)` to be divisible by the residue characteristic at every strict
dilated off-face grade. An upper-face version reverses the inequalities.

## 9. Sharp multi-factor warning

One scalar exposed face does not control a product of independent factorial
coordinates. The obstruction is not merely a missing proof. Consider two
neutral channel atoms with Wick vectors

```text
v1=(1,1,1),        v2=(0,0,3),                            (19)
```

which have the same total scalar grade `3`. At moment order `p`, the pure
`v1^p` channel has factorial valuation `3`, while the channel using `p-1`
copies of `v1` and one copy of `v2` has vector

```text
(p-1,p-1,p+2)
```

and total `p`-valuation

```text
v_p(binomial(p,1)) + v_p((p-1)!^2 (p+2)!) = 1+1=2.       (20)
```

Thus it undercuts the proposed scalar-face Frobenius channel for every large
`p`. With two factorial coordinates the analogous vectors `(1,1)` and
`(0,2)` already tie. This scalar-face method would therefore need additional
structure in higher dimension, for example a coordinatewise/orthant-exposed
vector face or a rank-one lock forcing all factorial coordinates to be
functions of one grade. A scalar Newton face alone is not sufficient.

THM-2095 supplies a project-internal effective proof of the constant-term
input: the logarithmic small-root identity converts initial moment vanishing
to contact with `ct`, while complementary-subset duality produces a nonzero
degree-`C` compound polynomial that bounds that contact.  Beyond the small-root
identity this proof is Galois-free.  THM-2067 is an alternate project-internal
Galois route to bare existence. A stronger alternate source remains J.J.
Duistermaat and W. van der Kallen, "Constant terms in powers of a Laurent polynomial,"
*Indagationes Mathematicae* (N.S.) 9(2) (1998), 221--231, Theorem 2 and
Remark 3.

## 10. Lean formalization ledger

The kernel-checked development is gathered by
`TournamentH7.GMC2Formalization`. The aggregator now includes a conditional
`NC2`/GMC(2) endpoint with its exact remaining internal interface visible.

1. `GMC2Reduction` and `GMC2ChargeGeometry` prove both strict charge branches,
   identify failure of one-sidedness with charge straddling, define full
   `NC2`, and prove `NC2 -> GMC(2)`.
2. `GMC2WickChannels`, `GMC2MomentRelations`, `GMC2MomentTransport`, and their
   integral copies prove the exact expansion (1) and its universal polynomial
   interpretation. `GMC2NormalizedMoment` proves cancellation of the common
   minimum factorial before reduction.
3. `GMC2LowestFaceExistence`, its exact-face package, the coordinate
   dictionary, and `GMC2FaceHeightFloor` prove rational exposing data, the
   balanced height floor, equality on the face, scaling to mass `p*m0`, and
   the strict integer off-face gap.
4. `GMC2DvdKInterface` states the one-variable input as an explicit
   proposition. Its mathematical content is now proved internally by
   THM-2095 (with THM-2067 as an alternate route), while the Lean proof remains
   to be implemented. `GMC2FaceSeed`,
   `GMC2FaceSeedChannel`, and the reference-channel bridge turn it into a
   nonzero exact face seed and an actual balanced multiplicity vector; no
   custom Lean axiom is declared.
5. The former descent gap is closed twice. `GMC2TorusDescent` and
   `GMC2FaceSeedDescent` preserve all moment relations and the nonzero seed in
   a number field. `GMC2IntegralSpecialization`,
   `GMC2IntegralTorusSpecialization`, and
   `GMC2IntegralFaceSeedDescent` instead specialize directly to a finite field
   of an output prime characteristic while universally preserving every
   integral zero relation.
6. `GMC2ChannelDilation` proves extension-by-zero and componentwise-dilation
   identities for mass, charge, bidegree, height, and coefficient products,
   together with injectivity, the exact antidiagonal image, and finite-sum
   reindexing.
7. `GMC2FrobeniusResidue` proves non-`p`-dilated multinomial isolation,
   multinomial Lucas, the whole-face identity `Qbar^p`, non-cancellation, and
   strict normalized-factorial divisibility. `GMC2ResidueAssembly` proves the
   abstract three-case sum theorem, while `GMC2NormalizedResidue` instantiates
   all three cases for the normalized Wick relation. `GMC2GoodReduction`
   certifies the longer number-field prime-reduction route independently.
8. `GMC2SupportFaceBridge` proves exact seed reindexing into support
   coordinates and packages the global floor, exact face height, and strict
   off-face gap. `GMC2NC2` proves the specialized normalized relation is both
   zero and nonzero. `GMC2HeightWitness` seals the coefficient lookup behind
   an opaque local definition, proves `heightWitnessSupplier_holds`, and
   derives the root-imported endpoints `nc2_of_dvdK1` and `gmc2_of_dvdK1`.

The first residue lemma uses multivariate exponent expansion in
characteristic `p`, not an explicit carry count. It is slightly stronger
than the manuscript narration: it does not need `p>m0`. Retaining `p>m0`
above is harmless and keeps the elementary two-digit Kummer explanation.

One explicit Lean interface remains: `DvdK1`, whose mathematical content is
proved by THM-2067 but whose root-factorization/Galois proof is not yet in
Lean. The former `HeightWitnessSupplier` interface is discharged by
`GMC2NC2.heightWitnessSupplier_holds`: sealing the `P.coeff` lookup avoids the
elaborator's earlier `whnf` explosion without new axioms or a heartbeat
increase.

Formalizing the stronger published DvdK theorem is unnecessary for this
paper proof: THM-2095 gives the required existence statement effectively and
without a Galois endgame. Formalizing its small-root identity plus compound
coefficient argument is still a separate project, so `DvdK1` remains visible
as a Lean theorem hypothesis rather than being hidden behind `axiom`, `sorry`,
or `native_decide`. Thus the paper proof is now internally closed, while the
Lean theorem remains honestly conditional on the one proved-but-unformalized
`DvdK1` proposition.

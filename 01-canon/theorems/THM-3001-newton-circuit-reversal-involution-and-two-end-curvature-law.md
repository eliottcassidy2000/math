---
id: THM-3001
title: "Newton-circuit reversal involution and the two-end curvature law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The exact reversal
  involution, class no-go, two-end expansion, scaled-reciprocal fixed-locus
  characterization, circuit antipalindrome, and parity wall are proved.  The
  former finite-degree curvature sign screen is repaired to its sharp
  O(1/d)-tolerant asymptotic form.  The two-number global classifier remains
  explicitly REFUTED by THM-3004.
source: klein-S428
audit: >
  Independent hostile audit ACCEPT after scope repair.  It rederived the
  reversal formulas and two-end signs, identified the missing O(1/d) tolerance
  in the necessary curvature screen, proved the converse fixed-locus theorem,
  and checked the odd/even chamber wall.  Its exact independent companion
  exhausts 3267 positive coefficient rows, 114 balanced two-cluster controls,
  seven fixed-locus parity controls, three equivariant path midpoints, and a
  constant-R=2 positive-real-rooted non-equality control (five roots counted
  by an exact Sturm chain).  Primary and independent normal,
  optimized, and stored transcripts all byte-match.
depends_on:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
related:
  - THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return
  - THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-2994-first-gap-hurwitz-hermite-biehler-prefix
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
  - THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation
script: 04-computation/gmc_newton_circuit_reversal_involution_and_two_end_curvature_law_thm3001.py
output: 05-knowledge/results/gmc_newton_circuit_reversal_involution_and_two_end_curvature_law_thm3001.out
script_sha256: 430e8a83795ab65bd3ab6df6baed427894c7917ae3939f3581991582338ef7e5
output_sha256: f72746c15b96e4f77df801e42dec908c6c2c4738783700f6212a3f60e94bf00f
independent_script: 04-computation/gmc_newton_reversal_fixed_locus_referee_thm3001.py
independent_output: 05-knowledge/results/gmc_newton_reversal_fixed_locus_referee_thm3001.out
independent_script_sha256: 1eae2e618d0607a7cf33db59ade085d3ec1f71dc7706f20332af956d4ba2b89c
independent_output_sha256: 7ba4dd7017b8dcf69e5f2150a1b8b2100a1f3c89029eec50e3c521f6f5542f5b
hash_basis: LF-normalized bytes
---

# THM-3001 -- Newton-circuit reversal involution and the two-end curvature law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Sections 1--3 are unconditional and exact.  Section 4 inherits THM-3000's
hypotheses.  Section 5 includes the exact fixed-locus strengthening.  Section 6
retains an explicitly refuted hypothesis and its misleading census only under
the displayed retraction.

## 1. The involution

Let `N(n)=sum_(i=0)^d a_i n^i` with every `a_i>0`, and let
`h_k=a_(d-k)/(a_d binom(d,k))`, `R_k=h_k^2/(h_(k-1)h_(k+1))` as in
THM-2997 (1).  Let `N*(n)=n^d N(1/n)`, i.e. `a*_i=a_(d-i)`.

**Theorem (exact, no asymptotics).**

    h*_k=(a_d/a_0) h_(d-k),   0<=k<=d,                  (1)
    R*_k=R_(d-k),             1<=k<=d-1.                (2)

*Proof.* `h*_k=a*_(d-k)/(a*_d binom(d,k))=a_k/(a_0 binom(d,k))` and
`h_(d-k)=a_k/(a_d binom(d,k))`, giving (1).  The constant `a_d/a_0` cancels
from the second-order ratio, giving (2).  QED.

Equivalently, in third-difference form,

    Delta^3(log h*)_j=-Delta^3(log h)_(d-3-j),           (3)

so **coefficient reversal negates every Newton circuit**.  Since `R_k` is also
invariant under root scaling `r->lambda r` (`h_k->lambda^k h_k`), reversal acts
on root multisets as `r -> 1/r` up to scale.

Verified exactly on real-rooted, rational, and non-real-rooted
positive-coefficient controls.

## 2. Class-level no-go for global no-return

Call `N` *ratio-monotone* if `R_k>=R_(k-1)` for every `2<=k<=d-1` (the global
no-return property of THM-2989/2991).

**Theorem.**  Let `H` be any class of positive-coefficient polynomials closed
under reversal.  If every member of `H` is ratio-monotone, then every member of
`H` has a **constant** ratio sequence.

*Proof.* For `N in H` also `N* in H`.  Monotonicity of `N*` reads
`R*_k>=R*_(k-1)`, i.e. by (2) `R_(d-k)>=R_(d-k+1)`, i.e. `R_j>=R_(j+1)` for
`1<=j<=d-2`.  Monotonicity of `N` gives the reverse inequality on the same
range.  Hence `R_j=R_(j+1)` throughout.  QED.

**Corollary.**  None of the following, in any combination, can imply global
no-return, because each is reversal-closed and each contains members with
non-constant ratio sequence:

- positivity of all coefficients (reversal permutes coefficients);
- all roots real and negative (`r>0 => 1/r>0`);
- PF-infinity (Aissen--Schoenberg--Whitney--Edrei: `prod(1+rt)`, `r>0`, and
  `1/r>0`);
- Hurwitz stability (`Re(rho)<0 <=> Re(1/rho)<0`);
- strict ULC (by (1) log-concavity of `h` and of `h*` are the same statement).

Non-constancy witnesses: `(n+a)^m(n+b)^m` for `a!=b,m>=2`, and geometric root
sets of length at least three with ratio different from one; the script prints
exact ratio sequences.  Within the positive
real-rooted/PF-infinity class, equality in every Newton inequality gives
`R_k=1` and forces all roots to coincide.  In the unrestricted
positive-coefficient universe, **and even on some positive-real-rooted
examples**, a constant ratio sequence need not be the Newton-equality sequence:
the independent referee includes an exact `R_k=2` control.  Thus the class
no-go above concludes only that the ratio sequence is constant; it does not
conclude `R_k=1` or coincident roots.

**Relation to THM-2991.**  THM-2991 is *strictly stronger in its own
direction*: it adds an arbitrarily long improving **leading prefix** -- a
one-sided, non-reversal-closed hypothesis -- and produces a return strictly
**below `R_1`**.  Section 2 does not reach that.  What section 2 adds is the
explanation of why THM-2991's PF-infinity / Hurwitz / strict-ULC decorations
are *automatically inert*: they are reversal-closed, so they could never have
helped, and the construction had to be one-sided for a structural reason.
Any future family-specific no-return proof must invoke a **reversal-breaking**
invoice.  THM-2989/2997's wall invoice is one-sided (it constrains only the
top coefficients), so it qualifies -- but section 5 makes the price explicit.

## 3. The reciprocal curvature

For a root measure `mu` (roots of `N`, i.e. `N=a_d prod(n+r_m)`) write
THM-3000's cumulant curvature

    C(mu)=(3kappa_2(mu)^2-2kappa_1(mu)kappa_3(mu))/kappa_1(mu)^4.  (4)

Let `mu*` be the push-forward of `mu` under `r->1/r`; this is the root measure
of `N*`.  Both `C(mu)` and `C(mu*)` are scale-invariant.

## 4. The two-end curvature law

**Theorem** (hypotheses of THM-3000 at both ends: bounded normalized log jets
for `N` and for `N*`).  For each fixed `k>=2`,

    log(R_k/R_(k-1))       = +C(mu)  d^(-2)+O(d^(-3)),   (5)
    log(R_(d-k+1)/R_(d-k)) = -C(mu*) d^(-2)+O(d^(-3)).   (6)

*Proof.* (5) is THM-3000 applied to `N`.  For (6), apply THM-3000 to `N*`,
whose root measure is `mu*`, to get `log(R*_k/R*_(k-1))=C(mu*)d^(-2)+...`, then
use (2) twice: `R*_k=R_(d-k)` and `R*_(k-1)=R_(d-k+1)`, so
`log(R*_k/R*_(k-1))=-log(R_(d-k+1)/R_(d-k))`.  QED.

So **the two ends of the Newton circuit are governed by two explicit cumulant
numbers**.  When the curvatures stay a fixed distance from zero, the bottom and
top signs are eventually `sign C(mu)` and `-sign C(mu*)`.

**PROVED NECESSARY CONDITION (quantitative form).**  For a sequence of degrees
`d->infinity` satisfying the two uniform bounded-jet hypotheses, asymptotic
global no-return requires

    C(mu_d)>=-O(d^-1),   C(mu_d*)<=O(d^-1).              (7)

Consequently `liminf C(mu_d)>=0` and `limsup C(mu_d*)<=0`.  If the two
curvatures are independent of `d`, or converge with a nonzero limiting margin,
this reduces to the formerly displayed sign screen `C(mu)>=0>=C(mu*)`.  The
`O(d^-1)` tolerance is load-bearing when a curvature approaches zero; the
`O(d^-3)` remainder in (5)--(6) does not determine its finite-degree sign.

Verified numerically at `d=120,240,480` on balanced and unbalanced two-cluster,
three-cluster, and geometric families; predicted and measured end values agree
to the expected `O(d^(-3))`.

## 5. Newton-ratio-palindromic families and the reversal fixed locus

**Theorem (exact characterization).**  The following are equivalent:

1. `R_k=R_(d-k)` for every `1<=k<=d-1`;
2. there are `A,B>0` such that `h_(d-k)=A B^k h_k` for every `k`;
3. `N*(x)=A^-1 N(x/B)`;
4. the root-parameter multiset is closed under `r->1/(B r)`.

For positive real root parameters, condition 4 says exactly that the empirical
measure of `log r` is symmetric about `-log(B)/2`.  Thus log-symmetry is not
merely sufficient: it is the complete reversal-fixed locus seen by the Newton
ratios.

*Proof.*  Put `y_k=log h_k` and `z_k=y_(d-k)-y_k`.  Palindromy of `R` says

    2z_k-z_(k-1)-z_(k+1)=0,

so `z_k` is affine and `h_(d-k)=A B^k h_k`.  Since
`h_(d-k)/h_k=a_k/a_(d-k)`, this is equivalent to

    a_k=A B^k a_(d-k),   N*(x)=A^-1 N(x/B).            (8)

Comparing root multisets gives condition 4.  Every implication reverses, and
(2) gives the ratio palindrome.  QED.

The circuit coordinates `c_k=log(R_k/R_(k-1))` therefore satisfy

    c_k=-c_(d+1-k).                                     (8a)

If `d=2m+1` is odd, the central coordinate is fixed by this involution, hence
`c_(m+1)=0` and `R_(m+1)=R_m`.  If `d=2m` is even, the two central coordinates
are opposites.  In particular the symmetric positive-coefficient path

    gamma(t)=(1-t)N+tN*,   gamma(1-t)=gamma(t)*,         (8b)

has a self-reciprocal midpoint; in odd degree that midpoint lies on the exact
central circuit wall.  This is the chamber-holotopy content of reversal: an
equivariant path cannot pass through the fixed locus while retaining a strict
central orientation.

**Instances, all exact in the output.**

- every balanced two-cluster `(n+a)^m(n+b)^m`: `log` roots are `log a`, `log b`
  with equal weights, symmetric about their midpoint, hence `R_k=R_(d-k)`
  exactly.  The five printed controls
  `(a,b,m)=(1,2,4),(1,3,6),(2,5,8),(1,10,10),(3,7,5)` additionally rise
  strictly through the midpoint and then fall; the general theorem here needs
  only exact palindromy and nonconstancy, not an unaudited global shape claim;
- every geometric root set `{q^i}_(i<d)`: `log` roots form an arithmetic
  progression, symmetric.  Exact palindromic `R` printed for
  `q=1/2,2/3,9/10,1/3`.

So the **simplest non-degenerate real-rooted positive family already refutes
global no-return**: its ratio path is palindromic and nonconstant, hence cannot
be monotone.  The listed controls have an exact midpoint maximum.

**Lineage.**  MISTAKE-335 already records "the exact reciprocal symmetry" of
the balanced two-cluster family `(x+1)^n(x+B)^n` and its consequence
`R_(n+1)=R_(n-1)`.  Section 5 is not new for that single family.  What is new
here is (i) the general involution (1)--(3) valid for **every**
positive-coefficient polynomial, of which the two-cluster symmetry is one
instance; (ii) the exact characterization of the palindromic class as
**log-symmetric root profiles**, its parity wall, and its equivariant path; and
(iii) the class-level no-go of section 2, which shows that the failure is
structural rather than a property of one lucky family.  MISTAKE-335's rule --
separate a directional turn from a baseline crossing -- is respected: sections
2 and 5 claim only the directional/monotonicity statement, never a crossing
below `R_1`.  Crossing remains THM-2991's.

## 6. RETRACTED: the two-number sign classifier is FALSE

> **RETRACTION (klein-S428, 2026-07-31).  The hypothesis stated in this section
> is REFUTED by [THM-3004](THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation.md).**
> Witness: `N(n)=(n+1)^2(n+3)^2(n+8)`, degree `5`, all roots real and positive,
> `R_1=256/215, R_2=1849/1600, R_3=10000/8643, R_4=4489/4000` -- down, up, down,
> so the circuit has **two** sign changes while `C(mu)=+0.1414` and
> `C(mu*)=+0.4129` are both positive.  The correct law is a **cluster count**:
> with `m` well-separated root clusters the circuit has up to `2m-3` sign
> changes, attained.  The classifier is true exactly for `m<=2`.
> The `42/42` census below is not evidence: it held cluster **sizes** equal in
> every three-cluster row, which is the axis the failure lives on (MISTAKE-337).
> Sections 1--5 and 7--9 of this file are unaffected; the corrected necessary
> condition (7) still stands, but it is **not** sufficient.

**HYP (REFUTED -- see the retraction above; retained for the correction
lineage).**  For a positive real-rooted `N` with bounded jets at both
ends, the pair `(sign C(mu), sign C(mu*))` determines the global shape:

    C(mu)>0>C(mu*)  =>  R strictly increasing,
    C(mu)<0<C(mu*)  =>  R strictly decreasing,
    C(mu)>0, C(mu*)>0  =>  interior maximum,
    C(mu)<0, C(mu*)<0  =>  interior minimum.                (9)

Only the quantitative necessary half (7) of the first line is proved.  The census in
the output tests `42` families (two-cluster over nine multiplicity splits and
four root ratios, three geometric, three three-cluster) at `d=60` in exact
rational arithmetic: **agreement 42/42**.  Illustrative rows:

    1^(1/6) x 2^rest : C=+0.03518, C*=-0.08538 -> INCREASING   (actual INCREASING)
    1^(1/4) x 5^rest : C=+0.29297, C*=+0.18750 -> INTERIOR-MAX (actual INTERIOR-MAX)
    1^(2/3) x 2^rest : C=-0.01388, C*=+0.05126 -> DECREASING   (actual DECREASING)
    balanced 1,2     : C=+0.03704, C*=+0.03704 -> INTERIOR-MAX (actual, exact palindrome)

The hypothesis is *not* implied by section 4, which controls only `O(1)`-many
circuits at each end.  Closing it needs an interior theorem (a single sign
change of `Delta^3 log h` in `k`), i.e. a **discrete unimodality** statement
that section 4 cannot supply.  A hostile would be a family with
`C(mu)>0>C(mu*)` and a strict interior dip.

## 7. Consequence for the first-gap lane

THM-2997 (25) gives, for the first-gap wall-stripped core,

    C(mu_M) -> 21630685837/71563480803 > 0,

so the **bottom** of the circuit is settled asymptotically.  Since its degree is
linear in `M`, (7) makes the remaining necessary screen

    C(mu*_M)<=O(M^-1),   hence limsup_M C(mu*_M)<=0.     (10)

the same cumulant functional of the **reciprocal** root measure, computable from
the **bottom** coefficients `a_0,a_1,a_2,a_3` of the core -- equivalently from
the first three log jets of `N*_M`.  This is a bounded new obligation of exactly
the shape THM-2989/2997 already discharge at the top, applied to the reversed
core.  A reciprocal curvature bounded below by one fixed positive constant
would **refute eventual** global no-return.  A single positive finite-`M` value,
or a positive value of order `M^-1`, need not: that was the scope defect repaired
by the hostile audit.

## 8. Boundaries and losses

- (1)--(3) and section 5 are exact for every `d`; sections 4 and 6 are
  asymptotic and inherit THM-3000's graded-jet hypothesis at **both** ends.
  Unbounded jets at the bottom break (6) exactly as THM-3000 section 6 breaks
  (5).
- Section 2 bounds hypothesis classes, not individual polynomials.  It does not
  say a particular reversal-closed class has no ratio-monotone members; it says
  the class cannot *characterize* them.
- The map `N -> N*` destroys nothing (it is an involution on positive-coefficient
  polynomials) but it exchanges "leading" and "trailing" data: a proof invoice
  stated only for `a_d,a_(d-1),...` says nothing about the top of the circuit.
  Restoration sidecar: the reversed core's own wall/jet package.
- Nothing here proves GMC(2), ULC, arbitrary-radial NC2, or removes any
  continuation wall hypothesis; MISTAKE-211 still applies.

## 9. Reproduction

    python3 04-computation/gmc_newton_circuit_reversal_involution_and_two_end_curvature_law_thm3001.py
    python3 -O 04-computation/gmc_newton_circuit_reversal_involution_and_two_end_curvature_law_thm3001.py
    python3 04-computation/gmc_newton_reversal_fixed_locus_referee_thm3001.py --output .scratch/thm3001.referee.normal.out
    python3 -O 04-computation/gmc_newton_reversal_fixed_locus_referee_thm3001.py --output .scratch/thm3001.referee.opt.out

Five checks (I involution, II reversal-closure of the standard classes with
exact PF-3 Toeplitz minors and strict-ULC controls, III two-end law, IV exact
palindromy, V classifier census); all report `True`.

The primary normal/optimized/stored transcript is 7,442 LF bytes; the
independent normal/optimized/stored transcript is 815 LF bytes.  Frozen hashes:

    primary script      430e8a83795ab65bd3ab6df6baed427894c7917ae3939f3581991582338ef7e5
    primary output      f72746c15b96e4f77df801e42dec908c6c2c4738783700f6212a3f60e94bf00f
    independent script  1eae2e618d0607a7cf33db59ade085d3ec1f71dc7706f20332af956d4ba2b89c
    independent output  7ba4dd7017b8dcf69e5f2150a1b8b2100a1f3c89029eec50e3c521f6f5542f5b

The hostile audit treats the `42/42` classifier census only as a preserved
negative lesson.  Its mathematical promotion covers sections 1--5 and the
quantitative screen (7), not the refuted implication (9).

**QED.**

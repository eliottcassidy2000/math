---
id: THM-3001
title: "Newton-circuit reversal involution and the two-end curvature law"
status: PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428
depends_on:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
related:
  - THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return
  - THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-2994-first-gap-hurwitz-hermite-biehler-prefix
script: 04-computation/gmc_newton_circuit_reversal_involution_and_two_end_curvature_law_thm3001.py
output: 05-knowledge/results/gmc_newton_circuit_reversal_involution_and_two_end_curvature_law_thm3001.out
script_sha256: 75411c55dbda4b223a9cbbf4758ed9fdc97e902103d703a2eca64ed79b51b231
output_sha256: 90bd4140a2a2ddc053c0ae0e43904493c5ed26fc43dd88a175ab3019fccd91f2
hash_basis: LF-normalized bytes
---

# THM-3001 -- Newton-circuit reversal involution and the two-end curvature law

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.**

Sections 1--3 are unconditional and exact.  Section 4 inherits THM-3000's
hypotheses.  Section 6 is an explicitly unproved HYPOTHESIS with a census.

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

Non-constancy witnesses: every `(n+a)^m(n+b)^m`, `a!=b`, and every geometric
root set; the script prints exact ratio sequences.  Newton equality forces
`R_k=1` for all `k` only when all roots coincide.

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
numbers**, and the bottom and top signs are `sign C(mu)` and `-sign C(mu*)`.

**PROVED NECESSARY CONDITION.**  Asymptotic global no-return requires

    C(mu)>=0>=C(mu*).                                    (7)

This is a cheap two-scalar screen that replaces a per-edge circuit
computation at both ends.

Verified numerically at `d=120,240,480` on balanced and unbalanced two-cluster,
three-cluster, and geometric families; predicted and measured end values agree
to the expected `O(d^(-3))`.

## 5. Newton-ratio-palindromic families: log-symmetry

**Theorem.**  If the root multiset satisfies `mu*=mu` up to scale -- equivalently
if the empirical measure of `log r` is **symmetric about a point** -- then

    R_k=R_(d-k) for all k,                               (8)

so the ratio sequence is an exact palindrome and cannot be monotone unless it
is constant.  By (5)--(6) its two ends then carry opposite signs `+C(mu)` and
`-C(mu)`.

*Proof.* Scale invariance of `R` plus (2).  QED.

**Instances, all exact in the output.**

- every balanced two-cluster `(n+a)^m(n+b)^m`: `log` roots are `log a`, `log b`
  with equal weights, symmetric about their midpoint.  Checked for
  `(a,b,m)=(1,2,4),(1,3,6),(2,5,8),(1,10,10),(3,7,5)`: `R_k=R_(d-k)` exactly,
  strictly rising for `k<=m` and strictly falling from `k=m+1`, with the turn
  exactly at the midpoint;
- every geometric root set `{q^i}_(i<d)`: `log` roots form an arithmetic
  progression, symmetric.  Exact palindromic `R` printed for
  `q=1/2,2/3,9/10,1/3`.

So the **simplest non-degenerate real-rooted positive family already refutes
global no-return**, with an exact interior maximum at the midpoint, and no
construction is required.

**Lineage.**  MISTAKE-335 already records "the exact reciprocal symmetry" of
the balanced two-cluster family `(x+1)^n(x+B)^n` and its consequence
`R_(n+1)=R_(n-1)`.  Section 5 is not new for that single family.  What is new
here is (i) the general involution (1)--(3) valid for **every**
positive-coefficient polynomial, of which the two-cluster symmetry is one
instance; (ii) the exact characterization of the palindromic class as
**log-symmetric root profiles**, which also covers every geometric root set;
and (iii) the class-level no-go of section 2, which shows that the failure is
structural rather than a property of one lucky family.  MISTAKE-335's rule --
separate a directional turn from a baseline crossing -- is respected: sections
2 and 5 claim only the directional/monotonicity statement, never a crossing
below `R_1`.  Crossing remains THM-2991's.

## 6. HYPOTHESIS: the two-number sign classifier

**HYP (unproved).**  For a positive real-rooted `N` with bounded jets at both
ends, the pair `(sign C(mu), sign C(mu*))` determines the global shape:

    C(mu)>0>C(mu*)  =>  R strictly increasing,
    C(mu)<0<C(mu*)  =>  R strictly decreasing,
    C(mu)>0, C(mu*)>0  =>  interior maximum,
    C(mu)<0, C(mu*)<0  =>  interior minimum.                (9)

Only (7) -- the necessary half of the first line -- is proved.  The census in
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

so the **bottom** of the circuit is settled asymptotically.  By (7), the
remaining necessary condition for that family is the single new scalar

    C(mu*_M) <= 0,                                       (10)

the same cumulant functional of the **reciprocal** root measure, computable from
the **bottom** coefficients `a_0,a_1,a_2,a_3` of the core -- equivalently from
the first three log jets of `N*_M`.  This is a bounded new obligation of exactly
the shape THM-2989/2997 already discharge at the top, applied to the reversed
core.  Note that a positive `C(mu*_M)` would **refute** global no-return for the
family outright, so (10) is a decisive test, not a convenience.

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

Five checks (I involution, II reversal-closure of the standard classes with
exact PF-3 Toeplitz minors and strict-ULC controls, III two-end law, IV exact
palindromy, V classifier census); all report `True`.

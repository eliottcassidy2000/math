---
id: HYP-2104
status: SCOPE-CORRECTED CONNECTION; elementary interval criterion PROVED in THM-398; unrestricted C-prime REFUTED; primitive candidate OPEN; no Vitali-cover theorem
source: opus-2026-06-03-S573
related:
  - THM-398
  - HYP-2103
  - HYP-2102
  - THM-369
---

# HYP-2104: strict-measure and boundary-witness handoff, scope corrected

> **Correction lineage, 2026-10-05.** The former claim that a fixed-radius
> periodic arc family is a genuine fine Vitali cover was false. Its interval
> containment test uses connectedness, not a Vitali covering theorem. The
> unrestricted C′ is also false: n=3, S={3,6} has M(S)=1/3. Only the primitive
> candidate gcd(S)=1 remains OPEN. This note does not promote that candidate
> or the historical sampling percentages into a theorem.

**The precise survivor.** For a fixed finite speed set S, let
M(S)=max_t min_v ||vt|| and G(S)={t: ||vt||>1/n for every v}.
Then M(S)>1/n iff the strict safe set G(S) has positive Lebesgue measure.
The weak safe set at >=1/n may be nonempty even when its measure is zero;
for S={1,2}, n=3, its two points are1/3 and2/3. This distinction concerns
measurable finite unions of intervals and their boundaries, not a
nonmeasurable Vitali set. A null weak safe set may also be empty, so bulk
measure alone does not certify the weak witness required by LRC.

**Gcd normalization comes first.** M(dS)=M(S), since t->dt is surjective
on the circle. For a general S, set S0=S/gcd(S). If no member of S0 is
divisible by n, t=1/n supplies a weak witness. Otherwise the unproved
primitive C′ candidate would give strict slack. This is a sufficient
reduction to a conjecture, not an equivalence with LRC or a proved
measure/construction classification of all speed sets.

**Elementary interval criterion, n>=3.** For fixed v, the danger set
D_v={t: ||vt||<1/n} has arcs of radius1/(nv), separated by positive gaps.
A connected open interval is contained in D_v iff it lies in one arc,
with the endpoint inequalities interpreted weakly. Length at most2/(nv)
is necessary for containment, not sufficient: the location matters.
Thus a component of G(S\{v}) longer than2/(nv) proves strict slack.
If no such component exists, no conclusion follows from this length test.
The exact common-centre endpoint test is LemmaD of the corrected THM-398.

**Why this is not a fine Vitali cover.** For fixed v the listed arcs have
one fixed positive radius, not arbitrarily small radii through each
covered point. Bounded eccentricity alone does not supply that missing
hypothesis. Adding all smaller subintervals would change the family and
would not prove the desired arithmetic non-covering statement. The
Vitali covering theorem and Lebesgue density theorem are not invoked.

**Historical finite evidence.** The old S573 experiment reported
72.4%,78.7%,88.9%,91.5%,96.8% successes for its sampled multiple-of-n
rows at n=6,8,10,12,14. These are sampled criterion outcomes, not proved
fractions of an infinite class. The uncut/un-normalized universal
multiple-of-n claim is refuted above; primitive alignment remains open.

**Connection to the current certificate work.** A point can be invisible
to a nonatomic bulk measure while still requiring a proof. By contrast,
a declared discrete full-support measure assigns positive mass to each
integer, so an exact global zero residual excludes every integer
exception. This is a change of measure and domain, not a transfer of the
Vitali covering lemma. A finite sample or small residual without an
explicit atom bound still does not certify all integer inputs.

See [THM-398 with exact scopes](../../01-canon/theorems/THM-398-lrc-reduction-to-Cprime-and-dominance-dodge.md),
[the selector/certificate audit](../results/vitali_selector_certificate_20261005.md),
and [the atomic-prefix certificate package](../results/atomic_prefix_certificate_20261005.md).
The earlier reflection and S573 output remain historical provenance.

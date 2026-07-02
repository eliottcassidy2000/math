# Two localizations, and what they leave behind

*klein-2026-07-01-S86. A reflection on HYP-3832/3833 — encoding the certificate on the apex complex and
localizing the moment LP, two "localization" moves that share a shape: each makes the object we care about
manifest, and each leaves the same global remainder.*

The session asked two questions that look unrelated — encode the flip-rank certificate as a co-cycle on the
PSL(2,7) complex and test its soundness; localize the covering-min moment LP with the Γ₀(14) congruence and see
if the minimum turns positive — and they came back the same shape. Both are localizations. Both work on the
thing I wanted. Both leave exactly the hard part behind, and in the same way. That symmetry is the keeper.

The cochain computation was the more sobering and the more honest. I built the GF(2) complex on the apex group
and computed its Betti numbers, and the first one is not zero: `b1 = 14` for the small generators, fifty for a
richer set. So the bare complex is not a coboundary expander — there are cocycle certificates the local square
test simply cannot see, and no amount of staring at the raw topology removes them. That is the anti-LTC
obstruction refusing to die on contact with the apex group, and I made myself write it down rather than bury
it. But the same computation carried a real positive: the specific certificate I care about, the Paley/QR sign,
is a coboundary — it lifts to a difference of vertex values, which is the definition of locally testable. And
the links of the complex are complete bipartite graphs, the best possible local expanders, which is exactly the
prerequisite the locally-testable-code construction needs. So the picture is not "it works" and not "it fails";
it is "the object is testable and the substrate is valid, and the residual `b1` cosystoles are what the base-
code layer is for." The apex group gives a good building; it does not, by itself, furnish the room.

The moment localization was cleaner and, I think, genuinely illuminating about why the analytic method stalls.
The unlocalized first-moment bound on the lonely measure is `1 - 2(n-1)M`, and at the covering-min that is
about `-1`: the runners' danger arcs total nearly two, they cover the circle twice over, and the union bound is
worse than useless. That is the sum-of-squares stall in one line. Then the identification that made the session:
localizing with the Γ₀(14) congruence is not a new device — it *is* the mult-by-`n` phase coordinate I have
used since the phase-residue work, because the level `N=14` multiplier is exactly mult-by-`n`. And in that
coordinate the runners are an arithmetic progression of step `n`, and an arithmetic progression of step `n`
leaves a gap of `n` at the observer, so the minimum is `n/Φ6`, positive, above `1/n`. The moment method's stall
was never a failure of the method; it was a failure of coordinates. In the raw circle the certificate is a
cancellation the union bound cannot see; in the congruence coordinate it is the trivial fact that a progression
tiles with its step. `min m0 > 0` under localization, and the reason is that localization is the change of
variables in which the positivity is manifest.

What both moves leave behind is the same object, and naming it is the honest close. The cochain localization
makes the *specific* certificate testable but leaves the *general* cosystoles for the base code. The moment
localization makes the *construction's* covering-min positive but leaves the *minimum over all configurations*
— the actual Lonely Runner — untouched. In both cases the localization handles the object I can name and
exhibit, and defers the object that is a quantifier over everything. That is not a coincidence of these two
computations; it is what localization *is*. You localize by fixing a symmetry or a coordinate, and fixing it
trivializes the instance while the universal statement — over all cocycles, over all configs — sits outside the
fixed frame. The lesson I want to carry is to expect this: a localization that makes your example manifest is
telling you precisely which quantifier you have not yet earned. Here the two unearned quantifiers are the same
kind — "for all cocycle classes" and "for all speed configurations" — and both are the step from a certified
construction to a universal theorem.

The through-line to the two-halves picture from last session is worth marking. The certificate side (the apex
group, `√p`, the cochain) and the measure side (the covering-min, `1/√(πn)`, the moment LP) were supposed to be
complementary, and this session they behaved identically under localization — same partial success, same
residual quantifier. That is a hint that the certificate and the measure are not merely complementary but
governed by one obstruction wearing two costumes: the passage from a fixed instance to a universal statement,
whether the instance is a QR co-cycle or an arithmetic-progression cloud. Whatever eventually closes the Lonely
Runner will, I suspect, be the thing that earns that universal quantifier — a coboundary/cosystolic expansion
theorem on the group side, or a min-over-configs discrepancy bound on the runner side — and the two may turn
out to be the same theorem read in the group and in the circle. For now I have two localizations that work, two
remainders that match, and a clearer name for the wall.

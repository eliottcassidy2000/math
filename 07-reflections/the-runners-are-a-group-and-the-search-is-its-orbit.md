# The runners are a group, and the witness search is its orbit

*mac-mini-2026-07-01-S79. Reflection on HYP-3795 (with opus-S13, kps-S7, klein-S68).*

The owner pointed at a metaphor — pushing Verblunsky coefficients to the unit circle is a recursive picture
of runners on a loop — and asked for a dictionary of functions between points on a loop, enough of them that
they might operate group-like. The dictionary came easily: rotations, the antipode, dilations, the affine
`vx+a`, Blaschke and Cayley and Möbius, the Gauss map and the Farey mediant, the runner flow, the three-gap
return, the Atkin–Lehner involution, the Ramanujan transform, the Dedekind cocycle, the sawtooth and the
distance. Two dozen maps, and they close up: rotation and reflection and dilation generate the affine group
of the circle; Gauss and the mediant generate `PSL(2,Z)`; Blaschke fills in `PSL(2,R)`; and the arithmetic
apex is `Gamma_0(14)` with its finite shadow `PSL(2,7)`, the Klein quartic the project has circled for
sessions. The dictionary is not a list. It is a group, seen through many windows.

But the sentence that reorganized everything was smaller. **The runners are the dilation subgroup.** A
runner of speed `v` is the map `x -> vx` on the loop. So the whole configuration of runners is a tuple of
group elements, and asking for a lonely time is asking whether some *other* group element — a dilation `a`,
a choice of rational clock `a/q` — carries the configuration into the safe region where no runner is near
the origin. The witness search is a group acting on itself. `M(S) = max_q (1/q) max_{a in (Z/q)^*} min_i
||a v_i||_q`: over all resolutions `q`, does the orbit of the residue configuration under the unit group
`(Z/q)^*` meet the safe box `[rq, q-rq]`? The Lonely Runner Conjecture, in the covering case, is the
statement that every covering orbit meets that box.

Once it is an orbit, the analytic engine is forced on you, and it turns out to be the one kps had already
found. Count the good clocks — the units `a` that land the configuration safely — and expand the box's
indicator in characters. The inner sum over units is a Ramanujan sum. So the orbit-count is
`phi(q)(6/7)^{13}` plus a signed sum of Ramanujan sums over the resonances `sum k_i v_i`. The main term is
just the volume of the box times the size of the group: fifteen good clocks, if the runners were
independent. They are not. For the construction the residues line up in an arithmetic progression, the most
resonant configuration there is, and the Ramanujan error `c_{183}(61) = -60` drags fifteen down to two. Two
good clocks: `14/183` and its mirror. And two is exactly the number of atoms in the extremal measure, the
rank of the Toeplitz moment matrix, the number of binding runners. The group count, the moment rank, and the
binders are one integer. Four sessions and four vocabularies — orbits, moments, residue bands, phase — were
all naming the same two points.

The pattern that transcends the theorem: **when a search has a symmetry, it is not a search over a set but
an orbit under a group, and the group hands you both the coordinates and the estimate.** We had been
treating the witness `t` as a free continuous parameter to be hunted with Fourier tails and equidistribution.
It was never free. It is a group element, the dilation clock, and the runners are group elements too; the
loneliness is one orbit meeting one box, and the error in the count is the group's own Fourier — the
Ramanujan sums. The metaphor the owner offered, runners as points on a loop, was literal. The loop is a
group, the runners live in it, and the whole difficulty — the resonant AP that pulls fifteen clocks down to
two — is the single fact that the extremal runners form a subgroup-like arithmetic progression, maximally in
tune with itself. The proof, if it comes, will be a bound on how far out of tune a covering configuration is
forced to be.

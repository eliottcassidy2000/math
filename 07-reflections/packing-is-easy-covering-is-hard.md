# Packing is easy, covering is hard: a symmetry-group asymmetry

*klein-2026-07-01-S72. A reflection on HYP-3804 — a safari of new tournament invariants, and the one that
surprised me: the rainbow number.*

Invited to roam the tournament world, I ran a safari of new invariants on the cube model, and most of them
did what safaris do — one was a clean dead-end (the skew-adjacency spectrum `T-T^T` distinguishes almost
nothing: six distinct spectra among fifty-six classes at `n=6`), one was a pretty confirmation (the
directed triangles of the quadratic-residue tournament form a `2-(q,3,(q+1)/4)` design, and only the
doubly-regular ones do). But one was a genuine surprise, and it is the dual of last session's flip-rank.

Last session I asked for the *covering* minimum: the fewest free arcs so that flipping them reaches every
isomorphism class. Call it `rho(n)`. The natural dual is the *packing* maximum: the largest subcube whose
completions are all in *distinct* classes — the rainbow number `R(n)`. Both are bounded by the information
count `log2 |G_n|`: packing from below (`R <= floor`), covering from above (`rho >= ceil`). The surprise is
that these two bounds behave completely differently. **Packing always achieves its floor** — `R(n) =
floor(log2 |G_n|)` exactly, for every `n` I could check. **Covering does not** — `rho(n)` exceeds its
ceiling starting at `n=6`, and the gap `rho - R` grows: `0, 0, 1, 2`. You can always find a perfect rainbow
subcube; you cannot always cover with the ceiling many bits.

Why should covering be harder than packing here? The two are usually thought of as a matched pair — packing
numbers and covering numbers of a code bracket each other around the sphere-packing bound. But those are
statements about an *unstructured* space. Here the space is a cube with a *symmetry group* acting, and the
objects we are packing or covering are its *orbits*. That changes the game asymmetrically. To pack — to find
a subcube of distinct classes — I have freedom: I only need to *avoid* the collisions the group creates, and
a cube has enough room to route around them. To cover — to hit every class with a subcube — I have no
freedom: I must include a representative of each orbit, and the group's folding forces some orbits to share
a face whether I like it or not. The `S_n` action *identifies* points, and identification is an obstruction
to covering (you can't separate what the group has glued) but an *aid* to packing (fewer distinct things to
fit). So the same folding that lets the classes pack tightly is what stops them from being covered cheaply.

That reframes last session's `n=6` transition. I had said the balanced-cut encoding "fails" at `n=6`
because the two triangle-halves are too symmetric. The rainbow number says the failure is not really about
that one shape — it is intrinsic: at `n=6` the *covering* problem first outgrows its information floor,
while the *packing* problem never does. The extra bit `rho` needs beyond the ceiling is the price the
symmetry group charges for covering, and it first comes due at `n=6` — the same size where so much in this
project turns hard.

The lesson generalizes past tournaments. When you have a group acting on a coded space and you want to
represent every orbit by flipping a few coordinates, expect packing and covering to diverge: the orbits
will pack at the information floor and cover above the ceiling, and the gap measures how tightly the group
folds the space. It is a cheap diagnostic for "how much symmetry is in the way" — run the packing and the
covering, and read the gap. For tournaments the gap opens at `n=6`; for another family it would open
wherever the group's identifications first outrun the cube's freedom. The safari's real catch was not a new
invariant but a new *contrast*: two quantities that agree in a structureless world and split apart exactly
in proportion to the symmetry, and the splitting is the signal.

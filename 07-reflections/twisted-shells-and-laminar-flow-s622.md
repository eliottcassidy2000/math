# Twisted involutions on flow shells — the n=14 dodge (S622)

The prompt was almost entirely metaphor: twisted involutions on flow shells, water and laminar
flow, the cross the lattice grid makes, and Kravitz's normalization. I want to record how literally
each metaphor turned into the mathematics, because the translation *was* the result.

**The flow shells.** All of LRC(14) is the residual C′(14): a config with a multiple of 14 must be
loose. The elementary witness `t=1/14` dies exactly on that runner (it sits on the observer). The
question is which *other* clock to look at. The clocks that beat `1/14` are `t=a/m` with `m ≤ 2n-1
= 27`, and `27 = 3³` — a tower of 3-adic shells. So "flow shells" are not a metaphor at all; they
are the divisor tower of the pair-sum-sieve modulus, the only place a new witness can live.

**The twist.** On a shell `m`, the witness `t=a/m` reads each runner as `a·v mod m`. The multiplier
`a` is a free perspective: it *rotates* every runner around the clock. A config is loose if some
rotation lifts all runners out of the central danger band `(m/14, 13m/14)`. For the shallow shells
`m ≤ 13` the band is just `{0}` and the twist does nothing — that is the old one-clock lemma. The
twist only matters on the deep shells `15 ≤ m ≤ 27`, where the band is three residues wide and you
have to *spin* the configuration to thread it through. The twist is the perspective key, wearing a
working hat.

**The involution and the cross.** On the deepest shell `m = 27`, threading the band means avoiding
`{0, ±1}` for every runner, and that happens iff some unit pair `{u, -u}` is left empty by the
config. The pairing `u ↦ -u` is the twisted involution; `(ℤ/27)*` splits into nine of its orbits.
The runners must block all nine pairs to trap the config. And the place the involution breaks — its
fixed locus — is exactly the multiples of 3, the inner shells, with the apex `0` at the center.
That is the cross the lattice grid makes: lay the residues on the 27-cycle and the `±` symmetry
folds it across the `0`–axis; the fold line passes through the 3-shells. A runner that lands on the
cross cannot be spun off it by any perspective; it has to be dodged on a coprime shell instead.

**Laminar versus turbulent.** Here is where the water picture earned its place. A near-tight config's
residues form a contiguous arc on the clock — a laminar, non-crossing flow. A contiguous arc of
thirteen almost always leaves an antipodal pair untouched (all but the two half-line arcs), so
laminar configs almost always have a free perspective to escape through. A config whose residues are
*spread* across the whole clock — turbulent — blocks every pair, but a spread config is a dominant
one, and dominance is the old measure dodge (Lemma B). So the dichotomy is physical: **laminar flow
escapes by a twist; turbulent flow escapes by dominance.** Nothing falls between. Over forty-six
thousand multiple-of-14 configs, the residual was empty.

**What is actually proved, and what is not.** The dodge itself is a theorem and a one-line proof —
exhibit the witness `t=a/m`, done (THM-412). What is conjecture is the *coverage*: that
twist-on-shells (`m ≤ 27`) together with dominance leaves nothing behind. I believed for an hour
that the laminar half had a clean lemma — "a thirteen-arc always frees a pair" — and the computer
found the two exceptions (`{1..13}`, `{14..26}`, the two halves of the line) before I wrote it into
canon. The honest statement is a union over shells, not a single deep-shell lemma, and the half-line
exception is precisely *why* one shell can never suffice. That correction is the shape of the open
problem: prove the dichotomy shell by shell.

The satisfying part is that the improvement is concrete. THM-398 reduced LRC(14) to C′ and then
left a 3.2% "all-short" residual its measure tool could not reach. The twisted-shell dodge reaches
exactly that residual — the 846 configs where the dominance dodge fails are all spun loose by a
shell twist. The frontier didn't move because someone integrated a harder bound; it moved because
the runners were read in the right perspective, on the right shell, and the laminar ones simply flow
around the obstruction.

# The two halves of the project meet at the atoms

*klein-2026-07-01-S70. A reflection on HYP-3802 — putting a tournament on the roots of -1, and finding the
repo's two research programs standing on the same seven points.*

This repository has always had two mandates that felt only loosely related: the **tournament-parity**
program (Rédei's theorem, the odd-cycle collection formula, the tiling model of orientations of `K_n`) and
the **lonely-runner** program (LRC(14), the covering-min `n/Phi6`, the lonely measure on the circle). The
owner's suggestion — put a tournament on the atoms of the lonely measure, use the dihedral group, add the
seventh point — turned out to be the place where the two programs are literally the same object.

Here is the coincidence, and it is not a coincidence. The tight core `{1..13}` binds at `M = 1/14` on a set
of **six atoms**, and those six atoms are exactly the units `(Z/14)*` — the **primitive 14th roots of
unity**. Six, because `phi(14) = 6`. Add the one odd point that is not a unit, `7/14 = 1/2` (the point the
antipode fixes), and you have **seven** points: the seventh roots of `-1`, spaced for a dihedral `D_7`. On
seven points there is a tournament — the quadratic-residue / Paley tournament — and it is the most
symmetric one, with the maximal odd-cycle count and an automorphism group of order 21. So the lonely
measure's support *is* a tournament's vertex set, and the prime that governs both is the same prime: `7 =
n/2`, the prime of the 7-vanishing (the Fourier coefficients `s(t)` die on multiples of 7), which is also
the number of tournament vertices, which is also the modulus `Z/7` whose quadratic residues orient the
arcs. The tournament program's `Z/7` and the runner program's 7-vanishing are one prime wearing two hats.

And the antipode closes it. In the runner picture, `iota: t -> 1-t` is the complement — the sign involution
of all my recent sign-theory. In the tournament picture, reversing every arc is the **converse** `T^op`. On
the seven roots of `-1`, these are the *same reflection* of the dihedral group, and the Paley tournament is
self-converse under it. So "prime 2" — the involution, the `+-` of the distance to the nearest integer —
is the tournament converse; "prime 7" is the vertex set. The whole edifice sits on the two prime factors of
`n = 2 * 7`, exactly as the S67 prime lens predicted: 2 is the sign, 7 is the arithmetic.

Cyclotomy is the glue, and it paid a concrete dividend. Reading the six-atom measure through orthogonal
polynomials, its moments are **Ramanujan sums** `c_14(k)/6` (independently mac-mini and kind-pasteur reached
the same Ramanujan-sum moments from the flat-extension side — HYP-3793 — so three roads meet here too), and
its Verblunsky coefficients are the clean ladder `1/6, 1/5, 1/4, 1/3, 1/2, 1`. The recursion terminates at
`|alpha_5| = 1` because there are exactly six atoms; the para-orthogonal polynomial is the cyclotomic
`Phi_14` itself. The lonely runner, the orthogonal polynomial, and the cyclotomic field are three
descriptions of one six-point set.

What this buys for the finish is a single refined invariant — the **binding fingerprint** — that both
programs can read: the binding modulus `q*` (the atom denominator), the gap `M = k/q*`, the atom count (=
Verblunsky termination = number of primitive roots), the Verblunsky ladder, and the atom-tournament's
parity invariants. In this language the covering-min theorem is one sentence: **covering forces `q*` up from
`n` to `Phi6`.** The tight AP lives on the `n`-lattice (six primitive `n`-th roots, `M = 1/n`) but is not a
covering set; a covering set must carry a multiple of every `q <= n`, and that constraint pushes its atoms
off the `n`-lattice onto the `Phi6`-lattice — two atoms at `+-t*`, `M = n/Phi6`. The depth `q*` jumps
exactly one continued-fraction level, `[0; n] -> [0; n-1, n]`, and that jump is the whole gap between the
floor and the covering-min.

The lesson is the one the owner has been pushing toward for several sessions: give the same structure enough
encodings — intervals, moments, Verblunsky coefficients, atoms, a tournament, a cyclotomic field — and the
encodings start to constrain each other. The tournament program knew about Paley and Rédei and the OCF; the
runner program knew about the covering-min and the lonely measure; neither alone saw that the tight
measure's atoms are a Paley tournament's vertices on the roots of `-1`. Standing where the two programs
overlap, the finish becomes a statement about a single arithmetic invariant, `q*`, that both sides can prove
things about — the tournament side through the dihedral symmetry of the atoms, the runner side through the
covering constraint on their denominators. The runner is a tournament seen on a loop; the atoms are where
they touch.

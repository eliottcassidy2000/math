# Late incoming voltage and character comparison

**Status: VERIFIED comparison and bounded proof review; FINITE-EXACT controls
only where specified below.** Read directly with `git show` from incoming commit
`20c9c7e0a500077380bc172d18f332d6c446caca`, before rebasing this worktree.
No frozen audit, producer file, or prior receipt was changed. The incoming
order-six census is consistent with our independent census; it was not rerun
in this late review. No external priority claim is made.

Incoming paths are under
`05-knowledge/results/session318_reservoir_descent_20261113/`:

| File | SHA256 of exact source bytes |
|---|---|
| `VOLTAGE_LIFTS.md` | `1de2a3c301cc51baf937625e40aa9d62a438f71b6374f84c28e2e4b32c036f16` |
| `CHARACTER_CHANNELS.md` | `b34d39d4287b42f83339f45a005c74acc3eeab91cc62051605dcde462c651ddf` |
| `check_wiring.py` | `8736ce0bf4bd9d80ac0a2f11ef5e5968ab5c798331034ad471761866833d553e` |
| `replay_characters.py` | `810162cb5d819bf3902d9acb1ecd214b00fee1cf213fe886f7c82336ddd870a1` |

Compared local snapshots: `pairings/SIGNED_LIFT.md` at
`b33b6554becc179a60879302750c14643386187664ac0c30444bc67303e284b8`,
and `README.md` at
`cf56124cfa82636d0f9a0a6d3f44254fecce9b2874c5d62ceb6515e70b0cedd7`.

## Exact overlap and the useful additional involution formula

The two-sheet block `[[a,b],[b,a]]`, signed difference `a-b`, gauge covariance,
zero/magnitude/diagonal guards, coboundary test, connected-support flip kernel,
and order-five/order-six collapse are the same mechanisms as our S1–S3/S5.
They are concurrent overlap, not additional novelty for our packet. Our
order-seven complete matrix, 184-vertex signed connectivity certificate, and
`x^88+x^456` result extend beyond the incoming wiring check's range.

The incoming note makes one useful specialization more explicit. For a
prescribed base involution `g`, append zero-labelled edges `i--g(i)` to the
ordinary signed-ratio constraints. These encode `e_i=e_g(i)`, the condition
that the chosen lift itself has order dividing two. If the combined graph is
balanced, and a baseline solution has `a_C` unflipped and `b_C` flipped
**g-fixed** chiral fibres in component C, its fixed-count factor is

    x^(number of g-fixed singleton fibres) product_C (x^(2a_C)+x^(2b_C)).

A component without fixed fibres contributes 2. This is correct: the two
solutions on a component exchange those two counts. Our S4 already imposes
`epsilon_i epsilon_g(i)=1` in its general cycle sum; the incoming expression
is a closed componentwise evaluation of that restricted sum. All lifts of a
base involution need not themselves be involutions.

A cheap independent hostile confirms why the added edges matter. Take four
chiral fibres, `g=(01)(23)`, and signed matrix

    D = [[0,0, 1,1], [0,0,-1,1], [1,-1,0,0], [1,1,0,0]].

Let merged weights be 2 on the support and zero elsewhere, so reconstruction
gives nonnegative integer ordinary blocks. Exhausting the sixteen flip-bit
vectors gives exactly `(0,1,0,1)` and `(1,0,1,0)` as wiring lifts. Both square
to the global sheet flip; neither is involutive. The appended zero constraints
therefore reject both. This is an abstract matrix control, not a claimed
realized tournament incidence example.

For an application requiring H, additional colors, or other observers,
restrict the allowed base permutations to preserve those labels as in our S3.
Ordinary dart degree alone is `H/Aut`, not H. The incoming theorem starts from
actual ordinary matrices; our parity/nonnegativity reconstruction guards remain
necessary when proposing formal `(A,D)` data instead.

## The character note retains genuinely additional sidecars

The incoming `A_v` counts all class incidences for every reflection-invariant
cube translation v. Its `A_kappa` is our ordinary C, while `A_0` is the
**presentation-mass** diagonal `diag(H/Aut)`. That mass is not our S4 diagonal
of converse-fibre sizes 1 or 2. Also its syndrome space named D is different
from our signed matrix D. These normalizations and types must remain explicit.

The new retained columns are indexed by a translation character and a
reflection syndrome. Expanding the two character sums proves the Gram identity
`sum_v chi(v) A_v = sum_d b_(chi,d) b_(chi,d)^T`. Reflection on a syndrome
coset translates by d, so each column has converse parity `chi(d)`. These
identities recover every translation matrix, not only its degree or the one
native signed block. The corresponding normalized operators are compressions
of commuting cube translations, but compression does not preserve products.
This strengthens the reason that an upstairs tensor decomposition need not
supply a downstairs class recurrence.

An independent literal cube/permutation calculation at orders three and four
verified every character Gram identity and column parity (16 identity gates).
At order four, ordinary codes `(0,2,4,5)` have masses `(1,1,5,1)`. The altitude
translation `{0,3}` and mirrored-pair translation `{{0,2},{1,3}}` give

    A_alt = [[0,0,1,0], [0,0,1,0], [1,1,2,1], [0,0,1,0]],
    A_pair= [[0,0,1,0], [0,0,0,1], [1,0,4,0], [0,1,0,0]].

The `(0,2)` entry, indexed by class codes, of
`A_alt M^-1 A_pair - A_pair M^-1 A_alt` is `-1/5`. Thus the two symmetrically
normalized operators do not commute. At order three the only nonidentity
translation gives no such failure. This reproduces the claimed first boundary
without importing either incoming checker.

The normalized odd block is entrywise dominated in absolute value by the even
block restricted to chiral fibres. If a connected component of the latter
has an edge to a fixed fibre, its corresponding substochastic restriction
leaks, giving spectral radius strictly below one. The domination and finite
leak argument are sound. This adds an analytic boundary sidecar beyond our
sign-support components; its hypothesis is not asserted at all orders.

The outer-ring formulas count characters and syndromes of the **presentation**
representation, and yield inertia upper bounds under compression. They do not
factor the class indicators. The separate cyclic-substitution formulas count
vertex automorphisms and anti-involutions, not Hamilton path marks. The stated
recurrences follow from three quotient rotations/reflections and recursive
block automorphisms; they supply no new scalar A000568 kernel by themselves.
No Clifford or Bott identification is needed for these real-character facts.

## Code scope and integration consequence

`check_wiring.py` enumerates all fixed-spine bit patterns and unrestricted
permutation minima at orders 4–6. It tests that weighted multiset refinement
is discrete and that the signed-support component sizes are `(1)`, `(1,1)`,
and `(22)`. Its explicit exception gate survives optimization. The code does
not independently audit every general lift statement, every color, or a
preferred path-level action; those require the corresponding proofs/sidecars.
`replay_characters.py` retains the full cube character decomposition and tests
its Gram identities and normalized noncommutation through six. In this late
review that script was read, not fully rerun; the independent control above
covers only three and four.

No mathematical conflict was found. Integration should credit the common
signed-lift result, retain the explicit base-involution component formula,
and identify all-translation character/syndrome data as additional information.
The current README's boundary remains correct: fixed-order walk closure and
presentation tensor growth do not themselves produce next-order classes or
the unrestricted Hamiltonian spectrum.

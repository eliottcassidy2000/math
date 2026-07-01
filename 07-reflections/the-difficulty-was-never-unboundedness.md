# The difficulty was never unboundedness

*mac-mini-2026-07-01-S78. Reflection on HYP-3794.*

For fifteen sessions the covering-min program carried one sentence as its load-bearing premise: *the
covering reduction does not bound the speeds, so unbounded far configurations are the real difficulty.*
Everything downstream inherited it. The proof was split into "bounded" (search, lazy-cut) and "unbounded"
(analysis). The analytic effort — equidistribution on the lonely set, the signed correction as a Fourier
tail, the `>= 7`-huge-speed cross-harmonic residual — was all aimed at the limit `w -> infinity`, and its
open piece was an *effective discrepancy bound*, the thing you need when speeds genuinely run to infinity.

They do not. Malikiosis, Santos, and Schymura proved in 2024 that to verify the Lonely Runner Conjecture
for `n+1` runners you need only check velocities up to `C(n+1,2)^{n-1}`. For our fourteen runners that is
`91^{12}`, about `3 x 10^{23}`. A counterexample, if one exists, lives entirely below that ceiling. There
is no infinity in the problem. There never was. The "unbounded far configurations" were a finite (enormous,
but finite) window that the reduction simply hadn't been told the size of.

This is worth sitting with, because it is a particular kind of error — not a wrong computation but a wrong
*horizon*. The repo grepped clean for Malikiosis, Santos, Schymura, Rosenfeld, and "product bound." The
frontier had moved — eight, nine, ten runners fell in 2025 and 2026 — and it moved by a method the project
had independently rebuilt and then walked past. Rosenfeld's proof for eight runners is elementary: a product
bound on top, prime filtering from below, and where the two meet, contradiction. No Fourier, no singular
series, no moment matrix. And his prime filtering — "if a prime `p` covers the first half of the residues
mod `(k+1)p`, it divides every counterexample" — is, line for line, the band-prime reduction the repo wrote
as HYP-3750 and then set aside because "counting cannot close Step 3." Counting couldn't; the *computational
per-prime check* could, and did, for the cases the world has actually closed.

So the reframing cuts two ways. It humbles the analytic program: the hardest-fought open piece, effective
equidistribution as `w -> infinity`, is answering a question the problem does not ask. And it re-elevates the
tool the program already had: prime filtering, the method of the actual frontier, sitting in the repo under
a pessimistic verdict. The right move is not to prove something about infinity but to *effectivize the
ceiling* — to push `3 x 10^{23}` down toward the searchable, one prime at a time, exactly as Rosenfeld did.

The pattern that transcends the theorem: **before pouring years into the hard part, check that the hard part
is real.** "Unbounded" and "finite-but-astronomical" feel the same from inside a proof — both defeat the
naive search, both invite analysis — but they are different problems with different tools, and the
literature may already have collapsed one into the other. A project deep in its own vocabulary is precisely
the project most likely to miss that the outside world has renamed its central obstacle. The cheapest,
highest-leverage act in a long research program is not another lemma; it is a literature search that asks
whether the mountain you are climbing is still there.

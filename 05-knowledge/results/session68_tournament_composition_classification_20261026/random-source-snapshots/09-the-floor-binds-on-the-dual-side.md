# The floor binds on the dual side

*mac-mini-2026-06-29-S31. The open question was whether the metagraph reproduces the cusp value 4cos²(3π/7). The answer made me move to the other metagraph — the one we keep calling first-class and keep not using.*

## Looking on the wrong side

I went in expecting the tournament metagraph to hand me `4cos^2(3\pi/7)`. It is the natural place to look: the cusp is the transitive limit, the binding mode is the 3-cycle (THM-588's Fiedler eigenfunction, the unique quadratic invariant), and the apex prime 7 should show up in the `Z_7` circulant tournaments. So I computed them — all eight circulant tournaments on `Z_7`, the size-three connection sets. Their autocorrelation gaps are `0.308` for the six generic ones and `2.0` for the two Paley/Fano sets `\{1,2,4\}, \{3,5,6\}`, which are also the `H`-maximizers (`189 = 27\cdot 7` against `175`). That is a clean, satisfying picture — the Paley tournament is the flat, optimal, most-Hamiltonian object, exactly as THM-586 says — and it has nothing to do with `0.198`. The binding value is *smaller* than the smallest tournament gap. It is sub-tournament. The tournament metagraph does not contain it.

The reason is one inequality: a tournament connection set on `Z_7` has size three (it must split `Z_7^*` with its negative), and the binding object is a **doublet**, size two. You cannot make a tournament out of two residues. The doublet is below the floor of what a tournament can express.

## The value was an even-graph eigenvalue all along

So I asked what a size-two set *is*. A doublet `\{a,b\}` depends only on its difference `d=b-a`; its autocorrelation lives on `\{0,\pm d\}`, which is `2I + A(C_7)` — twice the identity plus the adjacency of the seven-cycle `\mathrm{Cay}(Z_7,\{\pm d\})`. And then the value falls out exactly:

> `4\cos^2(3\pi/7) = 2 + 2\cos(6\pi/7) = 2 + \lambda_{\min}\!\big(A(C_7)\big).`

The binding value is the bottom of the seven-cycle's spectrum, shifted by two. The seven-cycle is the minimal connected `Z_7`-circulant **even graph** — the cusp of `E_7`, the even-graph dual metagraph. The number we were hunting was never a tournament invariant. It was an even-graph eigenvalue, and we had been looking for it in the wrong metagraph.

This is the kind of correction that should have been obvious in advance, because THM-588 told us. The tournament metagraph has *no* first-order invariant (`\mathrm{mult}(1)=0`): the cut space, the scores, the hierarchy — none of it survives the quotient. The only invariant is the cyclicity, the cycle-space part. The metagraph itself says the binding content is cyclic, not cut. The Lonely Runner floor is the `\mathrm{BOUNDED}` half, the second moment, the cycle space; the witness/score side is the cut space, and we already retired it as off-path. Of course the floor binds on the dual. The primal carries the optimum — the Paley tournament, flat spectrum, `H`-max; the dual carries the floor — the minimal cycle, the bottom eigenvalue. Optimum on `G_n`, floor on `E_n`. Cut and cycle. We have been writing "even graphs are first-class, `E_n` is the dual" at the top of every session and reaching for `G_n` every time.

## Dual minimal cycles, and the apex sets the length

S30 said the metagraph 3-cycle mirrors the LRC doublet, and that was almost right and on the wrong axis. They are not a same-side mirror; they are dual. The 3-cycle is the minimal relation in the cut-space metagraph `G_n` — a `Z_3` object, gap `1`, the Fiedler mode. The doublet is the minimal relation in the cycle-space metagraph `E_n` — a `Z_7` object, the seven-cycle, gap `2+\lambda_{\min}(C_7)`. Both are the minimal **cycle**; the apex prime sets its length. Three on the tournament side, seven on the even-graph side, and the seven is seven because `14 = 2\cdot 7`.

That last sentence is a formula. For the Lonely Runner at `n = 2p`, the apex cusp binds at the `Z_p` doublet, whose value is the `E_p` cusp eigenvalue:

> `2 + \lambda_{\min}(C_p) = 2 - 2\cos(\pi/p) = 4\sin^2\!\big(\tfrac{\pi}{2p}\big).`

I checked it down the family: `p=3` gives `1` (LRC6), `p=5` gives `0.382`, `p=7` gives `0.198` (LRC14), `p=11` gives `0.081`, `p=13` gives `0.058`. It decays like `\pi^2/p^2`. The floor is generous at small apex primes and tightens as the prime grows — larger apex, smaller floor, harder problem. There is a real asymptotic here: the per-level coupling the Lonely Runner can rely on shrinks quadratically in the apex prime, and `n=14` sits at `p=7`, where it is still a comfortable `0.198`. Whatever the full floor program needs, the atom it rests on is `4\sin^2(\pi/2p)`, the spectral gap of the minimal cycle in the dual metagraph.

## What I'll carry forward

The concrete next question is whether the descent's R-tail (THM-578) literally lands on `C_p` — whether the 2-adic descent's odd core, at the binding level, is exactly the minimal even graph, so the floor's last bound is verbatim the `E_p` cusp eigenvalue rather than a quantity that merely equals it. If so, the floor is closed by a single even-graph spectral fact, uniform in the apex prime. And the wider lesson is a habit, not a theorem: when the cut side goes quiet, the content is on the cycle side. The witness was cut and off-path; the floor is cycle and binding; the obstruction (S24) that cannot be a coboundary is the cyclicity; and now the binding value is a cycle eigenvalue in the dual. Every time the project has reduced, it has reduced to the cycle space. The even graphs were first-class the whole time.

See [[the-minimal-relation-binds-the-cusp]] (HYP-3585, the same-side version corrected here), [[even-graphs-as-first-class]] (the dual `E_n`), [[relations-not-things]] (the cyclicity is the only invariant), THM-588 (no cut invariant, only the quadratic), THM-586 (Paley = the optimum). New: HYP-3590.

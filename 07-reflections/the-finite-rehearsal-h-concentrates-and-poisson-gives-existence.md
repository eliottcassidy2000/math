# The finite rehearsal: H concentrates, and Poisson gives existence for free

*mac-mini-2026-06-29-S21. The owner asked to challenge assumptions, invent more approaches, and actually pursue one to a proof. I pursued the finite rehearsal and it gave a clean theorem and a new method. New: HYP-3560.*

## The rehearsal is provable, and it proves

The first of the owner's seven approaches said: the metagraph's exact moments are the same Siegel–Rogers machinery as the LRC, but finite and checkable — prove a variance inequality on the metagraph, then transfer. So I did. The number of Hamiltonian paths `H` in a random tournament has `E[H] = n!/2^{n-1}`, and its second moment, by relabeling symmetry, collapses to a single sum against a reference path. The coefficient of variation has a closed form that is a pure permutation statistic:

`CV(H)^2 = \frac1{n!}\sum_{\pi'} c(\pi')\,2^{j(\pi')} - 1`,

where `c(\pi')=1` unless `\pi'` reverses a consecutive pair, and `j` counts the consecutive pairs it keeps in order. The values — `1/3, 1/3, 19/60, 13/45, 131/504, 131/560` — decrease, and they decrease to zero, for a reason that is the whole point: the consecutive-integer adjacencies of a random permutation are asymptotically **Poisson(1)**, and `CV^2 + 1 = E[2^{\text{asc}}\,\mathbf 1(\text{desc}=0)] \to e^{1}\cdot e^{-1} = 1`. The two Poisson factors cancel exactly. **H concentrates**; its variance is Poisson-like, `Var(H)\sim E[H]`, the second moment diagonal-dominated with an off-diagonal `2^{\text{overlap}}` tail that `n!` outgrows.

This is not a toy. klein's THM-588, pushed the same day, proves the metagraph has *no* linear invariant and *exactly one* quadratic — the 3-cycle count. There is no first-order content to bound; the entire proof lives at the second moment. The rehearsal is not a warm-up for the real thing; it *is* the real thing, run on the one object where you can compute it to the last digit.

## Poisson does two jobs at once

The proof used a Poisson approximation, and that is the method worth keeping — sharper than I expected, because Poisson delivers two things the LRC needs from opposite ends. It gives **concentration**: the count is tightly peaked, `CV^2 \to 0`. And it gives **existence**: a Poisson(`\lambda`) count with `\lambda>0` is at least one with probability `1-e^{-\lambda}>0`. The second is the owner's approach 2 — forced existence, existence without construction — but reached through the probabilistic method and Chen–Stein rather than Ky Fan. If the count of lonely points is asymptotically Poisson with positive mean, a lonely point exists, and you never had to build it. Chen–Stein makes this quantitative: a total-variation bound between the true count and its Poisson limit, with an explicit error from the dependency neighbourhoods — exactly the resonance-overlap pairs that are the off-diagonal of the second moment.

So concentration and existence are not two problems. They are the two readings of one Poisson limit, and the dependency structure that controls the Poisson error is the same pair-overlap tail the metagraph rehearsal made me compute.

## The transfer, and why the congruence is what closes it

The LRC sheet count `N_R` has the identical shape: a diagonal that is its mean, plus an off-diagonal sum over pairs of sheets weighted by how much their resonances overlap. On the metagraph the overlap weight was `2^{j}`, `j` the shared consecutive arcs; on the runner it is a **congruence** condition — which speeds, mod 14, hit both sheets. That is precisely where Han–Lee enters (HYP-3553): the congruence Siegel second moment bounds the overlap tail with the covering built in as `\Gamma_0(N)`, and — the crucial word — *set-independently*. The metagraph rehearsal shows the diagonal-dominated Poisson structure is real and the tail is beatable; the congruence Siegel formula beats it uniformly over all covering sets, depending only on `N`. The Poisson limit then hands back both the concentration (THM-579's `CV(N_R)^2` small) and the existence (a lonely point, typically), with the worst case — the disproof — pushed into a large deviation that the set-independent bound forbids.

## The other doors I'd open

Pursuing one approach surfaced more. A **large-deviation rate function** for the H- or sheet-distribution would make the disproof an explicitly exponentially-unlikely event, and the rate is the worst-case the congruence bound must clear. The **Walsh/Fourier** basis that diagonalizes the metagraph (THM-584) is the same exponential-sum basis the LRC danger functions live in, so the metagraph's exact Walsh second moment is a literal template for `\sum|\hat c(14N)|^2`. And THM-588's vanishing linear invariant is a statement I had not appreciated: the **Eisenstein/score part carries no obstruction** — the LRC's binding content is purely the cusp form, the quadratic 3-cycle, matching the project's long-standing "the floor is R-even, the first moment vanishes." The score sequence, the thing that looks like it should matter most, is invisible to the obstruction; only the cyclicity survives the quotient. That is the kind of assumption-challenge the owner asked for, and it came free with the rehearsal.

See [[the-metagraph-is-a-finite-siegel-transform]] (HYP-3554, the second moment), [[the-covering-is-a-congruence-subgroup]] (HYP-3553, the set-independent transfer), klein THM-588 (no linear invariant), [[the-modular-tournament]]. New: HYP-3560.

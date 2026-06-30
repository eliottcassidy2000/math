# The metagraph variance was already a theorem (THM-219), seen from a new side

*mac-mini-2026-06-29-S22. The owner handed me a closed form for Var(H) and asked me to extend it. Verifying it walked straight into a theorem the project proved long ago — from a completely different direction. New: THM-589.*

## The closed form, and where it already lived

Last session I proved `CV(H)^2 = \frac1{n!}\sum_{\pi':\text{no unit descent}}2^{\#\text{unit ascents}}-1 \to 0`. The owner took that sum, called it `G(n)`, and handed back a closed form: `Var(H) = \frac{n!}{4^{n-1}}(G(n)-n!)`, with `G(n)` a count over compositions of `n` into odd parts, weight `k!\,2^{\#\{\text{parts}\ge3\}}`, OGF kernel `W(x)=x(1+x^2)/(1-x^2)`. I verified it three independent ways — successions, odd-compositions, the OGF — and they all returned `1, 2, 8, 32, 158, 928, 6350, 49752`.

That sequence is not new to the project. It is `W(n)` from THM-219 (kind-pasteur's NUD–Poisson–Euler connection), and it is the same `W(n)` from the S90/S112 simplicial-Rédei marathon — the one whose eighth term `49752` is the subject of MISTAKE-025 (the `+8` error). The metagraph variance is THM-219's `W(n)`, exactly. The owner's `G(n)` is THM-219's `W(n)`. My succession sum is THM-219's `W(n)`. Three sessions, three frameworks, one integer sequence, and it was sitting in canon the whole time wearing a different hat.

## Why this is a real unification, not a coincidence

THM-219 was stated abstractly: permutations with no unit descent, weighted by `2^{\#\text{unit ascents}}`, with `NUD(n)=A000255(n-1)` and a "Poisson–Euler" structure nobody had attached to a tournament quantity. THM-589 says what it *is*: the second moment of the Hamiltonian-path count over random tournaments. And the Poisson–Euler structure THM-219 noticed is *precisely* the Poisson(1) proof I gave last session. `CV^2+1 = E[2^{\text{asc}}\mathbf1(\text{no desc})]`; the no-unit-descent count is `NUD(n)\sim n!/e` (the **Euler** factor, `1/e`), and `E[2^{\text{asc}}]\to e` (the **Poisson(1)** generating function at `2`); their product is `e\cdot e^{-1}=1`. THM-219's two names — Poisson, Euler — are the two factors of the limit. The abstract permutation statistic and the tournament second moment and the probabilistic concentration are one thing, and the closed form is the bridge that made me see it.

## What the closed form buys: the exact rate

`\to 0` is concentration; the closed form gives the *rate*. The leading correction to `W(n)/n!=1` comes from the single odd-composition with one part equal to 3 — `(n-3)` ones and a single `3`, contributing `2(n-2)/(n(n-1))`. So

`CV(H)^2 = \frac{W(n)}{n!}-1 = \frac{2(n-2)}{n(n-1)}+O(n^{-2}) \sim \frac2n,\qquad n\,CV(H)^2\to 2,`

verified through `n=20`. The variance is `Var(H)\sim 2\,E[H]^2/n` — a sharp, simple concentration rate, and the `2` is literally the weight `2` that a part-of-size-`\ge3` carries in the OGF kernel. The combinatorics of the kernel *is* the asymptotics of the variance.

## The recurrence is the Hecke operator, now on the variance

`NUD` satisfies `NUD(n)=(n-1)NUD(n-1)+(n-2)NUD(n-2)`. That two-term recurrence is the vertex-addition operator — the modular `T`, the Hecke-like raising map `G_n\to G_{n+1}` I floated two sessions ago — now acting concretely on the *second moment*. The variance is not just a closed form; it is a linear recursion in `n`, the Hecke operator made arithmetic. The speculative "vertex-addition is Hecke" has a worked instance: it is how `Var(H)` is built from one `n` to the next.

## The transfer, sharpened by klein's counterexample

The day I was proving the metagraph side, klein (S4) was testing the LRC side, and found the thing that makes the rehearsal matter: the LRC's `CV(N_R)^2` is **set-dependent and unbounded** — push `R` dense and add speed 7 and it blows up. The metagraph variance never does this: `CV(H)^2\sim2/n`, set-free, classical. The contrast is the whole argument. The runner's second moment is unbounded *because it is computed set by set*; the cure is to stop computing it set by set, which is exactly the `\Gamma_0(N)` congruence move — make the moment depend only on the covering modulus `N`, and the unbounded corner klein found is absorbed into a subgroup index. The metagraph `W(n)` is the template that *is* bounded; Han–Lee's congruence Siegel formula is the lift that makes the runner's moment bounded too, by the same Riordan-classical second-moment structure rather than a bespoke one. The obstruction is a classical succession count. That is the most reassuring thing the closed form says: there is nothing exotic at the second moment, only `A000255`, Poisson, and Euler.

See [[the-finite-rehearsal-h-concentrates-and-poisson-gives-existence]] (HYP-3560, the Poisson(1) proof), THM-219 (NUD–Poisson–Euler, the same `W(n)`), [[the-covering-is-a-congruence-subgroup]] (HYP-3553, the `\Gamma_0(N)` lift), [[the-metagraph-is-a-finite-siegel-transform]] (HYP-3554). klein THM-588/S4. New: THM-589.

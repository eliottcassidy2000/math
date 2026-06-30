# Two properties of one relation

*mac-mini-2026-06-29-S26. The owner asked to use the recent ideas — the obstruction lens, relations-not-things, the moment method, Γ₀(N) — to reframe the remaining LRC targets. They collapse.*

## Five obligations were always two

For a long time the LRC(14) proof has carried a list. There is the witness density (THM-527, "the linchpin"): show `G2 > 0` forces `M ≥ 1/14`. There is the witness floor for `k = 8..13`. There is the gK8 concentration, `max_E ≤ 10·cap`. There is the doublet R-tail, `R(M) = O(1/M)`. And over all of them the gatekeeper: the covering floor `R' > 0`, equivalently `CV(N_R)^2` bounded. Five analytic obligations, each with its own machinery, each re-proved per configuration.

Seen through the last six sessions, they are two — and the two are not new obligations but two *properties of a single object*. The object is the danger relation `D(v,t) = [\lVert vt\rVert < 1/14]` composed with itself, `D \cdot D^{\top}`. That composition is simultaneously the second moment, the pair-correlation that klein's THM-588 proved is the *only* invariant, and the R-equivariant obstruction class whose measure (S23) is the floor. Everything left to prove is two facts about it:

> `D \cdot D^{\top}` is **essential** (it does not factor) and **bounded** (`CV(N_R)^2` is finite), under the right change of base.

Essential gives existence — a lonely point. Bounded gives uniformity — the gatekeeper. They are the σ-odd and σ-even halves, the two measures of one obstruction.

## Which obligation is which half

The witness side is essentiality. The witness density `G2 > 0` (A) and the witness floor (B) are, reframed, the statement that the obstruction is *nonzero* — its counting measure `\varphi(n)/2`, the saddle index, the R-odd part that survives at the extremal where the floor's Lebesgue measure vanishes (S23). And essentiality has a name now: the relation does not separate. `D(v,t) = f(v)g(t)` would be a coboundary, rank one, a danger that factors into a speed-part and a time-part and therefore covers — a disproof. The product `vt` inside `\lVert vt\rVert` is the one bilinear thing that refuses to factor (HYP-3564), so `D` is full rank, essential, and the lonely point exists because the relation cannot be written as a single product. The witness obligations stop being "prove an inequality for every config" and become "the class is essential" — a topological, set-independent statement.

The floor and cap side is boundedness. The gK8 concentration (C) is the second moment controlling the maximum — Chen–Stein, the Poisson tail (S21). The doublet R-tail (D) is literally the tail of `D \cdot D^{\top}`, the pair-correlation decaying. The gatekeeper `CV(N_R)^2` is the variance of that composition. These are not three problems; they are the boundedness of one second moment, and klein-S4 told us exactly what is wrong with the way we have been bounding it: *per configuration it is unbounded*. That is not a failure of estimation. It is the symptom of collapsing along `\mathbb{Z}_{14}` — the covering — which degenerates at the cusps, where the metagraph's collapse along `S_n` stays clean (klein-S5). The cure is not a sharper bound; it is a change of base. Pull the second moment back along the `\Gamma_0(N)` congruence (HYP-3553), and it depends only on `N` — `\psi(14)=24`, `\varphi(14)=6`, `J_2(14)=144` — and on nothing about the speeds. "The floor must manufacture the transitive symmetry it lacks," klein wrote; that is a relational instruction, and it discharges C, D, and the gatekeeper at once, with one congruence second moment instead of three analyses.

## The rehearsal is finished, on the metagraph

What makes this more than relabeling is that both halves are *already proved* on the finite mirror. Boundedness: `CV(H)^2 = W(n)/n! - 1 \sim 2/n` (THM-589) — the second moment of the Hamiltonian-path count is bounded, cleanly, under the `S_n` collapse, with the exact rate. That is the theorem the LRC's `CV(N_R)^2` must match, and the only difference is the base. Essentiality: `b_1^-/b_1 \to 1/2` (HYP-3565) — the R-odd obstruction is *half* the metagraph's cycle space, robustly nonzero, not a marginal coincidence; and `SC = \mathrm{trace}(R) > 0` (THM-587) forces self-complementary tournaments to exist by Lefschetz, without anyone building one. The metagraph has already shown both: a bounded second moment and a robustly essential obstruction. The LRC asks for the same two facts about the same kind of object, with the covering as the only complication, and the covering is a change of base.

## What is actually left

Stated this way, the proof has one concrete next computation and one structural fact. The computation: evaluate the `\Gamma_0(14)` congruence second moment (Han–Lee, arXiv:2507.05905) and check it bounds `CV(N_R)^2 < \mathrm{cap}_r/(1-\mathrm{cap}_r)` set-independently. The fact: the danger relation does not factor — essential — which the bilinearity of `vt` already gives in spirit and the Borsuk–Ulam counting measure certifies at the extremal. Five obligations become essential-times-bounded; the witness/floor dichotomy becomes the counting/Lebesgue measures of one obstruction; the covering becomes the symmetry you change base to; and the whole thing reads, finally, as a sentence about a relation: *it does not factor, and composed with itself it stays small, once you look at it in the right frame.*

See [[relations-not-things]] (HYP-3564), [[the-measure-of-the-obstruction]] (HYP-3562), [[the-covering-is-a-congruence-subgroup]] (HYP-3553), [[the-metagraph-variance-was-already-a-theorem]] (THM-589), [[the-finite-rehearsal-h-concentrates-and-poisson-gives-existence]] (HYP-3560). klein-S4/S5 (the reference-collapse), THM-527/579 (the obligations). New: HYP-3570.

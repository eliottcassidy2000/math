# The measure of the obstruction

*mac-mini-2026-06-29-S23. The owner asked to consider obstruction theory, the objects we have studied, and their measure. The three turn out to be one thing seen three ways.*

## One class, named at last

For two dozen sessions the project has produced existence theorems and measured them, without once naming what was being measured. Rédei: every tournament has a Hamiltonian path (`H \ge 1`). Eplett: self-complementary tournaments exist (`SC > 0`). The Lonely Runner: a lonely time exists (floor `> 0`). klein: the metagraph is connected (algebraic connectivity `4`). These look like four problems. They are one: each is the **nonvanishing of the equivariant Euler class of the complement involution `R`**, and each "measure" we have computed — the path count, the `SC` count, the lonely measure, the spectral gap — is that one class evaluated by a different measure.

The cleanest instance is the metagraph, because there obstruction theory is exact and finite. `SC = P_n(-1) = \mathrm{trace}(R)` is the **Lefschetz number** of the complement. The Lefschetz fixed-point theorem says: if `\mathrm{trace}(R) \ne 0`, then `R` has a fixed point. `SC = 2, 2, 8, 12, 88` — always nonzero — so a self-complementary tournament always exists, and *we never had to build one*. That is existence without construction, the topological kind, and it is the finite rehearsal of the lonely point existing because the danger cover has a hole.

## The obstruction has two measures, and that is the whole `σ`-story

Here is the thing the computation made me see. The obstruction does not have *a* measure; it has *two*, and the project's `\sigma`-even/`\sigma`-odd split — the split that runs through every object, the floor-versus-witness, the SOS-versus-Borsuk-Ulam, the `R`-even-versus-`R`-odd — *is* the split between them.

The `\sigma`-even side is the **Lebesgue measure**: the floor, `\mathrm{meas}(\text{lonely})`, the bulk, the part you bound with sums of squares. The `\sigma`-odd side is the **counting measure**: the number of lonely points, the units mod `n`, the Euler characteristic, the index. And they are not redundant, because the Lebesgue measure can vanish while the counting measure does not. I checked it: the extremal `\{1,2,3\}` has Lebesgue floor exactly zero — and a lonely set that is two isolated points, the units mod 4. The covering `\{2,3,4\}` has Lebesgue floor `1/8`. So at the extremal — exactly the hardest configuration, the would-be counterexample — the measure you have been bounding *disappears*, and the only thing left holding the obstruction up is the **count**. That is precisely why the Borsuk–Ulam index, the `\sigma`-odd witness, the odd degree `(p-1)/2`, has to exist as a separate object: it is the measure that survives when Lebesgue dies. The witness is not a second tool bolted onto the floor; it is the *other measure of the same class*, the one that detects the obstruction at the extremal.

`SC = (R\text{-even}) - (R\text{-odd})` is the alternating sum of the two — the Euler number is the Lebesgue measure minus the counting measure, the Lefschetz number that is nonzero precisely because the two do not cancel.

## Disproof is a coboundary

Phrasing it as a class sharpens what a disproof would be. The obstruction is a cohomology class; a disproof is the class being **exact** — a coboundary — every measure of it zero, the hole filled, the lonely set empty. A proof is the class being **essential**. And essentiality is a *topological* property: it depends on the homotopy type of the quotient, on the congruence structure of the covering, not on which speeds you chose. That is the deepest reason the `\Gamma_0(N)` move is right and the set-by-set bound is wrong. klein found that `CV(N_R)^2` is unbounded set by set; of course it is — measuring a topological invariant pointwise on a moduli space will diverge at the cusps. The invariant itself, the essentiality of the class, lives upstairs in the equivariant cohomology of the level-`N` quotient, indexed by `N` and nothing else. The floor's positivity is not a delicate analytic estimate you re-verify per speed set; it is the essentiality of an Euler class, and you read it off the topology once.

## And the moment method is how you weigh it

Obstruction theory says *when* the class is essential. The measure says *how much*. The moment method (THM-589) is the bridge: the first moment `E[\#\text{lonely}] = E[\#\text{fixed}]` is the Siegel–Burnside mass — the average measure, positive iff existence is typical; the second moment `W(n)`, the metagraph variance with `CV^2 \sim 2/n`, is the concentration that turns "average positive" into "always positive." So the equivariant homology (klein, the `R`-odd Betti numbers), the moment method (the mean and `W(n)`), and the Lefschetz number (`SC = \mathrm{trace}(R)`) are three calculations of one measure-valued obstruction — the homological, the probabilistic, the fixed-point. They had to agree, and they do.

Name the class. The floor stops being an estimate and becomes the essentiality of the Euler class of the complement; the witness stops being a second method and becomes the counting measure that outlives Lebesgue at the extremal; the covering stops being a constraint and becomes the level structure the class lives over; and the disproof becomes the one thing a nonzero cohomology class cannot do — bound a coboundary.

See [[the-metagraph-variance-was-already-a-theorem]] (THM-589, the second moment), [[the-covering-is-a-congruence-subgroup]] (HYP-3553, the set-independence = essentiality), [[the-one-involution-three-spectra]] (the `R`-spectrum), klein HYP-3544 (equivariant homology), THM-581 (the Borsuk–Ulam witness), HYP-3242 (`\chi_{\text{meas}}`). New: HYP-3562.

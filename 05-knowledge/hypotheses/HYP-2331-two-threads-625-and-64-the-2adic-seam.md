# HYP-2331 — Two parallel threads on the arc's 2-adic seam: Erdős 625 (σ, the order-2 involution) and Erdős 64 (2^k, the doubling tower)

**Session:** S653
**Status:** CONFIRMED (625 σ-bound formalized + verified; 64 Sidon/summand framing, with an honest
self-correction; both are open/near-open and I claim no resolution)
**Provenance forward:** math-lean `Math/Combinatorics/CochromaticBound.lean` (sorry-free)
**Prompt:** continue 625; integrate a second parallel thread (Erdős 64 = Erdős–Gyárfás) via Sidon
sequences / the "cauldron game" / the "summand graph"; work back and forth.

---

## Thread A — Erdős 625, cont.: the σ-asymmetry bound (formalized)

Closing the S652 handoff. `ζ` (cochromatic) is bounded by *both* chromatic numbers across `σ = complement`:

> **`ζ(G) ≤ min(χ(G), χ(Gᶜ))`**, hence **`χ(G) − ζ(G) ≥ χ(G) − χ(Gᶜ)`** — Erdős 625's gap is at least the
> *σ-asymmetry of the chromatic number*.

*Proof.* A proper colouring of `G` has independent classes (`ζ ≤ χ`); a proper colouring of `Gᶜ` cuts `G`
into cliques (`ζ ≤ χ(Gᶜ)`). **Formalized (math-lean, sorry-free):** `cochromColorable_of_colorable`,
`cochromColorable_of_compl_colorable`. Verified on small `G(n,½)` (`two_threads_625_and_64_s653.py`):
`ζ ≤ min(χ, χ̄)` always; the bound `χ−ζ ≥ χ−χ̄` holds.

---

## Thread B — Erdős 64 = the Erdős–Gyárfás conjecture (1995), via Sidon / summand / cauldron

**The problem (verified):** every finite graph with min degree ≥ 3 contains a cycle of length `2^k`
(a *power* of 2: 4, 8, 16, …). OPEN; `$100`/`$50`. (The literal "even cycle" version is the *solved*
THM-443; the open one is powers of 2.) Counterexample needs ≥ 17 vertices (cubic ≥ 30); Markström found
24-vertex graphs whose only power-of-2 cycle is `16`.

**The cycle = additive-relation dictionary.** A cycle of length `L` is an `L`-term `±` relation among the
edge steps. So a `2^k`-cycle is a `2^k`-term additive relation, and *avoiding* short cycles is *avoiding*
short additive relations — Sidon territory.

**Honest correction (caught this session).** I first claimed "Sidon ⟹ no 4-cycle in the Cayley graph." It
is **false**: an abelian Cayley graph `Cay(ℤ/n, {a,b,…})` *always* has the **parallelogram** 4-cycle
`(0, a, a+b, b)` (verified). The Cayley graph is the *wrong* model. The correct object is the repo's
**summand graph**:

> **`S` is Sidon `⟺` no two distinct pairs of `S` have equal sum `⟺` the summand graph on `S` is
> C₄-free.** (A 4-cycle there is exactly `a+b = c+d`, a Sidon violation.) Verified exactly: `Sidon ⟺
> summand-C₄-free` on all test sets. **So the summand graph (repo) IS the Sidon / C₄ detector.**

**The connection to Erdős–Gyárfás.** Sidon sets give the **extremal C₄-free graphs** (Erdős–Rényi/Brown)
— the densest graphs with no 4-cycle. *These are exactly the hard instances of Erdős–Gyárfás:* a C₄-free
min-degree-3 graph must (conjecturally) still have an 8- or 16-cycle, and Markström's near-counterexamples
(only a 16-cycle) live here. **To *disprove* 64 you would have to avoid `C₄` and `C₈` and `C₁₆` and …
simultaneously** — a Sidon-type condition at *every doubling level* (a "doubling-tower Sidon") while
keeping min degree ≥ 3. That stacked avoidance is why the conjecture is believed true; it is the additive
heart of the problem.

**The cauldron game** (repo, S618) is the same engine one level down: colour `[1,N]` avoiding a
monochromatic sum `a+b` (a Schur "boil" = a 3-term relation); Schur's theorem = you *cannot* avoid boils
forever. Erdős–Gyárfás is the `2^k`-cycle analogue: you cannot avoid power-of-2 additive relations forever
at min degree 3.

---

## The cross-connection: both threads live on the arc's 2-adic seam

This is the back-and-forth payoff:

| | the "2" | the object |
|---|---|---|
| **Erdős 625** | the **order-2 involution** `σ` (complement) | `ζ` = the σ-fixed chromatic number; `χ−ζ` = σ-asymmetry |
| **Erdős 64** | the **powers of 2** (the doubling tower) | a `2^k`-cycle = a doubling-closed additive relation |

> **625 is the order-2 (involution) face and 64 is the `2^k` (doubling) face of the *same prime 2*.** The
> arc's **2-adic seam** — the half-turn `σ` (S638/S643) and the doubling tower `⟨2⟩` (S640) — is exactly
> where both Erdős problems sit. The cube-root / odd machinery that dominated S637–651 is *orthogonal* to
> both: **these two are the project's pure 2-adic Erdős problems**, one for each meaning the arc gives the
> prime 2.

## Convergence with concurrent sessions (S709, S710) — independent agreement, not overturn

While this ran, the fleet was on the same problem from the other direction:
- **S709 / MISTAKE-064 / THM-445** flagged that Erdős 64 is the **power-of-2** (Erdős–Gyárfás) version,
  *not* the even-cycle version (THM-443/S708). Thread B here independently states the same distinction
  (the even version is "the *solved* THM-443"); we agree, no conflict.
- **S710 / HYP-2314 / THM-446** ("Sidon × cauldron × summand-graph → Erdős 64") reached, independently and
  in more depth, the *same* Thread-B core: Sidon ⟺ C₄-free summand graph; cauldron/Schur (3-term) ⊂ Sidon
  (4-term = first power of 2) ⊂ B_h (2h-term) as one additive-relation ladder; the hard core =
  "dyadically-Sidon" / high-girth graphs; and it adds the additive-energy quantification
  (`E(S)=2|S|²−|S|` floor, via S706/THM-441) and cage/random-cubic verification I did not run.

**This is convergence, not duplication to overturn.** S710 is the more complete Thread-B artifact; HYP-2331's
distinct contributions are (i) **Thread A** — the *formalized* σ-bound `ζ ≤ min(χ,χᶜ)` (math-lean d7570e3),
which S710 does not touch, and (ii) the **cross-connection**: framing 625 and 64 as the *two faces of the
single prime 2* (involution σ vs doubling tower) on one 2-adic seam — the "back-and-forth payoff" the prompt
asked for. Readers wanting the deep Sidon/cauldron development should use HYP-2314; this doc owns the
625-side bound and the two-problems-one-prime synthesis.

## New threads / handoffs
- **625:** formalize `χ − ζ ≥ χ(G) − χ(Gᶜ)` numerically (tie `CochromColorable`/`Colorable` to
  `chromaticNumber`); does the chromatic σ-asymmetry `χ(G) − χ(Gᶜ)` itself grow (the route to the
  Heckel–Steiner `n^{1/2}`)?
- **64:** formalize "Sidon ⟺ summand-graph C₄-free"; explore the doubling-tower-Sidon obstruction (avoid
  `C_{2^k}` ∀k) — can min degree ≥ 3 coexist with it on Cayley/incidence graphs? (the additive route to
  Erdős–Gyárfás).
- **Cross:** is there a single "2-adic forcing" lemma — min-degree/σ-structure forces *either* a σ-asymmetry
  *or* a `2^k` cycle — unifying the two seams?

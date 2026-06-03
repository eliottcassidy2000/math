---
source: claudebox-2026-06-03-S585
status: REFLECTION + RESULT — two-lens (functional ⟷ rep-theory) reading of Lemma A/B; the
  good-measure as a relation-lattice theta graded by character length; the AP-translation flip.
tags: [LRC, representation-theory, characters, fusion, theta-function, relation-lattice,
  additive-energy, circuit-free, Lemma-A, Lemma-B, currying, tournament-arcs]
---

# Currying produces arcs, fusion produces hardness: two lenses on the lonely runner

**Prompt (user):** think in terms of abstract functions (Haskell) and representation theory —
curried functions produce tournament arcs; bounce between the two whenever stuck; use the
particularities of one as noise to inspire the other. Lemma A (randomness, circuit-free ⇒ G≥δ)
and Lemma B (3-term relation = a literal fold).

I want to record the *path*, because the bouncing is where the result came from.

## Bounce 1 — rep theory: a 3-term relation is a fusion rule

The characters of `ℝ/ℤ` are `χ_v(t)=e^{2πivt}`, and they multiply by adding indices:
`χ_a·χ_b=χ_{a+b}`. So Lemma B's "3-term relation `v_c=v_a+v_b` is a fold" is, said in the dual,
a **fusion rule** `χ_a⊗χ_b=χ_c` — the speed set fails to be closed under tensor product. That
reframing told me what "circuit-free" should mean: the characters are *multiplicatively
independent*. And it suggested the right invariant of the whole set: the Laurent polynomial
`p_S(x)=Σ x^v`, a virtual character, with 3-term structure in `p_S²` and 4-term energy in `|p_S|⁴`.

## Stuck — why is circuit-free actually *safe*? → Bounce 2: Haskell

Knowing 3-term = fusion didn't tell me why no-fusion gives a lonely time. Here the functional lens
unstuck me. The all-lonely measure is `∫ Π_i (1 − 1_{B_i})` — a **fold over the runners** of the
single-runner good-indicators. `foldr (*) 1` of a list of functions is a *separable* object, and
the integral of a separable product is a constant-term extraction. Writing `1_{‖vt‖≥δ}=G(χ_v(t))`
with one fixed Laurent series `G`, the fold becomes `∫ Π_i G(z^{v_i}) dz/z = [z⁰]Π_i G(z^{v_i})`.

## Bounce 3 — back to rep theory: it's a theta over the relation lattice

That constant term has a name. `[z⁰]Π_i G(z^{v_i}) = Σ_{m: Σ m_i v_i=0} Π ĝ(m_i)` — a sum over
the **relation lattice** `Λ=ker(v)`. It is a theta function of `Λ`, equivalently the multiplicity
of the trivial character (a Molien constant term) in `Π G(χ_{v_i})`. And the `m=0` term is exactly
`(1−2δ)^k` — the independence baseline, the `e^{−2}>0` that *is* Lemma A's mechanism. Every other
term needs a genuine relation. So **Lemma A = "the relation-lattice theta sum is dominated by its
constant term."** The whole dichotomy is now one picture: circuit-free = no short lattice vectors
= constant-term domination; 3-term = a short vector = a real correction.

## Stuck — why is 4-term *richness* harmless? → Bounce 4: a particular as noise

The user flagged that high-energy (4-term-rich) circuit-free configs are also safe, and that this
had been a worry. I couldn't see why from the lattice alone — until a *particular* threw off noise.
Searching for the highest-energy 3-term-free set, the computer returned `{14,16,18,20,22,24}` =
`2·{7,…,12}`: a **translated arithmetic progression**. It has the *same* additive energy as
`{1,…,6}` (146 — energy is translation-invariant) but *zero* 3-term relations (the sums escape the
set) — and it is safe (margin +0.225), while `{1,…,6}` is tight. That one example was the whole
answer. It says energy and 3-term-count are *independent* knobs, and only the latter is the gap.

Back in the theta picture this is a clean grading. A length-`L` relation contributes `Π ĝ(m_i)`
with `L` factors of `ĝ(m)=O(sin 2πδ)=O(δ)=O(1/(k+1))`. Two-term relations are impossible (distinct
speeds), so the minimum is length 3. Therefore
`meas(good) = (1−2δ)^k + O(N₃·δ³) + O(δ⁴)`, and **each extra term in a relation costs a factor
`~(k+1)/2`** (verified: the 3-term/4-term correction ratio grows 1.98→4.75 across k=4..10). The
4-term energy lives one full order of `δ` down from the 3-term hardness — suppressed, and more so as
`k` grows. That is exactly why 4-term-rich configs stay on the safe `(1−2δ)^k` plateau.

## Bounce 5 — Haskell again: currying makes the tournament

The repo's tournament arcs finally fell out of the same `ev`. Curry `ev v t = {v·t}` not on `v`
(which gives the character) but evaluate at the *witness* time `t*` and compare runners: the arc
`a→b` is `sign Im χ_{v_a−v_b}(t*)`. The arcs are indexed by the **difference set `S−S`**, not the
speeds themselves — the tournament is a function of the *differences*, just as the energy is. For
`{1,2,3,5,8}` (Fibonacci, a Sidon set) the witness tournament has score sequence `[1,2,2,2,3]` and
all 15 differences distinct. The orientation lives on `S−S`; the hardness lives on the 3-term
sums in `S`. Difference structure vs sum structure — the same `±` that splits cut-space from
cycle-space in the tiling model splits the tournament (differences) from the loneliness obstruction
(sums).

## The synthesis

One object, three faces: `meas(good)` is (functionally) a fold over runners, (analytically) a
theta sum over the relation lattice `Λ=ker(v)`, (rep-theoretically) the trivial-character
multiplicity in `Π G(χ_{v_i})`. Lemma A is constant-term domination, quantified by a length grading
that makes the 3-term count `N₃` the unique leading obstruction and demotes energy to a `δ`-smaller
correction. Lemma B is what happens when `Λ` has a short vector: the fusion confines the orbit and
the theta sum tips into a resonance that can cancel the constant exactly (`{1,…,6}` → 0). The path
mattered: every step I was stuck in one lens, the other — or a particular it coughed up — moved me.

## The transcending pattern

Differences orient; sums obstruct. Currying the evaluation map one way gives the characters whose
*sums* (fusions) are the hardness; the other way gives the arcs indexed by *differences*. The
lonely-runner gap is the constant term of a representation built from a speed set, and the whole
difficulty of the conjecture is the few low-degree invariants — the 3-term fusions — that perturb
that constant. Keep pulling and the conjecture looks like asking when a Molien series stays positive.

**Artifacts:** `04-computation/lrc_character_sieve_exploration_s585.py` (+`.out`); new **HYP-2120**.
Builds on HYP-2064/2065 (sieve vs thin core = the same split), HYP-2115 (CRT/sieve thread).

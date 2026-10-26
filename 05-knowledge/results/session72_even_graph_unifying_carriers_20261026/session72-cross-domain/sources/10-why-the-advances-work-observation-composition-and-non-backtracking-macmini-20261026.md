# Why the advances work: observation, composition, signs, spectra, and the non-backtracking walk

*mac-mini-2026-10-25-pyth, 2026-10-26. A REFLECTION. Under the repository's precedence rule it records why things
worked, not what is true. Every mathematical claim below points to its canon or result file, and its status is the
status of that file. The owner asked for "a cohesive and unifying set of theories that understand how and why we have
been able to make all of the novel advances in this repo".*

**What already exists and is not repeated here.**
- kind-pasteur S128c140, *the presentation-first method*: every problem as (object, monoid, action).
- synthesis_20260905_overview: recovered missing coordinates.
- strong_atom_continuation_vitali HISTORY.md: observations that survive operations.
- deathstar S76: twelve meta-lenses.
- macmini S159: the sign-reversing involution engine.
- kps S135: the first-untested-point law.
- opus S399: the LRC triage law.
- The 60 cards of META-PATTERNS.

This note adds a single schema that organises those documents. It shows the schema working on today's results, where
four independent sessions converged on one theorem. It ends with the predictions the schema makes.

## 0. The schema

Each advance is a statement about a quadruple

    (X, M, O, ρ)

- **X**, an object or space: labelled tournaments, primitive triples, odd integers, runner speeds.
- **M**, a monoid of operations acting on `X`: arc flips, complement lines, re-markings, Berggren branches, Collatz
  steps, Farey moves. This is kind-pasteur's presentation.
- **O**, an observation `X → Y`: H, a Hamilton-path language, a residue class, a radical, a score sequence, a census.
- **ρ**, a level: a modulus, a vertex count `n`, a height, or a depth. The observation is computed exactly up to it.

The repository's theorems answer one of five questions about this quadruple.

| Question | Form | Canonical instances |
|---|---|---|
| **Reconstruction** | Is `O` injective? If not, what is the fibre? | THM-4807 and three concurrent proofs; the `e*h*e` degree lemma (HYP-9416 audit); THM-4804 module shattering (lines destroy modules) |
| **Descent / continuation** | Is `O` an `M`-congruence? If not, what sidecar repairs it? | THM-840 (sidecar rank = `dim N(ker R)`); vitali C1–C5; Session70 (missing linear rank `2^k − d_t`); the recovery principle (a target descends iff it is constant on fibres) |
| **Composition** | Does `O` factor through a composition law on atoms? | THM-1370 step 1, THM-1862, THM-1950, THM-2183 (ordered joins; H multiplicative); the Berggren free monoid (unique words); aliquot and F1 evaluations |
| **Cancellation** | Does an involution act on the fibres of `O` with a sign? | THM-062, THM-068, THM-589, THM-1805, S159; the THM-4808 cusp pairing |
| **Level** | How does `O` behave under finite quotients as `ρ` grows? | the twin CRT obstruction and THM-4808 (c) (fixed `ρ`); growing_sieve_path W1–W3 (exact level `√x`); THM-4628, THM-4663, THM-4780 (character sums); THM-4808 (G) (non-backtracking spectra) |

An advance happens when the right quadruple is named, and then either exactness is proved or the exact obstruction is
computed. Nearly all of this repository's "novel" content is of the second kind. It names the missing coordinate,
finds the hostile pair, and computes the rank, the fibre or the residual.

## 1. Observation: complete, lossy, and growing

The five questions are three regimes of one question: how much of `X` does `O` keep?

**Complete observations.** An observation can look lossy and still be injective. The mechanism is always the same: a
structural theorem pins the invisible part to a visible certificate.

- **Legal paths.** Hamilton paths miss some arcs. THM-4807 shows that every missed arc `u → v` is the chord of a legal
  path from `u` to `v`, and that reversing it would close a Hamiltonian cycle. Its orientation is therefore written in
  the data, as the endpoints of a legal word.
- **The concurrent proofs.** The other three use a different certificate, the Busch–Jacobson–Reid separator `z`.
  Session70 reads it through two traced two-step detours, guard_reconstruction through two-step counts, and
  growing_sieve_path through one-way reachability after deleting `z`.
- **The sandwich.** The audit lemma behind HYP-9416 (1) reconstructs a strong two-port state from `e*h*e` by degree
  recognition.

Lesson: before declaring an observation lossy, look for a forced witness. A structure theorem about the invisible
objects usually provides one (META-PATTERNS #19, "a nearly true theorem hides a second coordinate").

**Lossy finite observations.** A fixed modulus sees only predicates that congruences force.

- **Twins.** The CRT residue obstruction shows that every twin class modulo `M` also contains composite pairs.
- **The branch character.** THM-4808 (c) gives the family `(16M+2, 4M+1) ≡ (2, 1)` with opposite character. The
  determinant of a Berggren word is invisible to every fixed modulus, even though the integer determines it.
- **Shared mechanism.** Strong approximation makes every finite quotient surjective. Any predicate that is not
  periodic is then lost.
- **Positive use.** Fixed moduli prove impossibility. Mod 3 gives the single edge of `J`, and Hardy–Littlewood
  admissibility describes the patterns a fixed modulus kills.

**Growing observations.** Between these extremes is the observation whose level grows with the object.

- **Eratosthenes.** Level `√x` decides twin-ness exactly (growing_sieve_path W1–W3).
- **Density.** At that exact level the naive density is wrong by the Mertens factor. Our pre-registered ratios were
  1.1116, 1.1312 and 1.1448 at `10⁸`, `10⁹` and `10¹⁰`, rising toward `4e^(−2γ)`.
- **Beyond.** Sieve weights cannot pass the parity barrier. Brun's upper bound and Chen's P2 are the transported
  limits.

**Why this regime was productive here.**
- Exhaustive small-`n` enumeration turns "what does `O` lose" into a finite fibre computation: 191,536 classes at
  `n = 9`, and 9,355,949 strong classes at `n = 10` in the audit.
- A hostile pair turns a vague doubt into a falsifiable statement.
- The repository insists on typed status: PROVED, FINITE-EXACT, CITED, OPEN. That keeps the regimes from being blurred
  into each other.

## 2. Composition: atoms, joins and their crevasses

**The tournament factorisation.** `HP(A ⇒ B) = HP(A)·HP(B)`. Three consequences:
- H is multiplicative over strong components (THM-1370 step 1);
- the components can be read off the cut positions of the path language (THM-4807 (e));
- every reduction "to the strong case" in the repository rests on it: THM-1862, THM-1950, THM-2183, THM-4806 (R)(ii).

**Gaps come from the atoms' floor.** The gaps 7 and 21 live in the crevasses of the floor of strong tournaments: 3,
5, 9, 15, 25, 45, … combined with multiplicativity. S60 compares the aliquot gaps `{2, 5}` to this as a "crevasse
dictionary". That comparison is typed as ANALOGY.

**Sitting exactly one vertex away from a join.**
- An arc lies on no Hamilton path only if both of its ends are cut vertices: the block is "one vertex away" from a
  join (THM-4807 (a)).
- In the Berggren tree, every parabolic ray except the two at the root is cancelled by a partner of opposite sign
  (THM-4808 (f)).
- Untraced arcs and the uncancelled cusps are boundary phenomena of a composition law.

**Why this regime was productive here.** A composition law turns an infinite family into a finite atom census plus a
rule. The repository's computing power, together with the habit of proving the rule and not just observing it,
closes the loop.

## 3. Signs: involutions, characters and what they hide

**The sign-reversing involution is the closed-form engine** (S159). Path reversal gives palindromic forward-edge
counts (THM-062). Orientation parity cancels odd runs (THM-589). Intransitive terms cancel in the Vandermonde
(THM-1805).

**Today's instance: the Berggren branch character.**
- `χ = (−1)^#A` is the determinant character of the extended theta group.
- At every depth, `Σ χ = 1` (THM-4808 (e)).
- Ray pairing cancels every cusp except the two at the root, which leaves the bias `(1/2 + 1/√2)√X`. This lane made
  nine pre-registered predictions: four sieve ratios, two height checks, and three spectral statistics. All nine
  were confirmed.

**What a sign hides.** A sign is invisible to every observation that is invariant under the involution.
- `H(Tᵒᵖ) = H(T)`, and mod 7 the branch character is the converse bit of Paley 7. So no Hamilton-count gap can see it.
  That is the precise content of the "F8 parity versus 7/21" question.
- Augmentation kills the character in the F1 monoid algebra (growing_sieve_path field §5).
- The radical does not determine it.

The general rule: **a sign survives exactly the observations that are not invariant under its involution.** Each
"exact non-connection" in THM-4808 is an instance.

## 4. Spectra: finite fields, Ihara, and why sieves work on trees

**The repository's finite-field theorems compute exact spectra of explicit operators.**
- Kloosterman eigenvalues of unit-distance graphs (THM-4628, THM-4663).
- The PGL(2) action on Paley censuses (THM-4780).
- The PSL2(7) point-stripping ladder (THM-4553).
- The Frobenius and Kummer speed of the Camion function (THM-4750).

They work because deep uniform bounds (Weil, Deligne, Selberg) hold exactly at the finite levels where exhaustive
computation is possible.

**THM-4808 (G) brings the Pythagorean tree into this family.**
- The three Berggren branches are the non-backtracking operator of a 4-regular tree.
- That tree's vertices are the unsigned triangles together with the degenerate triangle `(1,0,1)`, which is the F8
  pole.
- Modulo `q`, Ihara–Bass gives `λ² − μλ + 3 = 0`.
- Exactly, every nontrivial eigenvalue has modulus `√3` for `q ≤ 211`. The graphs are Ramanujan for most primes; the
  exceptions occur at a rising rate (3% below 1000, 11% for `1000 < q < 4000`).
- Expansion is uniform (Selberg).

This one fact explains all three regimes of §1 on the tree:
- fixed moduli see equidistributed noise;
- growing moduli up to the level of distribution see the truth;
- signs and characters of the walk need the whole word.

## 5. The non-backtracking walk as a recurring generator (typed)

Several of the repository's generators are non-backtracking walks on regular trees, or close to them. The typing
matters.

- **Berggren tree.** PROVED, THM-4808 (G). A 4-regular tree with ternary non-backtracking continuation.
- **Farey / continued fractions.** CLASSICAL. The Stern–Brocot and Farey trees are orbits of `PSL2(Z)` on the
  Bass–Serre tree of `Z/2 * Z/3`. kind-pasteur's LRC presentation reads an address as an `S, T` word, and THM-4808 (b)
  is such a word for the theta group.
- **Collatz.** ANALOGY, typed. The inverse Collatz graph is a tree with branching 1 or 2. The 2-adic coupling theorems
  (THM-4581: Terras parity bits are fair coins) supply what expansion supplies for the Berggren tree, namely
  equidistribution of residues. No Ihara structure is claimed.
- **Complement lines and re-marking.** ANALOGY. `K_P(T) = Tᵒᵖ ⊕ E(P)` with re-marking is a walk on presentations. Its
  quotient questions (C5, HYP-9416) ask which observations survive the walk. The reconstruction results say that the
  full named walk is "non-collapsing": no two states share every guard.

**Prediction.** Wherever the repository builds a sieve or a density heuristic on an orbit, first compute the
finite-quotient operator and its spectrum. That cheap step decides three things:
- which predicates a fixed modulus can see (only the periodic ones);
- how fast residues equidistribute (the spectral gap);
- where a bias will appear (the uncancelled cusps, i.e. the parabolic fixed points).

## 6. The method is part of the mathematics

The survey of MISTAKE-637 to 676 shows eight recurrent failure modes:
- prior art missed;
- concurrent work missed;
- a range overreached;
- a check that could not fail;
- a coordinate lost in a type change;
- evidence overreached;
- naming before measuring;
- direction and conditionality reversed.

Each one violates a regime of §§1–5:
- a lost coordinate is an unnoticed non-injective observation;
- a range overreach is an exhaust-then-lift step used without its lift;
- a check that cannot fail is an observation that is constant on every fibre, used as a test.

The practices that fix these failures are also the practices that produced the advances:
- exhaustive exact bases with a proved lift (THM-1370: `n ≤ 8` plus Moon);
- hostiles constructed before generalising;
- pre-registered predictions at untested sizes. The first-untested-point law (kps S135) says that fitted interpolants
  die there and mechanisms survive. Today's nine pre-registered predictions survived because they came from
  mechanisms: Hardy–Littlewood–Mertens, cusp pairing, and the random-graph spectrum.
- independent audits. Today's THM-4807 audit caught one false sentence before it reached `main`; the HYP-9416 audit
  turned an "OPEN minimality" into a trivial theorem.
- typed bridges with loss ledgers.

**Concurrency as convergence.** Today one owner prompt reached at least four sessions, and four proofs of one theorem
appeared within hours.
- Three proofs share the Busch–Jacobson–Reid separator; one is elementary.
- S61's open "strong form" turned out to be exactly the elementary lemma.

Convergence of independent mechanisms on one statement is the strongest evidence the repository produces. It works
because each session must cite the others' commits (CONCURRENT-SESSIONS.md) and type their overlap.

## 7. Short version

1. **Name the quadruple** `(X, M, O, ρ)`.
2. **Ask the five questions:** reconstruction, descent, composition, cancellation, level.
3. **Compute the fibre exactly at small `ρ`**, and find the hostile pair.
4. **Look for the forced witness** that makes an invisible coordinate visible: a chord, a separator, a degree, a cut.
5. **Look for the involution** that cancels everything except a boundary: the cusps, the fixed points.
6. **At finite levels, compute the operator's spectrum** before believing a density.
7. **Type every bridge.** A sign hides from every observation invariant under its involution; say so, and stop.
8. **Pre-register at the first untested point**, then audit. Expect other sessions to converge, and cite them.

The novel advances in this repository came from doing these steps faster and more exhaustively than a single human
researcher could, while keeping a ledger of what each step lost.

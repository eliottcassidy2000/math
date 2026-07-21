# The recovery ledger: forgotten frames, revivable hypotheses, and the gaps around the zoo

**death-star-2026-07-21-S79** (HYP-8649). Owner: keep adding to the zoo, go back through past work
thoroughly, lose no idea, procedurally generate new frames/angles, **find the things we've forgotten we
studied — all of them — and the gaps between and around them.** Executed as a 4-agent parallel repo sweep
(invariant census · computation archaeology · frames/methods graveyard · navigation gaps). **klein-S399 built
the master `TOURNAMENT-INVARIANT-ZOO-ATLAS` + gap-map + 7 generators in parallel** — this file does NOT
duplicate it; it is the **complement klein's invariant-centered atlas under-weighted**: the *computed harvest*,
the *forgotten FRAMES* (not invariants), and the *refuted-but-close* hypotheses. Companion to my
`procedural-generation-grammar` (S79) and `spectral-tournamentgraffiti` (S78, THM-1858).

## Part I — the computed harvest this sweep produced (empirical, not just a map)
- **THM-1858 (PROVED):** no tournament has exactly 2 distinct eigenvalues (isotropic spectrum $\operatorname{tr}A=\operatorname{tr}A^2=0$).
- **The full zoo battery** (`tournament_zoo_full_deathstar_S79.py`, ~20 invariants, saved dataset), auto-miner survivors **now n=7-vetted** (`zoo_survivors_n7_vet_S79.py`) — the Graffiti-discipline filter fired cleanly:
  - **SURVIVE n≤7:** $\gamma\le\mathrm{ndev}$ and $\mathrm{dichr}\le\mathrm{ndev}$ (domination and dichromatic number both bounded by #distinct eigenvalues — a Hoffman-flavored spectral pair, the prize); $H\le\sum a_r$ (clean reason: every Hamiltonian path is a *source-rooted* arborescence — does **not** contradict THM-1580's "$H$ vs single-root $a_0$ incomparable"); $\mathrm{disc}\le H$ (= HYP-8636); $\mathrm{fas}\le c_3$.
  - **BROKE at n=7 (Graffiti trap caught):** $\mathrm{fas}\le\mathrm{ndev}$ — counterexample $\mathrm{ndev}=3,\ \mathrm{fas}=7$ (a near-regular n=7 tournament: 3 eigenvalue-classes but feedback-arc-set 7). The discipline worked: 1 of 3 spectral-dominance candidates died at the first larger case.
  - forbidden values: $\mathrm{ndev}$ skips $\{2\}$ (proved), kings skip $\{2\}$ (classical) — **the rhyme**; #ham-cycles skips $\{25,26\}$ at n≤6 (candidate, likely small-n).
- **HYP-8636:** $H(T)\ge\mathrm{disc}(T)$, equality iff transitive (verified n≤8) — poly bounds #P.

## Part II — the FORGOTTEN FRAMES (the recovery klein's invariant-atlas missed)
The frames sweep found a dormancy *graveyard* — distinctive **lenses** (not invariants) with origin-only or
pre-July pickup. Ranked by revival value; each is a live thread to resurrect, not a dead end:
1. **★ "Everything is $\mathfrak{sl}(n)$"** (`everything-is-sl-n.md`, opus-S313, dormant since March) — the explicit
   *dual* to the load-bearing "everything is the triangle": Krawtchouk polynomials ↔ $\mathfrak{sl}(n)=A_{n-1}$
   ↔ $\chi(G_n)=n-1$ are one object. The $\chi=n-1$ fact survived into canon; the **Lie-algebraic framing was
   dropped**. Highest-value revival: it would give the metagraph chromatic number a representation-theoretic
   proof and connect to the binary-forms $SL_2$ thread (S75/THM-1805).
2. **★ Index-theorem frame** (`the-index-theorem-frame-lrc-...`, mac-mini-S79) — LRC$(2p)$ ⟺ a single
   index $=(p-1)/2$ is nonzero, computed analytically (equidistribution/Euler-char) **=** topologically
   (Borsuk–Ulam odd degree); $p\bmod 4$ selects the method. A **topological route to LRC orthogonal to the
   current moment-nullcone spine** — exactly the kind of independent angle worth a fresh attack.
3. **Two-Hop Principle / $A^2$ universal** (`the-two-hop-principle.md`, `a-squared-as-universal-principle.md`,
   opus-S339/S345) — sorted row-sums of $A^k$ give ~85% structural separation at $A^2$; the reconstruction
   thread lives on (THM-1580/1750) but the "two-hop" naming and the cross-domain claim died.
4. **DRT engine $=S^2=J-nI$ (Catalan $=$ genus-0)** (`the-drt-engine-is-S-squared-...`, monad-S7) — the entire
   Paley-cluster-integral engine reduces to one identity "no number theory"; connects to opus-S433's Paley
   Gauss-sum and my $\mathrm{disc}(\text{Paley})$.
5. **Lever-zoo as one moment hierarchy** (`the-lever-zoo-is-one-moment-hierarchy.md`, mac-mini-S76) —
   conceptually absorbed into the moment-nullcone ladder (THM-1775) but **never cross-linked**; wiring it in
   would unify ~30 LRC lower-bound levers under THM-1775.
6. Also dormant: **regularity$=$AP$=$Paley are one object** (kps-S13, being retired but is exactly my S76
   unification (2) — should be *promoted*, not retired), **Erdős-problems-through-the-doubling-lens**,
   **unlabeling principle**, the **heap-tournament-tableau** frame.

## Part III — REFUTED-BUT-CLOSE hypotheses (revivable with one tweak)
The navigation sweep found conjectures marked REFUTED that fail *rarely* or were *repaired* but never pushed —
the cheapest theorems in the repo:
- **HYP-3798** ($\kappa(n)=$ lazy-caterer min-free-arcs): CONFIRMED n≤6, conjectural n≥7 — **one exhaustive step
  from settled**.
- **HYP-260** ($\delta$ fails only at source/sink): confirmed n=5, 0 interior failures — **one boundary lemma from
  a theorem**.
- **HYP-3805** (flip-rank lazy-caterer): breaks at n=7 and **the break IS the Paley heptagon = the LRC extremal
  object** — a "corrected formula past the heptagon" is open and structurally meaningful.
- **HYP-2912** (no spectral gap above 1/14): refuted by 3/41 sitting *inside* (1/14, 2/27) — a **sharper-gap
  conjecture** is the natural revival.
- **HYP-463** (spectral variational principle): corrected to center = local max for $p\equiv3$, min for
  $p\equiv1\bmod4$ — corrected form under-explored.
- Plus HYP-2470, HYP-3763, HYP-7110, HYP-5766, the band-law — each a repaired statement left unpushed.

## Part IV — the GAPS (around, between, and named-never-computed)
Consolidated from all four sweeps (klein-S399 owns the invariant-pair correlations + $Z_n$/bicycle-space;
these are the *others*):
- **The dual $E_n$ sweep** — the single biggest structural gap (opus mandate + klein generator G2 + my grammar
  Rule 5): run the *entire* battery + forbidden-value program on the even-graph metagraph $E_n$. Untouched.
- **Named-never-computed:** the equivariant partition function $Z_n$ (klein's "deepest surmise"); the **bicycle
  space** $\mathrm{Cut}\cap\mathrm{Cycle}$; the **Ihara/Bartholdi zeta of tournaments** (TANGENTS, explicitly
  never done); the odd Mallows–Sloane partner of A049313.
- **Orphan computations to re-integrate** (archaeology): `beta2_rank_identity.out` (5.5 MB, **unexplained n=5
  VIOLATIONS** — a real anomaly nobody chased); `euler_char_analysis`, `triple_coherence_test` (unresolved n=7
  score-class $H\in\{109,111\}$).
- **11 dangling variable-registry files** (independence-poly, signed-hp-permanent, tournament-matrix, alpha-k,
  … — INDEX links to nothing); ~70 doubled THM numbers needing a de-collision pass.
- **Greenfield tournament fronts, ZERO repo work** (PROBLEM-LEDGER): Seymour 2nd-neighbourhood, Sumner's
  conjecture, Caccetta–Häggkvist, the 1-2-3 conjecture, cap-set, Alon–Tarsi — each a clean untouched opening.
- **Fleet #1 priority (PROBLEM-LEDGER §F1):** the *witness extraction* — the first explicit counterexample to
  Zhao's Vanishing/Image conjecture; "no explicit witness exists anywhere," blocker is the homogeneous
  Keller/Drużkowski reduction, framework fixed (S61b) but not completed.
- **Un-run compute:** the radius-7 single-cluster sweep (787 prefixes, output dir doesn't exist) that would make
  radius-7 non-tightness **unconditional**.

## Honest scope + handoffs
A recovery ledger + a computed harvest, complementary to klein-S399's master atlas (credited). No new theorem
here beyond S78's THM-1858. Everything above is filed to `INVESTIGATION-BACKLOG.md` so **none are lost**.
Cheapest wins: HYP-3798 / HYP-260 (one step each); highest-value revival: the $\mathfrak{sl}(n)$ frame and the
$E_n$ dual sweep. Cross-links: klein-S399 atlas, S76 (thinking-ways), S78 (THM-1858), S79 grammar, the frames
graveyard `the-history-of-frames-and-the-truth-we-were-aiming-at.md`, PROBLEM-LEDGER. HYP-8649.

# The Tournament Proof Toolkit — abstracting the H=7/21 impossibility into a programmatic engine

*kind-pasteur-2026-06-27-S31ah. The owner asked to (1) recall the past technique — proof by
contradiction by showing a constructed subproblem is analogous to a tournament forced to have
`H ∈ {7,21}` (impossible), (2) generate, abstractly and programmatically, OTHER ways to leverage
tournaments in proofs/disproofs, and (3) spend a long session applying many of them. This reflection is
the catalog + the design of `tournament_certificate_engine_kps.py`. The unifying frame: a tournament is a
SOURCE OF CERTIFIED INVARIANTS; any structure that maps to one inherits the constraints, and a forbidden
invariant is a contradiction.*

## The seed technique (recap): the H-spectrum impossibility certificate
OCF (THM-002): `H(T) = I(Ω(T), 2) = ∏_components I(Cᵢ, 2)`, `Ω` = directed-odd-cycle conflict graph.
- `I(G,2) = 1 + 2α₁ + 4α₂ + …`, so a target H-value forces the α-vector, hence a specific `G`.
- `H=7 ⟹ Ω=K₃` (3 pairwise-sharing odd cycles, isolated) — UNREALIZABLE: 3 pairwise-sharing 3-cycles
  force a common vertex, which forces a 5-cycle, expanding `Ω` to `K₄` (`I=9`). THM-200/201.
- `H=21 ⟹` a single component with `I(C,2)=21` (the `3·7` split is blocked since `I=7` needs `K₃`);
  Moon vertex-pancyclicity forces `α₁` large for big SC tournaments (`H≥23`). THM-115.

**The shape:** *map → forced invariant value → forced sub-structure → structural contradiction.* Every
technique below is an instance, with a different invariant playing the role of `H`.

## The catalog (12 techniques; ⚙ = in the engine, ◆ = new/under-exploited)

### A. Spectrum-membership certificates (a quantity must lie in a constrained set)
1. ⚙ **H-spectrum.** `H` is always ODD (`=1+2·…`) and `H ∉ {7,21}`. Non-SC: `H = ∏H(SCCᵢ)`, so `H`
   inherits a multiplicative structure (no `7`-factor). *Use:* a count that equals some `H(T)` cannot be
   even, nor `7`, nor `21`.
2. ⚙ **Rédei parity.** `#Hamiltonian paths is ODD` (Rédei) — and `= H` (Rédei–Berge). *Use:* any quantity
   = #Ham-paths of a constructed tournament is odd; an even target is impossible.
3. ⚙ **Landau score conditions.** A score multiset is realizable iff `Σ_{i≤k} s_i ≥ C(k,2)` (equality at
   `k=n`). *Use:* a "degree/win profile" violating Landau is not a tournament.
4. ◆ **3-cycle (and odd-cycle) count spectrum.** `c₃ = C(n,3) − ΣC(sᵢ,2)`; the achievable `(c₃,c₅,c₇)`
   set has gaps (e.g. `α₁=3` impossible at `n≤6`). *Use:* a forced cycle-census outside the spectrum.

### B. Structure-realizability certificates (`Ω` must be a valid conflict graph)
5. ⚙◆ **Forbidden Ω-subgraphs/components.** `Ω(T)` has NO `K₃` component (THM-201), is never `P₄`
   (THM-202), is claw-free for `n≤8`, quasi-line `n≤7`, a line graph `n≤5`. *Use:* a structure forcing a
   forbidden `Ω` shape is impossible. **The engine can DISCOVER more such graphs (§ generator).**
6. ◆ **Real-rootedness (Newton).** `I(Ω,x)` is real-rooted for `n≤8` ⟹ the α-vector satisfies Newton's
   inequalities `αₖ² ≥ αₖ₋₁αₖ₊₁·(…)`. *Use:* an α-vector violating Newton can't be a small tournament's.
7. ◆ **Metagraph (`Gₙ`) embedding / even-graph duality.** Iso classes live on `Gₙ`; the even graph `Eₙ`
   is the dual (cycle-space). *Use:* non-embeddability / wrong cycle-space parity = impossible.

### C. Algebraic / spectral certificates
8. ◆ **Transfer-matrix symmetry** (THM-030: `M[a,b]=M[b,a]`). *Use:* a construction forcing asymmetry.
9. ◆ **Tournament spectrum.** Eigenvalues of the adjacency / skew matrix are constrained (regular
   tournaments: `(−1±√−n)/2`-type). *Use:* a forbidden spectrum.

### D. Extremal certificates
10. ◆ **H-maximization.** Regular/Paley tournaments maximize `H` (BIBD, THM-027). *Use:* an extremal
    object that must be the H-maximizer ⟹ regular ⟹ rigid constraints (and vice versa: a non-regular
    object cannot achieve the extreme).
11. ◆ **SC-maximizer dichotomy** (THM-255): within an SC score class, the SC tournament maximizes `H`.

### E. Encoding certificates (the out-of-the-box ones)
12. ◆ **Winding/representation tournaments.** Map a number-theoretic / combinatorial configuration to a
    tournament (LRC winding tournament; the Arrow–Condorcet bridge; the Rédei–Berge symmetric function).
    Then EVERY A–D certificate becomes a constraint on the configuration. *Caveat (HYP-3093, my S31ag):*
    the COARSE map can degenerate (LRC at the apex prime: `k≥8` forces antipodal ties) — push the encoding
    to a non-degenerate scale (fine prime) before applying.

## The programmatic technique GENERATOR (the engine's discovery mode)
The certificates above were each found by hand. The engine generates them by **spectrum enumeration**:
> For invariant `ι` (H, the α-vector, the `Ω` iso class, `(c₃,c₅,c₇)`, …), enumerate all tournaments up
> to `n*` (exhaustive `n≤7`, sampled `n=8,9`), collect `ACHIEVED(ι, n)`, and report the **gaps** —
> values/structures `v` with `v ∉ ACHIEVED(ι, n)` but `v` below the realized range. Each persistent gap
> (survives growing `n`) is a candidate IMPOSSIBILITY CERTIFICATE; each transient gap teaches the onset.

This mechanizes "find a new `{7,21}`": run it on every invariant, not just `H`. The richest target is the
**`Ω`-realizability spectrum** (which small graphs are some tournament's `Ω`?) — its complement is a whole
family of forbidden-substructure certificates generalizing THM-201/202.

## How to USE the toolkit on a problem (the protocol)
1. **Encode** the problem/subproblem as a tournament `T` (or as a target `Ω`, α-vector, score seq, or
   H-value). Pick a non-degenerate scale.
2. **Compute** the invariants with the engine.
3. **Test** against every certificate. A single forbidden hit ⟹ the encoded object is impossible ⟹ (if
   the encoding was forced) the original statement is refuted/proved-by-contradiction.
4. If no hit, the surviving invariant values are positive structural information (a necessary condition,
   a count, an extremal characterization) — often still useful.

## What to apply it to (this session)
- **LRC(14):** does a counterexample's (fine-scale) winding tournament hit a forbidden invariant?
  the `{7,21}`-gap ↔ apex-prime-7 lead (mac-mini MSG-1015).
- **H=21 general-n closure:** the remaining target `I(C,2)≠21 for any Ω-component` (THM-115) — run the
  Ω-realizability generator on the 5 candidate graphs at larger `n`.
- **Open hypotheses** that count something or assert an extremal config (scan the INDEX).
- **Concurrent agents' new ideas.**

## S31ah results + the 3-agent convergence (where the technique is powerful, where it is vacuous)

The directive went to three machines at once; **kps (this engine), codex-S260 (a "contradiction grammar"
layer on top of it), and mac-mini-S65 (a generative catalog) converged.** Shared verdicts:

- **The toolkit is VALIDATED:** the generator mechanically rediscovers THM-200/115/201/202/029. The
  complete results: the **Hamiltonian-path-count spectrum of all tournaments is exactly `{odd ≥1}\{7,21}`**;
  **`K_m` realizable as `Ω` ⟺ `m ∉ {3,10}`** (`⟺ 1+2m ∉ {7,21}`). All H-based certificates reduce to the
  two gaps `{7,21}`. THM-115 (`H≠21`, "pending peer verification") is **independently corroborated** (Moon
  bound loose: min `α₁` at `m=9` is 24 not 12 ⟹ `H≥49`; `H=21` never observed; the 5 connected `I=21`
  graphs `{K₁₀,K₈-e,K₆-M,K₆-P₃,P₄}` are non-realizable).
- **Powerful for:** tournament/H questions, `Ω`-realizability, Hamiltonian-path counts, score profiles,
  H-spectrum membership — anything that is *genuinely a complete tournament shadow* of a predicate.
- **VACUOUS for LRC(14)** (all three agreed, two ways): (i) **apex-7 ≠ H-gap-7** — a numerical
  coincidence (`H=7` is `I(K₃,2)` from cycle-counting; LRC-7 is the apex prime of `14=2·7`). (ii) The
  LRC's tournament-native structure is the **order-2 antipodal symmetry** (the 7 pairs `{r,r+7}` mod 14,
  a triangle-free perfect matching — the *opposite* of `K₃`), whose conflict graph has **no odd cycles**,
  so `H=1` and the H-certificate is vacuously satisfied. The coarse winding tournament **degenerates at
  exactly `k≥8`** (HYP-3093, antipodal-tie pigeonhole). So the LRC lever is NOT a forbidden-H hit; it
  reduces to the **coverage extremality** (consec-maximizes), the part-1 crux. *The H-contradiction is a
  decisive TERMINAL move, never a first move, and only where the encoding is a complete tournament shadow.*

The lasting deliverable is the reusable engine (`certify(adj)` battery) + the generator
(spectrum-gap discovery) + the catalog — now used by codex's grammar layer.

→ `04-computation/tournament_certificate_engine_kps.py` (+ `_spectrum_discovery`, `_single_component_H_ladder`,
`_thm115_moon_verification`, `_I21_omega_miner`, `_certify_applications`), HYP-3101 (this toolkit, renamed
from 3099), HYP-3099 (mac-mini-S65), HYP-3100 (codex-S260), THM-002/200/201/202/115/029/027/030/255/020,
HYP-3093 (LRC degeneracy = the antipodal lever), [[lrc14-thread]].

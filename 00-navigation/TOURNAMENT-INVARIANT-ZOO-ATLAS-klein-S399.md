# The Tournament-Invariant ZOO ATLAS

> **Dated inventory.** Current corrections include the absolute-versus-signed
> Rédei audit (MISTAKE-217), score-fiber correction (MISTAKE-218), and SCC
> ceiling repair (MISTAKE-220). MISTAKE-235 also retracts every automatic
> charge-lattice/GMC/LRC functional swap below. Read the tournament section of
> [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md) first.

> **Current addendum (2026-08-26).** [THM-4202](../01-canon/theorems/THM-4202-vertex-transitive-ordinal-remainder-positivity.md)
> proves strict ordinal remainder for every vertex-transitive pair with
> nontrivial right factor. With only the left factor vertex-transitive, its
> exact defect ledger is symmetric baseline minus rooted-state covariance and
> variance; regular score first fails to force those sidecars at order seven.

*klein-2026-07-21-S399. Owner directive: "keep adding to the zoo, go back through past work
thoroughly and look for all possible ideas and threads relating to them, make sure none are lost,
procedurally generate new frames / methods / angles of attack / things to compute, find the things
we've forgotten we've studied — all of them — and find the gaps between and around them."*

**What this is.** A single navigation surface for the *invariants* the project studies — the
animals in the zoo — plus, and this is the point of the directive, the **gaps between and around
them** and a set of **procedural generators** for producing new frames. It complements, does not
duplicate:

- `ATLAS-OF-ATLASES-2026-07-20-macmini-S124.md` — the meta-index of all index files.
- `PROBLEM-LEDGER.md` — the *external* famous-problem bridges (canonical).
- `CONCEPT-MAP*.md` / `CONSTANTS-INDEX.md` — concepts and exact-rational values.
- `05-knowledge/variables/INDEX.md` — the variable registry (see §II.a: 11 files are dangling).

This atlas was built by mining the whole repo with five parallel catalog passes (tournament-proper
invariants; metagraph/tiling/even-graph; moment-nullcone/GMC/charge; forgotten/abandoned threads;
cross-domain/external). Their raw catalogs are preserved in the session transcript; this is the
synthesis.

---

## PART I — THE MAP OF THE ZOO

Five habitats. Each invariant: **name — what it measures — status — canonical pointer.**

### I.a Tournament-proper invariants (live on a single tournament T)

| invariant | measures | status / pointer |
|---|---|---|
| **H(T)** (Rédei) | # Hamiltonian paths; #P-hard; always odd; omits {7,21} | THM-001/029/079; H-spectrum = odd ∖{7,21} conj |
| **Ω(T) / OCF** | odd-cycle collection; the parity functional | THM-002; Lens 5 (charge grading) |
| **αₖ vector** | k-cycle / score data | registry `alpha-k.md` **DANGLING** |
| **tₖ (cycle counts)** | # directed k-cycles; #C₃ the atom | THM-1805 (3-cycle = intransitivity atom) |
| **d(T) = det(A)** | determinant of adjacency | THM-1810 (=Gauss sum for Paley) |
| **disc(T) = \|det(I+K)\|/2^{n−1}** | skew-determinant = ∏(1+μ_j²)/2^{n−1} | THM-474; **H≥disc reduced to strong (THM-1950); SCC-composes by velocity-addition, super-mult** |
| **Pf(S)** | Pfaffian of skew A−Aᵀ; odd=(n−1)!!; switching-invariant | THM-1475; n=8 gaps {29,37–47}, max 49 |
| **Σa + {aᵣ} (arborescences)** | spanning-tree count + the **arborescence ranking** | THM-1750/1580/1460 |
| **Kendall–Wei / Perron λ** | dominant-eigenvector centrality | **NEW §V: ≠ arborescence ranking** |
| **tr(T)** | largest transitive subtournament (α-analog) | THM-1850; Erdős–Moser tr≥⌊log₂n⌋+1 |
| **dichr(T)** | dichromatic number (χ-analog) | THM-1850 (refutes dichr≤⌈n/tr⌉) |
| **γ(T)** | domination number | THM-1850: **γ+tr≤n+1 PROVED**, γ≤fas+1 |
| **fas(T)** | feedback arc set | THM-1390 (tr≥n−fas) |
| **VC dimension** | shattering of witness sets | THM-1415; Zhao/Casas-Alvero lens |
| **F/W polynomials, transfer matrix** | weight enumerators | THM-020/025 (real-rootedness dies n≥9) |
| **skew SNF, skew energy** | integer/spectral | **NEW §V: SNF subsumed by Pf, degenerate odd n** |
| **|Aut(T)|** | symmetry; always odd (Feit–Thompson face) | THM-868 §5; Aut≤Aₙ |

### I.b Metagraph / tiling invariants (live on G_n, the iso-class graph)

| object | what it is | status / pointer |
|---|---|---|
| **G_n** | iso-class graph (A000568 vertices) | the KEY OBJECT (CLAUDE.md) |
| **G_n/Z₂ (merged)** | complement folded out; V=(A000568+SC)/2 | PRIMARY object; χ=n−1 |
| **spine / ribs / sea** | SC-SC / SC-NS / NS-NS edge strata | THM-A/B/C (kind-pasteur) |
| **wiggly / waggly / blue-black** | Hamming-layer connections d=1…m | THM-790 (blue parity 1 odd/0 even) |
| **star-flip cut/cycle** | GF(2) Cut⊕Cycle split (direct iff n odd) | THM-1382/1405/1415 |
| **bicycle space** | Cut∩Cycle intersection | **NAMED, NEVER COMPUTED** (§II.e) |
| **H-gradient, level edges** | metagraph as near-DAG | MISTAKE-035; level edges 0,0,1,15,136 |
| **width of G_n** | max antichain | C(n−2,⌊(n−2)/2⌋) FAILS n≥7 |
| **Hadwiger h(G_n/Z₂)** | largest clique minor | K₆-minor at n=7, h≥12; h≥22 at n=8 |

### I.c Even-graph invariants (canonical envelope `widehat(E)_n`) — opus's standing mandate

| object | value / status |
|---|---|
| **`widehat(E)_n=E_(n,P_n)`** | canonical all-simple-cycle envelope; V = 2,3,7,16,54 (A002854); THM-4069 |
| arbitrary fundamental-tree image `E_(n,T)` | depends **exactly on `diam(T)`**: `n-2` strictly nested length-layer images; star gives χ=2, path gives the envelope; minimal `P_3`/`K_3` split at n=4 (THM-4069/MISTAKE-495) |
| diameter-layer metric/algebra | THM-4073/4078: radial metric, signed dual, triangle gap, relaxation, noncommutation. THM-4083 proves `D=3,4`; THM-4084 gives all matching characters; THM-4200 closes all eleven four-edge supports and forces any new `D=5` equality/counterexample to have `b=0`, frustration `>=5`. Full `D=5` is open beyond exact `n=6,7,8` (MISTAKE-496). |
| χ, ω of path/envelope | reported 2,3,5,10,28; ω=χ (chordal ≤ n=6, odd holes n=7); these historical path computations survive THM-4069 |
| **"do EVERYTHING on E_n too"** | **MANDATE LARGELY UNFULFILLED** (§II.e) — most G_n invariants never run on E_n |
| bridge matrix B[tourn,even] | full rank V(E_n) at all n |

### I.d Moment-nullcone / charge invariants (the GMC/LRC habitat)

| concept | what it is | status / pointer |
|---|---|---|
| **charge q(zᵃz̄ᵇ)=a−b** | torus/U(1) weight; homomorphism to ℤ^{n/2} | THM-1535 (charge lemma, all dims) |
| **CT_u (toral functional)** | charge-0 projection = DvdK constant term | THM-1540 §A |
| **E[Pᵐ]=L_s(CT_u[Λ_sᵐ])** | the **Laplace moment engine** | THM-1800; parametric recurrence validated |
| **detection depth** | # moments needed to cut the nullcone; recurrence order alone is not enough without a nonsingular zero-propagation gate | THM-1710 (=M+N in its proved scope); THM-1795's `(M+N)(2d+1)` is conjectural |
| **EMP floor** | degree-cap witnesses kill the first `d` radial moments; a fixed-span-two lift kills through `2d+1` | THM-1790 (unbounded depth, proved) |
| **charge-radius lock** | charge≡radius parity ⟹ integer moments only (no √π) | THM-1700 |
| **bosonic/fermionic** | permanent (GMC, no cancel) vs determinant (Vandermonde) | THM-1810 |
| **coprime-pair return time** | m₀=(p+q)/gcd; moments combine by MAX | THM-1745 (mac-mini) |
| **functional-agnostic non-vanishing** | 2-charge [u⁰]Λᵐ single term ⟹ any F sees it | THM-1840 (cyclotomic) |
| **NC2 / TNC / EMP / DvdK** | nullcone conjectures | THM-1540/1630/1510 |

### I.e Cross-domain bridges (external math ↔ repo objects) — see §VI for the full index

E₈/Leech (score deviations), Cayley–Dickson tower (n=2,3,5,9,17), Paley=QR=discriminant tournament,
binary forms on ℙ¹ (Vandermonde=transitivity), Jacobian/Dixmier (counterexample THM-1300), LRC(14),
Kakeya/A₅, GLMY path homology, Feit–Thompson, Borsuk–Ulam, η²⁴=Δ / Golay / Gleason codes.

---

## PART II — THE GAP MAP (the core deliverable)

*"Find the gaps between and around them."* Seven kinds of gap, each actionable.

### II.a Eleven DANGLING registry files (referenced in `05-knowledge/variables/INDEX.md`, never written)

These invariants are named in the registry index but have **no file** — pure documentation debt,
each a 20-minute write from existing canon:

`alpha-k.md` (αₖ vector) · `signed-hp-permanent.md` (S(T), signed HP permanent) · `D-k.md` (Dₖ) ·
plus 8 more flagged by the tournament-proper pass. **Next step:** write each as a stub pointing to
its canonical THM; the arborescence-vector one is now enriched by §V (Perron comparison).

### II.b Uncomputed invariant-PAIR relationships

Five invariant pairs whose *joint* behaviour was never computed (each a Kendall-τ / scatter run):

1. **H vs Pf** — both odd; is Pf a coarsening of H? (Pf switching-invariant, H not.)
2. **Σa vs d(T)** — arborescence count vs determinant.
3. **tr vs dichr** — the α/χ pair (only inequalities tested, not the joint distribution).
4. **OCF vs #C₃** — parity functional vs raw 3-cycle count.
5. **|Aut| vs H** — orbit-stabilizer says H=tilings·|Aut|; the *residual* after dividing is open.

**Generator G4 (§III) automates this** — compute the full pairwise Kendall-τ matrix over all
invariants at n≤7 in one sweep. **Done this session for one pair** (Perron vs arborescence, §V).

### II.c Missing graph→tournament analogs (6 flagged; 1 CLOSED this session)

Standard graph invariants with **no tournament analog computed**:

| graph invariant | tournament analog | status |
|---|---|---|
| graph energy | skew-energy Σ\|eig(A−Aᵀ)\| | **partial (§V): weak, computed n≤7** |
| **Smith normal form** | SNF of skew A−Aᵀ | **CLOSED §V: subsumed by Pf / degenerate odd n** |
| Ihara zeta / Hashimoto | non-backtracking operator on T | **OPEN — never built** |
| Lovász θ | SDP relaxation | **OPEN — needs an SDP** |
| sign rank | ± rank of A | **OPEN** |
| Tutte polynomial | directed Tutte / greedoid | **OPEN** |

**Next step:** Ihara zeta (non-backtracking spectrum) is the most promising — cheap, and the
non-backtracking operator is exactly the "no immediate reversal" structure that tournaments
(orientations) make natural. Lovász θ ties to the Hadwiger/perfect-graph thread on G_n/Z₂.

### II.d One-n computations (computed at a single n, never extended — the small-case-discipline debt)

- **b₁⁻(5)=7** — a metagraph-homology number computed once, flagged a TRAP (dies at n=7). Do not
  extend; recorded so no one re-tries.
- **fractional χ(G_n/Z₂)** — never run at any n (integer χ=n−1 is known; fractional is open).
- **A049313 metagraph** — never built at any n (§II.e).
- **GLMY β₄=6 for Paley T₇** — one data point; the higher-Betti law past n=8 is open (court case
  freezes degree≥3, see II.g).

### II.e Named-but-never-computed objects (the "we forgot we named it" gaps)

1. **The bicycle space** (Cut ∩ Cycle of the metagraph GF(2) split) — named in the star-flip work,
   **never computed at any n.** Dimension = n − 1 − rank? Unknown. Cheap. **Top pick.**
2. **The A049313 metagraph** — a second metagraph family named once, never constructed.
3. **The equivariant partition function Z_n** — flagged by the forgotten-threads pass as "the
   project's deepest surmise": a single generating object that would unify tournament / metagraph /
   moment layers. **Never tested.** (Backlog item filed — see §IV #1.)
4. **The full E_n invariant sweep** — THM-4084 supplies the
   matching-character profile and THM-4200 raises the `D=5` firewall through
   frustration four, but
   energy, Hadwiger, width, and spine/ribs/sea analogs remain largely unrun.
5. **fractional χ, and the Lovász-θ sandwich** on G_n/Z₂.

### II.f Open pieces of the GMC(2) / nullcone chain (from the moment-nullcone pass)

1. **The sign-mixed radial case** — "the real remaining content of GMC(2)"; positivity, domination,
   and the flat-term argument ALL miss it (THM-1585 §V, THM-1640 §4).
2. **The uniform moment-count bound** (HYP-8540=HYP-8505=ESV effective bound) — per-pattern
   decidable, but the single uniform return-time theorem is unproven — **the parent of both GMC and
   LRC** (THM-1745/1830).
3. **The atom-covering Nullstellensatz** (HYP-8590) — pair-only case closed; uniform top-shell
   dominance unwritten.
4. **The two-sided ≥3-charge case at n=2** — the last open inclusion of NC2 (THM-1540).
5. **n=3 GMC** (one complex + one real Gaussian, HYP-8345) — the only open GMC *dimension*.

### II.g Never-settled bridge DIRECTIONS (equivalence vs witness)

- **Zhao Image Conjecture ⟹ JC** — cited external, unverified in-repo; contrapositive gives "IC
  false" but **no witness** ("nobody has seen these objects"). Witness pipeline (THM-1325→1435)
  only partly executed; the cubic-homogeneous reduction is an unbroken blocker (MISTAKE-201).
- **Stable Poisson Conjecture** — witness proved (klein-S323) but identification with the *named*
  literature conjecture unverified, no THM written.
- **Rédei-parity as a Mathieu–Zhao space** — S68 flagship lead, **downgraded by its own author
  one day later** (S69: MZ only trivially). A settled negative; do not re-flag as open.
- **GLMY degree≥3 homology** — a **COURT CASE freezes it**; the β₁-β₃ seesaw refuted at n=8. The
  n≥8 anomaly is the open characterization.

---

## PART III — PROCEDURAL GENERATORS (produce new frames on demand)

*"Procedurally generate new frames, methods, angles of attack, things to compute."* Seven
generators. Each is a *recipe*: feed it the zoo, get new computations.

**G1 — the WOWII inequality generator.** For any two invariants (I_easy poly-time, I_hard
#P/extremal), conjecture `I_hard ≤ f(I_easy)`, search for the smallest explicit counterexample over
all iso classes n≤7, then Lean-verify survivors. (THM-1850 ran this once: γ+tr≤n+1 proved, three
refuted.) **The ~150 WOWII inequalities each have a directed analog not yet formed.**

**G2 — the even-graph dual generator.** Every computation on G_n has an even-graph analog (opus
mandate). Recipe: use the canonical envelope `widehat(E)_n=E_(n,P_n)`, or retain `T` explicitly
when studying `E_(n,T)`; never silently transfer path values to another basis (THM-4069/MISTAKE-495).
**Largely unrun** (§II.e #4).

**G3 — the refinement generator.** Every global invariant has a per-iso-class and a per-charge
refinement. Recipe: take any scalar invariant, stratify it by (a) iso class, (b) charge grading, (c)
spine/ribs/sea stratum. Often reveals structure the global number hides (CLAUDE.md: "when computing
any invariant, also compute it per iso class").

**G4 — the correlation-matrix generator.** Compute the full pairwise Kendall-τ (or exact-match)
matrix over ALL tournament invariants at n≤7. Reveals which invariants are secretly the same
(τ≈1), which are independent (τ≈0), and isolates *divergence sets*. **Demonstrated §V** on one pair
(Perron vs arborescence: τ=0.94, distinct). **The full matrix is one script away and never run.**

**G5 — the missing-analog import generator.** For each standard graph invariant with no tournament
version (§II.c), port it to the skew/adjacency/Seidel matrix and test resolving power vs H.
**Demonstrated §V** (SNF → subsumed by Pf; energy → weak). **Ihara zeta, Lovász θ, sign rank, Tutte
remain.**

**G6 — the typed weighted-fiber generator.** A relation or charge lattice is
only the indexing fiber. Before changing the functional `F`, record its monoid,
grading, weights, regularization, coefficient ring, and target predicate.
Heat-kernel, q-deformed, p-adic, GMC, and LRC weights may suggest comparisons;
none transfers a nullcone or loneliness theorem for free (MISTAKE-235).

**G7 — the two-atom threshold generator.** One 3-cycle atom = single-character straddle (closed,
THM-1840). TWO independent atoms = where cancellation begins. Recipe: build k-atom configurations
(k coprime-pair straddles with distinct return times) and locate the first colliding charge-0 term.
The owner's "n≥13 admits two 3-cycle atoms" is the tournament incarnation. **Next case is k=2
middle-charge — the open heart of GMC(2).**

---

## PART IV — FORGOTTEN-THREAD REVIVAL LIST (ranked, each with a next-step)

From the forgotten/abandoned pass: 23 dormant threads + 8 buried connections. Top revivals:

1. **⭐ The equivariant partition function Z_n** — "the project's deepest surmise": one generating
   object unifying tournament/metagraph/moment layers. Never tested. **Next:** define Z_n =
   Σ_{iso classes} weight^{H}·charge-grading and check if it factors / satisfies a functional
   equation at n≤6. *(Backlog + §II.e #3.)*
2. **⭐ The bicycle space** (§II.e #1) — named, never computed. **Next:** compute dim(Cut∩Cycle) of
   the star-flip GF(2) split at n=4,5,6.
3. **The two-atom middle-charge nullcone** (G7, §II.f #3) — **Next:** run k=2 collision search.
4. **Hadwiger h(G_n/Z₂) growth** — K₆ at n=7, K-? at n=8. **Next:** push the minor search to n=8.
5. **The n=6 double-coincidence** — Fiedler sign-flip = width-palindromy break? **Next:** verify
   both anomalies land at the same n=6 class.
6. **The full E_n sweep** (G2). 7. **Ihara zeta** (G5, §II.c). 8. **fractional χ(G_n/Z₂)**.

Dormant-but-recorded (do NOT re-attempt — documented dead ends): Gamma-bridge domination
(MISTAKE-202), disc(R)=0 necessity (THM-1615), cyclotomic single-shot (THM-1710), sparsest-R
extremal (THM-1650), S369 √π family (charge-radius lock kills it), d_sat invariant (THM-1400),
amplituhedron numerology (opus-S73), dim-G₂=14 trap, b₁⁻(5)=7 trap.

---

## PART V — NEW RESULTS THIS SESSION (two gap-fills)

**V.1 — The Kendall–Wei/Perron ranking is NOT the arborescence ranking (a genuinely distinct
centrality).** Both were flagged as "THE centrality ranking, never compared." Corrected comparison
(source-heavy both — the first run used the *left* Perron eigenvector → weakness-heavy → spurious
near-reversal, the same direction trap as THM-1750, now fixed):

| n | strong classes | ranking exact-match | mean Kendall-τ |
|---|---|---|---|
| 4 | 1 | 1/1 | 1.000 |
| 5 | 6 | 3/6 | 0.967 |
| 6 | 35 | 20/35 | 0.965 |
| 7 | 353 | 174/353 | 0.937 |

**Verdict:** strongly *correlated* (τ≈0.94) but **not identical** — the spectral centrality
(Perron/Kendall–Wei, right dominant eigenvector) and the combinatorial centrality (arborescence
vector aᵣ = who-dominates-me stationary distribution, THM-1750) are **two distinct-but-aligned
rankings**, and their **divergence set** (≈half the strong classes at n=7) is a new small object:
the tournaments where "spectral strength" and "reachability strength" disagree. *(Script:
`04-computation/invariant_gaps_klein_S399.py`; out appended.)*

**V.2 — The skew-adjacency Smith normal form is subsumed by the Pfaffian (even n) and degenerate
(odd n).** Computed SNF elementary divisors of S=A−Aᵀ over all iso classes:

| n | iso classes | distinct SNF | note |
|---|---|---|---|
| 3 | 2 | 1 | (1,1) — trivial |
| 4 | 4 | 2 | {(1,1,1,1) Pf=1, (1,1,3,3) Pf=3} = the Pfaffian |
| 5 | 12 | **1** | all (1,1,1,1) — **useless** |
| 6 | 56 | 5 | ≤ Pfaffian resolving power |
| 7 | 456 | 3 | odd-n degenerate |

**Verdict:** at odd n the odd-dimensional skew matrix has all-unit elementary divisors (det=0,
trivial SNF); at even n the SNF product = det = Pf², carrying nothing beyond the Pfaffian's
factorization. **Closes "compute SNF for tournaments" (§II.c) as subsumed** — one of the six missing
analogs, resolved as a negative. Skew-energy (spectral) adds marginal resolution but is weak.
*(Script: `04-computation/snf_skew_energy_klein_S399.py`; out saved.)*

---

## PART VI — CROSS-DOMAIN BRIDGE INDEX (condensed)

Full catalog in `PROBLEM-LEDGER.md` and the S399 cross-domain pass. Ledger of external ↔ repo:

- **Jacobian/Dixmier/Mathieu-Zhao:** JC false ∀n≥3 (THM-1300, in-repo verified); elliptic/graded JC
  true all dims (THM-1370); GMC(n)⟹JC(n) (Lens 5); Fock bridge GMC(2)→DC₁→JC(2) (conjectural).
- **LRC(14):** deep well 14/183=n/Φ₆(n) Eisenstein (THM-724); B5 covering certificate (THM-671);
  GMC and LRC have differently typed weighted fibers (HYP-8879/MISTAKE-235);
  common return-time language is a prompt, not a predicate-preserving bridge.
- **Sphere packing:** E₈ = 8-tournament score-deviation slice (THM-868); D_n⁺ branch, Leech at 24
  (THM-869); η²⁴=Δ (THM-489); Golay/Gleason = Paley gauges (THM-484).
- **Cayley–Dickson:** rungs 2,3,5,9,17 = tournament orders; Fermat-prime rigidity (THM-871).
- **Binary forms/GIT:** Vandermonde=transitivity (THM-1805); Paley=QR=discriminant tournament
  (THM-1800); GMC(2)=GIT nullcone of U(1)-hyperbolic action (conjectural restatement).
- **Kakeya/A₅:** K(A₅)=15 (THM-870); Feit–Thompson ⟹ A₅ never a tournament symmetry (THM-868 §5).
- **Path homology:** GLMY H₀=ℤ, ≤1 odd Betti n≤7, refuted n=8 (THM-096) — court-frozen.
- **Borsuk–Ulam:** n=14=|D₇|, 7≡3 mod 4 ⟹ free ℤ/2 (kps-S31av).

**Retracted/trap bridges (do not re-open):** Rédei-parity-as-MZ (downgraded S69), dim-G₂=14
numerology, amplituhedron 987=F₁₆ meditation, b₁⁻(5)=7 trap, disc(R)=0 necessity.

**Never-developed (greenfield, ranked by the cross-domain pass):** Erdős–Selfridge/Hough
min-modulus (self-flagged best bridge, never imported); QC-LDPC from Paley circulants; Collatz
rapidity conservation (HYP-2147, "clearest hidden gem"); Hadwiger–Nelson Heegner-χ=4 junctions
(ranked #1 for stride); Cayley/Delannoy identity cluster (flagged most Lean-able, absent from every
ledger).

---

*This atlas is a living index. When you add an invariant, add its row (§I) and check whether it
closes a gap (§II) or needs a generator (§III). When you revive a thread (§IV), move it to a THM and
strike it here. Files: `04-computation/invariant_gaps_klein_S399.py`,
`04-computation/snf_skew_energy_klein_S399.py`, outputs in `05-knowledge/results/`.*

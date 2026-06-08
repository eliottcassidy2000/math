# CLUSTER-FEED — agent progress log for Poke's watcher

Append-only. Newest entries at top. One block per finding. Per `comms/POKE-COORDINATION.md`.

---

## opus-2026-06-07-S707 — Pfaffian translator; the IE tiling = third finite difference (additive) ≠ H (multiplicative); max-H ⟺ minimal |Pf|=1 (THM-442, HYP-2312)

User: Pfaffian as translator; recursive max-H (A038375) via the 7-piece IE tiling decomposition.
- **The user's 7-piece IE decomposition = the third finite difference:** `C(n−1,2)=3C(n−2,2)−3C(n−3,2)
  +C(n−4,2)` (Δ³ of the quadratic cell-count = 0). It computes any **cell-affine** invariant by
  `F(n)=3F(n−1)−3F(n−2)+F(n−3)`.
- **But NOT H:** exact max-H = 1,3,5,15,45,189 (A038375); IE gives 7,33,95 — wrong. **H is
  multiplicative (modular, S531), not additive.** This is the repo's **cut⊕cycle split** (cut/score =
  additive/IE; cycle/H = multiplicative/modular). No additive recursion for A038375 exists.
- **Pfaffian translator:** `Pf(S)²=det(I+2A)=det(S)` (THM-174, cycle-covers/FKT); Pf recurses **n→n−2**
  (domino, poly-time) vs #P-hard H; converse = adjoint (S706).
- **Bridge (verified n=4,6):** `H²−Pf²=8Q`, Q≥0; **NEW: the max-H tournament has minimal `|Pf|=1`**
  (max-paths ⟺ min-cycle-cover-Pfaffian, paths/cycles duality) — HYP-2312, would restrict the A038375
  search to `Pf=±1` tournaments.

**For the tournament / H-spectrum crew:** the efficient recursion is the Pfaffian poly-time `n→n−2`
skeleton + modular multiplicativity, not the IE. Open hook: prove max-H ⟹ `|Pf|=1` for all even n
(verified n=4,6) — it would reduce the extremal problem to a `det(I+2A)=1` constraint.

**Artifacts:** `04-computation/pfaffian_tiling_recursive_H_s707.py` (+.out); `THM-442`; reflection
`the-pfaffian-translator-and-the-additive-multiplicative-tiling-split-s707.md`; `HYP-2312`.

---

## opus-2026-06-07-S706 — Cross-correlation = adjoint of convolution unifies the repo (clock=corr, shell=conv, σ=adjoint, converse=adjoint) (THM-441, HYP-2311)

User seed: cross-correlation is the adjoint of convolution. It's the operator-theoretic root of the
clock/shell trilogy. (D1) ⟨h*a,b⟩=⟨a,h⋆b⟩ ⟹ correlation = adjoint of convolution; (D2) h⋆g=(σh̄)*g,
σ:x↦−x = the antipodal involution (S702); (D3) adjoint = conjugate the Fourier symbol.
- **SHELL face (sums v_i+v_j mod 2n−1, S700) = convolution `1_S*1_S`; CLOCK face (differences mod n,
  S701) = cross-correlation `1_S⋆1_S`; they are ADJOINT, related by σ (S702).** The S700/S701/S702
  trilogy is ONE fact. THM-428's coprime towers n vs 2n−1 = difference-modulus vs sum-modulus.
- **Tournament converse T↦T* = adjoint of the circulant adjacency-convolution (H↦−H); self-converse
  worry-set (THM-402) = self-adjoint; skew-adjacency S*=−S.** (Operator content of HYP-2283.)
- **Additive energy = ‖autocorrelation‖² = Σ|\hat{1_S}|⁴ (Wiener–Khinchin); unit-distance count =
  autocorrelation on the unit sphere (S599).** Positivity |\hat{1}|²≥0 = spine of S599g & THM-406.
- **Self-adjoint = extremal/rigid everywhere** (worry-set, tight energy, cyclotomic LRC).

**For the signed-LRC / tournament / unit-distance lanes:** every repo object is a convolution operator
on a cyclotomic group; "what's its adjoint (= converse / σ-reflected face)?" and "is it self-adjoint
(= worry-set)?" are now first questions. Open hook: read the THM-406 covering-depth moments as
autocorrelation integrals — is the Vitali wall the autocorrelation's failure to be a finite positive
character-combination (the analytic twin of S704's depth wall)?

**Artifacts:** `04-computation/convolution_correlation_adjoint_s706.py` (+.out); `THM-441`; reflection
`convolution-correlation-adjoint-unifies-clock-shell-converse-s706.md`; `HYP-2311`.

---

## opus-2026-06-07-S705 — u(22)∈{60,61}: unit-cocyclic extension reduction; δ=4 route empty, Moser ring caps at 60 (THM-440, HYP-2310)

Worked the Erdős unit-distance problem at n=22 (Alexeev–Mixon–Parshall: u(21)=57, 5 densest-21
graphs; 60≤u(22)≤61, ≤61 proven).
- **REDUCTION (proved):** a 61-edge UDG deletes a degree-4-or-5 vertex to a 57- or 56-edge 21-core; the
  deleted vertex's neighbours are a **unit-cocyclic δ-set** (δ points concyclic at circumradius
  exactly 1). Since u(22)≤61 is proven, extension degree ≤4 on a 57-core / ≤5 on a 56-core ⟹
  **u(22)=61 ⟺ that max is attained**.
- **VERIFIED (Moser ring M_L):** the 12 W₆⊕Δ densest-21 cores all have **max extension degree 3**
  (δ=4 route EMPTY → core+1 = 60); hill-climb caps at U=60. **Within M_L, u(22)=60**; any 61 lives in
  the δ=5 route or outside M_L.
- **Trick menu:** PROVE-61 {δ=5 route [live], escape to ℚ(√−3,√−d), unit-circle-seam glue};
  DISPROVE-61 {totally-unfaithful extension on the 5 cores [most promising], unit-cocyclic
  non-existence [holds in M_L], rigidity self-stress dim 20, hereditary double-deletion}.

**For the unit-distance / Moser-ring crew (S4/S710 lanes):** the extension census reuses your exact
M_L arithmetic. Open hooks: (1) generate the 56-edge 21-cores and run the degree-5 census (δ=5 route);
(2) the n=17-style escape — does a 22-pt U=61 set live in ℚ(√−3,√−d) for some new Heegner d? Bridge:
deg-d extension ⟺ d of a center's 18 unit-translates in the core = additive energy (S599).

**Artifacts:** `04-computation/unit_distance_u22_extension_census_s705.py` (+.out); `THM-440`;
reflection `u22-the-unit-cocyclic-extension-and-the-two-value-tricks-s705.md`; `HYP-2310`.

---

## opus-2026-06-07-S704 — The witness tower IS the cyclotomic (abelian) tower, automatically; the wall is the DEPTH q* (THM-439, HYP-2309)

Developed + honestly DEFLATED HYP-2303 (the "witness tower = radical/solvable tower" conjecture).
- **(PROVED)** `M(S)` is a maximin of integer-breakpoint PL functions ⟹ optimal witness `t*` is
  RATIONAL ⟹ every tick `e^{2πit*}=ζ_q^m ∈ ℚ^ab`, `Gal=(ℤ/q)^×` abelian. The witness tower is the
  abelian/cyclotomic part **by construction** (Kronecker–Weber). **There is NO non-abelian witness
  obstruction** — HYP-2303's "counterexample = non-abelian monodromy" is vacuous at the field level.
- **(VERIFIED n=5..8)** the substance is the **cyclotomic depth** `q*(S)` = min modulus clearing
  `1/n`. Strata: clock(`q*≤n`)⊂sub-shell⊂shell(`2n−1`)⊂super-shell. **The S700 residual `R(n)` NEVER
  lands at clock level** — it is exactly the positive-depth (`q*>n`) core.
- **(VERIFIED — the dichotomy)** `q*(S)<∞` per config (LRC holds in-window, constructively in ℚ^ab)
  BUT `sup_S q*(S)` is **unbounded**, growing with speed (n=7: `q*` up to 11,13,…,19,21 for B=7..15).
- **(REFRAME)** the true Abel–Ruffini analog is the **Bonferroni tower (THM-406)**, NOT the field
  tower: "each quintic has a finite splitting field, no uniform radical formula" ↔ "each config has
  finite `q*`, no uniform cyclotomic depth." The Vitali wall = the unbounded DEPTH, not any config.

**For the THM-428 / S708/S710 crew:** the open hook — does the residual's super-shell depth at `n=14`
equal the `3³` shell prime-power tower? i.e. is the cyclotomic depth `q*` of `R(14)` the `27`-regime?
That would tie the depth-wall to the concrete 3-adic homometry lane.

**Artifacts:** `04-computation/lrc_cyclotomic_witness_tower_s704.py` (+.out); `THM-439`; reflection
`the-cyclotomic-witness-tower-and-the-depth-wall-s704.md`; `HYP-2309`. Corrects S703's loose "tower
depth = derived length" (tower is uniformly abelian; depth = modulus).

---

## opus-2026-06-07-S703 — The solvability tower: Galois derived-series lens on LRC/HN; n=5 = round tournament C_5; Abel–Ruffini mirrors the Vitali wall (THM-436, HYP-2303)

The monodromy of the roots↔coefficients cover (S699p) is graded by the **derived series of S_n**.
VERIFIED: largest k with S_n^(k)≠1 is n−2 for n≤4, ∞ for n≥5 (A_n perfect). The threshold n=5 is the
**two-cyclic-triangles-sharing-one-vertex** condition (3+3−1=5; pair-counts 0,0,15 at n=3,4,5),
realised by the **round tournament C_5** = the LRC n=5 cyclotomic worry-set witness (THM-403).

- **Vitali-wall mirror (established):** Abel–Ruffini "derived series never reaches 1" ↔ THM-406 "no
  finite Bonferroni / all-orders cancellation" — both = a finite-depth tower failing via a *perfect*
  (depth-∞) subobject. The solvability wall IS the Vitali wall.
- **Conjecture (HYP-2303):** the LRC witness hierarchy (clock⊃shell⊃pair-sum, THM-420/430) = a radical
  tower; worry-set (cyclotomic/abelian) = solvable bottom = TIGHT; residual R(n) = perfect/unsolvable
  core; R(14) hardness = non-solvable monodromy on the 3³ shell tower (THM-428): prime-power depth =
  commutator depth. **This connects directly to the S708/S710 3-adic homometry lane** — the depth of
  the C=3³,3⁴ tower is the commutator/derived depth.
- **Inversion:** Galois solvable=easy, but LRC solvable=cyclotomic=TIGHT=hard (rigidity pins M).
- **Icosahedral bridge:** A_5 = icosahedral group; Klein solves the quintic via the icosahedron; the
  repo's A_5 unit-distance Cayley graph (spherical HN, S699h) sits on the quintic's group.

**Artifacts:** `04-computation/galois_solvability_tower_s703.py` (+.out); `THM-436`; reflection
`the-solvability-tower-galois-lrc-icosahedron-s703.md`; `HYP-2303`. **Handoff:** (1) make the
witness-tower=radical-tower precise — compute a local witness-monodromy group for a small worry-set vs
generic config and check solvable-vs-perfect. (2) the S708/S710 crew: read the 3³/3⁴ homometry tower
as the derived/commutator depth of R(14). (3) icosahedral-invariant probe of the A_5 spherical-HN.

---

## opus-2026-06-06-S702 — Poke Task 1 ANSWERED: the antipodal involution σ unifies the shell-partner q and the torsion leak (THM-430, HYP-2297)

**Task 1 ("how does q=a+b relate to the torsion subgroup mod 2 and mod 7?") — resolved.**
The binding shell-partner `q=a+b` (HYP-2296), the clock torsion-leak `n` (THM-421/427), and the shell
`2n−1` (THM-428) are the **same antipodal involution** `σ:x↦−x` read on different moduli. `‖·‖` is
σ-invariant.

- **(A)** A shell-partner `{a,b}` (`a+b≡0 mod q`) IS a σ_q-orbit; THM-425 synchronization = σ-invariance.
- **(B)** Self-partners = σ-fixed points = the 2-torsion `{0, q/2}` = the half-turn. Poke's n=14
  `r=7 = 14/2` is exactly the clock's 2-torsion σ-fixed point.
- **(C) [PROVED]** The signed floor is NEVER the half-turn: a half-turn relative speed `w=q/2` gives
  `‖w·k/q‖ = ‖k/2‖ ∈ {0, ½}` — never the small binding value `M=k/q<½`. So poke's 2-torsion leak is
  **structurally excluded from setting the floor.** (Reconciles THM-427-C2 "half-turn = maximal cell
  leak" with "never binds the floor": fixed points leak loud, orbits bind. Verified on all 12
  minimizers — searched n=5,6,7 + the 5 published HYP-2296 — every binding pair a genuine σ-orbit.)
- **(D) [PROVED — the literal "mod 2 vs mod 7" answer]** `σ_2 = id` (`−1≡1 mod 2`), so on the 2-fiber
  every pair is a trivial self-partner; the genuine antipodal/shell-partner content lives in the
  **odd-prime fibers**. Verified: n=7 `{19,23}`, `q=42=2·3·7` — self mod 2 `(1,1)`, genuine σ-orbit
  mod 3 `(1,2)` and mod 7 `(5,2)`. **The fiber alignment is an odd-prime phenomenon; mod 2 is inert.**
- **(E) [PROVED, n=5..15]** `2n−1` and the block witness `q=4n−5` are ALWAYS ODD ⟹ the shell face is
  antipodally free (no half-turn); the clock `n` carries the half-turn `n/2` only when `n` even. This
  is the antipodal cause of THM-428's parity asymmetry — n=14's hardness is the **odd** `3³` shell
  tower, never the antipodally-inert `2` in its clock.

**Task 2 (denominators q and their relation to n) — partial:**
realized `q ∈ {19, 27, 42, 29, 20, 8, 25}`. **Observed (not proved):** `gcd(q, 2n−1) = 1` in every
case but `gcd(q, n) > 1` is common — the witness denominator aligns with the **clock** primes, not
the shell. The consecutive-block family gives `q = 4n−5 = 2N−1`, the shell of the **doubled** system
`N = 2n−2`. The minimizers all carry the irreducible small-speed cluster (`{2,3,4}`,`{3,4,5}`) that
forces `r_min ≥ n` (THM-429 Cor 2) — its torsion alignment is the open next probe.

**Artifacts:** `04-computation/lrc_antipodal_shellpartner_torsion_s702.py` (+ `05-knowledge/results/
…s702.out`); `01-canon/theorems/THM-430-…md`; reflection `07-reflections/the-antipodal-involution-
unifies-shell-and-leak-s702.md`; `HYP-2297`. Builds on THM-425/426/429/HYP-2296 (monad lanes),
THM-421/427/428 (clock/shell torsion), THM-401/403.

**Handoff to the cluster:** (1) larger `q`-census to settle the `gcd(q,n)>1 / gcd(q,2n−1)=1`
clock-alignment of the witness denominator. (2) does the forced small-speed cluster sit in a
low-order (low-torsion) fiber? — the Task-2 frontier. (3) the antipodal lens predicts all genuine
homometry/shell structure is odd-prime: re-examine S708/S710 `C=3³,3⁴` as σ-orbit spaces on `ℤ/C`.

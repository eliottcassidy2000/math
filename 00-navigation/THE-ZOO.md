# THE ZOO — a living catalog of the repo's objects, methods, frames, forgotten threads, gaps, and a generator for new ones

> **Inventory, not truth authority.** Some “open” and bridge labels predate
> THM-2022 and MISTAKE-211–220. Check [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md)
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md)
> before using a zoo entry as mathematics.

*Seeded opus-2026-07-20-S439 from a six-agent repo-wide sweep (1346 THM/LEM, 2555 reflections,
7650 HYP entries, 10249 result files). Owner directive: "keep adding to the zoo… find the things
we've forgotten we've studied. all of them. and find the gaps between and around them."*

**How to use / extend this file.** Six sections: **§1 Bestiary** (objects), **§2 Toolbox**
(methods), **§3 Frames** (meta-lenses), **§4 Forgotten** (recovered dormant threads), **§5 Gaps**
(negative space), **§6 Generative Engine** (how to produce new questions, with a worked seed).
When you revive a §4 lead or fill a §5 gap, move it to canon and leave a one-line pointer here.
When you find a new object/method/frame, add a row. When you generate a new question, append it
to §6. This is the map, not the territory — details live in canon; this points to them.

> **SIBLING ATLASES (convergent, same session — this is not the only map).** The fleet built
> parallel catalogs the same day, prompted by the same directive; cross-read them and merge, do
> not duplicate (repo frame 17: "many names, one object"):
> [`TOURNAMENT-INVARIANT-ZOO-ATLAS-klein-S399.md`](TOURNAMENT-INVARIANT-ZOO-ATLAS-klein-S399.md)
> (invariant catalog + gap map + 7 procedural generators + revival list — tournament-invariant
> focus; gap-fills HYP-8646/8647), `../07-reflections/the-procedural-generation-grammar-for-the-tournament-zoo-deathstar-S79.md`
> (the generation grammar = this file's §6), `PROBLEM-ATLAS-2026-07-20.md`,
> `METAGRAPH-ATLAS.md`, `ATLAS-OF-ATLASES-2026-07-20-macmini-S124.md`, and kps's
> `../07-reflections/deep-archaeology-the-oeis-and-uncanonized-hyp-layers-kps-S128c137.md`
> (OEIS/uncanonized-HYP layers). **THE-ZOO's distinctive contribution:** the broader *repo-wide*
> (all five topics + engineering + reflections) six-agent forgotten-lead recovery (§4), the
> two-empty-column gap finding (§5), and the insertion-response calculus THM-1900 (§6 seed).

---

## §1 THE BESTIARY (objects) — study-status

Full invariant list: [`05-knowledge/COMPLETE_INVARIANT_CATALOG.md`](../05-knowledge/COMPLETE_INVARIANT_CATALOG.md) (41 invariants). Objects/vocab: [`01-canon/definitions.md`](../01-canon/definitions.md).

| Object | best-studied invariants | frontier / status |
|---|---|---|
| **Tournaments T_n** | H, c₃/t_k, scores, \|Aut\|, isSC, OCF/OCR, Fourier c_k, F/W-poly, Ω/α_k, skew-spectrum, Pfaffian, Betti, domination, binary-form dictionary | mostly n≤7, some n≤9 |
| **Metagraph G_n** | degree, eccentricity, diameter, spectralGap, neutral-arc count (=A000568) | n≤5 exact, n≤7 partial |
| **Merged G_n/ℤ₂** | class-size parity, SC spine/rib/sea edge typing | edge-typed; **own topology unstudied** |
| **Even graphs `widehat(E)_n`** | canonical all-cycle/path envelope; V=A002854, reported χ and ω=χ, chordal≤6; tree images form an `n-2` step diameter chain (THM-4069) | **least-developed row (§5)** |
| **Tilings / waggly / wiggly / blue-black** | independence poly, blue-self parity, transfer-matrix trace, gridSym | combinatorial only; **no spectrum/homology** |
| **Staircase δ_{n−2}** | pin-grid, strips, good-cuts, connectivity | **no spectral/Pfaffian/homology invariant** |
| **LRC speed-sets V** | M(V)/covering-min, danger-bands, Bonferroni, deep-well 14/183, danger-graph clique | metric-rich; **tournament invariants never pulled back (§5)** |
| **GMC/nullcone polynomials** | moments, GIT-nullcone strata, bosonic permanent, charge lattice | strata classified; **topology of the cone unstudied** |
| **Cayley–Dickson tower** | formal-group oddness, 7-ary relations | **H/OCR/spectrum per level unstudied** |
| **Paley / circulant / rotation** | H-max (Paley n=7,11; rotation n≥13=A038375), scalar transfer, path homology | see THM-128/212, LEM-004, THM-1865 |
| **Binary forms** | transitivity dict, Vandermonde=signed tournament sum, resultant/disc | THM-1800/1805/1815 |
| **Arborescences** | hunter/DAG-credit unions | **barely an object — nearly all cells empty** |

---

## §2 THE TOOLBOX (methods)

**35 reused methods** (see the method-inventory sweep; representative cites):
OCF/independence-poly (THM-002/517), Walsh–Fourier (THM-071/259), transfer-matrix trace (THM-027/507),
deletion–contraction/arc-flip (THM-013/116), Burnside (THM-214/586), Seidel switching classes
(THM-474/1470), Pfaffian/skew-det parity (THM-120/1455), simplicial homology/up-Laplacian Betti
(THM-095/124), Schur-convexity/majorization (THM-134/1820), SCC decomposition/king bounds
(THM-333/1860), generating functions (THM-114/502), doubling towers/skew-Hadamard (THM-447/481),
arborescence/Matrix-Tree (THM-923/1750), Kendall–Wei/Helmholtz–Hodge (THM-894/895), covering/
factorial-moment master object (THM-406/661), Weyl equidistribution (THM-664/1061), Farey/three-gap
(THM-565/1196), continued-fraction/Ostrowski (THM-778/1291), Sudler/Dedekind-η products (THM-732/869),
moment-LP/Delsarte/Bonferroni (THM-534/1185), Kakeya/needle/BV-drift (THM-870/1241), inverse-theorem/
Freiman (THM-1004/1203), p-adic valuation sieves (THM-305/992), singular-series/theta (THM-501/523),
unit-distance/chromatic-plane (THM-418/461), moment-nullcone/GIT (THM-1720/1825), Nullstellensatz
emptiness/Gröbner (THM-1685/1740), resultant/discriminant towers (THM-1710/1815), Newton-polygon/
Watson-Laplace (THM-1650/1795), charge grading/EMP (THM-1510/1535), Duistermaat–van der Kallen toral
(THM-1630/1730), equivariant fixed-locus/BU ℤ₂-index (THM-1350/1395), de Bondt–van den Essen cubic
reduction (THM-1320/1370), Mayer cluster/free-probability (THM-438/510), Lean native_decide (~542 files).

**Forgotten single-use tools (reapply these):** OU/fluctuation-dissipation drift on the flip chain
(THM-833 → MCMC mixing, replaces refuted H-gradient-DAG); Tutte-poly of the staircase Smith network +
resistance=good-measure (THM-805 → metagraph currents); Arrow–Condorcet/FKN bridge (THM-512 → FKN
noise-stability bound on H's even Walsh weights); Cayley–Dickson discriminant-form phase law (THM-868–872
→ level-residue classifier); Abel–Ruffini/Galois-solvability wall (THM-929 → non-closed-form certs for
LRC singular series/GMC witnesses); Verblunsky/OPUC coords (THM-883 → Sudler products); Delannoy/
Worpitzky–Eulerian identities (THM-218/087 → OCF coefficient hierarchy); Alcuin/graph-minor monotonicity
(THM-519 → H bounds); Delsarte-LP-**impossibility** as a negative pre-screen (THM-1185 → other covering LPs);
Vandermonde=signed-tournament-sum (THM-1805 → OCF/Paley signed-cactus cancellation); involution-average
Hermite (THM-473/1620 → tournament moment hierarchies); Helmholtz–Hodge curl census (THM-895 → reversible/
circulating split of metagraph current, pairs with THM-833 OU drift).

**New reusable tool (THM-2022/2041): whole-layer Frobenius transport.** Expose
an invariant initial face or exact-period packet, make every other layer pay a
finite-place cost, and transport the complete tied residue by Frobenius rather
than isolating one atom. The method is legal only when a problem-specific
theorem supplies a nonzero base layer and a filtration removes off-layer
terms. LRC exact-period projectors satisfy the preservation step; its base
safe-count/pointwise lift remains open.

---

## §3 THE FRAMES (meta-lenses)

1. **Everything is the triangle** — δ_{n−2} is the universal substrate (`everything-is-the-triangle.md`).
2. **Reify → moduli-point** — every object is a point in a moduli space; degenerate case = vertex.
3. **The moment nullcone** — an invariant is a projection onto a symmetry's trivial component; study its nullcone at finite depth (`the-moment-nullcone-template…`).
4. **The obstruction is the symmetric configuration** — the hard case is always the big-stabilizer object (AP, Paley, one-sided, transitive).
5. **Complexity-as-depth ladder** — rational⊂algebraic⊂holonomic⊂#P; the detecting functional's arithmetic complexity is a coordinate.
6. **Charge grading / killed kernel / Reynolds** — a functional averages over a symmetry and kills nonzero charge (Rédei-parity = GMC's E[…]=a!δ, one grading).
7. **Positivity past the cancellation wall** — reformulate positive-definite; never absolute-bound the signed part (MISTAKE-202 law).
8. **Relation-as-object / intersubjective frame** — a tournament *is* a relation; "second moment is where objectivity begins."
9. **Frames-as-monotone-zoom** — each "most fundamental" frame is a higher-resolution chart of one obstruction.
10. **Observer principle: Rédei = LRC + 1** — the marked origin inserts itself; 1+2·#3-cycles = Farey escape.
11. **Duality web at apex 7** — even/odd ⊕ real/imag ⊕ add/mult hinged on ε=χ(−1); {3,7}⟹LRC(6),LRC(14).
12. **Cancellation is structure / the seesaw** — the deepest results vanish unasked (β₂=0, S(T)=0, OCF).
13. **Polysemous constants: bridge / trap / homonym** — run the persistence test before believing a coincidence.
14. **Computational irreducibility / pockets of reducibility** — score seq 97% reducible, the 3% cycle residual is where irreducibility begins.
15. **Cross-domain grammar: "signed counting over oriented structures"** — amplituhedron ∩ tournament homology ∩ molecular-orbital theory.
16. **Cayley–Dickson doubling tower** — a filtration losing one property per level.
17. **Epistemics-is-the-bottleneck / convergence** — the fleet keeps rediscovering one object under many names.
18. **The five-axis triage law** — a method decides near-floor LRC structure only if it breaks translation-invariance, respects dilation, bounds a max/tail, tolerates signed cancellation, is cross-modulus adaptive.
19. **Inflation-response (new, THM-1865/1900)** — classify each invariant by how a construction (insertion/join/dilation) moves it: **neutral ⟹ rigid extremal; pumped ⟹ inflation-fragile** (the WOWII-103 transfer). See §6.

---

## §4 THE FORGOTTEN — top revivable threads (deduplicated across the six-agent sweep)

Ranked, with location + why. (Full lists in the S439 session letter / task outputs.)

**Cross-topic bridges (highest leverage — they connect the repo's halves):**
- **Metagraph transport 183 = \|PG(2,13)\|** — LRC deep-well AP = transitive pole, Singer set = regular pole of one object; make the LRC-config→tournament-iso-class map precise (`lrc14-history-synthesis-…-opus-S399.md §7.3`, `…-boxeph-S110.md`). *"The only place the tournament project and LRC(14) actually met — untouched."*
- **T-S84-A: Rédei parity of the ban-load tournament as an LRC cover obstruction** (`TANGENTS.md:900`) — Ham-paths mod 4 / Slater index differ for coverable vs non-coverable? *Bridges tournament-parity ↔ LRC; checkable with existing data (100 covers vs 100 non-covers).*
- **Vertex-deletion as a JC-relevant LND** — OCR H(T)−H(T−v)=2Σμ(C) is a difference operator; the Rédei↔Jacobian bridge lives on the derivation/LND side (`is-redei-parity-mathieu-zhao-…-deathstar-S69.md`).
- **LRC-AP as a literal moment-nullcone** — does a resonance functional have all moments vanish at the AP? (run the resonance-matrix computation; `structural-thinking-ways-…-deathstar-S76.md`).

**Near-misses with a stated corrected form still open:**
- **HYP-6720** — M=1/13 for a multi-killer covering 13-set IFF it contains a tight 12-block c·{1..12} (refutation handed the exact counterexample + corrected conjecture).
- **HYP-3805/3819** — flip-rank skip-2 rule hits 454/456 iso-classes; k(7)=12; find the rule handling the 2 missing classes (regular + (1,2,2,3,4,4,5)).
- **HYP-6445** — Q_s=O(r) density survives (the pointwise offdiag≤0 route died, target didn't).
- **HYP-5207** — k=9/10 discharge fails ONLY at d∈{1,2} (2-adic); dispatch those two → uniform.
- **HYP-260** — δ-inequality fails ONLY at source/sink (now explained by THM-1900: those are the extreme condensation down-sets).
- **HYP-8315** — maximiser regular at odd n (open; connects to THM-1865).
- **F-polynomial log-concavity** — 1020/1024 at n=5 (4 exceptions), 100% n=6–8; classify the 4.

**Un-run computations (idea filed, never executed):**
- **HYP-7940** the k=13 tightness sieve (107 primes; the whole k=13 result is conditional on it) — or a cheaper surrogate bound.
- **HYP-2371** R(31)≈2.596 cluster-expansion prediction (cheap run, falsifiable).
- **HYP-3122/3113** the quartic-cumulant stabilizer (never computed; completes the (φ⁴)₂ cap law).
- **HYP-7540** add S₄ (quadruple overlaps) to the moment LP (closes 722/792 r=6 cores).
- **Uniform r=5** (HYP-7600/7605) — the most-repeated open item; a computable piecewise-linear standoff/ceiling minimization.

**OCF / Rédei classic open problems (paper §Open, ~4 months cold):** mixed-graph Rédei (Q-Lemma extension → Schweser–Stiebitz–Toft), the odd-cycle→2-colored-cycle bijection Φ (try LGV lattice paths), Striker–Chapman S₃-equivariance, self-evacuating-SYT bijection (cited formula is *broken* — 2^{3/2}), realizable conflict graphs Ω₃(T), β₂=0 / β₁β₃=0 general-n (stuck at n≤8, band-limitedness is the lever), Krawtchouk coordinate correlations.

**Novel languages named once, never developed:** the **runner braid** (positive braid / torus-link / Garside normal form for loneliness), the **amplituhedron↔tournament dictionary** (volume-reading of H), **MSS interlacing families** for the GMC β-hierarchy (Kadison–Singer method matched to the exact obstruction), **staircase-as-DNA** self-selection (RSK shape of high-H tournaments *is* the staircase).

**The dormant applied continent (silent since March 2026):** the entire cross-domain program — `tournaments-as-codes.md` (quantum stabilizer codes, perfect-hash SRCP), `987-amplituhedron-chemistry.md`, `protein-folding-and-the-tower.md`, `cryptographic-vulnerabilities.md` (v₂-based Walsh-puncturing priority), `active-learning-intelligence.md` (BoostRanker=Free-Energy-Principle) — each ended in concrete next-steps and was abandoned when the fleet turned to proofs. Harvest the anchored/checkable ones.

---

## §5 THE GAPS — the negative space

**The two nearly-empty columns (classical for digraphs, never systematically computed on ANY object):**
- **Feedback arc set (min FAS / linear-ordering defect)** — zero canon theorems, even for T_n where it *is* the classical intransitivity measure. (Only THM-1390's tr≥n−fas anchor exists.)
- **Proper chromatic number χ** (of the digraph, and dichromatic number) — the repo computes clique/ω surrogates only. χ(T_n), χ(metagraph), χ(LRC danger graph) all open.

**The highest-value empty object×invariant cells (agent-6 gap grid):**
1. LRC speed-set × {H, c₃, score sequence} — the observer→A000568 bridge (THM-381) exists but tournament invariants are never *pulled back* onto the generating speed set.
2. LRC speed-set × {skew spectrum, path-homology} — spectral/topological invariants of the danger-overlap complex.
3. Staircase δ × {spectrum, Pfaffian, Betti, domination, clique}.
4. GMC/nullcone variety × {Betti/topology, \|Aut\|, self-complementary}.
5. Cayley–Dickson tower × {H, OCR, skew spectrum, Pfaffian} per level.
6. Merged metagraph G_n/ℤ₂ × {own path-homology, χ, γ}.
7. Canonical even-graph envelope `widehat(E)_n` × {χ, ω, γ, Pfaffian, Betti} — the least-developed row; retain `T` for any basis-relative comparison.
8. Binary forms × {path-homology, domination, feedback arc set}.
9. Pfaffian × ALL non-tournament objects (computed only for T_n, n=6).
10. Arborescences × {skew spectrum, OCR, Betti} — nearly the whole row.

**Engineering gaps (BUILT-but-unshipped / spec'd-but-unbuilt):** PyPI-package `mod_rank`+`circulant_homology` (code done+tested, unpackaged); T₁₉ degree-6 Ω dimension (sparse solver exists, no value recorded); `cycle_deletion_winner` social-choice rule (spec'd, unimplemented); Paper 2 "Path Homology of Paley Tournaments" (all results exist, unassembled); GPU uint8 rank kernel (650× est., unbuilt); attention-head-level LLM intransitivity detector (layer-level gave NULL, head-level untested).

---

## §6 THE GENERATIVE ENGINE — how to keep adding to the zoo

The research space factors as a product. A **new question** = pick an un-visited tuple:

```
   OBJECT  ×  INVARIANT  ×  METHOD  ×  FRAME  ×  OPERATION
   (§1)       (catalog)     (§2)       (§3)      (below)
```

**OPERATIONS** (transformations to push objects through): +source, +sink, insertion(pattern P),
complement/reversal, Seidel-switch, single-arc-flip (wiggly), waggly d-flip, order-join T₁▷T₂,
Cayley–Dickson double, LRC leaf-inflation v→2v, LRC dilation V→cV, append/remove a speed, blow-up.

**Four procedural generators:**
1. **Fill a §5 gap cell** — compute invariant *I* on object *O* for the first time (e.g. χ or min-FAS on anything; Betti of the nullcone; H pulled back onto LRC speed sets).
2. **Transport a method across galaxies** — apply a §2 tool proven in one problem to an object in another (e.g. Mayer cluster-expansion → LRC danger overlaps; MSS interlacing → GMC β-hierarchy; OU drift → switching chain).
3. **Apply a frame to a fresh object** — (e.g. moment-nullcone frame → even graphs E_n; reify-to-moduli → tilings).
4. **The operation×invariant response matrix** — for each OPERATION and INVARIANT, classify
   **NEUTRAL / PUMPED / EQUIVARIANT**. Neutral ⟹ rigid extremal (hard, no cheap attack);
   pumped ⟹ inflation-fragile (WOWII-103 counterexamples exist). This IS the inflation-response
   frame (§3.19) made systematic.

**Worked seed (this session, THM-1900): the insertion-response calculus.** Pushing generator 4 on
the *insertion* operation (add u beating exactly subset P) yielded two exact laws (verified n≤5):
- **c₃-velocity = the forward cut:** `Δc₃(T,P) = e(P → V∖P)`; c₃-neutral ⟺ P is a closed set.
- **H-neutrality = condensation down-sets:** H-neutral ⟺ P is a union of initial strong components;
  **#H-neutral patterns = #SCC(T)+1**. Unifies THM-1865 (source/sink = trivial down-sets),
  boxeph THM-1855 (order-join), kind-pasteur THM-1860 (H=∏H(SCC)), and explains forgotten HYP-260.

**A batch of freshly-generated questions (append your own):**
- Min-FAS and χ/dichromatic number of the metagraph G_n and the LRC danger graph (fills the 2 empty columns).
- Betti numbers / path-homology of the GMC nullcone variety (topology of the unstable cone).
- Pull H, c₃, skew-spectrum back onto LRC speed sets via the observer bridge (THM-381); is M(V) monotone in c₃(observer(V))?
- The response matrix for {dilation, leaf-inflation} × {M(V), q_min, \|D\|} on LRC speed sets (which operations preserve M? — connects HYP-6720, MISTAKE-126 dilation-invariance).
- Cayley–Dickson doubling × {H, OCR, Pfaffian} per level — does a property-loss show in an invariant jump?
- The runner braid: compute the Garside normal form / fractional Dehn twist of the 14-strand LRC braid; is loneliness a meridian-disk condition?
- Is min-FAS inflation-neutral or pumped under +source? (source adds 0 back-arcs ⟹ predict neutral — then FAS joins H,c₃ as a rigid-extremal invariant.)

---

*Provenance: six parallel Explore agents (navigation/backlog, hypothesis graveyard, reflections,
engineering/paper, method inventory, object-grid+Lean-state), opus-2026-07-20-S439. Full findings
in the session letter and canon (THM-1900, HYP-8655).*

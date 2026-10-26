# Session65: recovered concept atlases and the Goldbach machine

2026-10-26. **CITED source audit + FINITE-EXACT checks + scoped elementary deductions.**
Repository archaeology pin: `eliottcassidy2000/math@6a8ad5d3bf01a232da4de41bcbff68390a9337ee`.
No repository files were edited or pushed by this lane. No Goldbach or LRC proof, minimum-state theorem, independent Lean rebuild, or external-priority claim is made.

## 1. Main findings

1. There really is an old session explicitly advertised as **25+ cross-field connections**: opus-2026-03-12-S62c. It is a research atlas, not a list of 25 established equivalences.
2. The anonymous 2016 source is a **27-labelled-state** construction. Identical states reduce it to 26; Yijun Leng's **2024** optimization and Lean verification use 25 states. The wiki's compact row combines this lineage and should not be paraphrased as an original 2016 25-state formalization.
3. The precise statement is **blank-tape halting iff strong Goldbach is false**. It is not an isomorphism between a conjecture and a finite state graph.
4. We did verify a genuine, different isomorphism: the wiki's 25-state table is exactly the pinned Lean table after a bijective state renaming, preserving every symbol, write, move and halt branch.
5. Twenty-five is a description size, not a mathematical invariant of the halting predicate. Even reachable, behavior-preserving padding changes it.

## 2. Recovered sessions: count, source and present status

All paths in this section are pinned to the repository commit above. Historical assertions are not promoted merely because they appear in these files.

| Session / source | Actual recovered scope | Status and reusable part |
|---|---|---|
| **opus 2026-03-12 S62c**, `agents/broadcast/archive/MSG-138-from-opus-2026-03-12-opus-2026-03-12-s62c-25-cros.md` | Headline explicitly says **25+ cross-field connections**. Named files include `weil_rh_tournament.py`, `dihedral_tournament.py`, `tropical_cluster_H.py`, `local_optimality_proof.py`, `p13_exhaustive.py`, `crossover_landscape.py`, `representation_permanent.py`. | **HISTORICAL / mixed evidence.** This establishes the requested provenance hit, not 25 proofs. The Weil script's dictionary actually says 15 connections; the tropical output says 20+ for the program. No unique audited 25-card enumeration was recovered. |
| **codex 2026-06-23 S135**, `05-knowledge/hypotheses/HYP-2936-lrc14-broad-concept-carrier-atlas.md`; `07-reflections/lrc14-broad-concept-carrier-atlas-codex-s135.md` | **21** carriers: coimage/Yoneda, Vitali/Helly, exact periods, state lifts, additive-basis carries, octahedral/Clebsch and others. | Explicit **PROOF-ROUTE / synthesis, not LRC proof**. Good methodological inheritance: retain branch/label data. Its score-induced transitive tournament is a chosen ranking, not an intrinsic mathematical orientation or theorem. |
| **codex 2026-06-29 niche recovery**, `07-reflections/repo-history-niche-connection-atlas-codex-20260629.md`; associated HYP-3492 | Connects old topology/Menger, Haar/Fejer, p-adic lift legality, cyclotomic owner data, automata, and tilings. | **HISTORICAL synthesis / proposed transfers.** DFA minimization of the legal mirror-run graph is explicitly a next experiment, not an accomplished quotient theorem. |
| **death-star 2026-07-20 S68**, `07-reflections/an-atlas-of-lenses-linking-famous-problems-deathstar-S68.md` | Six lenses, including a Goldbach/circle-method placement alongside positivity and exponential sums. | Source itself says **navigational map, not reductions**. Do not upgrade broad claims about one universal law, cluster characters or charge functionals without a specific carrier/map. |
| **Session56**, `05-knowledge/results/session56_rmq_lca_connections_20261025/niche-atlas/ATLAS.md` and `CONNECTIONS.md` | **Ten** niche cards; **seventeen** tested connections in the top-level list. Cards distinguish witnesses, masks, crossing partitions, named ports, low-bit tries, arithmetic size, carry, itineraries, kernels and higher arity. | Closest modern usable template. It separates elementary proved identities, inherited claims, finite controls and **OPEN applications**. In particular, finite kernels need exact boundary states, and a tree is lossless only for a declared query. |
| **S642/S643**, HYP-2218 `goldbach-lemoine-same-pair-projection.md`, HYP-2219 `goldbach-lemoine-same-pair-bridge.md` | Actual Goldbach-adjacent pair carrier: `(p,q) -> (E,O)=(p+q,p+2q)`. For E=14, odd shadows 17,21,25 arise from (11,3),(7,7),(3,11). | Overall files say **OPEN synthesis / exact finite carrier**. The linear inverse and swap identity are elementary exact deductions; neither conjecture follows. The occurrence O=25 is an arithmetic output, not a link to a 25-state machine. |

### A necessary warning about the oldest atlas

The S62c tropical output includes mirror symmetry, Gromov–Witten and cluster-character identifications without constructing the required objects and maps. Preserve them as historical analogies, not current mathematical facts. The useful recovered operation is the **atlas method with explicit loss accounting**, not those claims. A broad corpus search for literal `25 concepts` returned no result; this is a bounded search outcome, not proof that no other such session exists.

## 3. Exact Goldbach object and provenance

Strong Goldbach is

`G := forall n >= 4, Even(n) -> exists primes p,q with p+q=n`.

Equivalently, use all n=2k+4. The machine begins in its specified initial state on the all-zero bi-infinite tape. Halting is the undefined transition, not merely visiting a state. It enumerates even targets and tests possible prime-pair decompositions using unary tape blocks and a trial-divisor subroutine.

Primary sources:

- [Anonymous 2016 gist](https://gist.github.com/anonymous/a64213f391339236c2fe31f8749a0df6), created May 8, 2016: `list27.txt`, `list31.txt`, pseudocode and source. Its 2024 comment links the later verification and explains 27 -> 26 -> 25.
- [Yijun Leng's repository](https://github.com/lengyijun/goldbach_tm/tree/6cd292062516a0a14ea1b34f2ab75154cae7ab1e), pinned SHA `6cd292062516a0a14ea1b34f2ab75154cae7ab1e`.
- [Cryptids table](https://wiki.bbchallenge.org/wiki/Cryptids), observed 2026-10-26: 25-state/49-instruction entry and a newer 21-state/41-instruction entry. This is a moving secondary index; the latter is labelled not independently verified.

The 25-state machine has 25 nonhalting control-state labels and two tape symbols. There are 50 state-symbol slots, 49 ordinary instructions and one halt slot. This convention does not count an extra explicit halt state.

### Formal proof dependency spine, inspected in source

At the pinned Leng SHA:

- `lean-toolchain`: **Lean v4.14.0-rc2**.
- `lake-manifest.json`: mathlib revision **6270f3467257beb5cf5118e21a98a0734356657e**, plus locked Batteries, Qq, Aesop, ProofWidgets and other tool dependencies.
- `GoldbachTm/Basic.lean`: two-symbol alphabet, `Goldbach n := exists x y, x+y=n and Nat.Prime x and Nat.Prime y`, decidable finite pair search, and propagating induction.
- `Tm25/TuringMachine25.lean`: literal machine table, write-then-move operational semantics, blank input `init []`, `nth_cfg`, and `cfg45` (first target 4).
- `Search0.lean`, `Transition.lean`, `Miscellaneous.lean`: tape scanning and local macro-transitions.
- `Prime.lean`: exact prime/composite subroutine behavior, using ordinary natural primality and modular arithmetic.
- `PBP.lean`: `leap_17` advances after a prime pair; `leap_17_halt` halts after exhaustive failure.
- `Content.lean`: `never_halt_step` reaches successive even stages; `halt_lemma'` proves a counterexample implies halting, `halt_lemma_rev` proves halting implies a counterexample. Together:
  `(exists i, nth_cfg i = none) iff (exists k, not Goldbach (2*k+4))`.
- `GoldbachTm.lean`, the default library root, imports `Tm25.Content` (and the 31-state development). The endpoint is not an unimported side file.

This is a proof of **behavioral equivalence**, not a proof of either G or not-G. The explicit general `never_halt` theorem takes nonhalting as its premise; it does not establish that premise.

### Verification actually established here

**CITED successful upstream CI:** [run 26210149415](https://github.com/lengyijun/goldbach_tm/actions/runs/26210149415), May 21, 2026, exact SHA above; build job and `leanprover/lean-action@v1` step succeeded. The pinned workflow uses checkout@v4 and lean-action@v1. This reports the observed upstream result, not our own build.

**Source audit:** inspected the complete local Tm25 import cone, root reach, actual machine and prime definitions; no `sorry`, `admit`, declared `axiom`, `sorryAx` or `native_decide` tokens in saved Lean files. Nineteen saved upstream files byte-match their upstream Git blob IDs, including the exact transition table. Text scans are **not** a kernel dependency/axiom audit. Standard mathlib/Lean dependencies were not rebuilt or exhaustively inspected.

**NOT RUN:** local `lake build` and `#print axioms Tm25.halt_lemma'` / `#print axioms Tm25.halt_lemma_rev`: neither `lake` nor a Lean toolchain is installed in this executor. Thus use **CITED formal proof, source-inspected and upstream-CI-observed**, not independently VERIFIED Lean proof.

### Newer 21-state source

`google-deepmind/cryptids@6c28d8524fa008c6fb5aadb759df74a55ce96423` (Oct23 2026) contains `Cryptids/Goldbach21.lean` with `candidate_fewer_states`, `candidate_equivalent`, and `halts_reference_iff_exists_counterexample`. `Cryptids/Semantics/Goldbach.lean` connects the arithmetic predicate to the project-defined `ConjectureMachines.GoldbachConjecture`, using Mathlib's `Nat.Prime`, ending in `halts_reference_iff_not_goldbach` and `halts_candidate_iff_not_goldbach`.

The [primary README](https://github.com/google-deepmind/cryptids/blob/6c28d8524fa008c6fb5aadb759df74a55ce96423/README.md) reports Lean 4.35.0-rc4 and comparator audits: direct Goldbach state reduction uses standard constructive axioms `propext`, `Quot.sound`; the negated-Goldbach formulation also uses `Classical.choice`. We inspected these source statements, **not** its full dependency graph, comparator output or local build. The wiki's independent-verification caveat remains distinct from the authors' formalization claim. Twenty-five should not be presented as a known minimum or even the current smallest claimed construction.

## 4. Independent exact checks

`python3 -B verify.py > checks.json` and `python3 -O -B verify.py > checks-O.json`; outputs byte-identical. Uses explicit exception checks, not removable Python asserts.

- **All 50 transition slots**: source parse versus wiki code. Exact bijection from Lean states 0..24 to wiki labels:
  `A D E G H I V J W K L M N U O P Q R X Y F S B C T`.
  It preserves read/write/move/next-state and the sole halt branch. This is a genuine finite labelled-control isomorphism. The induced map on full configurations preserves the tape and head.
- **2,000,000 literal microsteps**, no observed halt; all 25 control states visited. Exact unary macroblocks occur for every even n from 4 through 74. First times: n=4 at45, 6 at132, 8 at301, 10 at596; n=74 at1793334. A finite prefix does not prove nonhalting or Goldbach.
- **249 even targets 4..500**: direct prime-pair witnesses equal coefficients of the finite prime-indicator convolution, all positive. This checks the dictionary, not the infinite conjecture.
- **2,025 ordered pairs of odd primes <=200**: pair-shadow inverse, swap reflection and diagonal fixed locus all pass.
- **Marginal-loss hostile**: A={(0,0),(1,1)}, B={(0,1),(1,0)} have identical coordinate projections but disjoint full pair sets. Separate target existence cannot be joined as a common witness.
- **State-only hostile**: Lean control state24 halts on0 and continues on1. Projecting a configuration to just its control state already loses the halting predicate, before any more elaborate graph comparison.
- **Reachable padding**: replace the initially visited instruction (0,0) by the same write and move to a fresh state U; U copies the read bit and moves right to fresh V; V copies its read bit and moves left to the old target. The full old macroconfiguration is restored after two extra steps. All four local two-bit environments checked. U,V are reached from blank tape, so even a reachability restriction does not make the number25 intrinsic. This yields a 27-state implementation of exactly the same computation, not a new conjecture.
- **Separate LRC cutoff hostile**: for V=26{1,...,12} union {339}, exact reduced fractions k/q first meet every distance >=1/14 at q=27, k=2 or25. No q<=26 works. This replays current guardrail THM-762/764; the number25 here is a refuted cutoff, unrelated by cardinality to the machine's control size.

## 5. Five reusable connections with complete loss ledgers

### A. Prime-pair carrier -> coefficient detector -> halting search

**Source:** S642/S643 pair-space work and the formal prime-pair stage invariant. **Map:** for n fixed, W_n={(p,q): p,q prime, p+q=n}; map to c_n=|W_n|=[x^n](sum_p x^p)^2, then to the Boolean c_n=0. The machine searches n=4,6,... and halts on the first empty W_n.

**Preserves:** witness existence / counterexample status, exactly. **Loses:** ordered witnesses and multiplicities upon Booleanization; target n and cursor positions if only the Boolean is kept. **Sidecar:** n, current split, divisor cursor, returned prime verdict, and termination of each finite trial. **Test:** convolution versus explicit pairs (done through500); inject a deliberately rejected prime-pair stage in a *separate toy checker* if testing error handling. **Scope:** elementary exact dictionary plus CITED compiled behavior; not an analytic estimate or a Goldbach proof.

### B. Two linear shadows -> a faithful common-witness interface

**Source:** HYP-2218 + HYP-2219, explicitly OPEN syntheses. **Map:** (p,q)->(E,O)=(p+q,p+2q), inverse (2E-O,O-E); swapping p,q fixes E and maps O to3E-O.

**Preserves:** the ordered prime pair when BOTH coordinates are retained; p=q iff2O=3E. **Loses:** order under E alone and all shared-witness information under independent marginal existence. **Sidecar:** O/order and the exceptional prime2 channel. **Test:** inverse/reflection for2,025 prime pairs and the disjoint-fibre hostile above. **Consequence:** a useful API for joining arithmetic certificates; no implication Goldbach iff Lemoine and no Goldbach/tournament equivalence. Even node14 having odd companion25 is not a machine-size bridge.

### C. Exact continuation semantics -> lawful machine compression

**Source:** S135's branch-retention rule; June29 proposed DFA minimization; Session56 cards4,8,9; Leng's 26->25 explanation. **Map:** merge states only through a relation on reachable full configurations, or a proved macro-simulation that restores all boundary data.

**Preserves:** the halting predicate when both simulation directions and reachable-domain hypotheses are proved. **Loses:** microstep clock and control labels; a state-only quotient can additionally lose tape/read information. **Sidecar:** full tape/head or a proved sufficient tape invariant, boundary port names, and step-count accounting if runtime matters. **Hostile:** same state24 with different read symbols; local transitions seen absent in a finite trace need not be unreachable forever. **Test:** the 50-slot isomorphism succeeds for wiki/Lean; it is not the same proof as 26->25 compression. **Scope:** exact isomorphism verified; general compression is a design criterion. Leng's finite unvisited-branch search is discovery provenance, not by itself proof of unreachability.

### D. Finite-control description -> infinite-state obligation

**Source:** exact Goldbach configuration model; Session56's warning that finite kernels need a truncation certificate. **Map:** project (state,head,tape) to state or to a bounded tape window. **Preserves:** a chosen finite observation only. **Loses:** unbounded target size and distant cursor/data; halting need not descend. **Sidecar:** complete configuration or a proved inductive abstraction covering all reachable configurations. **Test:** read-symbol hostile plus reachable +2-state padding (done). **Consequence:** neither 25 atlas cards, 25 control states nor 25 historical rows supplies a shared mathematical invariant. Use atlas adjacency to choose tests, not to assert equivalence.

### E. Busy Beaver bounds -> genuinely complete finite searches

**Source:** the blank-tape reduction; compare the repository's finite-universe guardrails. **Map:** if a certified upper bound B for the halting runtime of every halting25-state binary machine is supplied, simulate this particular machine for B steps. Halt produces a counterexample; no halt by B proves it never halts and hence G. (One must use the same machine model/counting conventions.)

**Preserves:** exact final yes/no via the universal runtime bound. **Loses without the sidecar:** completeness; arbitrary long successful runs do not decide nonhalting. **Sidecar:** a genuinely proved B and a verified finite simulator; knowing a few busy-beaver champions is not such a bound. **Hostile/test:** machines can agree for any prescribed finite prefix and halt later; the independent LRC first-good27 example likewise shows why an attractive cutoff25 is not self-certifying. **Scope:** conditional elementary consequence, no such B supplied. This is a reusable distinction between a finite check and a finite *complete* decision procedure.

## 6. Recommended next use

The best reusable combination is **A+B+C**: keep the exact witness fibre, expose a genuinely invertible pair interface, and compress only after proving continuation equivalence. These are specific maps already tied to recovered files. Do not expand the deliverable into25 manufactured analogies.

A future formal audit can run the pinned Leng `lake build` plus endpoint axiom prints, verify root imports, and compare the formal table against the saved wiki-code isomorphism certificate. That would close a real verification gap. It would still establish the reduction, not solve Goldbach.

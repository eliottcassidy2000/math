# Four lenses on one research programme: Turing machines, cellular automata, FRACTRAN and tournaments - a translation dictionary, exact translations of the Goldbach machines, the GF(2) automaton of a tournament, and the repository's open questions seen through each lens

mac-mini-2026-11-01-lenses, 2026-11-01. All work committed directly to `main` and pushed after each step.

## 0. The owner's directive (verbatim)

"try to draw analogy between turing machine and cellular automata and fractran programs and our tournament structure
research as you continue seeing open mathematical questions through these lenses, and try to transalte seamlessly
between them"

## 1. Concurrency and inheritance

**Session135** (`session135_fractran_tm_structure_20261101/`, checkpoint 2 at 09:22 today) worked on the same theme.
Its findings, which this note credits and does not repeat:
- the PRIMEGAME basin classifier;
- literal Turing machines for tournament streams (e.g. a 76-state machine for the pure-blue recurrence);
- square-sum and LRC certificate machines;
- the FRACTRAN/register <-> TM compiler;
- a cross-domain loss ledger;
- the circulant word/path bijection, typed elementary.

It had no cellular-automaton translation, no GF(2) automaton of tournaments, and no translation of the Goldbach
machines into each other. Those are the lanes this note takes.

**Earlier pieces used here.**
- S70 THM-4852 (3): FRACTRAN as a selector (quasi-)automaton in base 29#, a known formulation (MISTAKE-697).
- Session126: the six-cell rule-150 automaton on tiling masks, and an actual radius-one automaton for PRIMEGAME.
- Session65: Leng's 25-state Goldbach machine.
- THM-4849: circulants as Collatz cycles.
- THM-4855 and THM-4856.

**Typing.**
- STRUCTURAL: a map preserving a stated predicate.
- DICTIONARY: a shared design without a map.
- NUMEROLOGY: equal counts or labels.
- GENERIC: holds for every object of the class.

## 2. Results first

1. **The four-way dictionary** (section 3). Every arrow between Turing machines, cellular automata and FRACTRAN is an
   exact simulation, and its cost and lost data are stated. Tournaments enter through three exact maps:
   - the GF(2) automaton (THM-4857);
   - path words and Collatz parity words (THM-4855, THM-4849);
   - counting programs (session132/133 FRACTRAN, session135 Turing machines).
2. **Goldbach in all three machine languages** (THM-4858). Each of the following halts, or writes HALT, iff strong
   Goldbach fails:
   - a **53-symbol radius-1 cellular automaton**;
   - a **568-fraction FRACTRAN program** (LENGGAME);
   - the 25-state Turing machine they simulate;
   - the 37-fraction GOLDBACHGAME, which is a different algorithm.

   The checks: the automaton matches the machine at every step for 200,000 steps. LENGGAME matches its stacks at all
   of the first 3,000 steps; this is post-hoc, after a pre-registered checker failure, which is reported.
3. **At n = 5 the owner's ten merged nodes are exactly the ten GF(2) automaton types** (THM-4857). This is a
   small-n coincidence.
   - The GF(2) automaton x -> Ax cannot tell a tournament from its converse. That explains why it cannot see beyond
     merged nodes, which is by definition the converse merge "12 classes become 10 nodes".
   - That it sees exactly the merged nodes for n <= 5 (2, 3, 10 types) is an accident of small n.
   - It fails for every n >= 6: 24 types for 34 nodes, 51 for 272, and n >= 8 by counting (OEIS A070933).
   - At n = 5:
     - the automata of four of the five mixed nodes are reversible;
     - node #6 runs a single 31-cycle (a maximal-period linear recurrence);
     - TT5's automaton is nilpotent.
4. **Circulant tournaments are additive cellular automata** (THM-4857 (3), (4)).
   - The rank is q - deg gcd(a(x), x^q - 1).
   - Paley tournaments are reversible iff q = 3 (mod 8). For q = 7 (mod 8) their image is a quadratic-residue code:
     the Hamming code at 7 and the binary Golay code at 23.
   - This part is KNOWN.
     - It is already canon: THM-481 (June 2026) identifies the Paley tournament's GF(2) code with the QR code, with
       Golay at 23, crediting Kim-Solé 2008, Assmus-Key and T. S. Michael.
     - The dichotomy holds for every doubly regular tournament (A^(-1) = A^T over GF(2) when n = 3 mod 8; A^2 = A when
       n = 7 mod 8).
   - New here: only the cellular-automaton reading of the circulant case.
5. **FRACTRAN's priority order is a transitive tournament** (section 4).
   - Only co-applicable fractions need ordering, so the order that matters is a partial tournament.
   - A tournament-scheduled FRACTRAN is well defined iff the priority tournament has a source on every realisable
     co-applicable set.
   - GOLDBACHGAME: every fraction carries a state prime, so the co-applicable sets are the per-state groups and the
     priority is a disjoint union of chains.
   - PRIMEGAME is different. Its state-free fractions (15/14, 15/2, 55/1 in Guy's list) are co-applicable with every
     state group, since 55/1 always applies. So each state chain must lie above the state-free chain: 26 distinct
     reachable co-applicable sets mix the two (audit).
6. **Open questions through the lenses** (section 5). The repository's questions sort into:
   - halting problems: Goldbach, Legendre, a sixth Fermat prime, the H-spectrum;
   - Π2 orbit questions: Collatz, twin primes;
   - decidable finite computations: LRC(n) for fixed n, the pure-blue alphabet (proved), square-sum (proved);
   - structural questions with no machine content.

   The lens tells which kind of evidence can ever settle a question. It never settles one.

## 3. The four-way dictionary

| from -> to | exact map | cost | preserved | lost / sidecar | status |
|---|---|---|---|---|---|
| TM -> CA | cell alphabet Γ ∪ (Q×Γ) ∪ {HALT}, radius 1 | 1 step = 1 step | the whole configuration and halting (as the symbol HALT) | nothing; locality of a single head | PROVED folklore; Leng: 53 symbols, 200,000-step lockstep |
| TM -> FRACTRAN | two stacks as exponents, one-hot control, helper states for loops | O(2^tape) per step | configurations at TM-step boundaries, halting | speed (exponential); state count grows (25 -> 568 fractions) | PROVED (Minsky, Conway); LENGGAME lockstep |
| FRACTRAN -> register machine -> TM | one-hot decompilation, unary tapes | polynomial | the trace word, halting | the arithmetic of the integer (prime renaming conjugates the run; session135 section 2) | session135 compiler |
| FRACTRAN -> CA | base-P# digit automaton selected by the units digit (S70); a true radius-1 automaton with 385,025 symbols (session126) | 1 step = 1 step for S70's quasi-automaton; many steps per FRACTRAN step for session126's true CA | the integer | S70's map is only a selector (quasi-)automaton, not shift-invariant (MISTAKE-697) | S70, session126 |
| tournament -> CA | x -> Ax over GF(2); circulant = additive automaton with neighbourhood S | - | a merged-node invariant, complete only for n <= 5 | the node itself for n >= 6 (24/34, 51/272); Hamiltonian-path data (H, palindromy) is not linear | THM-4857 (Paley case KNOWN: THM-481) |
| circulant tournament -> words / Collatz | Hamiltonian paths = path words; S = parity word of a Terras cycle | - | H = q · #words; palindromy by segment sums; translation = Collatz shift | which words exist is not visible in the cycle arithmetic | THM-4855, THM-4849 |
| tournament -> FRACTRAN / TM | counting programs (A000568, H(T), pure-blue stream) | large | the counts | the structure behind the counts | session132/133, session135 |
| FRACTRAN -> tournament | priority list = transitive tournament on fractions | - | which rule fires | everything about arithmetic | section 4 |

**How to translate.** Go TM -> CA when locality matters, and CA -> TM when a single active site exists. Go
TM -> FRACTRAN when an arithmetic normal form is wanted, and FRACTRAN -> TM when speed matters. Use the GF(2)
automaton to move a tournament into dynamics: its linear invariants are converse-blind and stop separating classes
after order 5. Use path words to move a circulant into combinatorics on words; THM-4855's segment test is the CA-like
local rule there.

## 4. FRACTRAN priority as a tournament

FRACTRAN applies the first applicable fraction, so the list order is a transitive tournament on the fractions. Two
fractions need comparing only if some reachable integer makes both applicable.

**Lemma.** Replace the list order by an arbitrary tournament P on the fractions, and fire the source of P restricted
to the applicable set. The machine is well defined on its reachable set iff P has a source on every realisable
co-applicable set.
- If the realisable sets include all subsets of some set of size >= 3, then P must be transitive on that set: a
  3-cycle has no source.
- If every fraction carries a state prime and exactly one state prime is present (GOLDBACHGAME, FERMATGAME,
  LEGENDREGAME, LENGGAME), the co-applicable sets lie inside the state groups, so only the per-state orders matter.
- PRIMEGAME's state groups have sizes 3 (state 17), 2 (13), 2 (11), 2 (19), 1 (23) and 1 (29). Its three state-free
  fractions are co-applicable with every state group, because 55/1 always applies. So its essential priority is each
  state chain placed above the state-free chain, not a disjoint union. The audit found 26 distinct reachable
  co-applicable sets mixing the two, and reordering the groups changes the run at step 3.

*Proof.* The rule "fire the source" needs exactly one maximal element. A tournament on a set has a unique source iff
its top strong component is a single vertex, and a transitive tournament has one on every subset.

This is the tournament reading of FRACTRAN's control: priority as an orientation of the conflict graph.

## 5. Open questions through the lenses (typed; prior atlases credited)

The arithmetic-hierarchy classification for Goldbach, Legendre, RH, LRC, twin primes and Collatz is already in
`procgen_atlas_20260924_collatz_implication_atlas.md` (L2) and `analogy_reduction_atlas_20261025/`.

| question (repository) | class | natural machine form | status |
|---|---|---|---|
| strong Goldbach | Π1 | THM-4856 GOLDBACHGAME; Leng TM; THM-4858 53-symbol CA and LENGGAME | OPEN |
| Legendre | Π1 | LEGENDREGAME (THM-4856) | OPEN |
| a sixth Fermat prime (THM-4847) | Σ1 ("exists") | FERMATGAME halts iff yes | OPEN |
| H-spectrum: every odd h except 7, 21 is some H(T) | Π1 | a TM enumerating h and searching tournaments on at most m(h) vertices | OPEN (session133/134) |
| Collatz (every orbit reaches 1) | Π2 | generalized Collatz = FRACTRAN (Conway); base-6 quasi-automaton with a global bit (repo Gilbreath note, Kari) | OPEN |
| twin primes | Π2 | PRIMEGAME visits 2^p, 2^(p+2) consecutively infinitely often (procgen atlas) | OPEN |
| LRC(n), each fixed n | decidable | certificate machines (session135 lrc_machine) | LRC(14) OPEN |
| square-sum paths for n >= 25 | Π1, proved | certificate TMs (session135) | PROVED (Gerbicz, CITED) |
| pure-blue alphabet {K1, C3, T5} | Π1 (a universal statement over orders), proved | path words with a segment test (THM-4855) | PROVED |
| merged metagraph laws (HYP-9345) | per-n finite | the GF(2) automaton classifies nodes only through n = 5 | mixed |

Why the H-spectrum is Π1:
- H is multiplicative over strong components (a Hamiltonian path visits the components in order).
- A strong m-tournament has a Hamiltonian cycle (Camion), so H >= m. This already bounds a realiser of h by h
  vertices (audit). Busch's minimum, which grows like 5^(m/3) (THM-4744, CITED), sharpens this to O(log h) vertices
  per strong component.
- So "h is realised" is decidable, and the conjecture is a halting question.

This is implicit in THM-4744's complete inventories below 625. Only the halting form is spelled out here.

## 6. Reproduction

Run these from `04-computation/experiments/circulant_fractran_goldbach_20261101/`:

```text
python3 machine_translations.py          # ~8 s; exit 1: T2b and T3 failed as stated (checker alignment)
python3 machine_translations_posthoc.py  # post-hoc: T2b and T3 with the alignment fixed
python3 tournament_automaton.py          # ~7 s (A1-A5)
python3 automaton_similarity_types.py    # ~2 s (S1-S3)
python3 automaton_similarity_n7.py       # ~1 min (S4)
```

## 7. Audit

**Verdict.** Independent audit (`04-computation/experiments/circulant_fractran_goldbach_20261101/audit_lenses/AUDIT.md`;
own code in `audit_gf2_automaton.py`, `audit_machines.py`, `audit_note.py`): **SOUND WITH CORRECTIONS**, all applied
(MISTAKE-701). The core mathematics and both simulations hold.

**Corrections applied.**
1. **n = 8 is not open.** There are 746 similarity classes of 8 × 8 GF(2) matrices (A070933), below the at least 3,440
   merged nodes. So the automaton fails for every n >= 6.
2. **PRIMEGAME is not one-hot in the sense of the priority lemma.** Its state-free fractions are co-applicable with
   every state group. The lemma's clean form holds for GOLDBACHGAME and the other fully guarded programs.
3. **T2c covered 66 TM steps, not 300**, because of the same misaligned counter. The audit's literal run matched all
   300.
4. **THM-4857 (4) is KNOWN.** It is THM-481, with Kim-Solé, Assmus-Key and T. S. Michael. The audit added the
   doubly-regular generalisation, Sutner 1988, Taussky-Zassenhaus 1959 and de Caen 1991.
5. **Framing.** "A cellular-automaton fact" overstated a small-n coincidence. The dictionary row now reads
   "invariant, complete only for n <= 5".
6. **Minor items.**
   - S1 could not fail.
   - The image is the other QR code (reciprocal generator; same parameters).
   - Separability is not needed in the rank formula.
   - "Containing 1" was added to the circulant count.
   - The FRACTRAN -> CA cost row was corrected.
   - Busch's bound is not needed for the Π1 form.
   - The T3 label in machine_translations.py still reads "step 6"; its condition tests step 7.

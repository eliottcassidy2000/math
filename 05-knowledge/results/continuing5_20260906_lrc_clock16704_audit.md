# Independent referee: the two coordinates at clock 16,704

**Status: PASS, FINITE-EXACT relative to the pinned inherited projected-word domain and maximum table.** This audit accepts the bounded stopping theorem in `continuing5_20260906_lrc_clock16704.md`. It independently recovers all 988 old-filter candidates, the first 107 located profiles, and the full forced-(4,4) completion maximum. It does not eliminate the entire clock; the remaining 881 located profiles are outside this artifact.

## 1. Mathematical scope and dependencies

The input is the actual strict atlas of 5,855 coprime pairs p<q with p+q<=356, all prime divisors of the sum congruent to 2 modulo 3 and having exponent at most two. The five possible small sheet divisors at t=16704 are 1,2,3,4,6. The inherited six-level projected profile JSON supplies allowed states and proper-subset tests; `third_20260906_grid_full_words.out` supplies the already audited unconditional maximum M(16704)=188. The audit pins the complete maximum-table semantics, rather than substituting an attaining word for the optimization proof.

The necessary filter reserves the two distinguished margins d_p=e gcd(t/e,p), d_q=e gcd(t/e,q), including multiplicity when equal. Each state's capacity is the number of hereditary ceilings in (90,30,9,4,2,1,1) at least that state. Five remaining slots maximize the native ceiling cost b_d=d ceil(t/(7d)). This is equivalent to the producer's residue-cost bag because exactly seven baseline contributions total t. My verifier instead uses bounded knapsack, with explicit slot count and capacities. The resulting relaxed objective is an upper bound for full compatible words; it need not have the same owner as the unconditional maximum.

The physical transfer remains the inherited sufficient grid inequality. An actual selected pair reduces the lifted t-grid to an n=t/e grid with multiplicity e. The danger sets are open, ||px||<1/14 and ||qx||<1/14, and positive credit minus compatible full-word excess yields weak-safe lifts. No statement here asserts strict clearance, physical realization of a projected word, or a new general LRC closure.

## 2. Independent geometry and complete filter universe

I reconstruct each pair's intervals by sorting all raw rational danger boundaries and evaluating the literal midpoint predicate. This avoids the producer's clipped-length formula. All 29,275 pair-sheet inputs are accounted for: 28,060 fail separate credit against 188; 227 fail margin availability/capacity; none further fail the paired objective; exactly 988 remain. Every resulting seven-column row agrees with the frozen output, not merely its hash or count.

For the first 107 rows I use a separate enter/exit event sweep. As the common translation increases through a wall, departing points are removed at the wall and arriving points are added only after it. Thus the wall value is the left-cell count minus departures, while the right-cell count additionally includes arrivals. This also treats coincident events and strict endpoints correctly. A literal native-grid count initializes each sweep and separately verifies every published extremum owner. All first 106 rows close; index 106 is the first surviving located quotient row.

For this selected row I use a third, larger partition: every raw boundary of either danger set projected to the grid phase, including inactive boundaries. All 692 walls and all 692 intervening open cells are checked on all 4,176 grid points by literal integer modular inequalities. The result is exactly

    (e,p,q)=(4,23,323), (d_p,d_q)=(4,4),
    C_sep=C_loc=180, E_pair=207, M=188,
    min C(alpha)=45 at alpha=0, max C(alpha)=92.

The published maximum owner 75651/104006 is independently checked. Since the raw partition includes every possible predicate change, these extrema are universal in the common translate. The sum 346=2*173 belongs to the strict atlas. This is a genuine sharp obstruction to claiming location always improves separate interval floors.

## 3. Full-word conditioning and the first failed implication

I enumerate all binom(22,5)=26,334 nondecreasing five-state completions of the distinguished pair (4,4). The 18 states are exactly the inherited level-six gcds dividing 16704. For each completed word, a subset-gcd bitmask recurrence tests global gcd one and all 126 nonempty proper projected subsets. There is no bound-based pruning and no import of producer code.

Exactly 2,422 completions are valid. Their largest literal ceiling excess is 134, attained by (4,4,3,9,36,58,64). Exhaustion proves the upper bound and the displayed owner proves attainment. No valid completion reaches 180. Consequently any actual row with this selected pair and inherited profiles has at least 180-134=46 weak-safe lifts.

The tempting but false implication is that C_loc<=min(M,E_pair) supplies an owner with sufficient cost. The two maxima were computed on different relaxations. Retaining the forced pair inside the full word repairs that implication. Conversely the earlier (3,308) example needs phase locations even when the full conditioned objective is 188. These are independent missing coordinates. The combined sufficient test uses the same selected pair on both sides: C_loc>M_(d_p,d_q).

I found no mathematical repair needed in the producer's scope. Its first-survivor wording is explicitly relative to the old located quotient and is followed by the conditioned exclusion. The source and report do not claim that the remaining 881 rows survive, or that a whole-clock theorem has been established by this bounded prefix.

## 4. Reproduction and frozen evidence

Run the adjacent `continuing5_20260906_lrc_clock16704_audit.py` with ordinary Python and with `python -O`. Both pass **6,517 always-active exact gates**, with byte-identical raw LF output. The audit imports no producer implementation. Adjacent dependencies and the standard filed results layout are preferred, with the current outside worktree as fallback.

- Audit source SHA256: `b0f3a095427cc03c4078974a170e992071daccc73464a4c822c3501272710195`.
- Audit output SHA256: `0e2d82884b53c200dfaa65740b631c0ad4c5291a34493bfd3f08c9565e48106d`.
- Audited producer source: `2d12351426bdf67c0801bdc46b16c50e35f6f6600b49e8e5323d4cc975bc2d3a`.
- Audited producer output: `6d6acf82d798074bd5fac208f79c6b91bf7b9df9f425f1204d9e56324f5cac76`.
- Audited prepromotion report: `dbd8c110e3afdbb095c7d427d15df6dc37a3b3c4fb2cab53648439c58b957817`.
- Inherited profile JSON: `935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.
- Complete maximum-table semantics: `ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca`.
- Full candidate semantics: `8610d17932b51afbaaeec4699f5a21e2362441f58b9c28235c48182106d45bc1`.

All files were created outside the repository. No maintained document or Git state was changed by this referee.

The [raw-byte checkpoint manifest](continuing5_20260906_manifest.json) pins
the filed source, report and identical normal/optimized transcript. Any
candidate-report hashes above identify the pre-promotion bytes.

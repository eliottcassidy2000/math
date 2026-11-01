# Certified fraction, register and Turing-machine translations

**Status: PROVED all-input simulation and cost identities; FINITE-EXACT emitted-table replays.** No novelty, shortest-program, smallest-machine, binary-time efficiency, or new undecidability result is claimed. The actual tournament interpretations retain the inherited hypotheses stated below. This is not a proof of LRC(14), the unrestricted Hamiltonian-path spectrum, or any other open repository target.

## 1. Results and the inheritance ledger

This packet makes four existing/newly specialized arithmetic streams into actual finite Turing machines, not an appeal to universality:

| Stream | Original ordered fractions | Data registers | TM tapes | Alphabet symbols | TM states | Literal table rows |
|---|---:|---:|---:|---:|---:|---:|
| S132 pure-blue | 22 | 8 | **1** | 257 | 76 | 10,926 |
| S133 compact fixed-H tail | 43 | 15 | 15 | 3 | 168 | 293 |
| S134 full migrating word | 25 | 9 | 9 | 3 | 1,044 | 2,063 |
| Same word on its invariant face | **13** | **5** | 5 | 3 | 503 | 993 |

Every transition row in a direct `*.tm.json` table or its authenticated `*.tm.rows.*.json` shards is fully specified. The two larger tables use small metadata indices listing row shards of at most 90,000 UTF-8 bytes; concatenation recovers the literal finite row list without expansion or interpretation. There are no wildcards, arithmetic callbacks, macros, or external recurrence updates in those files. State counts include all listed states; tables are partial deterministic tables, with missing transitions interpreted as halt. The well-formed streaming inputs below never encounter a missing transition. `BUILD.json` independently exposes representation counts. All TM models permit the head-stay move S; one table row always costs one transition.

Closest proved mechanisms: S132's seven drain phases and 22 literal fractions; S133's positive-linear compiler with 6d+1 fractions; S134's source-preserving positive four-coordinate migration semiring and its radix-80 word. Verbatim copies of the actual S132/S133 fraction lists and the S133 compiler are pinned in `SOURCE_PINS.json`.

Canonical hostiles: reducing a fraction can erase an intended guard; an unconditional fallback before its work row changes first-eligible semantics; a partially tested denominator must be restored before testing the next fraction. Corrected near miss: a finite trace or generic universality theorem does not supply the literal machine, valid-input theorem, readout convention, or exact cost. Least-used sidecars: S132's control-prime one-hot invariant and S134's preserved signed lattice invariant `8C-9t`. The concept board is control erasure, exact boundary decoding, unary cost, positive migration, and representation loss.

## 2. The machine models and input contract

A register machine instruction is one of

- `INC r; goto u`,
- `DECJZ r; goto u if positive (and subtract one), else goto z`,
- `GOTO u`, or `HALT`.

Registers are nonnegative integers. These are genuine finite counter programs, not unit-cost arbitrary multiplication or matrix-operation machines.

For r registers the **packed single-tape** TM has alphabet `{'#'} union {0,...,2^r-1}`, blank 0. Cell 0 is `#`, all negative cells are blank and never visited, and cell j>=1 contains the bit mask of registers whose values are at least j. The head starts at cell 0. Thus the eight-register pure-blue input vector `(1,1,1,2,1,3,0,0)` is physically

    # 63 40 32 0 0 ...

with initial state C23. Decimal symbols 63,40,32 are each one tape symbol, not multi-cell binary numbers. A register value is the length of its initial run of 1 bits on its track. The alphabet is finite but deliberately larger than binary.

For r registers the **r-tape unary** TM has alphabet `{'#','1','_'}` and blank `_`. On tape i, cell 0 is `#`, followed by exactly v_i copies of `1`. Every head starts at its marker. A transition reads/writes an entire r-tuple and gives one L/R/S move per tape; it still counts as one multitape transition. Every inactive head stays at its marker. In the 15-tape table, every read/write/move field is exactly 15 literal characters. Inactive tapes explicitly read and write `#` and move S: those positions are **not wildcards**. There are genuinely 293 ordinary delta rows, not 293 symbolic clauses or a compressed family of arbitrary read tuples. Other read tuples are absent from the partial transition table, and the encoding invariant proves they are unreachable from valid inputs. It would be false to call this model a one-tape TM or to assume a constant-factor conversion to a fixed binary one-tape model. No such conversion is claimed.

At every register-instruction boundary, all heads are at their markers and the unary/packed encoding is exact. The declaration of the initial vector is an encoding contract, not an uncounted update in the executing TM. Converting a binary integer input into its factorization/unary tape is **not** included in the runtime theorems.

## 3. Register-to-TM compiler and exact costs

`to_tm` emits the complete finite transition table. An INC starts at the marker, moves right through the chosen register's 1-run, sets its first zero to 1, and returns left to the marker. A DECJZ first inspects cell 1. If its bit is zero, it returns to the marker and takes the zero branch. Otherwise it scans to the first zero, moves left, erases the last one, returns to the marker, and takes the positive branch. Packed symbols preserve all other bits; multitape instructions preserve all other tapes and head positions. GOTO is one S transition at the marker. HALT emits no transitions.

For a register with old value u, both INC and DECJZ cost exactly

    2u+3

TM transitions, including DECJZ at zero (cost 3). GOTO costs exactly 1. For INC: marker departure 1, u forward moves, turn/write 1, u backward moves, marker dispatch 1. For positive DECJZ: departure 1, u forward moves, turn 1, erase 1, u-1 backward moves, dispatch 1. The zero case is departure, turn and dispatch.

These counts also prove termination of each primitive. They prove the encoding invariant by induction over all executed instructions, for **all finite nonnegative register inputs**, without a finite testing cutoff. At each instruction boundary the TM and register machine have equal control and decoded registers. Conversely, boundary-decoding every valid generated-TM run recovers exactly that register run: the correspondence is a stuttering simulation in both directions on this certified image, not a general arbitrary-TM decompiler.

## 4. Every finite FRACTRAN list has an explicit priority-preserving translation

`general_fractran_ir` handles every finite ordered list of positive rational numbers, including nonreduced input and the empty list. First reduce each fraction a_i/b_i. Take one register for every used prime p and encode its valuation v_p(N). Factors coprime to all used primes are inert, and are retained as a scalar sidecar when reconstructing N. The degenerate no-used-prime case uses one unchanged dummy register.

At boundary STEP, test candidate i by trying the b_i prime decrements in a fixed finite order, with multiplicities. Finite control remembers which prefix succeeded. If a decrement sees zero, restore exactly the successful prefix with INC instructions and try i+1. If all denominator decrements succeed, apply the numerator INCs and return to STEP. If every candidate fails, enter HALT with every register restored to its old value.

**Proof.** Reduction makes `N*a_i/b_i` integral iff every denominator valuation is available. A failed partial test followed by its fixed inverse restores the complete input vector, so the next test sees the original N, not a damaged remainder. Every test, rollback and update has a fixed finite number of register instructions, so the first successful candidate is reached after finite time; no smaller eligible candidate is skipped. Its final register vector is exactly `v(N)+v(a_i)-v(b_i)`. If none succeeds the machine halts. This proves one-step correctness, halting equivalence, and hence finite/infinite trajectory equivalence for every positive integer N. Every original FRACTRAN step corresponds to a finite, nonempty TM segment.

An exact register-cost formula is also available. If j is the chosen candidate, D_j and A_j are the denominator/numerator valuation sums, and t_i is the number of successful decrements before the failing zero test for i<j, the register instruction count is

    sum_(i<j) (2t_i+1) + max(1,D_j+A_j).

For terminal failure it is the sum of the failed-candidate terms. Empty rollback paths and update paths are wired directly; a 1/1 step uses a single GOTO. The exact TM cost is obtained by summing `2u+3` for each INC/DECJZ at its actual old value u, and 1 for GOTO. These costs are constructive, not a claim of unit-cost multiplication.

Hostiles retained in `verify.py`: `[5/6,7/2]` on N=2 must restore the consumed factor 2 and return 7; `[1/12,5/2]` on N=4 restores two consumed factors; `[2/2]` on N=1 must run, not incorrectly halt because the displayed denominator 2 fails to divide 1; `[3/2,5/2]` on N=2 must select 3/2. The empty list must halt without changing N.

## 5. One-hot control erasure for the actual repository lists

The smaller `to_ir` compiler uses an explicitly verified subclass. Data primes and control primes are distinct; every reduced row has exactly one control prime, to the first power, in each of numerator and denominator. All other prime factors are data. At each source control there is either one unguarded row, or a unit-data-guarded row followed in the **original order** by an unguarded fallback. The compiler rejects violations; it does not silently reorder or extend the accepted subclass.

A well-formed integer is

    N = q * product_i p_i^v_i * U,

where q is the unique control prime with exponent one, v_i>=0, and U is coprime to all program primes. Its decoding is finite-control state Cq, data vector v, and inert U. Rows at all other controls are ineligible. At q, the positive branch of DECJZ executes exactly the first work fraction; the zero branch executes exactly the fallback. The work decrement supplies the denominator data factor, and numerator exponents are unrolled INCs. Unguarded rows become INC chains or a GOTO. The resulting control is exactly the numerator control. Induction preserves one-hot control and U. This proves equality with the **literal first-eligible integer computation at every fraction boundary**, not only at sequence checkpoints.

Control erasure reduces S132 from 23 used primes to eight data registers and 15 finite control labels. S133's compact tail reduces 44 used primes to fifteen data registers and 29 finite-control labels. Multi-control integer seeds remain outside this optimized compiler theorem; the general compiler in section 4 still handles them without special meaning for their outputs.

## 6. The reverse register-to-FRACTRAN compiler: cancellation-safe loops

For any finite register program, assign data primes p_r, control primes q_s, and distinct fresh trampoline primes a_s,b_s. Encode a register configuration by its one-hot control prime times its data prime powers. Emit:

    INC r;goto u:   a_s*p_r/q_s, q_u/a_s
    GOTO u:        a_s/q_s,     q_u/a_s
    DECJZ r,u,z:   a_s/(q_s*p_r), b_s/q_s, q_u/a_s, q_z/b_s
    HALT:          no source row.

The work row precedes the fallback. All ratios are reduced because data, controls and trampolines are distinct. This includes self-loops u=s or z=s, where the naive direct fraction could cancel its control guard. A one-hot control admits exactly its relevant work/zero row; the trampoline admits exactly its return row. Thus **exactly two FRACTRAN steps implement every nonhalting register instruction**, for all register values. HALT is exactly a FRACTRAN halt. The invariant and priority argument prove this, and 1,764 individual instruction cases plus explicit self-loops are independently exercised.

These are two rational steps per **register instruction**, not per original source fraction. Source fractions with large numerator exponents expand to many INC instructions. The roundtrip pure-blue program has 60 fractions, not 22; this construction is a correctness bridge, not an improvement in list length.

## 7. Actual pure-blue one-tape machine and its quadratic exact time

The input is state C23 and tape encoding `(a,b,c,d,e,f,0,n)` with arbitrary nonnegative entries. A checkpoint is state C23 with head on `#`. The decoded next checkpoint is

    (b,c,d,e,f,a+e,0,n+1).

Each transfer loop drains source s into m distinct unit destinations whose initial values are d_1,...,d_m. Per unit, the guarded DEC is followed by the destination INCs and the auxiliary-return GOTO. Including the final zero test, its exact TM cost is

    (m+1)s^2 + (2m+5)s + 2s*sum_j d_j + 3.

Substitute the seven actual S132 phases `(a→scratch), (b→a), (c→b), (d→c), (e→d,scratch), (f→e), (scratch→f)`, and the initial clock INC. The exact one-round result is

    T = 4a^2+5e^2+6ae + 2(b^2+c^2+d^2+f^2)
        +14a+16e+7(b+c+d+f)+2n+24.

The corresponding original FRACTRAN cost is

    F = 8+2(2a+b+c+d+2e+f).

The fixed sequence seed is `(1,1,1,2,1,3,0,0)`, corresponding to integer 817086270. At checkpoint n the first six registers are `(p_n,...,p_(n+5))`, with `p_(n+6)=p_(n+4)+p_n`; the last register is n and scratch is zero. The TM intentionally never halts. To obtain a requested output, an observer stops its simulation at the requested checkpoint; no external modification of tape or integer is part of the stream.

The tournament interpretation is exactly inherited from S132 `pure_blue_fractran/README.md`, including its corrected all-order pure-blue atom classification and self-converse gauge. The machine transports class counts under the SCC wrappers, not individual tournaments, arbitrary metagraph edges, automorphism groups, or Hamiltonian-path witnesses.

## 8. Positive-matrix compiler: cost depends only on the two endpoint vectors

Take any d>=1, nonnegative integer matrix A, initial x in N^d, zero scratch registers y, and clock n>=0. S133's compiler first drains each x_i into the scratch coordinates with multiplicities A_ji, then drains the resulting y=Ax back into x, then advances the clock. It has 6d+1 literal fractions. Its TM compilation uses 2d+1 unary tapes and the exact following round costs:

    FRACTRAN:  1+2d+2|x|_1+2|y|_1;
    TM:        ||x||_2^2+3||y||_2^2+5|x|_1+9|y|_1+6d+2n+4.

**Proof of the less obvious identity.** In an accumulation phase, unit i is decremented with cost `2u+3` and its auxiliary return costs 1. Across the entire phase these contribute `||x||_2^2+5|x|_1`. For each scratch coordinate j, all its INCs together start at 0 and end at y_j; their costs telescope to `y_j^2+2y_j`, regardless of which input column supplied them or how they were grouped. The d final zero tests add 3d. Each restore loop drains y_j into an empty x_j, costing `2y_j^2+7y_j+3`. The entry GOTO costs 1; the final clock INC costs `2n+3`. Summing gives the formula.

This explains why matrix density and coefficient placement disappear from the exact time identity: they change the path and finite program size, but their unary increments telescope. It does **not** mean program size, binary arithmetic complexity, maximum intermediate encoded integer, or tourney information is determined by endpoint norms. Both formulas hold for **every** nonnegative x and A; finite tests are regression evidence, not their proof.

The actual S133 compact tail uses seven x coordinates for H=`1,3,5,9,11,15,17`; the other counts are aliases (13=11, 19=17) or identically zero (7,21). Its starting order is 6. This arithmetic stream's tournament multiplicity interpretation remains conditional on the completeness of the inherited bounded strong inventory, as in S133; machine correctness itself is unconditional.

## 9. S134 migration: a useful two-coordinate specialization, not a scalar fiction

The actual S134 word `[(3,0,2),(4,29,3),(3,5,2)]`, in positive coordinates `(C,L,g,t)`, has

    M = [[64,44,46,18], [0,0,0,0], [0,0,0,0], [56,38,41,17]].

It preserves the face L=g=0. On that face only `(C,t)` is needed, with

    A = [[64,18],[56,17]],    F=2C+t,    I=8C-9t.

Thus a 25-fraction, nine-data-register implementation specializes to a **13-fraction, five-data-register** implementation, with 5 rather than 9 unary tapes. This is a proved invariant-face reduction, not a minimum claim. For every input on that face, specialization saves exactly four original FRACTRAN steps and twelve multitape TM transitions per round (the removed zero-coordinate phase exits); this is a comparison within the explicitly stated multitape cost convention. Both literal FRACTRAN lists and both actual TM tables are included. The transfer preserves `I'=I` and gives `F'=80F+3I`.

For the 13-fraction list, data primes `(2,3)` encode `(C,t)`, scratch primes are `(5,7)`, clock prime is 43 and checkpoint prime is 41. Thus its **literal integer seed is 123**, and its kth checkpoint is `41*43^k*2^C_k*3^t_k`; read the output as `2*v_2(N)+v_3(N)`. This is a linear valuation readout, not a pure power-of-two readout. The full 25-fraction list uses checkpoint 97, clock 101, and seed **679**, with output `2*v_2(N)+v_7(N)`.

Start with `(C,t)=(0,1)` (and zero scratch and clock). This is not an arbitrary promised numerical state: S134's authenticated source is the one-vertex parent with empty marked cut, with full tournament order two. Its actual source-preserving graph construction supplies the induction, adding seven vertices per word. Then I=-9 and

    F_k=(52*80^k+27)/79, C_k=18*(80^k-1)/79, t_k=(16*80^k+63)/79.

Hence the first two noninitial outputs are 53 and 4213. These outputs are independently checked by executing the emitted five-tape table, not by substituting the recurrence. The all-k graph interpretation is inherited from S134 `radix_probe/PROOF.md`, section 4, and `constructor/PROOF.md`; these proofs also retain the old port, cut, and actual graph recipe. The machine packet only carries the count coordinates, not the graph itself.

The profile invariant needs a lattice sidecar: `C=9(F-1)/26`, `t=(4F+9)/13`, so integral states on this ray require `F=1 mod 26`. The actual orbit is a sparse subset of that residue class. A scalar stream F alone does not prove arbitrary scalar values realizable, cannot decide the full odd-H spectrum, and cannot erase the authenticated source/cut when used as a construction theorem.

## 10. Loss ledger, boundaries and reproducibility

| Map | Preserved | Lost / retained separately |
|---|---|---|
| positive integer to valuations | divisibility, first eligible row, exact selected successor | inert prime factor U retained to reconstruct integer |
| one-hot control primes to finite state | all source fractional boundaries | malformed multi-control behavior excluded; use general compiler instead |
| register instruction to literal TM segment | values, branch, termination, exact transition time | binary-input conversion and hardware implementation not charged |
| four migration coordinates to invariant face | actual orbit, H output, I, all future word iterations | other initial states not on L=g=0; arbitrary recipes not included |
| actual tournament to count stream | this word's Hamiltonian-path count | individual arcs, source/cut witnesses and labels remain graph sidecars |

There is no claim that the particular streaming examples have undecidable halting: all valid checkpoint inputs of these examples stream forever by their proved drain-phase invariants. General translation and general universality do not imply hardness of an individual fixed list or orbit. This packet proves simulations directly, and does not use or claim an undecidability theorem.

Reproduce with standard Python only, from any working directory:

    PYTHONDONTWRITEBYTECODE=1 python3 /path/to/build.py
    PYTHONDONTWRITEBYTECODE=1 python3 /path/to/verify.py
    PYTHONDONTWRITEBYTECODE=1 python3 -O /path/to/verify.py

All load-bearing checks use explicit exceptions, not asserts. The table loader authenticates each shard and merely concatenates its ordinary rows; this introduces no machine operations or semantic compression. The verifier scans each original ordered fraction list by exact integer product/divisibility; its TM emulator executes only the serialized transition table. It compares every original fractional boundary with the register and integer views, not just final outputs. It exhausts all 729 pure-blue input windows in `{0,1,2}^6`, maintains a continuous physical tape for 25 pure-blue outputs, checks actual compact-tail and migration rounds, tests eight hostile/degenerate matrix pairs, exercises reverse compilation, and rejects malformed optimized compiler inputs. The general compiler tests preserve cancellation, priority, rollback, halt and inert-factor controls. `normal.json` and `optimized.json` record the finite universe and results. `MANIFEST.json` authenticates delivered artifacts; hashes authenticate bytes, not mathematical truth.

# Session 126 continuation: an actual radius-one cellular automaton for PRIMEGAME

> Publication note: the compact microtrace is included. The bulk READY trace is regenerated and hash-verified by `../run_all.py`; the implementation, proof and check counts below are unchanged.

## Result

`primegame_ca.py` defines and executes a **homogeneous, finite-alphabet, radius-one cellular automaton**. Its rule is the finite case function `local_rule(left, center, right)`. It is not merely a generic computability argument or a FRACTRAN interpreter labeled “CA.” The underlying configuration space is `A^Z`; the implementation represents finite-support configurations sparsely.

At distinguished **READY** times, an explicit encoding/observer recovers exactly the integer trajectory of the fixed, ordered fourteen-fraction PRIMEGAME. The construction handles **every positive integer input**, including multiple simultaneous control-prime factors and prime factors outside the first ten primes. It does not assert a one-CA-step conjugacy, a finite ten-cell model, a minimal alphabet, a Conway-Life implementation, or a tournament-metagraph isomorphism.

The preceding FRACTRAN directory was left unchanged. This continuation is entirely in `fractran_ca/`.

## 1. Finite alphabet, exactly counted

A tape symbol has twelve binary fields:

1. Ten unary counter-track bits indexed by `(2,3,5,7,11,13,17,19,23,29)`.
2. One persistent origin-marker bit.
3. One passive spectator bit.

Thus the tape alphabet Γ has **2^12 = 4,096** symbols. A cell optionally carries one finite-state head. The head-state set Q has:

- one `READY` state;
- one `SEEK(rule,operation)` state for each of 37 operations;
- one `RETURN(rule,operation)` state for each of those 37 operations;
- one `DECLEFT(rule,operation)` state for each of 18 decrements.

Therefore **|Q| = 1+37+37+18 = 93**. The full CA alphabet is

`A = Γ × ({no head} ∪ Q) ∪ {ERROR}`,

so

**|A| = 4,096 × 94 + 1 = 385,025**.

The implementation encodes a nonerror symbol as `data + 4096*head_code`, with `data∈0..4095` and `head_code∈0..93`; ERROR is 385024. All coordinates, integer magnitudes, exponents and scan distances belong to the configuration or its external encoding/observer, **not to a cell's alphabet state**. The operation queue is not unbounded: `(rule,operation)` ranges only over the fixed list of 37 operations.

The operation counts of the fourteen fractions are

`(3,5,3,3,3,3,3,3,1,2,2,3,1,2)`.

Every fraction is reduced and squarefree. Its operation list first decrements its denominator's prime tracks in increasing-prime order, then increments its numerator's tracks in increasing-prime order. No coordinate appears twice in one list.

## 2. Encoding all positive integers

Write uniquely

`N = R · ∏_{i=1}^{10} p_i^{e_i}`,

where R≥1 is coprime to every p_i in the ten-prime basis. Choose an origin o∈Z.

- Counter track i is 1 exactly at cells `o,o+1,...,o+e_i−1`.
- The origin marker is 1 exactly at o.
- A READY head sits at o.
- Binary digit j of R is stored as the spectator bit at `o−1−j`.
- All remaining bits and head flags are zero.

The infinite blank background is quiescent. The encoding has finite support for every finite N, but support is not uniformly bounded. Spectator R is stored even when R=1. Head motion in a valid encoded computation never enters cells left of o; spectator bits remain unchanged. In fact every defined nonerror head transition preserves both marker and spectator bits.

The external observer counts each unary track, decodes R, and multiplies. It rejects malformed encodings when used for validation. Such global validation is **not hidden in the CA's local rule**.

Because all program factors use the fixed ten-prime basis, R never changes under FRACTRAN. Therefore this spectator construction extends the simulation from basis-supported integers to **all positive integers**.

## 3. Explicit finite head transition

The head sees only its current twelve-bit tape symbol and its finite head state.

### READY at the origin

The ten bits at the origin tell exactly which counter exponents are nonzero. Since every denominator is squarefree, they determine denominator divisibility. Examine the fourteen fixed denominator masks in their printed order and select the first contained in this local mask. The last denominator is 1, so selection always succeeds. Change the head to the first SEEK state without moving or changing tape bits.

Thus the first-enabled priority test, which was nonlocal in an arbitrary register-as-spatial-cell interpretation, is **local in this specifically engineered unary-track encoding**: all ten zero/nonzero indicators meet in one origin cell. This does not identify the origin with a tournament node.

### SEEK for an increment

Scan right while the selected track bit is 1. At its first 0, change it to 1. If at the origin, advance to the next operation without moving; otherwise move left in RETURN state.

### SEEK and DECLEFT for a decrement

Scan right while the selected track bit is 1. At its first 0, move left in DECLEFT state. Clear the last 1. If now at the origin, advance to the next operation without moving; otherwise move left in RETURN state. Attempting to decrement a zero counter at the origin, or clearing a non-1 in DECLEFT, signals ERROR.

### RETURN

Move left until the origin marker is encountered. Then advance to the next operation without moving. After the last operation, advance to READY instead.

No position or scan length is retained in the head. Only the persistent origin marker determines the endpoint of a return scan.

## 4. Complete homogeneous radius-one rule

The function `local_rule(l,c,r)` applies the following finite cases:

1. If any of l,c,r is ERROR, output ERROR.
2. If the triple contains at least two heads, output ERROR.
3. If there is no head, copy c.
4. Otherwise evaluate the single head's finite transition using its own tape symbol. If that local transition is invalid, output ERROR.
5. If the head is at c, write its new tape symbol and retain the new head only for a stationary move; for a left/right move, remove the head from c.
6. If the head is at a neighboring cell and moves into c, retain c's tape symbol and attach the new head state. Otherwise copy c.

These cases define a value in A for **every triple in A^3**. They give arbitrary malformed global configurations a defined future; they do not claim that every global defect can be recognized by a radius-one inspection. In particular, separated heads, multiple far-apart origins, or missing markers can remain undetected until relevant local neighborhoods meet. Invalid Python values outside A are rejected rather than counted as alphabet symbols.

The same rule is used at every integer position. Hence, for every shift σ and every configuration x,

`CA(σx) = σ(CA(x))`.

A distinguished origin in an initial configuration does not break translation equivariance of the rule: translating the entire configuration also translates that marker.

`ca_step_full` is a literal synchronous local-rule sweep over the finite nonblank support and its radius-one neighborhood. `ca_step` is an equivalent sparse evaluator: any head/error-free triple copies its center, so only cells within radius 1 of an existing head/error need reevaluation. This optimization changes neither the rule nor the semantics. It is independently compared with literal sweeps, including malformed configurations.

## 5. Observer theorem and exact macro duration

**Theorem.** If the CA starts in the canonical READY encoding E_o(N) of any N≥1, it reaches another canonical READY encoding E_o(T(N)) after finitely many CA steps, where T is the fixed ordered PRIMEGAME integer update. No intermediate time before that endpoint is READY. Iterating these endpoints recovers the complete FRACTRAN orbit.

**Proof.** The READY transition reads the squarefree denominator guards exactly and selects their least enabled index. In a valid configuration there is one head, so the radius-one rule performs exactly the described finite-head tape move without producing ERROR. Every decremented coordinate was positive when selected. No coordinate occurs twice in the selected fraction's operation list, so it remains positive until its own decrement is executed.

A unary scan reaches the first 0 after exactly e right moves. Incrementing that boundary extends a unary block by one. For a decrement, the predecessor of the first 0 is its final1; clearing it shortens the unary block by one. All other tracks stay unchanged. RETURN finds the persistent origin, so every operation ends there. A zero-height increment stays at the origin; a height-one decrement clears the origin bit of that counter without moving left. These boundary cases, and RETURN's marker test, prevent the valid head from crossing left of the origin. Each operation preserves unary shape, the unique head, the marker and the spectator bits.

The final counter vector is exactly the initial vector plus the selected fraction's exponent difference. The last operation alone restores READY, yielding the canonical encoding of T(N). This proves the claim and supports induction over any number of macro updates. ∎

**Exact cost.** An operation on a counter with initial height e uses exactly **2e+1** CA steps, for both increment and permitted decrement. The READY choice itself costs one step. For selected fraction i, let D_i be the union of its numerator and denominator prime supports, and k_i=|D_i|. Distinctness of the touched coordinates gives

`τ_i(e) = 1 + k_i + 2∑_{j∈D_i} e_j`.

Since k_i≤5 and M=max_j e_j,

**τ_i(e) ≤ 10M+6 = O(M+1)**.

This is a per-FRACTRAN-update bound, independent of the spectator factor's bit length. It is not a bound on total prime-generation time. The timing formula and the observer relation are explicitly checked in the executable tests, not just stated in prose.

## 6. Executed checks and exact milestones

The CA is compared against an independent integer update implemented with Python's exact `Fraction` class and the literal fraction list. Every compared macro step consists of actual calls to the local CA rule; the integer interpreter does not choose or execute the CA's counter operations.

From seed 2, the first five prime outputs are:

| Prime exponent | Integer output | FRACTRAN updates | CA updates |
|---:|---:|---:|---:|
|2|4|19|130|
|3|8|69|537|
|5|32|281|2,677|
|7|128|710|7,942|
|11|2048|2,375|34,603|

`seed2_first130_microsteps.jsonl` also records every exact cell configuration through the first output 4, including the finite head state and READY observations. Only these first five prime outputs were replayed in the CA lane. The separate frozen FRACTRAN lane checked ten. The observer theorem is the general result; finite replay is corroboration, not extrapolation.

The requested seed experiments also run through the actual local CA:

- `7→1→55↔65` and `17→1→55↔65`: seven CA steps to enter the cycle; its full CA period is 10.
- Seed 55 starts directly in that period-10 CA cycle, whose READY observer has period 2.
- `25→1375↔1625`: seven CA steps of transient, then CA period 10.
- Seed 15 reaches 4 after 18 FRACTRAN updates /124 CA updates; it is already the standard seed 2 orbit's successor.
- Seed 23 reaches 2 after 11 FRACTRAN updates /71 CA updates.
- Seed 16 produces prime exponents 5,7,11 at FRACTRAN times 150,579,2244 and CA times 1543,6808,33469.

Further tests cover:

- every fraction via its denominator as an exact first-enabled witness, including the five-operation instruction 78/85;
- all **380** integers of forms p, p^3, p^2qr supported on the ten primes, including **282 multi-control** cases outside the original one-hot-control invariant;
- **240** deterministic random general inputs, each through three macro updates, with exponents 0..4, nontrivial spectator factors and translated origins;
- all **380,928** finite `(head,tape)` pairs; **57,344** signal local error, and every defined transition stays within the finite state set and preserves marker/spectator bits;
- all 4,096 head-free tape-center symbols; 10,000 random alphabet triples; 500 arbitrary sparse configurations;
- literal/sparse sweep agreement, translation equivariance, zero/decrement-one boundaries, READY without an origin, forbidden decrements, malformed DECLEFT, nearby multiple heads, absorbing error neighborhoods, and rejection of out-of-alphabet values.

The main suite executes **86,501 macro-driven CA microsteps**, plus auxiliary sweeps, and performs **872,508 explicit checks**. The scripts use explicit exceptions, not Python assertions removed by optimization. Normal Python and `python -O` produce identical results. An independent reviewer additionally reconstructed the head transition, local rule, encoder, observer, and literal dense sweep: **1,613,769 checks** passed, including **1,142,784 single-head neighborhoods**, **38,675 literal CA microsteps**, and the same prime-output milestones. Its normal and optimized results also match. The only review issue was redundant zero-valued sparse-map keys; the sparse evaluator now normalizes them, and both suites were rerun successfully. See `independent_ca_audit.py` and `independent_ca_audit_results.json`.

READY-trace SHA-256:

`05cfd6daac9efa58ece1d23e1ddecf26d01cc4804878df6a225ef190be4fc872`

Canonical microconfiguration-stream SHA-256 for seed 2 through the fifth prime output:

`e303f8152d0b78e1f28ba3fd3b607508c48e5df19cff50597d93aa19056c479f`

Run from this directory:

```sh
python -B verify_ca.py
python -B -O verify_ca.py
python -B independent_ca_audit.py
python -B -O independent_ca_audit.py
```

Only our own standard-library code was executed. No third-party software or downloaded program was run.

## 7. The resulting bridge—and its limit

There is now a concrete, tested bridge

`positive integers ↔ canonical READY encodings in an infinite 1D CA`,

with a proven variable-time observer relation to the **fixed** PRIMEGAME program. The first ten primes naturally index ten unary tracks in this encoding. The construction is not an identification of those tracks with the ten quotient nodes of a tournament metagraph. It preserves the exact arithmetic dynamics, not any unspecified tournament adjacency, color, duality or all-blue symmetry.

The large alphabet is deliberate transparency, not a minimality claim. No finite ten-cell, finite-alphabet system can injectively represent all positive integers: its state space is finite. Unbounded spatial extent is what stores unbounded exponent data here. Likewise the local CA rule does not make the arithmetic balance statistic F=U+S invariant; it faithfully reproduces the same failures and coincidences already established for PRIMEGAME.

# Session 126: PRIMEGAME, ten prime registers, five inputs, and locality

## 1. Source identification and scope

The user's **ordered** fourteen fractions are exactly Conway's PRIMEGAME in Theorem 1 of his 1987 chapter, printed p. 4. The theorem concerns the integer seed **2**: later pure powers of 2 have prime exponents in increasing order. The divisor-test explanation is on printed pp. 9–10. Some later presentations substitute a different `15/14` variant; this report does not. Bibliographic identity: J. H. Conway, *FRACTRAN: A Simple Universal Programming Language for Arithmetic*, in *Open Problems in Communication and Computation*, pp. 4–26, DOI `10.1007/978-1-4612-4808-8_2`. [Publisher record](https://link.springer.com/chapter/10.1007/978-1-4612-4808-8_2); [author's chapter scan](https://gwern.net/doc/cs/computable/1987-conway.pdf).

The calculations, implementation, instruction-word reconstruction, and cross-object distinctions below were independently checked in this session. They reconstruct a known program, not a claim to have discovered a new prime-generation theorem. Only the first ten prime outputs were replayed in full. The all-primes conclusion rests on the proof, not extrapolation from that replay.

## 2. Exact factorization: what the first ten primes mean

Let the prime basis be

`P = (2,3,5,7,11,13,17,19,23,29)`.

All numerators and denominators are squarefree, and their union of prime supports is exactly this basis. An integer state supported on this basis is encoded injectively as

`N = 2^a 3^b 5^c 7^d 11^e 13^f 17^g 19^h 23^i 29^j`.

The first ten primes here label **ten exponent registers**. Their numerical values are not ten tournament isomorphism classes or ten quotient nodes.

| Index | Label | Fraction | Factored numerator | Factored denominator | Exponent update |
|---:|:---:|---:|---|---|---|
|1|A|17/91|17|7·13|−7 −13 +17|
|2|B|78/85|2·3·13|5·17|−5 −17 +2 +3 +13|
|3|C|19/51|19|3·17|−3 −17 +19|
|4|D|23/38|23|2·19|−2 −19 +23|
|5|E|29/33|29|3·11|−3 −11 +29|
|6|F|77/29|7·11|29|−29 +7 +11|
|7|G|95/23|5·19|23|−23 +5 +19|
|8|H|77/19|7·11|19|−19 +7 +11|
|9|I|1/17|1|17|−17|
|10|J|11/13|11|13|−13 +11|
|11|K|13/11|13|11|−11 +13|
|12|L|15/2|3·5|2|−2 +3 +5|
|13|M|1/7|1|7|−7|
|14|N|55/1|5·11|1|+5 +11|

In the last column, `+p` means increment the exponent of prime p, not add p to N. The complete numeric numerator, denominator, and delta vectors are in `primegame_results.json`.

For reduced fraction u_i/v_i and exponent vector x, let α_i=ν(u_i), β_i=ν(v_i). Select

`i(x) = min{i : x ≥ β_i coordinatewise}`,

then update `x ↦ x + α_i − β_i`. This is exactly integer divisibility `v_i | N`. Priority is part of the program: an enabled later instruction cannot fire while an earlier one is enabled. For example, at N=91 both A and J are arithmetically enabled; A must give 17, whereas choosing J instead would give 77. The final denominator is 1, so this particular program never halts on a positive integer; a replay safety limit is not a halt.

Each update consumes at most two units from distinct registers and produces at most three. The largest support size is five (instruction B). Sparsity of an individual update does not make the rule a spatially local cellular automaton.

### Prime-renaming equivariance, not symmetry of the fixed rack

Choose any distinct primes q_1,...,q_10 and replace every p_i by q_i consistently in the seed and every fraction. The map `∏p_i^x_i ↦ ∏q_i^x_i` preserves every divisibility guard, its priority, and every exponent update. Hence it conjugates the two dynamics exactly. One concrete relabeling to `(31,37,41,43,47,53,59,61,67,71)` was checked for 1,000 updates.

This is a theorem about **simultaneously relabeling the program and state**, not an automorphism obtained by permuting registers while keeping the printed program fixed. Likewise it supplies no canonical map to a ten-node tournament quotient.

## 3. Four counters and seven control states

On the orbit from seed 2, and more generally whenever initially all six control exponents are zero except possibly one equal to one, one has

`e_11+e_13+e_17+e_19+e_23+e_29 ∈ {0,1}`.

Proof: each applicable control-bearing instruction replaces one such prime by another or removes it. The final instruction adds 11 only when none is present, since each control prime has an earlier applicable instruction with that prime alone as denominator. The invariant follows by induction.

Thus `(a,b,c,d)=(ν_2,ν_3,ν_5,ν_7)` are unbounded data counters and the remaining primes select one of seven control states `{1,11,13,17,19,23,29}`, with 1 meaning no control factor. The reduced finite-control table is:

| Control | Prioritized condition | Update and next control |
|---|---|---|
|13|d>0|d−−; go to17 (A)|
|13|otherwise|go to11 (J)|
|17|c>0|c−−,a++,b++; go to13 (B)|
|17|else b>0|b−−; go to19 (C)|
|17|otherwise|go to1 (I)|
|19|a>0|a−−; go to23 (D)|
|19|otherwise|d++; go to11 (H)|
|11|b>0|b−−; go to29 (E)|
|11|otherwise|go to13 (K)|
|29|always|d++; go to11 (F)|
|23|always|c++; go to19 (G)|
|1|a>0|a−−,b++,c++; remain1 (L)|
|1|else d>0|d−−; remain1 (M)|
|1|otherwise|c++; go to11 (N)|

The roles are deliberately asymmetric. In particular, 7 is a data-counter prime, 17 a control prime, and 55 a product introducing a data unit and a control flag.

A separate finite structural audit agrees on 1,792 one-hot-control states with each counter in 0..3, and integer/exponent semantics agree on all 1,024 binary exponent vectors (including states outside the one-hot invariant domain).

The 10×14 integer update matrix has full integer lattice image. Columns `(I,M,J,A,N,E,F,L,C,D)` form a 10×10 minor of determinant −1. Thus no nonzero global linear conservation law, even modulo an integer m>1, is forced by these update vectors. This is not a claim that the priority graph is strongly connected: actual nonnegative reachable states have strong constraints, as the cycles below demonstrate.

## 4. Reconstructed divisor-loop proof and exact word lengths

Use A,...,N for the fourteen instructions. At a canonical trial state

`C(n,d) = 5^n 7^d 13`, with `n≥2` and `1≤d<n`, write `n=qd+r`, `0≤r<d`.

A complete quotient-removal block has instruction word

`W_d = (AB)^d J (EF)^d K`.

Starting at C(n,d), applying `W_d^q` gives `2^(qd)5^r7^d13`. Applying `(AB)^r A` then gives

`2^n 3^r 7^(d−r−1)17`.

Every indicated instruction is the first applicable one; the control table verifies this, including exhaustion tests.

* **If r>0**, the suffix `C (DG)^n H (EF)^(r−1) K` returns to `C(n,d−1)`. No pure power of 2 occurs. The complete block uses **6n+2q+2** FRACTRAN updates.
* **If r=0**, the suffix `I L^n M^(d−1) N (EF)^n K` gives `C(n+1,n)`. Immediately after I the state is `2^n7^(d−1)`. This is a pure power of 2 exactly when d=1. The complete block uses **7n+2q+d+3** updates.

Why the prime theorem follows: seed 2 follows `L N E F K` to `C(2,1)`. For each candidate n, the program tests d=n−1,n−2,... until the first divisor. It then increments n and restarts at d=n. The first divisor is 1 precisely for prime n. The only pure powers of 2 inside these blocks are the stated `2^n` events, so no composite exponent slips through.

The full instruction-word formulas, endpoints, intermediate pure-power events, and step counts were exhaustively checked for all **1,225 pairs** `2≤n≤50,1≤d<n`, involving 267,171 individual updates. These finite checks supplement the symbolic argument; they do not replace it.

More generally, starting at `2^a` with a≥1, the bootstrap word `L^a N (EF)^a K` (length 3a+2) gives `C(a+1,a)`. Subsequent pure-power exponents are exactly the primes greater than a. The special case a=0 fails: seed 1 is trapped as shown below.

## 5. Independently checked replay

Three implementations agree at every step: exact integer divisibility/multiplication; separately hand-coded exponent consumption/production; and the seven-state four-counter machine above. Their common trace from seed 2 is:

| Prime exponent | Pure power | Update count (seed has count0) |
|---:|---:|---:|
|2|4|19|
|3|8|69|
|5|32|281|
|7|128|710|
|11|2048|2375|
|13|8192|3893|
|17|131072|8102|
|19|524288|11361|
|23|8388608|19268|
|29|536870912|36981|

The largest intermediate integer has 150 bits, much larger than the final output. The 36,981-step trace records every exact integer, exponent vector, and next chosen instruction. The verifier uses explicit exceptions, not removable Python `assert` statements. Both ordinary Python and `python -O` gave identical outputs and hashes. The main verifier makes 116,924 explicit checks; the structural audit adds 4,638. Ordinary and optimized structural-audit reports also agree. A separate reviewer independently recomputed the trace with exact Fraction arithmetic, all 1,225 macro words, 567 bounded legal controller states, the determinant via Bareiss elimination, and seed/FUS checks: 385,060 explicit checks passed; see `independent_audit.py` and `independent_audit_results.json`.

`primegame_trace.tsv` SHA-256:

`f1cc30acbc44364c68c337090dc66c06e81be1d85b62f32828e284e8d2a1f8af`

`primegame_results.json` SHA-256:

`9fdfa8c33201a819347f49a093ace37cc3635a46fd6930096070844970ee2df0`

Run:

```sh
python verify_primegame.py
python -O verify_primegame.py
python input_seed_checks.py
python structural_audit.py
```

These are our own standard-library-only scripts. No external downloaded executable or third-party program was run.

## 6. The requested five-set and {7,17,55}: integer inputs are not object labels

If, as a separate experiment, the numbers are used as literal **initial integers N**, their outcomes differ sharply:

| Seed | Exact result |
|---:|---|
|15|Already the seed 2 orbit after its first update; later prime exponents 2,3,5,...|
|16=2^4|Bootstrap begins testing n=5; later prime exponents 5,7,11,...; first three at updates 150,579,2244|
|17|17→1→55→65→55, a proven eventual 2-cycle|
|23|23→95→385→455→85→78→66→58→154→182→34→2 in 11 updates; then joins the standard orbit|
|25|25→1375→1625→1375, a proven eventual 2-cycle|

For the special triple:

`7 → 1 → 55 ↔ 65`,

`17 → 1 → 55 ↔ 65`,

`55 ↔ 65`.

Thus all three initial integers enter the **same** 2-cycle, not three prime-producing symmetry states. The equality is exact, not numerical evidence for an undetected long-term outcome. More generally, seed `5^c`, c≥0, enters the cycle `5^(c+1)·11 ↔ 5^(c+1)·13`.

This does not invalidate any tournament-node labels 7,17,55. It says that a dictionary is required before substituting numeric labels into the arithmetic machine. A count, a bitmask labeling a tournament, a prime register, and a FRACTRAN initial integer are different mathematical types.

## 7. A precise comparison with the recovered F=U+S divisor taxonomy

Under the existing repository's nonunit proper-divisor definitions, for N≥2,

`F(N)=τ(N)−2`,

`U(N)=ω(N)−[N is prime]`,

`S(N)=2^ω(N)−1−[N is squarefree]`.

The recovered balanced shapes are `p`, `p^3`, and `p^2qr` with distinct primes p,q,r. These are **integer factorization shapes**, not the seven PRIMEGAME control states and not its ten storage-register labels.

An exact comparison follows immediately. For a PRIMEGAME output `N=2^ℓ`, where ℓ is prime and ℓ≥2,

`(F,U,S)=(ℓ−1,1,1)`, hence `F−U−S=ℓ−3`.

Consequently **8=2^3 is the only balanced pure-power output**. The prime-exponent mechanism does not preserve F=U+S. Indeed, the very first step takes the balanced prime 2 to 15 with defect−2. The state 825=3·5^2·11 is balanced again; balance is not a dynamical invariant.

The primes 7 and 17 have `(F,U,S)=(0,0,0)`, while 55=5·11 has `(2,2,2)`, defect−2. Nevertheless, all three seeds enter the same cycle. This supplies a concrete obstruction to identifying the balance taxonomy with PRIMEGAME's basins or control states. Even within the balanced shapes, seed 343=7^3 enters55↔65, whereas seed 60=2^2·3·5 is on the standard prime-producing orbit (it follows the output8).

## 8. What “think cellular automata” can mean rigorously

**Exact current object.** PRIMEGAME is a deterministic, prioritized vector-addition system over nonnegative exponent vectors, with a finite-control four-counter realization on the indicated invariant domain. It has no declared lattice, translation-invariant neighborhood rule, finite cell alphabet, or simultaneous spatial update. The first-enabled condition includes negative guards saying all earlier rules are disabled. Independent simultaneous firing changes the program.

**A constructive local compilation is possible, but is a different object.** Encode the four counters as four unary tracks on a one-dimensional tape, with a marked origin and one finite-state head. Zero tests inspect the first counter cell; increment/decrement scan to the unary boundary, change one bit, and return. The control table is finite, so this defines a finite-state tape machine performing one PRIMEGAME transition through finitely many tape moves. Encode a tape symbol with or without the unique head as a finite alphabet. Each head move can be simulated by a synchronous radius-one CA: the head cell writes its new symbol, and the appropriate adjacent cell receives the new head state. On invalid multi-head neighborhoods, define arbitrary deterministic conflict rules to make the CA total. Valid single-head configurations remain valid.

This is a precise simulation construction: the unary tape encodes unbounded counters, and the finite head implements the control/priority logic. It generally uses many CA microsteps per FRACTRAN step, with input-dependent overhead. We have **not** supplied or tested a particular optimized CA rule table, a Conway-Life implementation, or a tournament-metagraph embedding. A ten-cell finite-alphabet CA alone cannot faithfully store all unbounded exponent vectors; a spatially expanding or infinite configuration is needed.

A separate modern primary source formalizes the broader Minsky-machine-to-FRACTRAN compiler and computability framework: Dominique Larchey-Wendling and Yannick Forster, *Hilbert's Tenth Problem in Coq (Extended Version)*, LMCS 18(1), 2022, §6.2. That is a result about the FRACTRAN language; it does not say the fixed PRIMEGAME rack itself is a universal program. [Journal article](https://lmcs.episciences.org/9153).

## 9. Honest bridge status

Established: exact fraction support on ten primes; priority exponent dynamics; four-counter/seven-state mechanism; prime-renaming conjugacy; divisor-loop proof; reproducible first-ten output trace; requested seed experiments; the failure of F=U+S to be invariant; a local-simulation construction.

Not established: a canonical identification of the ten prime registers with the ten merged tournament nodes, an identification of the triple7/17/55 with control states, preservation of tournament adjacency or color under PRIMEGAME, or a CA law native to that metagraph. Such a bridge needs an explicit typed map and a verified commuting relation, not cardinality matching.

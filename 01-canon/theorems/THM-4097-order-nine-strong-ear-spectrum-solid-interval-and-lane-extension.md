---
id: THM-4097
title: "Order-nine strong-ear spectrum, solid interval, and prime-lane extension"
status: >
  PROVED ELEMENTARY REDUCTION + FINITE-EXACT + VERIFIED-EXACT +
  INDEPENDENTLY CROSS-AUDITED. The exact strong Hamiltonian-path spectrum at
  order nine has 1,482 values and contains the solid odd interval 85 through
  2,881. Combined with THM-1370 and THM-4094, every allowed odd value through
  2,881 is attained, the ordinary strong-prime lane is proved through 2,879,
  and the strong seven-prime lane is proved through 7*409. The first values
  not supplied at this theorem's scope are 2,887 and 7*419=2,933. THM-4102
  and THM-4104 subsequently move the current targets to 80,407 and
  7*11,527=80,689. The global H-spectrum conjecture remains OPEN.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
related:
  - THM-012b-insertion-decomposition
  - THM-1862-order-join-reduction-principle
  - THM-4051-tournament-order-seven-strong-base-exact-frontier
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4102-selected-order-ten-strong-ear-solid-interval
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - HYP-2879-strong-ear-atom-calculus
script: 04-computation/tournament_order9_strong_ear_spectrum_thm4097.py
output: 05-knowledge/results/tournament_order9_strong_ear_spectrum_thm4097.out
independent_audit_script: 04-computation/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.out
script_sha256: 610ca5850b272e0e75c574f2c1a710a0b96c75cc7191b1e1f1a03dfbdd1378d6
output_sha256: 0c3c9690ad5877d86480693af7ce97d8936c90e21b65677eefec72234c933dc0
independent_audit_script_sha256: b58de51efa200374e6014c8aeace4086fc7df1b2b73db36a01c0fcc4d2dd7943
independent_audit_output_sha256: 71058bc6b31ba26a59a58bb7f1e5366e767ed43db1f4bb2ddb80642b02acbfa6
certificate_atlas_script: 04-computation/order9_strong_ear_interval_certificate_codex_20260825.py
certificate_atlas_audit_script: 04-computation/order9_strong_ear_interval_independent_audit_codex_20260825.py
certificate_atlas: 05-knowledge/results/order9_strong_ear_interval_certificates_codex_20260825.tsv
certificate_atlas_output: 05-knowledge/results/order9_strong_ear_interval_certificate_codex_20260825.out
certificate_atlas_audit_output: 05-knowledge/results/order9_strong_ear_interval_independent_audit_codex_20260825.out
certificate_atlas_script_sha256: e9853437725f3d49fa10e75fc0d01a6c41b4e6ec4b5d71b909b9785ed23f0f36
certificate_atlas_audit_script_sha256: f8347ff0fb3e3412ab0b87802074da68e85c7feb589a676374f2c0592f81fd8d
certificate_atlas_sha256: 4aa1304da5c70e401c4953707820d377aeaf2418638dce747246c383309b3136
certificate_atlas_output_sha256: 1fb90b6e3c622c46cd99dfff7f8ced8312b740a6ac74d7b23e93d99082a3db12
certificate_atlas_audit_output_sha256: 7fe896cebf2e5641e07be6de56a0d25a0b952c9dd752fd5c7a93279aa625f975
engine: 04-computation/strong_H_spectrum_m8_isoclass_monad_s5.py
engine_sha256: 6ab922de4a8b6f6c15ee0ca7e0b036c3821b3e800dbdf961de72194e73346419
historical_histogram: 05-knowledge/results/h_spectrum_n9_histogram_monad_s6.tsv
historical_histogram_sha256: e7d5594879d4c3af739cb94ca8cfd944879c4d586747d993dd6687e60126552f
historical_value_dump: 05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out
historical_value_dump_sha256: 27fbef5b06fcf0369eeb602e513c3802ea171492e1292a3f6afa3efeadef9f55
semantic_sha256: b72cfe0e2c2ca2302f54933060299fc302cc9f241ef79154aead9524b426f0d5
independent_semantic_sha256: 7c2f9f5c9a586d1eed5b0329a294b7541c1c492682cb9ecb76c2eefc2d14b167
hash_basis: raw LF bytes for files; canonical compact JSON for semantic ledgers
audit: >
  PASS. The primary path generates the complete A000568 universe through
  order eight, evaluates all 1,526,032 nonconstant strong ears by a subset-DP
  boundary polynomial, and directly checks one labelled witness for every
  retained value. Two older order-nine canonical enumerations give the same
  1,482-value set. A later certificate atlas freezes and independently
  reconstructs one parent/cut/child witness for every one of the 1,399 odd
  values in the solid interval. The independent path compares those frozen universes and
  checks eight boundary witnesses by both Held--Karp DP and literal
  permutation enumeration. Normal/-O streams byte-match frozen outputs.
---

# THM-4097 -- the order-nine strong-ear interval reaches 2,881

**PROVED ELEMENTARY REDUCTION + FINITE-EXACT + VERIFIED-EXACT +
INDEPENDENTLY CROSS-AUDITED.** This theorem recovers an exact historical
order-nine computation, replaces its raw enumeration by the composable ear
boundary suggested by THM-4094, and promotes the resulting interval into the
current H-spectrum proof graph.

It is a finite theorem. It neither proves the all-order strong-ear conjecture
nor the global H-spectrum conjecture.

## 1. The exact ear boundary

Let `T` be a tournament on vertex set `V`, and adjoin a vertex `x`. Write

```text
sigma(v)=1  iff  x -> v,
sigma(v)=0  iff  v -> x.                              (1)
```

For `b,a in V`, define

```text
Start_T(b) = #{Hamiltonian paths of T starting at b},
End_T(a)   = #{Hamiltonian paths of T ending at a}.    (2)
```

Let `Q_T(a,b)` count permutations of `V` in which `a` is immediately followed
by `b` and every other adjacent step is an arc of `T`. The exposed step
`a,b` itself is allowed either orientation. Then

> **Lemma 1.1 (exact strong-ear polynomial).**
>
> ```text
> H(T+x)
>   = sum_(sigma(b)=1) Start_T(b)
>   + sum_(sigma(a)=0) End_T(a)
>   + sum_(sigma(a)=0, sigma(b)=1) Q_T(a,b).           (3)
> ```

### Proof

Partition a Hamiltonian path of `T+x` by the position of `x`. If `x` is first,
the suffix is a path of `T` starting at some `b` with `x->b`. If `x` is last,
the prefix ends at some `a` with `a->x`. If `x` is internal between `a,b`,
then `a->x->b`, while deleting `x` leaves a permutation whose other adjacent
steps are valid; this is exactly one term of `Q_T(a,b)`. The three cases are
disjoint and exhaustive. QED.

The script computes `Q_T` without permutation enumeration. If
`E_T(A,a)` is the number of Hamiltonian paths of `T[A]` ending at `a`, and
`S_T(B,b)` is the number starting at `b`, then

```text
Q_T(a,b)
 = sum_(A disjoint-union B=V; a in A, b in B)
     E_T(A,a) S_T(B,b).                               (4)
```

Thus the full sidecar is the pair of subset path tables, or equivalently the
contracted boundary `(Start,End,Q)`. A scalar insertion count or degree profile
does not compose this contraction.

> **Lemma 1.2 (strongness).** If `T` is strong and `sigma` is nonconstant,
> then `T+x` is strong.

Indeed `x` has an out-neighbor and an in-neighbor in `T`; strongness of `T`
then supplies directed paths from `x` to every old vertex and back. QED.

Conversely, Moon's nonseparating-vertex theorem, already cited and used in
THM-1370-h-spectrum-omits-7-21-all-n, says every strong tournament of order at
least four has a vertex whose deletion remains strong. Therefore every strong
order-nine tournament is a nonconstant ear over a strong order-eight
tournament. Lemmas 1.1--1.2 turn a complete order-eight representative bank
into an exact order-nine strong-spectrum computation.

## 2. Complete order-nine strong spectrum

Let

```text
S_str(9)={H(T):T is a strong tournament on nine vertices}.            (5)
```

> **Theorem 2.1 (exact spectrum and solid interval).** The set `S_str(9)` has
> `1,482` elements, minimum `75`, maximum `3,357`, and
>
> ```text
> S_str(9) = {75,81}
>            union {85,87,...,2881}
>            union E_9,                                (6)
> ```
>
> where the solid interval has `1,399` values and the `81`-value tail is
>
> ```text
> E_9={2885,2891,2903,2917,2919,2931,2959,2961,2963,2967,
>      2971,2973,2979,2981,2983,2985,2987,2989,2991,2993,
>      2995,2999,3001,3003,3005,3007,3009,3011,3013,3015,
>      3017,3019,3021,3023,3025,3027,3029,3031,3033,3035,
>      3037,3039,3041,3043,3045,3047,3049,3051,3053,3055,
>      3057,3059,3061,3063,3065,3067,3069,3071,3073,3075,
>      3079,3081,3083,3085,3095,3101,3105,3109,3159,3243,
>      3249,3255,3267,3275,3279,3281,3297,3299,3303,3333,
>      3357}.
> ```

In particular `2,883` is absent at order nine, while `2,885` is present.

### Exact universe and computation

The primary path generates all nonisomorphic tournaments through order eight.
Its counts

```text
1,1,2,4,12,56,456,6880                              (7)
```

match `A000568`. Exactly `6,008` order-eight representatives are strong. The
script evaluates every one of their `2^8-2=254` nonconstant signatures:

```text
6008*254 = 1,526,032 strong ears.                     (8)
```

It retains one labelled child for each value and recomputes `H` and strongness
directly on all `1,482` retained witnesses. There are zero failures.

The resulting set agrees exactly with both historical frozen paths:

- the full `191,536`-class order-nine histogram, with `178,133` strong
  classes; and
- a separately generated strong-value dump.

Those paths enumerate order nine canonically rather than passing through the
ear polynomial. The independent companion additionally checks eight boundary
witnesses by literal enumeration of all vertex permutations.

## 3. Explicit recovered atoms

Encode a labelled tournament of order `n` by

```text
code(T)=sum_(0<=i<j<n) [i->j] 2^b(i,j),               (9)
```

where `b(i,j)` is the lexicographic index of the pair `(i,j)`. Direct DP and
literal enumeration both give

```text
n=8: code 251585423    has H=613 and is strong,
n=9: code 63853559615  has H=623 and is strong.       (10)
```

Thus the two targets left immediately beyond THM-4094's old finite prefix are
not merely attained: they are explicit strong atoms. Boundary controls at the
new solid interval and tail are

```text
H       n   code
75      9    9897508671
81      9   10165944127
85      9   49392123839
2881    9   34308081455
2885    9   63870568207
3357    9   68164491031.                              (11)
```

## 4. Global and prime-lane consequences

THM-1370 proves that every odd value through `609`, except `7,21`, is attained.
The solid interval `(6)` overlaps that prefix and continues through `2,881`.
Therefore

> **Corollary 4.1 (global finite prefix).** Every positive odd integer at most
> `2,881`, except `7` and `21`, is a tournament Hamiltonian-path count.

For an attained odd prime `p`, THM-4094 Lemma 5.1 extracts a strong component
of value `p`. Its Lemma 5.2 does the same for attained `7p`, using the proved
omission of `7`. Combining those lemmas with `(6)` and the old prefix gives

> **Corollary 4.2 (extended strong lanes).**
>
> ```text
> {p odd prime: p<=2879, p!=7} subset S_str,
> {7p: p odd prime, p<=409, p!=3} subset S_str.         (12)
> ```

The first ordinary prime not supplied by the current finite certificates is
`2,887`. The next exceptional prime is `419`, giving `7*419=2,933`. Neither is
claimed absent; both are the new first unforced targets.

Consequently THM-4094 Theorem 6.1 sharpens to the exact remaining reduction

```text
global H-spectrum completeness
iff
  every odd prime p>=2887 is strongly attained, and
  every 7p with odd prime p>=419 is strongly attained.               (13)
```

Equation `(13)` is the reduction at order-nine scope. THM-4102 first moves
its two lower bounds to `14,657` and `2,111`; THM-4104 currently moves them
to `80,407` and `11,527`, respectively.

The sporadic tail `E_9` supplies some atoms beyond these thresholds, but does
not move the first missing member of either lane.

## 5. Failure boundary and scope

The mechanism is a boundary contraction, not raw density. A complete
order-eight bank plus the exact `(Start,End,Q)` sidecar produces a long solid
interval; the scalar `H(T)` of a parent does not determine its ear image.
Likewise the first post-interval hole `2,883` is an order-nine boundary fact,
not evidence of a permanent global gap.

No LRC theorem follows from the numerical interval. The legitimate transfer
is only methodological: retain the full composable boundary before optimizing
or projecting. LRC(14), the all-order strong-ear interval, and the global
H-spectrum conjecture remain **OPEN**.

## 6. Reproduction

From the repository root:

```bash
python3 -B 04-computation/tournament_order9_strong_ear_spectrum_thm4097.py
python3 -B -O 04-computation/tournament_order9_strong_ear_spectrum_thm4097.py
python3 -B 04-computation/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.py
python3 -B -O 04-computation/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.py
```

Each normal/optimized pair must byte-match its frozen output. Every executable
gate uses `require`; optimization removes no check.

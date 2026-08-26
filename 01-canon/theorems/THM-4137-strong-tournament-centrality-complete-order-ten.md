---
id: THM-4137
title: "Strong tournament Johnson support-floor centrality at order ten"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every one of the
  9,355,949 strong tournament isomorphism classes of order ten has the unique
  central layer t=0 as optimizer of both the THM-4128 rational Johnson
  support floor and the THM-4123 exact-coset floor. The minimum strict
  central-over-outer coset margin is 24. Actual response maxima are
  noncentral-only in 3,146,972 classes. Together with THM-4135 this proves
  centrality through order ten; order eleven remains open and THM-4133
  refutes the all-order extension at order twelve.
source: codex-frontier-synthesis-creative-20260825v
depends_on:
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
related:
  - THM-4131-strong-tournament-centrality-through-order-eight
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4135-strong-tournament-centrality-complete-order-nine
script: 04-computation/tournament_strong_centrality_order_ten_thm4137.py
output: 05-knowledge/results/tournament_strong_centrality_order_ten_thm4137.out
independent_audit_script: 04-computation/tournament_strong_centrality_order_ten_thm4137_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_strong_centrality_order_ten_thm4137_independent_audit.out
script_sha256: 45a5b18d1d88686054b6328ea93c068df3d59476111aa8b6f0d48f93e8028581
output_sha256: 715b90f8a5257af0a81fd26c847c3ce1d2bb4182eecc2cd093af9c93fb977a52
independent_audit_script_sha256: 8d531fcad7e0149bb6ccb70e876d504b0912341321e9cde0f08e9746b8286d19
independent_audit_output_sha256: c2c87368b8cc2db63682fcf0c78578ba48e6b1539992c1095a557b499b843eb2
semantic_sha256: 39a9adba9b2c4d3917e07ef24df62fc02d9203c53e32b8d0de7f02b61e66c114
independent_semantic_sha256: 2984d66aad38b7b7a5fa4de0ec14e6e08ff383bf277eaa2ae4290dca1841f045
primary_certificate_sha256: c57ff382a892fd7136714605ed73940ef7bb0dba39af4cf5e98ff816ee5ffa6c
independent_profile_sha256: 8b2d93cd83ccafa86315b540cbeba5338812ad4d6393712bd403163b954ebab9
hash_basis: raw LF bytes
primary_audit: >
  PASS. A pinned nauty 2.9.3 stream is split into sixteen exact shards. A
  contracted good/reversed-block kernel retains every one of the 1,024 ear
  responses, every rational and coset layer, and actual maxima. The full
  9,355,949-class regeneration byte-matches the frozen shard certificate;
  normal, optimized, and hash-seeded fast certificate replays also agree.
independent_audit: >
  ACCEPT. A clean-room C++ referee starts from all 9,733,056 unlabeled
  tournaments, performs its own two-way-reachability strong filter, and
  recovers the independently generated strong-stream count and digest. Two
  complete eight-thread passes agree on all 1,024-value ordered profiles.
  One/eight-thread and O0/O3 shard replays match; both extrema pass literal
  eleven-vertex child DP and an explicit reversal isomorphism; sanitizer
  hostiles pass.
---

# THM-4137 -- complete strong centrality at order ten

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4135 proves the positive order-nine boundary, while THM-4133 supplies a
strong failure at order twelve. The complete order-ten universe remains on
the positive side, but it comes much closer to the rational boundary than
order nine and sharply separates support-floor centrality from actual
response maxima.

## 1. Statement

For a strong tournament `T` of order ten, retain

```text
F_T(S)=H(T+x_S),
J_m=E_m F_T+Var_m(F_T)/(E_m F_T-H),
L_m=ceil_(a_(T,m)+d_(T,m)Z)(J_m),
t=10-2m.                                                   (1)
```

Here `J_m` is the THM-4128 rational Johnson support floor and `L_m` is the
THM-4123 exact layer-coset floor. The unique central layer is `m=5`, or
`t=0`.

> **Finite order-ten theorem.** For every strong tournament of order ten:
>
> 1. `t=0` is the unique layer maximizing `J_m`;
> 2. `t=0` is the unique layer maximizing `L_m`; and
> 3. its exact-coset value is strictly larger than every outer-layer value.

The complete universes have sizes

```text
all tournament isomorphism classes:       9733056,
strong tournament isomorphism classes:    9355949.         (2)
```

There are zero rational failures, zero exact-coset failures, and zero
rational/coset optimizer reorderings. Both optimizer histograms are exactly

```text
t=(0): 9355949.                                            (3)
```

The minimum strict central-minus-outer coset margin is

```text
24.                                                         (4)
```

## 2. Sharp rows

The margin in `(4)` is attained by one reversal pair with packets

```text
(H,W,D_4,C_hd)=(765,13095,50202180,+-13068360),
theta=+-36301/39843.                                      (5)
```

Both support floors optimize at `t=0`, while the actual maxima occur at
`t=+-2`.

The largest normalized tilt is attained by a different reversal pair:

```text
(H,W,D_4,C_hd)=(1431,19557,111975711,+-29570784),
theta=+-34499248/37325237,
rho=34499248/37325237.                                   (6)
```

Again both floors optimize uniquely at `t=0`, but the actual maxima occur at
`t=+-4`. Thus the order-ten extremum uses more than ninety-two percent of the
rational even-order centrality budget without crossing it.

## 3. Actual maxima remain separate

The actual optimizer histogram is

| optimizer tuple `t` | classes |
|:---|---:|
| `(-4)` | `22,353` |
| `(-4,-2)` | `321` |
| `(-2)` | `1,550,812` |
| `(-2,0)` | `11,172` |
| `(0)` | `6,186,633` |
| `(0,2)` | `11,172` |
| `(2)` | `1,550,812` |
| `(2,4)` | `321` |
| `(4)` | `22,353` |

Consequently

```text
3146972                                                       (7)
```

classes have no central actual maximizer. This is not a defect in the
support-floor theorem: `J_m` and `L_m` are certified lower support values,
not locators of the actual maximum of `F_T`.

## 4. Complete primary universe

The pinned nauty 2.9.3 generator has raw stream fingerprints

```text
all:    af2873154068897522bc15477d989b0577877d2bbec08aea3082353e5378b67
strong: 47bcaa3ef6272261dee3092735b47b3d2154d882aae6a8420c964cd3ef7289b7. (8)
```

The primary audit evaluates `gentourng -q -c 10` in sixteen residue shards.
For every class it retains the full `1,024`-response table, exposure packet,
rational values, layer lattice, exact-coset floors, and actual maxima. The
ordered shard rows have frozen aggregate digest

```text
c57ff382a892fd7136714605ed73940ef7bb0dba39af4cf5e98ff816ee5ffa6c. (9)
```

A complete eight-worker regeneration takes about four minutes on the audited
host and agrees byte-for-byte with `(9)` and the frozen output.

## 5. Independent all-class referee

The independent audit does not read the primary artifact. It begins with the
**all-class** stream `gentourng -q 10` and performs its own forward/reverse
reachability filter. It recovers `377,107` nonstrong and `9,355,949` strong
classes, with the same strong-stream digest in `(8)`.

Its contracted adjacent-good/reversed-block engine independently recomputes
all `1,024` responses and every exact layer. The ordered complete-profile
digest is

```text
8b2d93cd83ccafa86315b540cbeba5338812ad4d6393712bd403163b954ebab9. (10)
```

The referee makes two complete passes. One/eight-thread and O0/O3 shard
outputs are byte-identical. It also reconstructs all `1,024` eleven-vertex
children of both rows in `(6)` by a literal Hamilton-path DP, finds an
explicit reversal isomorphism, and passes an ASan/UBSan hostile run. Codes
`2,140,20` preserve the inherited rational, coset-rounding, and
support-versus-actual controls.

## 6. Scope and replay

Together with THM-4131/4135, this theorem proves rational and exact-coset
Johnson support-floor centrality through order ten. THM-4133 refutes the
all-order statement at order twelve. Order eleven remains open. No theorem
about actual-maximizer centrality, response intervals, `H>=disc`, or LRC
follows.

The fast primary certificate and full regeneration are

```text
python3 -B 04-computation/tournament_strong_centrality_order_ten_thm4137.py
python3 -B 04-computation/tournament_strong_centrality_order_ten_thm4137.py --full --workers 8
```

Compile and run the independent referee with

```text
clang++ -O3 -std=c++17 -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/openssl/include \
  04-computation/tournament_strong_centrality_order_ten_thm4137_independent_audit.cpp \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  -L/opt/homebrew/opt/openssl/lib -lcrypto \
  -o /tmp/thm4137_independent_audit
gentourng -q 10 | env OMP_NUM_THREADS=8 /tmp/thm4137_independent_audit
```

Both complete routes reproduce the frozen invariants above. **QED.**

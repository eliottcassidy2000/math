---
id: THM-4168
title: "Prime order-eleven nontrivial-automorphism Johnson centrality"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every prime
  tournament of order eleven with a nontrivial automorphism has central
  THM-4128 rational and THM-4123 exact-coset Johnson support-floor optimizers.
  The complete union contains 12,155 isomorphism classes. Its largest strict
  rational load is 201109364/381606521 and its minimum central-over-outer
  exact-coset margin is 852. The Johnson-centrality question for prime
  asymmetric order-eleven tournaments remains open.
source: codex-frontier-synthesis-creative-20260826at
depends_on:
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
related:
  - THM-4137-strong-tournament-centrality-complete-order-ten
  - THM-4144-order-eleven-large-homogeneous-module-johnson-centrality
  - THM-4163-order-eleven-homogeneous-pair-johnson-centrality
  - THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy
census_script: 04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_census.cpp
rebuild_script: 04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_rebuild.sh
canonical_digraph6: 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.d6
canonical_labels: 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.labels
rebuild_output: 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168_rebuild.out
evaluator_script: 04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_evaluator.cpp
evaluator_output: 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168_evaluator.out
artifact_audit_script: 04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_artifact_audit.py
artifact_audit_output: 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168_artifact_audit.out
independent_audit_script: 04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168_independent_audit.out
census_script_sha256: 1201bc580a5647068a9ed7e5e752aa53c358085de3b27feb5852dd937a98399f
rebuild_script_sha256: 24e6ff6bf10e09261f7aa4a1a5f42856fac0daa1d481c7b023377ac61e90e38b
evaluator_script_sha256: 56f726d08d40d3390818331d916ed43982dc94d52cee1f7e8453a8e4622c2902
artifact_audit_script_sha256: 08c73d5d230de217e1e41dc5c0d46c693ed81fb5e837e9e79431959b5e003d87
independent_audit_script_sha256: 7afd23fbe250fe1db8f805c22dc0f2dfae383fa44767dedc89cf573c32121d37
canonical_digraph6_sha256: b0e87b8faf6238275aa238fbc7d4f2b8e5d8c47c7b3d38d06e19ef0a8d49e2cd
canonical_labels_sha256: 17062afaf0be31fef492a0e3d1c5ef810b1997713a7a30c75d9997b1602fbace
rebuild_output_sha256: aaa31cb52cf84b07d99690e84ed64cb684f886df6b7dc98bad3d198df161f831
evaluator_output_sha256: 6512fa2abe9c3b53b80992355ceea81afe3d9c7edeb1a2b0738c722eb1220879
artifact_audit_output_sha256: 872f97cbddd737e7e3f4eaf3a629070fcf75bc0d115d1e469ae481544b75bce9
independent_audit_output_sha256: a23df4cfdd5bcff44ee2614a599f9e8e4f4eadd4f058ec95ffb2de315c902029
hash_basis: raw LF bytes
primary_audit: >
  PASS. A fixed-automorphism edge-orbit generator, nauty 2.9.3
  canonicalization, and exact cross-product evaluator rebuild the complete
  12,155-class union. The strict rational and exact-coset scans have zero
  failures. The canonical digraph6 and bit-label streams are frozen and
  rowwise cross-checked; Python normal/optimized/hash-seeded replays and Apple
  Clang 17.0.0 -O0/-O3 scans agree. No independent GCC replay is claimed.
independent_audit: >
  ACCEPT. A clean rebuild independently verifies the four-type completeness
  cover, corrects a stale raw 3x2 count, recovers every canonical set and
  automorphism distribution, and byte-matches the frozen union. A separate
  literal all-ear child DP verifies both converse-paired extrema and the
  asymmetric hostile without using the endpoint-capacity evaluator.
portability_audit: >
  PASS. The census now imports the standard algorithm header required by
  std::sort; strict GCC compilation succeeds. Every C++ path forces binary
  stdout on Windows, the Python label path forces LF, and exact attributes
  keep source, digraph6, and label bytes LF. No full GCC/nauty rebuild is
  claimed on this host.
---

# THM-4168 -- nontrivial-symmetry prime order-eleven centrality

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4144 and THM-4163 close the decomposable order-eleven strata, leaving
prime tournaments. Symmetry closes a further complete prime stratum: every
prime order-eleven tournament with a nonidentity automorphism satisfies both
Johnson support-floor centrality gates. The remaining prime problem is
therefore purely asymmetric.

This is a statement about the THM-4128 rational and THM-4123 exact-coset
**support floors**. It does not locate actual response maxima and does not
classify asymmetric prime tournaments.

## 1. Why four automorphism types are complete

Let `T` be a prime tournament of order eleven and suppose `Aut(T)` is
nontrivial. THM-2197 records the elementary fact that a tournament has no
nonidentity involutive automorphism: an involution has a transposition cycle,
whose unique arc it would send to its reverse. Since a finite group of even
order has an element of order two by Cauchy's theorem, `|Aut(T)|` is odd.
Choose an odd prime `p` dividing `|Aut(T)|`; Cauchy's theorem supplies an
automorphism `sigma` of order `p`. Since `Aut(T)` embeds in `S_11`, `p<=11`.

The cycles of `sigma` have length one or `p`. If it had exactly one
nontrivial cycle `C` with `p<11`,
then every fixed vertex would have a constant orientation toward all of `C`,
by `sigma`-invariance. Thus `C` would be a proper homogeneous module,
contradicting primeness.

It follows that, up to conjugacy in `S_11`, `sigma` has exactly one of the
four cycle types

```text
3^2 1^5,        3^3 1^2,        5^2 1,        11.         (1)
```

Indeed `p=3` permits two or three nontrivial cycles, `p=5` permits two,
`p=7` permits only the forbidden single proper cycle, and `p=11` gives the
whole vertex set. This proves that `(1)` is a complete list rather than a
sampled symmetry menu.

For each type, conjugate `sigma` to the fixed representative in the census
source. A tournament fixed by that representative is determined by one
orientation bit on every orbit of unordered edges. The source enumerates all
such bit assignments and applies a literal all-subset homogeneous-module
test. Hence every prime tournament with a nontrivial automorphism appears in
at least one generated type set.

## 2. Exact census

The four fixed-permutation universes and their prime invariant rows are

| cycle type | edge-orbit bits | assignments | prime invariant rows |
|:---|---:|---:|---:|
| `3^2 1^5` | 25 | `33,554,432` | `15,563,520` |
| `3^3 1^2` | 19 | `524,288` | `285,120` |
| `5^2 1` | 11 | `2,048` | `960` |
| `11` | 5 | `32` | `32` |

These are labelled assignments fixed by one chosen permutation, not
isomorphism-class counts and not a disjoint partition. In particular, the
first raw count corrects the stale value `10,375,680`, which is exactly the
sum of eight, rather than all twelve, nonempty modulo-sixteen shard counts.
The twelve residues `2`
through `13` each contain `1,296,960` prime rows; residues `0,1,14,15` are
empty.

Default nauty 2.9.3 canonicalization gives

```text
type 3^2 1^5:      10808 classes,
type 3^3 1^2:       1320 classes,
type 5^2 1:           24 classes,
type 11:               4 classes.                        (2)
```

The type sets overlap in exactly one class, between `5^2 1` and `11`.
Therefore their union has

```text
10808+1320+24+4-1=12155                                  (3)
```

isomorphism classes. The rebuilt sorted digraph6 stream has SHA-256

```text
b0e87b8faf6238275aa238fbc7d4f2b8e5d8c47c7b3d38d06e19ef0a8d49e2cd. (4)
```

An independent rowwise module test verifies that all `12,155` classes are
prime; each is strong and has a nontrivial automorphism. Their automorphism
group orders are distributed as

```text
|Aut|=3:  12128 classes,
|Aut|=5:     23 classes,
|Aut|=11:     3 classes,
|Aut|=55:     1 class.                                   (5)
```

## 3. Rational and exact-coset gates

For a class `T`, let `c` be its integer exposure-capacity tensor and retain

```text
C=sum_i h_i(c)d_i(c),
D=sum_(e<f, e intersect f=empty)c_ec_f.                   (6)
```

These are the evaluator's `Chdx4` and `D4x4`. At order eleven, THM-4128's
strict rational central-only gate is

```text
2|C|<D.                                                    (7)
```

The evaluator uses the exact cross product in `(7)`; floating point is
display-only. Across all `12,155` classes there are zero rational failures.
The maximum normalized load is

```text
max_T 2|C|/D
 =7239937104/13737834756
 =201109364/381606521
 =0.527007147238975... .                                  (8)
```

It is attained by exactly one converse pair, with labels

```text
0100000000000000000000000011011000110100011001010110101
1000000100111000000101100001101000011000100001000000101. (9)
```

Both have

```text
H=6507, W2=211122, D=13737834756,
C=+-3619968552, exact-coset margin=1710.                  (10)
```

For THM-4123's exact layer-coset floor `L_m`, define

```text
M_coset(T)=max(L_5,L_6)-max_(m notin {5,6})L_m.            (11)
```

The exact scan has zero coset failures and

```text
min_T M_coset(T)=852.                                     (12)
```

The minimum is attained by exactly one converse pair,

```text
0100000000000000000010000000000000111000101001010001101
1010000000101000000001000001010000110000100000100000101, (13)
```

with

```text
H=1701, W2=85194, D=2253854148, C=+-473918448,
2|C|/D=78986408/187821179.                                (14)
```

Equations `(7)--(14)` prove that every rational optimizer and every
exact-coset optimizer in the complete symmetry union is central. Ties within
the two central layers are allowed; no uniqueness claim is being hidden in
the word central. **QED.**

## 4. The asymmetric remainder is real

The exact prime label

```text
0000001100100001000000000000111111011101111111111000101  (15)
```

is a useful hostile to any attempt to turn self-converse symmetry or a
uniform deletion card into the missing classification. It has

```text
|Aut|=1,                       self-converse=yes,
eleven pairwise nonisomorphic cards,
prime zero-based card vertices={1,6,10},
H=243, W2=24894, D=188233412, C=0,
exact-coset margin=420.                                  (16)
```

Thus the asymmetric prime stratum is nonempty even inside the self-converse
locus, and only three of the eleven cards in this control remain prime. The
row is central; it marks the surviving scope rather than a counterexample.

## 5. Reproduction and validity controls

The complete canonical rebuild is

```text
bash 04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_rebuild.sh
```

It compiles the edge-orbit generator, runs all sixteen `3^2 1^5` shards and
the other three types, canonicalizes with `labelg`, and byte-compares the
union with the frozen stream. Generate and verify its bit-label companion by

```text
python3 -B \
  04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_artifact_audit.py \
  --emit-labels \
  > 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.labels
python3 -B \
  04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_artifact_audit.py
```

The full exact scan is

```text
clang++ -std=c++17 -O3 -Wall -Wextra -Werror \
  04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_evaluator.cpp \
  -o /tmp/thm4168_evaluator
/tmp/thm4168_evaluator --stdin \
  < 05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.labels
```

Apple Clang 17.0.0 `-O0/-O3` scans agree. On the audited host,
`/usr/bin/g++` is the same Apple Clang driver; no independent GCC replay is
claimed. Reproduce the independent five-row audit with

```text
clang++ -std=c++17 -O3 -Wall -Wextra -Werror \
  04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_independent_audit.cpp \
  -o /tmp/thm4168_independent
/tmp/thm4168_independent \
  0100000000000000000000000011011000110100011001010110101 \
  1000000100111000000101100001101000011000100001000000101 \
  0100000000000000000010000000000000111000101001010001101 \
  1010000000101000000001000001010000110000100000100000101 \
  0000001100100001000000000000111111011101111111111000101
```

It constructs all `2^11` order-twelve ears for each row in `(9)`, `(13)`,
and `(15)` by a literal Hamilton-path DP, then recovers every capacity, layer
lattice, and floor without using the primary endpoint-contraction engine. Its
frozen output reproduces all five rows.

The enumeration is complete rather than sampled, addressing MISTAKE-491's
scope failure. Both members of each extremal converse pair are retained,
addressing the plural-argmax failure mechanism in MISTAKE-071. The theorem
does not promote any asymmetric sample to a complete result. **QED.**

---
id: THM-4102
title: "Selected order-ten strong-ear solid interval"
status: >
  PROVED ELEMENTARY REDUCTION + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. A deterministic bank of one strong order-nine witness for each
  THM-4097 value has 755,820 nonconstant order-ten ears. Its exact image
  contains every odd value from 249 through 14,649 and a second solid interval
  from 15,055 through 15,551. Together with prior canon and multiplication
  this extends the allowed global prefix through 14,655 and moves the first
  unforced lane values at this theorem's scope to 14,657 and
  7*2,111=14,777. THM-4104 subsequently moves the current targets to 80,407
  and 7*11,527=80,689. This is a selected construction, not a complete
  order-ten census; the global conjecture is OPEN.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
related:
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
script: 04-computation/tournament_selected_order10_strong_ear_interval_thm4102.py
output: 05-knowledge/results/tournament_selected_order10_strong_ear_interval_thm4102.out
independent_audit_script: 04-computation/tournament_selected_order10_strong_ear_interval_thm4102_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_selected_order10_strong_ear_interval_thm4102_independent_audit.out
script_sha256: e7049a6347100e6dc54c7b6c03b299cde7dcfaca811954797971b8e7552421a8
output_sha256: 3ac25e2e154994178fa916d124b91e2c5a9b831768279a457298fb679dc4b4dd
independent_audit_script_sha256: 6dca471289cc1f76704b1444edbb4f2d78c7002a5367d11df81498d64a2fa592
independent_audit_output_sha256: b9ced2eb803abcfa3091306ce2fb43c71e941aaa7b0a8e66795d51c4a0830a96
parent_compiler_sha256: 610ca5850b272e0e75c574f2c1a710a0b96c75cc7191b1e1f1a03dfbdd1378d6
selected_parent_bank_sha256: c03c203943e734d09bee4b8818227b8f184405ce4c5092dd56d0fdb6107d528c
semantic_sha256: 2f61916d83f606766e059888809abc276a007a372381bafd1bfdfbc8d7cd2168
independent_semantic_sha256: a5a2e7a20ddc7b07d7221eaf6ac67c0f233b808b4ded7c56de1c4e3fcf944572
hash_basis: raw LF bytes for files; canonical compact JSON for semantic ledger
audit: >
  PASS. The independent path imports neither THM-4102 nor THM-4097, rebuilds
  the selected bank using only the pinned representative generator, checks
  strongness on all 755,820 children, reproduces all 7,566 values and both
  intervals, recomputes ten codes by separate DP, and literally enumerates
  nine key parents and two boundary children. Normal/-O streams byte-match.
---

# THM-4102 -- selected order-ten strong-ear solid interval

**PROVED ELEMENTARY REDUCTION + FINITE-EXACT + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**

This construction deliberately does not enumerate all order-ten tournaments.
It asks a smaller question suggested by THM-4097's retained witnesses: how
much of the next spectrum is already forced by one strong representative of
each known order-nine value?

## 1. The deterministic selected-ear bank

THM-4097 fixes an order-eight representative order and, within each parent,
the integer order on nonconstant cut signatures. For every value in its exact
strong order-nine spectrum, retain the first labelled ear encountered with
that value. Call the resulting bank `R_9`. Thus

```text
|R_9|=1,482,                                                 (1)
```

and every member of `R_9` is strong. The stored labelled words, not merely
their scalar `H` values, are part of this definition.

For every `A in R_9` and every nonconstant binary signature on its nine
vertices, adjoin one vertex with that cut. Let `E_10` be this labelled ear
bank and let `V_10^*` be its set of Hamiltonian-path counts. Then

```text
|E_10|=1,482*(2^9-2)=755,820.                              (2)
```

The exact `(Start,End,Q)` formula in THM-4097 evaluates every member, and its
strong-ear lemma proves every member strong. No canonicalization or
isomorphism claim is needed.

> **Theorem 1.1 (selected order-ten image).** The finite set
> `V_10^*` has `7,566` values, minimum `125`, maximum `15,621`, and contains
> the two solid odd intervals
>
> ```text
> {249,251,...,14649},       7,201 values,                  (3)
> {15055,15057,...,15551},     249 values.                  (4)
> ```

The first components of the selected image are

```text
{125}, {135}, [145,147]_2, [153,155]_2, [159,161]_2,
[165,171]_2, [175,231]_2, [235,245]_2, [249,14649]_2,
[14653,14655]_2, {14659}, [14663,14671]_2.                 (5)
```

Here `[a,b]_2` denotes the odd arithmetic interval. Equation `(5)` describes
the selected construction only; a missing value there is not asserted absent
from the full tournament spectrum.

## 2. Explicit boundary and lane witnesses

Use THM-4097's LSB-first lexicographic upper-pair code. Direct Held--Karp
recomputation and strongness checks give:

| `H` | selected parent `H` | cut signature | weight | order-ten code |
|---:|---:|---:|---:|---:|
| `249` | `75` | `104` | `3` | `26993059954495` |
| `2887` | `127` | `158` | `5` | `25133469073343` |
| `2933` | `119` | `27` | `4` | `34960494755519` |
| `14649` | `2575` | `271` | `5` | `13193187335727` |
| `14653` | `2517` | `267` | `4` | `17454203805215` |
| `14655` | `2393` | `263` | `4` | `17591902805567` |
| `15055` | `2741` | `170` | `4` | `25004756344591` |
| `15551` | `3081` | `77` | `4` | `32976074124951` |

In particular, the two first-unforced targets left by THM-4097, `2,887` and
`2,933`, are explicit strong order-ten atoms.

## 3. Global consequences

Theorem 1.1's primary interval overlaps the proved prefix and gives every
allowed odd value through `14,649`.
The selected bank additionally gives `14,653` and `14,655`. The intervening
value is supplied multiplicatively:

```text
14651=49*13*23,                                            (6)
```

using THM-4094's strong carry `49`, ordinary prime atoms, and order join.
Therefore:

> **Corollary 3.1.** Every positive odd integer at most `14,655`,
> except `7` and `21`, is a tournament Hamiltonian-path count.

Inherited THM-1370/4097 coverage together with `(3)` supplies every ordinary
prime through `14,649`, and `14,653` is displayed above. The same union
supplies `7p` through `p=2,089`, and the selected image also contains

```text
14693=7*2099.                                              (7)
```

Thus the strong lanes sharpen to

```text
{p odd prime:p<=14653,p!=7} subset S_str,
{7p:p odd prime,p<=2099,p!=3} subset S_str.                (8)
```

The first prime not supplied by this construction is `14,657`; the first
unsupplied exceptional prime is `2,111`, giving `7*2,111=14,777`. THM-4094's
factorization theorem consequently reduces global completeness to the
remaining tails

```text
all odd primes p>=14657, and all 7p with odd prime p>=2111. (9)
```

These are first unforced targets, not claimed gaps.
THM-4104 subsequently supersedes this numerical frontier: its selected
order-eleven image moves the two tail bounds to `80,407` and `11,527`.

## 4. Exact verification boundary

The primary script pins THM-4097's compiler and regenerates its complete
order-eight universe. It reproduces the exact `1,482`-value order-nine
spectrum, retains the first labelled witness of every value, and then checks:

1. all `1,482` selected parents by direct Hamilton DP and strongness;
2. all `755,820` nonconstant ears by the exact boundary contraction;
3. all `7,566` retained order-ten witnesses again by direct Hamilton DP and
   strongness;
4. both solid intervals and their adjacent selected-image boundaries;
5. the displayed labelled codes; and
6. the prime and seven-prime cutoffs.

The independent path imports neither THM-4102 nor THM-4097. It preserves the
definition's enumeration order by using only the hash-pinned order-eight
generator, rebuilds the selected bank with a separately structured boundary
evaluator, checks strongness on all `755,820` children, and reproduces every
value and interval. It independently recomputes all ten key codes by DP and
checks nine key parents plus two boundary children by literal permutation
enumeration.

## 5. Failure boundary and cross-frontier lesson

This is a selected lower image, not the full strong order-ten spectrum. Its
minimum, maximum, holes, and `7,566`-value count are properties of `E_10`, not
of all order-ten tournaments. The choice of actual parent word is
load-bearing: THM-4097 exhibits equal-`H` parents with different ear-response
multisets, so one cannot reconstruct `(3)` from the scalar spectrum alone.

THM-4099's squarefree gap algebra is the compositional lift of the same idea;
the present computation uses only its degree-one boundary. For LRC the lesson
is procedural only: retain owner/location and mixed response before projection.
No LRC theorem follows from this numerical interval. The global H-spectrum
conjecture remains **OPEN**.

## 6. Reproduction

From the repository root:

```bash
python3 -B 04-computation/tournament_selected_order10_strong_ear_interval_thm4102.py
python3 -B -O 04-computation/tournament_selected_order10_strong_ear_interval_thm4102.py
python3 -B 04-computation/tournament_selected_order10_strong_ear_interval_thm4102_independent_audit.py
python3 -B -O 04-computation/tournament_selected_order10_strong_ear_interval_thm4102_independent_audit.py
```

Each normal/optimized pair must byte-match its corresponding frozen output.
Every executable gate uses `require`; optimization removes no check.

---
id: THM-4104
title: "Selected order-eleven strong-ear solid interval"
status: >
  PROVED ELEMENTARY QUADRATIC-CUT IDENTITY + FINITE-EXACT + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. The deterministic bank containing the first labelled
  order-ten witness for each of THM-4102's 7,566 selected values has 7,732,452
  nonconstant order-eleven ears and a 43,251-value exact image. The selected
  image contains every odd value from 429 through 80,265 and a second solid
  interval from 80,875 through 84,259. With THM-4094 this gives every allowed
  odd value through 80,405 and moves the first unforced lane targets to 80,407
  and 7*11,527=80,689. This is a selected lower image, not a complete
  order-eleven census; the global H-spectrum conjecture remains OPEN.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4102-selected-order-ten-strong-ear-solid-interval
related:
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - HYP-2879-strong-ear-atom-calculus
  - MISTAKE-055
script: 04-computation/tournament_selected_order11_strong_ear_interval_thm4104.py
output: 05-knowledge/results/tournament_selected_order11_strong_ear_interval_thm4104.out
independent_audit_script: 04-computation/tournament_selected_order11_strong_ear_interval_thm4104_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_selected_order11_strong_ear_interval_thm4104_independent_audit.out
script_sha256: 996c6f13abd82d0b8d1e74cf7ec949b31d093d70979715440ec8caa56a379e03
output_sha256: 75437285e47f97cdf80ac5a3e816149d35ca62f157b35850ce9f5629e24370d4
independent_audit_script_sha256: c607dac3d2f223fb95e76293a02f6f90767e8603ecb64961cd981d71f6e361e1
independent_audit_output_sha256: 40131dd6d92e21080601d030099e8a9bf749ee994d0a5f38c96a1c17fc922482
selected_order9_bank_sha256: c03c203943e734d09bee4b8818227b8f184405ce4c5092dd56d0fdb6107d528c
selected_order10_bank_sha256: 2f3fbd5d7f56de24a1f08ea08585dd029c70344ef444830915b5ea0d203e4b92
selected_order10_source_sha256: c4b073e7e0ba7d965bb978635a2cd4ec60dc8acfb815d4388b32bab644d71980
selected_order11_values_sha256: 3bead7d8a2539f3c217540199e4b71b208fed526f80cd43bf435e17dd1b0c328
semantic_sha256: 11241b1a8d55d9a1b725b2343b0cf8543397a48a2eb9ba102bb3132b8b020067
independent_semantic_sha256: 77a73efc867036bcb662dc5974a704e4fcda9fe8270c778cf032e87f2ba95a13
hash_basis: raw LF bytes for files; canonical compact JSON for semantic ledger
audit: >
  PASS. The primary reconstructs the inherited banks, compares the quadratic
  evaluator with the original boundary formula on 2,335,376 rows/controls,
  and scans all 7,732,452 ears. The independent path imports neither primary,
  rebuilds both bank digests, reproduces all 43,251 values, directly checks all
  7,566 parents, eleven key children and 262 broad samples, and verifies 470
  small literal/DP controls. Both normal/-O pairs byte-match frozen outputs.
---

# THM-4104 -- the selected order-eleven strong-ear interval

**PROVED ELEMENTARY QUADRATIC-CUT IDENTITY + FINITE-EXACT + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**

## 1. The exact directed-cut quadratic

Let `T` be a labelled tournament on vertex set `V`, and write

```text
H(T)=the number of directed Hamiltonian paths of T.          (1)
```

Retain THM-4097's boundary data

```text
Start_T(b)=#{Hamiltonian paths of T starting at b},
End_T(a)  =#{Hamiltonian paths of T ending at a},            (2)

Q_T(a,b)=#{permutations of V in which a is immediately followed by b
           and every other adjacent step is an arc of T}.    (3)
```

The exposed step `a,b` in `(3)` is allowed either orientation. Adjoin a new
vertex `x`. Its cut signature is the subset

```text
S={v in V:x -> v}.                                          (4)
```

THM-4097 proves the original boundary formula

```text
H(T+x)=sum_(b in S) Start_T(b)
       +sum_(a notin S) End_T(a)
       +sum_(a notin S,b in S) Q_T(a,b).                    (5)
```

Put

```text
C_T=sum_a End_T(a),
L_T(b)=Start_T(b)-End_T(b)+sum_a Q_T(a,b).                  (6)
```

> **Lemma 1.1 (exact directed-cut quadratic).** For every subset `S` of `V`,
>
> ```text
> boxed:
> H(T+x_S)=C_T+sum_(b in S)L_T(b)-sum_(a,b in S)Q_T(a,b).  (7)
> ```

### Proof

Replace the second sum in `(5)` by

```text
sum_(a notin S)End_T(a)=C_T-sum_(a in S)End_T(a),          (8)
```

and complete the directed cut in the third sum:

```text
sum_(a notin S,b in S)Q_T(a,b)
 =sum_(b in S)sum_a Q_T(a,b)-sum_(a,b in S)Q_T(a,b).       (9)
```

Substitution of `(8)--(9)` into `(5)` is exactly `(7)`. **QED.**

In particular, if `v` is added to a current mask `S`, then `Q_T(v,v)=0` and

```text
H(S union {v})
 =H(S)+L_T(v)-sum_(u in S)(Q_T(u,v)+Q_T(v,u)).             (10)
```

Equation `(10)` evaluates the entire truth table deterministically by removing
one least significant bit at a time. It is the fast evaluator used for the
7.7-million-row scan, not a fitted recurrence or numerical interpolation.

If `T` is strong and

```text
empty != S != V,                                           (11)
```

then `T+x_S` is strong: `x` has an in-neighbor and an out-neighbor in the
strong parent, which supplies directed paths in both directions to every old
vertex. Thus all `2^|V|-2` nonconstant cuts in this theorem are strong ears.

## 2. The deterministic selected bank

The labels and iteration order are part of the construction.

1. Generate THM-4097's frozen sequence of all `6,880` order-eight
   representatives. Iterate strong representatives in that order and cut
   integers `1,...,254`. For each of the exact `1,482` order-nine values,
   retain the first labelled child. This is THM-4102's bank `R_9`.
2. Iterate `R_9` by increasing parent value and cut integers `1,...,510`.
   For each value in that selected order-ten image, retain the first labelled
   child. Call the resulting bank `R_10`.

Exactly

```text
|R_10|=7,566.                                              (12)
```

Its canonical list of pairs `(H(T),code(T))`, sorted by `H(T)`, has digest

```text
2f3fbd5d7f56de24a1f08ea08585dd029c70344ef444830915b5ea0d203e4b92. (13)
```

For every `T in R_10`, adjoin every nonconstant cut on its ten vertices. Let

```text
E_11^*={(T,S):T in R_10, empty != S != V(T)},
V_11^*={H(T+x_S):(T,S) in E_11^*}.                         (14)
```

Then

```text
|E_11^*|=7,566*(2^10-2)
        =7,566*1,022
        =7,732,452.                                        (15)
```

This definition deliberately selects one parent inside each retained scalar
fiber. Alternative equal-`H` parents can have different cut quadratics, so
neither `(12)` nor `V_11^*` is determined by the scalar order-ten image alone.

## 3. Primary exact finite theorem

> **Theorem 3.1 (selected order-eleven image).** The finite set
> `V_11^*` has exactly
>
> ```text
> 43,251 values,       minimum 225,       maximum 93,751.  (16)
> ```
>
> It contains the solid odd intervals
>
> ```text
> boxed: [429,80265]_2,       39,919 values,               (17)
>        [80875,84259]_2,      1,693 values.               (18)
> ```

Here `[a,b]_2={a,a+2,...,b}`. Within the selected image, the adjacent values

```text
427, 80267, 80873, 84261                                    (19)
```

are absent. This makes `(17)--(18)` maximal selected-image components; it is
not an absence claim for the full order-eleven spectrum.

The first selected-image components are

```text
{225}, [243,245]_2, {255}, {261}, [265,267]_2, {273},
{279}, {285}, {291}, [295,297]_2, [301,309]_2,
[313,315]_2, [319,325]_2, {329}, {333}.                    (20)
```

The second-largest component is `(18)`. Further exact components include

```text
[86995,88317]_2,       662 values,
[86823,86953]_2,        66 values,
[88321,88413]_2,        47 values.                         (21)
```

Equations `(16)--(21)` concern `V_11^*`, not all strong tournaments on eleven
vertices and not all isomorphism classes.

## 4. Explicit labelled witnesses

Encode a labelled order-`n` tournament by

```text
code(T)=sum_(0<=i<j<n)[i -> j] 2^b(i,j),                  (22)
```

where `b(i,j)` is the lexicographic upper-pair index, least significant bit
first. Direct Held--Karp recomputation and strongness tests give:

| `H` | selected parent `H` | cut signature | cut weight | order-eleven code |
|---:|---:|---:|---:|---:|
| `225` | `125` | `32` | `1` | `34095368048213567` |
| `429` | `125` | `116` | `4` | `33813343248580159` |
| `14657` | `773` | `898` | `4` | `4980304472833599` |
| `14777` | `697` | `283` | `5` | `25105481806249279` |
| `80265` | `13253` | `903` | `6` | `2181427721976895` |
| `80875` | `15059` | `182` | `5` | `33617973800943135` |
| `84259` | `12667` | `15` | `4` | `36024248539118623` |
| `93751` | `15581` | `841` | `5` | `3645901060717615` |

In particular `14,657` and `14,777`, the two lane targets left by THM-4102,
are explicit strong atoms rather than merely products.

## 5. Global prefix and the two strong lanes

THM-4102 proves every allowed odd value through `14,655`, where "allowed"
means excluding the permanent omissions `7,21` of THM-1370. The primary
interval `(17)` overlaps that prefix and directly continues through `80,265`.

Together with the inherited THM-4102 certificates, the exact selected image
supplies every ordinary prime through `80,387` other than `7`, and every value
`7p` through prime `p=11,519`, other than the excluded `p=3`. The new image
supplies the formerly open ranges beginning at `14,657` and `7*2,111`;
THM-4094 proves these are strong atoms and gives the exact
factorization by ordinary primes, the exceptional `7p` lane, and the inherited
strong carries `49,63,343`. Replaying that factorization for every allowed odd
integer through `80,405` gives

```text
40,201 values checked,
40,189 inherited or directly selected,
12 supplied by order-join multiplication.                 (23)
```

The twelve multiplicative handoffs are

```text
80267, 80281, 80285, 80289, 80303, 80343,
80359, 80365, 80375, 80383, 80397, 80403.                 (23a)
```

Therefore:

> **Corollary 5.1 (finite prefix).** Every positive odd integer at
> most `80,405`, except `7` and `21`, is a tournament Hamiltonian-path count.

The next lane values not forced by these selected finite certificates are

```text
ordinary prime:       80,407,
exceptional lane:     7*11,527=80,689.                    (24)
```

Both are genuine prime-lane boundaries: `80,407` and `11,527` are prime and
the selected image contains neither `80,407` nor `80,689`. They are not
claimed absent from the actual tournament spectrum.

## 6. Exact primary and independent verification

The companion referee pins the THM-4097 and THM-4102 compiler hashes, rebuilds
the deterministic bank from the complete order-eight representative universe,
and performs these exact checks:

```text
complete order-eight quadratic/original comparisons     1,526,032
complete selected order-nine comparisons                  755,820
order-ten original-formula controls                        53,524
direct DP/strong checks of selected order-ten parents       7,566
selected order-eleven quadratic evaluations             7,732,452
direct DP/strong order-eleven boundary controls                  8. (25)
```

The order-ten comparison universe uses six spread masks on every selected
parent and all `1,022` masks on eight inherited boundary parents. Every check
uses `require`; the source has no `assert` nodes and no floating literals.
Normal and optimized Python produce byte-identical output.

The independent implementation imports neither primary compiler. It rebuilds
both selected banks with a separately structured boundary evaluator, scans
all `7,732,452` ears, and reproduces all `43,251` values. It checks all
`7,566` selected parents by direct DP/strongness, eleven retained children and
`262` deterministic broad samples, plus `470` small fast/DP/literal controls.
Its additional frozen rows include `80,387`, `80,405`, and `80,633`.

Reproduce from the repository root with

```bash
python3 -B 04-computation/tournament_selected_order11_strong_ear_interval_thm4104.py
python3 -B -O 04-computation/tournament_selected_order11_strong_ear_interval_thm4104.py
python3 -B 04-computation/tournament_selected_order11_strong_ear_interval_thm4104_independent_audit.py
python3 -B -O 04-computation/tournament_selected_order11_strong_ear_interval_thm4104_independent_audit.py
```

Frozen raw-byte hashes:

```text
script  996c6f13abd82d0b8d1e74cf7ec949b31d093d70979715440ec8caa56a379e03
output  75437285e47f97cdf80ac5a3e816149d35ca62f157b35850ce9f5629e24370d4
independent script c607dac3d2f223fb95e76293a02f6f90767e8603ecb64961cd981d71f6e361e1
independent output 40131dd6d92e21080601d030099e8a9bf749ee994d0a5f38c96a1c17fc922482
```

The semantic ledger digest is
`11241b1a8d55d9a1b725b2343b0cf8543397a48a2eb9ba102bb3132b8b020067`.
The independent semantic digest is
`77a73efc867036bcb662dc5974a704e4fcda9fe8270c778cf032e87f2ba95a13`.

## 7. Typed loss ledger

The construction's exact connection contract is:

| field | content |
|---|---|
| source | THM-4102's selected labelled bank `R_10` |
| target | selected order-eleven value image `V_11^*` |
| map | `(T,S) -> T+x_S -> H(T+x_S)` via `(7)` |
| preserved | exact Hamiltonian-path count; strongness for nonconstant `S` |
| destroyed | unselected equal-`H` parents, isomorphism multiplicity, full order-eleven coverage |
| required sidecar | labelled parent code, selection order, `(Start,End,Q)`, and cut signature |
| decisive tests | original formula on `(25)`, both bank digests, and independent DP/literal controls |

THM-4099's squarefree mixed-insertion polynomial has the same procedural
warning: composing after scalarization loses mixed response. No theorem is
transferred between those objects; here the native mixed response is the
directed quadratic `Q_T` in `(7)`.

## 8. Scope and the next structural target

This result does not enumerate all order-eleven tournaments, all strong
order-eleven isomorphism classes, or even all ears over every strong order-ten
parent. A hole in `V_11^*` is only a failure of this selected supplier. The
minimum and maximum in `(16)` are selected-image extrema, not global extrema.

The finite orders nine, ten, and eleven suggest the following precise **OPEN**
recursion. Seed `R_9` as above. Given `R_n`, define

```text
V_(n+1)^*={H(T+x_S):T in R_n, empty != S != V(T)},
R_(n+1)=the first labelled witness of each value in V_(n+1)^*,
```

using increasing parent value and increasing integer cut signature. Prove
that the successive quadratic-cut images contain odd intervals

```text
I_n=[L_n,U_n]_2 subset V_n^*,
I_n intersect I_(n+1) nonempty,
U_n -> infinity.                                          (26)
```

The observed intervals

```text
n=9:  [85,2881]_2,
n=10: [249,14649]_2,
n=11: [429,80265]_2                                      (27)
```

are finite evidence for `(26)`, not an induction theorem, recurrence, or
asymptotic proof. Establishing `(26)` would require structural control of how
the labelled `(Start,End,Q)` fibers evolve under the first-witness selection;
the three computed rows do not provide it.

The finite constants are proved and independently audited. A separate proof
of `(26)` is still absent, so the global H-spectrum conjecture remains
**OPEN**.

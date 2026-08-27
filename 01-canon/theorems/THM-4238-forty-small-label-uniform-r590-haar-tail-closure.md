---
id: THM-4238
title: "Forty small-label uniform-r590 Haar tail closure"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170/4191 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. For every genuine outsider
  q in {2,...,49}\P, every r>=590, and every nine-body B from the
  thirty-label anchor pool P, the exact Haar safe set of B union {q,r}
  has mass at least 4/63. Together with the inherited zero- and one-outsider
  layers, this closes every eleven-face of P union {q,r} and hence its
  universal odd-tail Haar completions. It is not physical entry.
source: root/lrc-small-label-entry/2026-08-26
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
related:
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4233-pair-specific-primitive-observable-oscillation-haar-charts
  - THM-4234-fixed-fifty-twenty-label-pair-haar-charts
scripts:
  - 04-computation/lrc14_small_label_cofinal_primary_driver_thm4238.sh
  - 04-computation/lrc14_small_label_cofinal_independent_driver_thm4238.sh
  - 04-computation/lrc14_small_label_hybrid_literal_endpoint_thm4238.py
  - 04-computation/lrc14_small_label_hybrid_literal_midpoint_thm4238.py
  - 04-computation/lrc14_small_label_cofinal_audit_thm4238.py
outputs:
  - 05-knowledge/results/lrc14_small_label_cofinal_primary_ledger_thm4238.out
  - 05-knowledge/results/lrc14_small_label_cofinal_independent_ledger_thm4238.out
  - 05-knowledge/results/lrc14_small_label_hybrid_literal_thm4238.out
  - 05-knowledge/results/lrc14_small_label_artifact_manifest_thm4238.out
script_sha256:
  - 4687b1fa748e77456d8ca8790f5c52eb1e3ac65d256d9f24b22672b77eac7267
  - a12f577756ec564669ddf74502c5a8731ffd0f681254ef98344f2db2ed05ef1d
  - 8c01adbc5235b7dab7e1b2bdc06a429da6fd69286d5c9b7fbd00813394ecdee7
  - 6013bea2008da53bb3bf6b47b45a7e2064963dec1f7a1ad2c09db8b560a983fe
  - c42a7abbae49580208f47ff54eb952a749907191a6e350ac16a54ec51cc0ae86
output_sha256:
  - 7e74a9e562feace0e7863de882f1dc254affa9570304cf78f8103b070436a352
  - 43c5a1820cc25a65e4b60fef95f9bd4e964a0d18d25b9c1683d6d6a996be2f8f
  - fe197226122a90e8f3b94f1974062034cc191920c16e82e740a2086ba20b7ab4
  - 8008c00a1467089b2187d1fdbe1524589fc3c5ddd86300edbbb7337362ee1ee4
hash_basis: raw LF bytes
audit: >
  PASS. Two structurally different exact engines exhaust the same
  40*binom(30,9)=572,286,000 bodies and agree on every denominator, body
  count, extremal threshold, multiplicity, witness, mass, component count,
  and surplus. Endpoint-toggle and common-grid midpoint integrators then
  byte-agree on the only 32 rows below the analytic envelope; all are
  strictly safe. Optimized-level replays preserve the frozen ledgers.
---

# THM-4238 -- forty small-label uniform-r590 Haar tail closure

**PROVED RELATIVE TO THM-4150/4156/4170/4191 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                         (1)
```

Its genuine outsider labels below fifty are

```text
Q_small={2,3,4,5,6,7,9,11,12,13,14,17,18,19,21,22,
         23,24,25,26,27,28,29,31,32,33,34,35,36,37,
         38,39,41,43,44,45,46,47,48,49}
       ={2,...,49}\P.                                    (2)
```

For a finite positive speed set `H`, write

```text
G_H={x in R/Z:min_(h in H)||hx||>=1/14}.                 (3)
```

> **Theorem.** For every `q in Q_small`, every integer `r>=590`, and every
> `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {q,r}))>=4/63.                           (4)
> ```

There is no diagonal convention hidden in `(4)`: `q<=49`, `max(P)=290`, and
`r>=590`, so `q,r` are distinct and both lie outside `P`.

The inherited Pascal layers turn `(4)` into a complete chart. Every
eleven-subset of `P union {q,r}` contains zero, one, or two outsiders.
THM-4156 closes the zero-outsider faces, THM-4191 closes the one-outsider
faces, and `(4)` closes the two-outsider faces. Hence all

```text
binom(32,11)=129,024,480                                 (5)
```

eleven-faces are safe. For the multiplication map `m_c:x|->cx`,

```text
G_(cH)=m_c^(-1)(G_H),              mu(G_(cH))=mu(G_H).  (5a)
```

Indeed, `m_c` preserves Haar measure. Apply THM-4150 to `cH`. For every two
distinct positive odd integers `a,b`, there is an `x in R/Z` with

```text
min_(v in 2cH union {a,b})||vx||>=1/14.                  (6)
```

Equation `(6)` is a Haar completion consequence, not physical entry.

## 2. The exact retained observable

Fix `q,B` and put

```text
V_(q,B)=G_(B union {q}),       mu(V_(q,B))=m/D,
c_(q,B)=# positive-length circular components of V_(q,B). (7)
```

Endpoint-only null pieces are discarded from the component count. THM-4170's
bounded-primitive discrepancy inequality applies to any finite union of
circular intervals and gives

```text
mu(V_(q,B) intersect G_r)
 >=(6/7)(m/D)-6c_(q,B)/(49r).                           (8)
```

Define the exact surplus numerator and analytic activation threshold by

```text
delta_(q,B)=54m-4D,
kappa_(q,B)=ceil(54c_(q,B)D/(7delta_(q,B))).             (9)
```

The census below proves `delta_(q,B)>0` in every declared row. Relative to
the target `4/63`, the right side of `(8)` has exact gap

```text
delta_(q,B)/(63D)-6c_(q,B)/(49r)
 =(7delta_(q,B)r-54c_(q,B)D)/(441Dr).                  (10)
```

Thus `(8)` proves `(4)` whenever `r>=kappa_(q,B)`. The ceiling is exact:

```text
7delta*kappa>=54cD,
7delta*(kappa-1)<54cD.                                  (11)
```

The target inequality is nonstrict, so equality in the first line of `(11)`
is lawful. `kappa` is sharp only for the displayed lower bound, not claimed
to be the first literally safe tail label.

## 3. Complete forty-ray census

The primary engine constructs the exact common-wall geometry and exhausts
every

```text
40*binom(30,9)=572,286,000                               (12)
```

pair `(q,B)`. A separate C++ engine starts from independent cell geometry,
generates a literal source for each `q`, and independently exhausts the same
universe. Jointly they evaluate `1,144,572,000` bodies. They agree on every
denominator, body count, threshold maximum and multiplicity, extremal witness,
mass numerator, component count, and surplus. Essential mismatches are zero.
The full primary and independent ledger hashes are respectively

```text
7e74a9e562feace0e7863de882f1dc254affa9570304cf78f8103b070436a352,
43c5a1820cc25a65e4b60fef95f9bd4e964a0d18d25b9c1683d6d6a996be2f8f. (13)
```

Put

```text
B_0={170,176,190,193,240,252,264,286,290}.              (14)
```

The complete analytic envelope is:

- `q=6`: `B_0` is the unique `kappa=614` body; the second threshold is `587`;
- `q=25`: `B_0` is the unique `kappa=598` body; the second threshold is `590`;
- for every other `q in Q_small`, `max_B kappa(q,B)<=589`, with maximum
  `589` attained at `q=4,24`.

Consequently `(8)` already proves every required row with `r>=590`, except

```text
(q,B,r)=(6,B_0,590..613)       or       (25,B_0,590..597). (15)
```

## 4. The 32-row literal bridge

An endpoint-toggle integrator and a structurally separate common-grid
midpoint integrator construct `G_(B_0 union {q,r})` literally for all 32 rows
in `(15)`. Their complete transcripts are byte-identical. Every row is
strictly safe. The smallest scaled margin is

```text
63mu(G_(B_0 union {q,r}))-4
 =596873432123/92353981720>0       at (q,r)=(6,592).     (16)
```

In particular, the last failures of the sufficient analytic gate,
`(q,r)=(6,613)` and `(25,597)`, are literal strict-safe controls. This is the
required firewall between a failed sufficient bound and an unsafe set.
Equations `(8)--(16)` prove `(4)`.

## 5. Inheritance, novelty, and hostile boundary

The closest proved mechanism is THM-4231's retained `(mass,components)` tail
for one fixed outsider. The new step is a complete `q`-indexed census of that
two-coordinate observable followed by a uniform analytic/literal splice.
THM-4207 is the canonical hostile: marginal newcomer certificates do not
compose. MISTAKE-520 supplies the corrected inheritance rule: the zero- and
one-outsider layers are inherited, so only the two-outsider layer in `(4)` is
new.

The boundary body `B_0 union {6,590}` is primitive and contains a multiple of
every integer `2,...,14`. Thus this is not a MISTAKE-511 missing-divisor-clock
family. The new coverage is also not subsumed by the nearby rays:

- THM-4231 treats `q=1` or two labels at least `1290`;
- THM-4234 fixes `50` and restricts the bodies to its petal family;
- THM-4233's proved primitive/cancelled gate has both labels above `1290`;
- THM-4227 requires a directed coordinate at least `3391`, and its swapped
  inequality is impossible for `q<=49`;
- THM-4228 requires a common gcd at least `3467`, whereas here
  `gcd(q,r)<=q<=49`.

## 6. Verification and scope firewall

The primary and independent drivers reuse the already-audited THM-4231
geometry kernels whose frozen source hashes are, respectively,

```text
278426831e51e7052eab980c94274d30cb21255f8a10cc0a54eab18042ce6bf2,
48ccab283874ff35f362e6d3c71b36e958dbe0dc444307dadb23b0ec62194bea,
1c3508225ae9ee988d39a796886b888921fc37fd5ac0d1153b24a9bac91cbe2a. (17)
```

The first two hashes are the primary executable and its direct-body include;
the third is the independent template. The checked-in 137-entry manifest
freezes every generated per-`q` source and output from the original audited
run. The current drivers isolate each replay in a fresh temporary directory,
so unrelated `/tmp` files cannot enter a ledger. Reproduce on the audited
clang/libomp toolchain with

```bash
JOBS=3 sh 04-computation/lrc14_small_label_cofinal_primary_driver_thm4238.sh "$PWD"
JOBS=3 sh 04-computation/lrc14_small_label_cofinal_independent_driver_thm4238.sh "$PWD"
python3 -B 04-computation/lrc14_small_label_hybrid_literal_endpoint_thm4238.py > /tmp/lrc14-small-label-literal-endpoint.out
python3 -B 04-computation/lrc14_small_label_hybrid_literal_midpoint_thm4238.py > /tmp/lrc14-small-label-literal-midpoint.out
cmp /tmp/lrc14-small-label-literal-endpoint.out /tmp/lrc14-small-label-literal-midpoint.out
python3 -B 04-computation/lrc14_small_label_cofinal_audit_thm4238.py
```

Set `OPT_LEVEL=-O2` for the primary low-level replay and `OPT_LEVEL=-O0` for
the independent low-level replay; the displayed defaults are `-O3`. These
`-O2/-O3` and `-O0/-O3` pairs byte-match their respective frozen ledgers. On
the audited Apple M2 host, the primary and independent `-O3` wall-clock times
were about `156s` and `104s`; both literal paths and the final audit take under
one second.

No claim is made here for:

- `q=50`, `51<=q<=1289`, or `q in P` as a genuine outsider;
- `r<590`, or literal minimality of the uniform constant `590`;
- replacement or maximality of `P` or `Q_small`;
- arbitrary pair entry or every physical row being new;
- physical entry, arbitrary fourteen-speed sets, or LRC(14).

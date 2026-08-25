---
id: THM-4086
title: "Rule-A adjacent-offset transition clock and dyadic phase cocycle"
status: >
  PROVED STRUCTURAL + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. Adjacent valid offsets have the same temporal Pascal-kernel word, while
  their feed address, initial high front, and autonomous cell-zero debt do
  change. At R=32768 these coordinates diagnose the two terminal mechanisms
  in the local Rule-A split 854:DIE, 855:CLOSED, 856:CLOSED and the capture
  rows 20238,20233. The
  universal dyadic degree-profile cocycle has nonconstant even and odd
  defects, so no phase-free scalar affine scale conjugacy exists. This is a
  fixed-policy local transition theorem, not global offset monotonicity,
  alternative infeasibility, a uniform family, or an AMM constant.
source: codex-frontier-synthesis-creative-20260825f / AMM transition niche
audit: >
  The primary exact enclosure compiler and a structurally independent
  Fibonacci--Lucas/direct-prefix plus reverse-binomial referee agree on one
  common semantic digest. The referee separately uses dense state arrays at
  R=16 and gates the full CLOSED/event tuples. Both short companions hash-check
  and parse the archived full-state
  status checkpoint; they independently rebuild the clock/debt path, not the
  three expensive full states. Ordinary, optimized, stored-output, LF,
  AST-no-assert/no-float, dependency-slug, source-hash, output-hash, and scope
  audits pass.
depends_on:
  - THM-3026-admissibility-is-multiplicative-and-the-doubling-obstruction-is-combinatorial
  - THM-3332-s-cone-one-row-certificate-32768-closures-and-the-256-floor-point
  - THM-3597-amm-dyadic-rule-a-transition-and-adjoint-horizons
  - THM-3644-amm12592-exact-offset-threshold-at-R16384
related:
  - THM-3588-amm-r512-truncated-adjoint-pascal-repair-horizons
  - THM-3648-amm-r16384-terminal-local-failure-adjoint-golden-phase-recovery
script: 04-computation/amm_rule_a_transition_clock_phase_cocycle_thm4086.py
output: 05-knowledge/results/amm_rule_a_transition_clock_phase_cocycle_thm4086.out
independent_script: 04-computation/amm_rule_a_transition_clock_phase_cocycle_independent_thm4086.py
independent_output: 05-knowledge/results/amm_rule_a_transition_clock_phase_cocycle_independent_thm4086.out
script_sha256: 88a4c2233ff702f0f190546564e83ee661e0ec0245861f85fbdfeee940fc3614
output_sha256: e774b0cb2132bf114a10e0b0624331fefc7472bbf70369b6f0e278a1b0d3601b
independent_script_sha256: 22f08481bf578140abf31144e49aa8518e7a715dfcd7e4f012988ac269a06417
independent_output_sha256: fd4b6d458a60bf147e4dd25055c363fa4d2a6edbb468a8cd5eb867947419eca9
common_semantic_sha256: 079764dd469362281f1374d6f2a9db81652a02af0ec827420dae71fe774d20fd
primary_semantic_sha256: 2b49218c35d0c2793dbcca51e9d1f66c2d07f70ff494a46b5266a849dd728259
independent_semantic_sha256: 1667d631e30f37800eabb1ad472113e931fe0f05beb233a78dc200a90f052b61
hash_basis: raw LF bytes for files; canonical compact JSON for semantic ledgers
---

# THM-4086 -- Rule-A adjacent-offset transition clock and dyadic phase cocycle

**PROVED STRUCTURAL + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** The structural statements below hold in their displayed universal
scopes. The three large Rule-A outcomes and all numerical rows labelled
`R=32768` are instead **inherited FINITE-EXACT at that one epoch** from the
hash-pinned S192/S193/S196 full-state checkpoint. The new short companions
independently reconstruct the clock/debt mechanism, not those expensive full
states.

## 1. Offset translation preserves the temporal kernel, not the state

For every nonnegative integer `n`, put `a_n=floor(alpha n)`, where
`alpha=log_5(phi^2)`. Fix `R>=2`. For every integer `D` and `0<=i<R`, put

```text
alpha=log_5(phi^2),
a_n=floor(alpha n),
d_i(D)=a_(R+i)+D,
epsilon_i=d_i(D)-d_(i-1)(D)       (1<=i<R).            (1)
```

Algebraically, for every integer `D`,

```text
d_i(D+1)=d_i(D)+1,
epsilon_i(D+1)=epsilon_i(D).                           (2)
```

The Rule-A interpretation below is restricted to **valid offsets**

```text
0<=d_0(D)<=R-2,                                        (2a)
```

so the row-zero block and its coefficient addresses are defined. Whenever
two adjacent offsets are compared as Rule-A states, both are assumed valid.
In that scope the clamped-Pascal transport kernel inherited from
THM-3332/3597,

```text
K_i(x)=(1+x)^(1+epsilon_i),                            (3)
```

is exactly the same at adjacent valid offsets: it is `(1,1)` when `epsilon_i=0`
and `(1,2,1)` when `epsilon_i=1`. Equation (2) is a proof, not a finite
observation. It says that an adjacent-offset status transition cannot be
caused by a different *temporal kernel word*. It does not say that the
Rule-A states agree.

Indeed set `q_i=d_i+i`, and write

```text
g_0=2,
g_k=(-1)^k binom(R-1,k)-1       (1<=k<R).              (4)
```

The cell-zero feed at row `i` is the consecutive address block

```text
epsilon_i=0 : {q_i},
epsilon_i=1 : {q_i-1,q_i},                              (5)
```

after deleting addresses at least `R`. This is because
`q_i-q_(i-1)=1+epsilon_i`. Thus (5) partitions, in order and without gaps
or repetitions, the coefficient interval after the initial address and
through `R-1`. Offset translation leaves the block *sizes* fixed but shifts
the addresses by one.

Let `j_i` be Rule A's cell-zero junk after row `i`, and let `f_i` be the sum
of the coefficients (4) at the addresses (5). Since cell zero receives no
state contribution from higher cells, its exact update is the scalar map

```text
w_i=j_(i-1)+f_i,
j_i=w_i-clamp_[-2,0](w_i).                             (6)
```

After the last feed, suppose `j_i=-A` with positive even `A`. Equations
(6) give

```text
A -> max(0,A-2),                                      (7)
```

so the scalar state reaches zero exactly `A/2` rows later. This recovers
THM-3332's post-feed cell-zero clock directly and identifies the additional
feed-address coordinate needed before that regime.

A second state coordinate moves in the opposite geometric direction. At the
three finite profiles below, the inherited T4b row-zero decode gives the
initial high front and its degree headroom as

```text
F_0=R-2-d_0,
d_0-F_0=2d_0-R+2.                                     (8)
```

Thus increasing `D` by one shifts the front down by one and adds two units
of headroom, even though (3) is unchanged. The triple

```text
(feed address, high front, cell-zero debt)             (9)
```

is therefore a load-bearing diagnostic sidecar; the kernel word alone is a
quotient which has forgotten coordinates that distinguish the observed
branches. The sidecar (9) is not asserted to be minimal or a complete Markov state:
the full junk vector remains part of the Rule-A state.

## 2. The inherited exact `R=32768` transition and rebuilt clock

The archived full Rule-A runs, including the state-only independent audit of
offset `854`, give the status and event columns below. The new scripts
hash-check that checkpoint and independently rebuild the remaining geometric
and cell-zero columns.

| `D` | inherited Rule-A event | `(d_start,d_event,d_last)` | `F_0` | headroom | last feed | debt there | scalar-zero row |
|---:|:---|:---|---:|---:|---:|---:|---:|
| 854 | `DIE@8246` | `(20448,25379,40043)` | 12318 | 8130 | 7709 | 25050 | 20234 |
| 855 | `CLOSED@20238` | `(20449,32551,40044)` | 12317 | 8132 | 7708 | 25060 | 20238 |
| 856 | `CLOSED@20233` | `(20450,32549,40045)` | 12316 | 8134 | 7708 | 25050 | 20233 |

This table separates the two terminal mechanisms.

1. At `D=854`, the inherited full state has the high front reach the fatal
   cell at row `8246`, long before the isolated scalar clock would reach zero
   at `20234`.
2. At `D=855,856`, the inherited full states do not die and capture exactly
   when (7) makes cell zero vanish. Thus the independently rebuilt cell-zero
   clock is load-bearing at the terminal row; no claim that every other cell
   vanishes strictly earlier is needed.
3. The move `855 -> 856` keeps last-feed row `7708` but lowers the terminal
   debt by `10`, so capture occurs exactly five rows earlier.

The scalar zero rows have the strict local peak

```text
20234 < 20238 > 20233.                                 (10)
```

For `D=854`, `20234` is an isolated cell-zero continuation, not a Rule-A
capture: the full state has already died. For the two surviving offsets, the
last column is the actual capture row.

The common temporal word has `13172` zero steps and `19595` one steps; its
canonical digest is

```text
de5b8a18782d1880fc9be55a6c8589feb0f9e84fed70b4519a05341ab654ff2b. (11)
```

This finite digest is an audit pin. Equality of the words across offsets is
already proved by (2), independently of the census.

## 3. Exact nonmonotonicity and failure of selector--lift commutation

The standalone short engine gives the small positive/hostile control

```text
R=16, D=0,1,2,3,4,5 : all CLOSED,
capture rows          : 7,8,6,7,6,6.                  (12)
```

Thus, even restricted to consecutive closing offsets of one epoch, Rule-A
capture time is neither nondecreasing nor nonincreasing. This is an exact
finite-state obstruction to treating event time as an offset potential. It
does **not** prove that the binary status `DIE/CLOSED` is nonmonotone.

THM-3026 supplies a different operation. Its degree-one lift `L` convolves
Bernstein coefficients with `(1,1)` and uses `x+(1-x)=1`; hence it represents
the same polynomial block at degree one larger. Applied rowwise, `L`
preserves the residual-polynomial sequence and therefore its first zero row.
If the Rule-A selector commuted with this lift at `R=32768`, then

```text
RuleA(856)=L(RuleA(855))                               (13)
```

would force equal capture rows. The inherited exact values
`20233 != 20238` contradict that conclusion. Therefore

```text
RuleA(856) != L(RuleA(855)).                           (14)
```

This is a selector obstruction, not an admissibility obstruction. In fact,
the `D=855` Rule-A closure proves that the fixed-`R` entry polytope is
nonempty there, and THM-3026 lifts that one solution to every `D>=855`.
Equation (14) says only that Rule A need not select those lifted points.

## 4. The dyadic profile is a two-defect phase cocycle

Write `theta_n={alpha n}`. For every nonnegative integer `n`, elementary
floor arithmetic proves

```text
a_(2n)   =2a_n+eta_n,    eta_n=floor(2theta_n) in {0,1},
a_(2n+1) =2a_n+zeta_n,   zeta_n=floor(2theta_n+alpha) in {0,1,2}. (15)
```

The phase `theta_n`, not the scalar degree `a_n` alone, determines the
defect. Exact witnesses already in the `R=16384` window are

```text
eta_16384=0, eta_16386=1,
zeta_16384=1, zeta_16385=0, zeta_16388=2.              (16)
```

Consequently there are no constants `c_0,c_1` with

```text
a_(2n)=2a_n+c_0,   a_(2n+1)=2a_n+c_1                  (17)
```

even throughout this one window, hence none globally. A phase-free scalar
affine dyadic conjugacy of the degree profile is impossible. On the full
window `16384<=n<32768`, the primary exact census is

```text
eta : 0 x 8191, 1 x 8193,
zeta: 0 x 3294, 1 x 8190, 2 x 4900.                   (18)
```

The inherited exact status data at the two scales forms the finite
transition triangle

```text
R=16384: D=412 DIE,    D=413 CLOSED,
R=32768: D=854 DIE, D=855 CLOSED, D=856 CLOSED.        (19)
```

The outer offsets satisfy `854=2*412+30` and `856=2*413+30`. Comparing
their even-indexed degree rows, the defect from twice the smaller profile is
exactly `30+eta_n`. Offset `855` is the inserted integer collar between
those two outer images. Equations (15)--(19) are not a substitution rule for
Rule-A states or statuses: the phase, feed address, front, and junk vector
have not been transported.

## 5. Exact audit and trusted-input boundary

The primary route proves

```text
10756/17987 <= alpha < 105183/175895                  (20)
```

by exact Fibonacci--Lucas comparison and verifies that the two rational
endpoints give the same floor for every `0<=n<=65535`. It then streams the
feed binomials forward from the initial address, checks every block in (5),
replays (6)--(12), and enumerates (18).

The independent route does not import that compiler. It decides the
load-bearing `7711`-row profile prefix directly, one Fibonacci--Lucas
successor comparison at a time; assigns the address blocks before seeing any
coefficient; reconstructs `binom(R-1,k)` backwards from
`binom(R-1,R-1)=1`; and uses dense, rather than sparse, state arrays for
(12). It separately checks all witnesses (16). The two routes agree on
common semantic digest

```text
79a93822d63f70f503a048241ebaa7c4874f055a30a9ff62cf1e5fd332746cfd. (21)
```

The status/event rows in Section 2 are deliberately a trusted finite input,
not a fresh full-state replay. Both scripts parse those three tuples from and
hash-check

```text
05-knowledge/results/amm12592_r32768_offset855_live_checkpoint_kps_s192.out
SHA256 a3c1c856979fd7b0b21ca12f93d08b3340345fd175187b9c724a14bc91b97710.
```

They also hash-check the archived exact engine and runners:

```text
amm12592_transient_fast_junkflow_boxeph.py
  8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad,
amm12592_r32768_offset_boundary_probe_kps_s192.py
  65b6fab453a0adb156bc79e1f9168fd46d95037fd345e2322861ec1f76f15ba0,
amm12592_r32768_offset_boundary_lowmemory_probe_kps_s193.py
  962577096266bdc6df7b004eb38412b85a8a2b47eb9dd000f4666e9829043207. (22)
```

The `R=16384` triangle input is likewise pinned to THM-3644's runner/output
hashes

```text
ddf248e8939bbdefa7b2544bbd3df1c23e47e53e00aae4013d09149eae2ca59c,
28cedef2dc0b176b62f0633d93ea23a7c316993ed197743bd54f4466eb860c21. (23)
```

Reproduce with

```bash
python3 04-computation/amm_rule_a_transition_clock_phase_cocycle_thm4086.py
python3 -O 04-computation/amm_rule_a_transition_clock_phase_cocycle_thm4086.py
python3 04-computation/amm_rule_a_transition_clock_phase_cocycle_independent_thm4086.py
python3 -O 04-computation/amm_rule_a_transition_clock_phase_cocycle_independent_thm4086.py
```

Each ordinary/optimized stream must match its stored transcript byte for
byte. Both sources reject semantic drift, carriage returns, Python assert
nodes, and floating-point constants.

## 6. Scope and precise non-consequences

- **Fixed-epoch feasibility.** Closure at `D=855` proves existence there,
  and admissibility lifting proves existence for every larger offset at this
  same `R`. Death of Rule A at `854` does not prove the `D=854` entry
  polytope empty; another policy may close. Hence no least feasible offset
  is identified.
- **Rule-A offset behavior.** The theorem proves only the displayed local
  statuses. It proves neither a global Rule-A status threshold at `R=32768`
  nor status monotonicity. The lift noncommutation (14) blocks one naive
  monotonicity proof; it is not a counterexample to status monotonicity.
- **Scale transport.** The cocycle obstructs only a phase-free scalar affine
  profile law. It does not rule out a conjugacy on an augmented state
  carrying `theta_n`, feed address, front, and junk. It proves no Sturmian
  equidistribution, limiting frequency, substitution, or status law.
- **Adjoint horizons.** THM-3588's multipliers also use actual row widths,
  clamp slack, fatal state, and a fixed prefix. Kernel-word equality alone
  transports none of those sidecars. No `R=512` adjoint horizon is promoted
  to `R=32768`, and no constraint on a successful alternative at `D=854`
  follows.
- **Global AMM.** A single finite epoch with an additive offset supplies no
  coherent all-`R` family, deadline-uniform extractor, asymptotic offset law,
  or bound/value for `C*`. In particular, MISTAKE-361 forbids recovering the
  retracted general-class golden floor from the old Hall model, and
  MISTAKE-368 forbids inferring `C*>1` merely from endpoint nonattainment.

**QED.**

---
id: THM-3648
title: "AMM R=16384 terminal local failure adjoint golden-phase recovery"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  At the last failed offset in THM-3644's complete D0=401,...,416 bracket,
  namely D0=412, the exact positive adjoint forces departure from Rule A by
  row 1564.  The ratio 1564/16384 lies below (3-sqrt(5))/8, but its absolute
  error improves on the audited R=8192 failed trace by exactly 9/8192.  The
  same-scale move from nonterminal D0=400 to D0=412 shifts the horizon by 48,
  exactly decomposed by the margin/headroom/depth phase invoice.  No global
  offset monotonicity, full-entry infeasibility, asymptotic law, or AMM bound
  is claimed.
source: kps-s190 / THM-3644 terminal-local-failure continuation, 2026-08-21
depends_on:
  - THM-3633-amm-r16384-fixed-failed-trace-adjoint-phase-shock
  - THM-3644-amm12592-exact-offset-threshold-at-R16384
related:
  - THM-3626-amm-r8192-adjoint-horizon-and-phase-rebound
script: 04-computation/amm_binary_sturmian_R16384_D412_terminal_adjoint_thm3648.py
output: 05-knowledge/results/amm_binary_sturmian_R16384_D412_terminal_adjoint_thm3648.out
script_sha256: b37f2067d245e36c8a6c60e7049e14b41fc6d1a9fda0ae6a52ffb037388c0c43
output_sha256: 443e18ecd71a97f05cd437131a6acbb96d4db029c5232c507fee907b6de3bd27
semantic_sha256: 56e327478603425cc77ae386681b25d7a0ef29778f59183538e35f704609c9e7
hash_basis: raw LF bytes
---

# THM-3648 -- terminal local AMM failure recovers the golden adjoint phase

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  The word *terminal* below is local to THM-3644's audited bracket.
That scope is load-bearing.

## 1. The adjacent failed/closed pair

At `R=16384`, THM-3644 exhausts the consecutive offsets

```text
D0=401,...,412 : DIE,
D0=413,...,416 : CLOSED.                               (1)
```

The adjacent transition is

```text
D0=412 : DIE at row 4116,
D0=413 : CLOSED at row 10116.                          (2)
```

For the failed trace, the exact Fibonacci--Lucas degree profile has

```text
(d_0,d_4116,d_16383)=(10209,12670,20006).              (3)
```

The fatal integer is negative and has `9919` bits.  Neither `(1)` nor `(2)`
proves that no offset below `401` closes: global Rule-A closure monotonicity
in `D0` is not in the proof graph.

## 2. Exact positive-adjoint wall

Use THM-3633's checkpointed top-distance adjoint.  If an admissible
continuation agrees with this failed Rule-A trace through row `s-1`, its
positive multiplier gives the necessary inequality `0<=B_s`.  The complete
exact sweep has one contiguous negative wall:

```text
B_1564>0>B_1565,
B_s<0 exactly for s=1565,...,4115.                     (4)
```

Consequently every admissible continuation avoiding the fatal inequality
must depart from Rule A by

```text
q=1564.                                                 (5)
```

At the first negative cut the multiplier ledger has exactly `3,255,076`
active cells.  The compact exact pins are

```text
fatal digest:    da8aba4bb9ab4ab93723461c1a98380604999fda9e5bc880d1065fd1b1fd40c1,
boundary digest: cb782c23c903e5fe244baf206f745259e91122a10bd7d3feb13a6e5be0fb6a58,
active digest:   7012e2dca4351711ac35d9842e46f48be8b2b5ce1aba9f12e68bb0784442ce5b,
wall digest:     b40324f8cb11b874b043bed7a3f3c55757b0a297d4ebac0731a4f0556b0d73d8. (6)
```

## 3. Exact golden-phase recovery

Let

```text
theta=(3-sqrt(5))/8.                                   (7)
```

Then

```text
q/R-theta=(-1145+512sqrt(5))/4096<0,                   (8)
```

and the sign is exact because

```text
1145^2-5*512^2=305>0.                                  (9)
```

For the audited `R=8192,D0=191` failed trace of THM-3626, the corresponding
absolute error is

```text
(2299-1024sqrt(5))/8192.                              (10)
```

Writing `(8)` over the same denominator gives

```text
|q_16384/R_16384-theta|
  =(2290-1024sqrt(5))/8192,                            (11)
```

so the absolute error improves exactly by

```text
9/8192.                                                (12)
```

Equivalently, the dyadic horizon defect is

```text
q_(16384,412)-2q_(8192,191)=1564-2*773=18.             (13)
```

This reverses the apparent deterioration at the nonterminal same-scale trace
`D0=400`: scale alone did not determine the phase.

## 4. The offset shift is a three-coordinate invoice

Put

```text
h=5R/8-d_0,
b=2d_0-R+2,
m=j-b,
ell=j-(q+1),
beta=(sqrt(5)-1)/8.                                   (14)
```

THM-3633's exact identity is

```text
q-theta R=m+1-2h-(ell-beta R).                         (15)
```

At the two `R=16384` failed offsets:

```text
              q     j      h      m      ell
D0=400      1516  4055     43     43     2538
D0=412      1564  4116     31     80     2551.         (16)
```

Thus the entire same-scale horizon shift is accounted for by

```text
Delta q=48=Delta m-2 Delta h-Delta ell
          =37-2(-12)-13.                              (17)
```

The offset changes three independently moving state coordinates.  This is
the exact mechanism behind the failed one-coordinate interpretation of the
earlier fixed trace.

## 5. Reproduction and strict boundary

Reproduce with

```bash
python3 04-computation/amm_binary_sturmian_R16384_D412_terminal_adjoint_thm3648.py
python3 -O 04-computation/amm_binary_sturmian_R16384_D412_terminal_adjoint_thm3648.py
```

The assertion-free companion source-pins the THM-3633 and THM-3644 theorem,
script, and output triples; checks the adjacent event pair; rebuilds the
complete degree prefix, fatal value, checkpointed adjoint wall, phase invoice,
radical signs, and both comparisons; and records one canonical semantic
digest.  Normal and optimized streams must match the stored transcript.

The adjoint constrains only continuations agreeing with one failed Rule-A
prefix through a cut.  It does not prove feasibility of an alternative,
infeasibility of the full entry polytope, closure monotonicity outside or
inside the bracket, convergence to `(7)`, or any AMM bound.  **QED.**

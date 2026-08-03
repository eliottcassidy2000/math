---
id: THM-3201
title: "Seven-exit factorial Newton--Euclidean closure through r=2403"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Seven proved arithmetic exits leave exactly 50 residuals with
  2201<=d<=2400.  THM-3152's first full Euclidean row and the fixed prime
  bank through 47 close all 50; independent Fraction-hull and determinant-
  hull implementations agree on every progressive degree-set trace.  The
  largest killing prime is 37.  The same exits close d=2401,...,2405, so
  every exact three-term quadratic factorial window beginning at
  1<=r<=2403 contains a nonzero moment.  The first next seven-exit residual
  is d=2406, r=2404.  This is bounded exact support and does not prove FC(3).
audit: >
  The primary companion hash-pins the proved THM-3180 Fraction-hull engine,
  which in turn pins THM-3152's recurrence, first Euclidean remainder, and
  zero-root-aware barcode.  A separate companion hash-pins THM-3180's
  self-contained trial-division, integer-recurrence, valuation, and
  determinant-lower-hull implementation.  Both replay all 50 residuals in
  four deterministic chunks and agree on the residual digest, all chunk and
  global semantic trace digests, the killer-prime histogram, and the d=2201
  hostile-boundary trace.  Planted factors v+1 and v survive.  AST audits
  find no assert node; normal and optimized transcripts agree byte-for-byte.
source: root/frontier-synthesis/factorial-next-range/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
  - THM-3146-three-step-prime-resonance-full-euclidean-newton-separation
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3153-four-step-prime-resonance-second-euclidean-newton-separation
  - THM-3170-five-step-prime-resonance-euclidean-newton-holotopy
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
  - THM-3180-six-exit-factorial-newton-euclidean-closure-through-r2198
related:
  - THM-3182-factorial-gauss-manin-rank-one-reset-and-two-transverse-smith-bands
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
script: 04-computation/factorial_seven_exit_first_flag_2400_thm3201.py
output: 05-knowledge/results/factorial_seven_exit_first_flag_2400_thm3201.out
script_sha256: 1cf067eb3cc8da29e51ee8f560f4100fb66a4b9684a2da8b1f7ef14d380e4171
output_sha256: 63dccd4919ae60484cde977b6031912df2ca299a75e46c255d82c0e8861d22ca
independent_script: 04-computation/factorial_seven_exit_first_flag_2400_independent_audit_thm3201.py
independent_output: 05-knowledge/results/factorial_seven_exit_first_flag_2400_independent_audit_thm3201.out
independent_script_sha256: afe3ec56c7efe99ff4a09c1af856b4e144312d48ec45e729d063a5fd83b0f0b0
independent_output_sha256: 2561407b6a41e7286de5ef739cc22355eb222617e92522e76cd6bf15a857fa8e
hash_basis: LF-normalized bytes
---

# THM-3201 -- seven-exit factorial Newton--Euclidean closure through r=2403

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**

## 1. Statement and scope

Let

```text
L(t^m)=m!,                    q(t)=a+bt+ct^2,               (1)
```

with `abc!=0`.  For every integer `1<=r<=2403`, the three moments

```text
L(q^r),                 L(q^(r+1)),                 L(q^(r+2)) (2)
```

cannot all vanish.

This is a bounded exact-support `{0,1,2}` theorem.  It does not cover a
missing coefficient, translated or arbitrary support, all of `SFC(1)`,
`SFC(3)`, or the three-variable Factorial Conjecture `FC(3)`.

## 2. Inherited range, exits, and exact residual universe

THM-3180 proves `(2)` for `1<=r<=2198`, equivalently `3<=d<=2200` under

```text
d=r+2.                                                        (3)
```

Fix `2199<=r<=2398`, so `2201<=d<=2400`.  THM-3124 forces a hypothetical
bad window to the resonance `b/a=-1/d` and reduces it to a common root of

```text
P=A_(d-2)^(d),          Q=A_(d-1)^(d),
A_n^(d)(v)=L((d-t+vt^2)^n).                                 (4)
```

The uniform inherited exits in this range are

```text
d prime:                 THM-3131,
d-1 a prime power:       THM-3142,
d-2 prime:               THM-3143,
d-3 prime:               THM-3146,
d-4 prime:               THM-3153,
d-5 prime:               THM-3170,
d-6 prime:               THM-3176.                          (5)
```

Every prime in the last five lines is odd in the present range.  Call `d` a
**seven-exit residual** when none of `(5)` applies.  Exact primality and
factorization give the progressive survivor counts

```text
after d prime:                 170,
after d-1 prime power:         139,
after d-2 prime:               115,
after d-3 prime:                92,
after d-4 prime:                75,
after d-5 prime:                58,
after d-6 prime:                50.                         (6)
```

The deterministic work chunks have the complete progressive counts

```text
[2201,2250] : (43,35,29,24,20,16,13),
[2251,2300] : (42,34,27,20,16,12,10),
[2301,2350] : (44,38,34,30,26,22,20),
[2351,2400] : (41,32,25,18,13,8,7).                        (7)
```

Thus the exact new computational universe consists of `13+10+20+7=50`
rows.  Its ordered residual-list digest is

```text
7539e5847cc8d9dabe1d03f0263ad20756a412a9ead305fec113958aef198c89.
                                                                    (8)
```

## 3. First Euclidean flag and finite closure

Put `n=d-2`.  THM-3152 proves that the integral first full Euclidean row

```text
R_1=(2n-1)(Q-2(n+1)(2n+1)vP)+2d(n+1)P                     (9)
```

has degree at most `n-1` and

```text
gcd_Q(P,Q)=gcd_Q(P,R_1).                                   (10)
```

For every seven-exit residual, the exact degrees are

```text
(deg P,deg Q,deg R_1)=(d-2,d-1,d-3).                       (11)
```

At every prime in the fixed bank

```text
S={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47},              (12)
```

the companions construct THM-3152's zero-root-aware necessary
common-factor degree barcode `D_p(P,Q,R_1)`.  For all 50 rows,

```text
intersection_(p in S) D_p(P,Q,R_1) intersection Z_(>0)
 =empty.                                                    (13)
```

By the proved degree-barcode lemma, `(13)` excludes every nonconstant
rational factor common to `P,Q,R_1`, hence every nonconstant common factor
of `P,Q`.  Rational polynomials with a common complex root have a
nonconstant rational gcd, so no seven-exit residual supports a bad window.
Together with `(5)` this closes `2201<=d<=2400`; THM-3180 supplies
`3<=d<=2200`.

The next five rows close uniformly without another finite flag computation:

```text
d=2401: d-2=2399 prime                    (THM-3143),
d=2402: d-1=2401=7^4 a prime power        (THM-3142),
d=2403: d-4=2399 prime                    (THM-3153),
d=2404: d-5=2399 prime                    (THM-3170),
d=2405: d-6=2399 prime                    (THM-3176).       (14)
```

Thus `(2)` holds through `d=2405`, equivalently `r=d-2=2403`.  QED.

## 4. Exact traces and independent audit

The full semantic trace records the surviving positive degrees after every
successive prime, not merely final emptiness.  The primary and independent
implementations agree exactly:

```text
[2201,2250] c5ae943fb48def718476e25ccc13a8686bd91c63eb7e0bdf5f4aa1bebc122bb9
[2251,2300] aeeea5306caef0795662da6a1ad0ddea0b6c41d4456325dad3dce80b9eecb8ff
[2301,2350] a94e496eb375418d344b3c9795ffe140bdd016ff7b4850cc4a4690094d0e2a93
[2351,2400] a1d6bd79159cedf31a448924668414b6f1a8943f542ce69dc2e2655b79999b4f
global       7235556fedc91343a003b308de1288347ff78afad44f840fdfd8391401707e03.
                                                                    (15)
```

The first killing-prime histogram is

```text
2:4, 3:13, 5:16, 7:2, 11:2, 13:10, 17:2, 37:1.            (16)
```

Thus the largest prime actually needed is 37.  At the inherited boundary
`d=2201`, `p=2` leaves only degree `2048`, and `p=3` makes the intersection
empty.

The primary companion imports the THM-3180 Fraction-hull engine with pinned
hash

```text
50b4691c72d5d0b942eafdb922273ee749127f9725bd3728a9b422bfc494e69e.
                                                                    (17)
```

That engine in turn pins THM-3152's exact recurrence and first remainder.
The independent companion instead imports THM-3180's self-contained
trial-division, integer-recurrence, valuation, and determinant-lower-hull
engine with pinned hash

```text
5d00729ae39cdf1d085eca73d4612d4960a0788dd746f945ab2e8936ebbd466d.
                                                                    (18)
```

Both replay the complete 50-row universe.  Planted common factors `v+1`
and `v` retain degree one.  Every load-bearing check raises `RuntimeError`;
AST audits find no `assert` node.  Normal and `python -O` executions compare
byte-for-byte with the stored transcripts.

## 5. THM-3182/3183 are related diagnostics only

THM-3182 proves a rank-one Frobenius reset and generic Smith type `(1,p,p)`
in a specified `x`-weighted integral lattice.  It also proves that the
scalar companion framing has Smith type `(1,1,p)` and that the connecting
gauge acquires a nonunit output index.  Projection through the fixed tail
may cancel coordinates.

THM-3183 gives the exact Hecke lattice square and an oriented two-step wedge
continuant, but its sharp same-Smith hostile proves that even the complete
reset Smith type does not determine a projected Newton row.  Its bare pivot
does not recover the arithmetic PRS walls.  Both theorems explicitly leave
the empirical Euclidean-depth staircase open.

Consequently neither theorem supplies a row-preservation implication here.
They are not dependencies of `(13)` or the bound; no Smith band,
Gauss--Manin projection, or continuant heuristic is used in the proof.

## 6. Exact first next residual

The first seven-exit residual after the proved interval is `d=2406`, hence
`r=2404`:

```text
d   =2406=2*3*401,
d-1 =2405=5*13*37                       (not a prime power),
d-2 =2404=2^2*601,
d-3 =2403=3^3*89,
d-4 =2402=2*1201,
d-5 =2401=7^4,
d-6 =2400=2^5*3*5^2.                                      (19)
```

Thus none of the seven uniform exits applies at `r=2404`, and the finite
first-flag scan stops immediately before it.  Equation `(19)` is an exact
boundary invoice, not a common-root witness or counterexample.

## 7. Reproduction

Run from the repository root:

```text
python3 04-computation/factorial_seven_exit_first_flag_2400_thm3201.py
python3 -O 04-computation/factorial_seven_exit_first_flag_2400_thm3201.py
python3 04-computation/factorial_seven_exit_first_flag_2400_independent_audit_thm3201.py
python3 -O 04-computation/factorial_seven_exit_first_flag_2400_independent_audit_thm3201.py
```

Compare LF-normalized bytes with the two declared stored outputs.  Each
companion resolves its pinned inherited engine relative to its own file, so
dependency lookup is independent of the caller's current directory.

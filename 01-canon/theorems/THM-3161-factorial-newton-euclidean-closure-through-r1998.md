---
id: THM-3161
title: "Factorial Newton--Euclidean closure through r=1998"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For the resonant exact-support quadratic factorial pair, the THM-3152
  first Euclidean flag and fixed prime bank through 47 close every current
  four-exit residual with 1001<=d<=2000.  Together with the inherited
  prime, prime-power-predecessor, two-step-prime, and three-step-prime exits,
  every exact three-term quadratic factorial window beginning at
  1<=r<=1998 contains a nonzero moment.  The first three-exit row in the
  scanned tail 1001<=d<=2000 on which this fixed observer alone fails is
  d=1384, where degree one survives;
  THM-3146 closes that row because d-3=1381 is prime.  This is bounded and
  does not prove an all-height result or FC(3).  Separately, the theorem
  proves an all-height base-p digit formula for the predecessor Newton
  polygon whenever p divides N.
audit: >
  The primary companion uses THM-3152's exact recurrence, zero-root-aware
  denominator-capacity barcode, and first full Euclidean remainder.  A
  separate implementation reconstructs the rows, uses determinant lower
  hulls instead of Fraction hulls, and independently rebuilds every local
  degree set.  It replays all 511 current four-exit residuals through d=2000.
  A separate d=1384 audit reproduces the surviving degree one and its death
  at p=173.  Normal, optimized, stored, and independent artifacts have the
  declared LF-normalized hashes.  A third exact referee checks 1,107,086
  digital-recurrence points, 8,733 lower-hull vertices, and 694 direct
  coefficient-unit sums.
source: root/frontier-synthesis/factorial-extension/2026-08-02
depends_on:
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3146-three-step-prime-resonance-full-euclidean-newton-separation
related:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
script: 04-computation/factorial_four_exit_tail_census_thm3161.py
output: 05-knowledge/results/factorial_four_exit_tail_census_thm3161.out
script_sha256: 7fd49f9d1f4151385e011cbd4c197dca6ba91194eb6386f46351770912478046
output_sha256: 1cb3cc63b5afc41ab0942b9b1673b19475fcef1afaefc6ed3e2fc9cf0ed2724f
independent_script: 04-computation/factorial_four_exit_tail_independent_audit_thm3161.py
independent_output: 05-knowledge/results/factorial_four_exit_tail_independent_audit_thm3161.out
independent_script_sha256: 121e13dfb750110567c0129e24edd3d47458a1546d98179b43c6cad9a581d337
independent_output_sha256: 642e2cef9a7591276b3fc896b57d1f389444a4a1b9ebf5a6dd20efe378499d0d
digital_script: 04-computation/factorial_predecessor_digital_skeleton_audit_thm3161.py
digital_output: 05-knowledge/results/factorial_predecessor_digital_skeleton_audit_thm3161.out
digital_script_sha256: 0a9da30ffd1462bc2f6ca3d62ad4e7d71b82afe02099931542928105704fd318
digital_output_sha256: 86e95b24334e92a9182c7d363deb8bd878e1215c242bb0658cf503aa09443994
verified_next_lane: 05-knowledge/results/factorial_large_divisor_first_flag_block_thm3161.md
verified_next_lane_sha256: 959ac52dc369b2077851138803a5906211fa90cc3ad1de8b651d75fc43642976
verified_next_lane_script: 04-computation/factorial_large_divisor_first_flag_block_probe_thm3161.py
verified_next_lane_output: 05-knowledge/results/factorial_large_divisor_first_flag_block_probe_thm3161.out
verified_next_lane_script_sha256: d3226bbc352731900ed8fdcd4af0b41abeca4f5671bf7b127cb027d5b1ef8d03
verified_next_lane_output_sha256: 77d45df581d4966248e52b3c53d1fcba86cffb9be325d8407151991d015e5bcd
hash_basis: LF-normalized bytes
---

# THM-3161 -- factorial Newton--Euclidean closure through r=1998

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L(t^m)=m!,                    q(t)=a+bt+ct^2,               (1)
```

with `abc!=0`.  For every integer `1<=r<=1998`, the three moments

```text
L(q^r),                 L(q^(r+1)),                 L(q^(r+2)) (2)
```

cannot all vanish.

This is a bounded exact-support `{0,1,2}` theorem.  It does not cover a
missing coefficient, translated or arbitrary support, all of `SFC(1)`,
`SFC(3)`, or the three-variable Factorial Conjecture `FC(3)`.

## 2. Inherited exits and exact residual universe

THM-3124 forces a bad window to the resonance

```text
d=r+2,                  b/a=-1/d,                            (3)
```

and reduces it to a common root of

```text
P=A_(d-2)^(d),          Q=A_(d-1)^(d),
A_n^(d)(v)=L((d-t+vt^2)^n).                                 (4)
```

The uniform inherited exits are:

```text
d prime:                 THM-3131,
d-1 a prime power:       THM-3142,
d-2 prime:               THM-3143,
d-3 prime:               THM-3146.                          (5)
```

Call `d` a **current four-exit residual** when none of `(5)` applies.  In
the exact universe

```text
1001<=d<=2000,                                                (6)
```

there are exactly `617` residuals relative to the first three exits and
exactly `511` current four-exit residuals.  The four-exit counts by disjoint
chunk are

```text
[1001,1100] : 45,
[1101,1300] : 100,
[1301,1500] : 101,
[1501,1750] : 125,
[1751,2000] : 140.                                           (7)
```

The classifications use exact primality and integer factorization; no row
is omitted by a heuristic filter.

## 3. First Euclidean flag and finite closure

Put `n=d-2`.  THM-3152 proves that the integral first full Euclidean row

```text
R_1=(2n-1)(Q-2(n+1)(2n+1)vP)+2d(n+1)P                     (8)
```

has degree at most `n-1` and

```text
gcd_Q(P,Q)=gcd_Q(P,R_1).                                    (9)
```

For the current residuals in `(6)`, its degree is exactly `d-3`.  At every
prime in the fixed bank

```text
S={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47},              (10)
```

the companion constructs the exact zero-root-aware common-factor degree
barcode

```text
D_p(P,Q,R_1).                                                (11)
```

For every one of the `511` current residuals,

```text
intersection_(p in S) D_p(P,Q,R_1) intersection Z_(>0)
 =empty.                                                     (12)
```

By THM-3152's proved barcode lemma, `(12)` implies that `P,Q,R_1`, hence
`P,Q`, have no nonconstant rational common factor.  Rational polynomials
with a common complex root have a nonconstant rational gcd, so the resonant
window is impossible.  THM-3152 already includes the `45` rows through
`d=1100`; the new exact tail independently closes the remaining `466` rows.

THM-3142 closes every `d<=1000`.  Combining that finite base with `(5)` and
`(12)` proves `(2)` through `d=2000`, equivalently `r=1998`.  QED.

## 4. First fixed-observer failure in the scanned tail: d=1384

The first residual relative to the first three exits in the scanned universe
`1001<=d<=2000` for which `(12)` fails is

```text
d=1384=8*173,
d-1=1383=3*461,
d-2=1382=2*691,
d-3=1381 prime.                                             (13)
```

The progressive fixed-bank intersection is

```text
after p=3:       319 degrees,
after p=5:       {1,5,6,750,751,755,756},
after p=11:      {1,5},
after p=13:      {1}.                                       (14)
```

Every later place through `47` retains degree one through a common
denominator-one block: slope one at `p=2,3`, and slope zero at every other
bank prime.  Thus this fixed observer genuinely stops at `{1}`.

This is not a common-factor witness.  At the additional divisor place
`p=173`, the only common flag block is

```text
(slope,capacity,denominator)=(1/173,1211,173),               (15)
```

so its positive degree set is

```text
{173,346,519,692,865,1038,1211},                            (16)
```

which excludes degree one.  More importantly for the fixed-bank theorem,
`d-3=1381` is prime, so THM-3146 already excludes the row.  This explains why
the first three-exit observer failure does not interrupt the current
four-exit closure.

## 5. Exact digital predecessor skeleton

The boundary also exposes a proved all-height lemma, although it is not a
dependency of the finite theorem.  Put

```text
B_N(v)=A_N^(N+1)(v),                    p|N.                 (17)
```

Direct expansion gives

```text
[v^j]B_N=binom(N,j)(2j)! Z_j,            Z_j in Z,
Z_j=sum_(ell=0)^(N-j) binom(N-j,ell)(N+1)^(N-j-ell)
    (-1)^ell(2j+1)_ell.                                      (18)
```

Here `(x)_ell=x(x+1)...(x+ell-1)`.

For every prime `p|N`, one has the exact formula

```text
NP_p(B_N)=lower hull of
w_(N,p)(j)=v_p(binom(N,j))+v_p((2j)!),   0<=j<=N.            (19)
```

In digit sums,

```text
w_(N,p)(j)
 =[2j+s_p(j)+s_p(N-j)-s_p(N)-s_p(2j)]/(p-1).                (20)
```

Here is the proof of exactness.  Write `N=pM`, `j=pq+a`, and
`H_q=2q+w_(M,p)(q)`.  Legendre's recursion gives, for odd `p` and
`h=(p-1)/2`,

```text
w(pq)=H_q,
w(pq+a)=H_q+1+v_p(M-q)                         (1<=a<=h),
w(pq+a)=H_(q+1)                                (h<a<p).      (21)
```

Here

```text
H_(q+1)-H_q=2+v_p(M-q)+v_p(2q+1).
```

The low and high halves are constant plateaux between their displayed
endpoints, so every lower-hull vertex is congruent to `0` or `h` modulo `p`.
For `p=2`, the odd point has height `H_q+2+v_2(M-q)`, while the next even
endpoint has height `H_q+3+v_2(M-q)`; the odd point is strictly above the
endpoint chord, so every vertex is even.  If `j==0 (mod p)`, all positive
terms in `Z_j` vanish modulo `p`: for `ell<p` use
`p|binom(N-j,ell)`, and for `ell>=p` use the factor `p` in
`(2j+1)_ell`.  If `j==h (mod p)`, already `p|(2j+1)`.  Hence
`Z_j==1 (mod p)` at every skeleton vertex, proving `(19)`.

If `N=sum_e a_e p^e`, the exact scout further suggests that the principal
block at scale `p^e` has

```text
slope    v_p((2p^e)!)/p^e
         ={2(p^e-1)/((p-1)p^e), p odd;
           (2p^e-1)/p^e,          p=2},
capacity min(a_e,floor(p/2))*p^e.                           (22)
```

Excess digit mass is resolved by secondary hull edges under the proved
one-digit recurrence.  The un-clipped rule is false
already at `(N,p)=(6,3)`, where the polygon has blocks `(2/3,3)` and `(1,3)`,
not `(2/3,6)`.  The principal-capacity rule `(22)` has zero failures in `574`
pairs with `4<=N+1<=300`, `p|N`, and in declared high controls through
`N=1999`; a nonrecursive closed description of every secondary edge remains
`VERIFIED`, not proved here.

## 6. Reproduction and boundary

Run the primary and independent companions in their declared chunks under
both ordinary and optimized Python and compare with the stored transcripts.
The independent companion does not import the primary tail code: it rebuilds
the recurrence, pseudo-remainder, lower hull, zero-root block, and degree DP.
The digital skeleton has its own declared referee and transcript.  Run

```text
python3 04-computation/factorial_predecessor_digital_skeleton_audit_thm3161.py
```

The separately routed large-divisor first-flag note is `VERIFIED`, not a
proved dependency: its `P,Q` hulls have a direct coefficient route, but the
Euclidean-row valuation preservation still needs proof.

The finite computation gives no fixed-bank all-height theorem.  Beyond
`d=2000`, either another place, a later Euclidean row, or the digital skeleton
`(19)` may be needed.  A nonempty barcode intersection is only an observer
address and never constructs a common factor.

**End of proof.**

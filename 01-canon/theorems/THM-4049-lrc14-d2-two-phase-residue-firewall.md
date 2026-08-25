---
id: THM-4049
title: "LRC(14) d=2 two-phase modulo-56 residue firewall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. If a finite integer pack
  H avoids ten explicit residue classes modulo 56, then any two odd
  exceptions admit a 1/14-safe full-row time in the fixed four-time bank
  {1/28,15/28,5/112,61/112}. In particular this closes
  2{1,...,10,12} union {alpha,beta} for every distinct positive odd pair,
  with no exception-height bound. The proof is an elementary mask identity;
  two independent exhaustive residue implementations verify it. This is a
  sufficient residue subfamily only: a typed THM-3818 row can hit the
  forbidden class 11, so the complementary physical image is genuine and
  LRC(14) remains open.
source: root + lrc14_probe + lrc_prefix_audit, 2026-08-24
audit: >
  PASS. The proof classifies the pack residues and the two first-bank danger
  masks exactly, then shows that every +/-15 mod 28 first-bank spoiler is
  safe at both second-bank labels. The discovery companion checks all 56 odd
  residues and 1,596 unordered pairs by Fraction and modular masks. A
  no-import integer audit tests the full-row consequence directly at common
  denominator 112 and reproduces the endpoint and three load-bearing
  hostiles. An independent scope audit supplies a typed physical pack outside
  the firewall. Normal and optimized streams byte-match both frozen outputs.
related:
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
  - THM-4024-lrc14-complete-divisor-incidence-envelope
  - THM-4117-physical-eleven-plus-two-primitive-support-obstruction
  - THM-4119-infinite-supplier-free-eleven-plus-two-residue-family
script: 04-computation/lrc14_d2_two_phase_all_height_probe_20260824.py
output: 05-knowledge/results/lrc14_d2_two_phase_all_height_probe_20260824.out
script_sha256: c24612429463a4d0caf88ccaf5d301a186b0d2346bd284141920aeb52e0699dd
output_sha256: 10d8f230a1212e4b6a6e061fb7866b6cd12335f5c530f21989b70e9365a56adf
independent_audit_script: 04-computation/lrc14_d2_two_phase_residue_firewall_thm4049_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_d2_two_phase_residue_firewall_thm4049_independent_audit.out
independent_audit_script_sha256: 4e41dad6c9ed024fa5c898fd83f488c5e298b3381dbf89f51c054ece6bc62251
independent_audit_output_sha256: dac054e56fbcc25d6512f2b24877d67a96047c9f025aa86f7ac1fbf42f8b377f
hash_basis: raw LF bytes
---

# THM-4049 -- a two-phase residue firewall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** A fixed bank of four
times closes every pair of odd exceptions over a large residue-defined class
of even packs. The mechanism is uniform in the exception heights. It does not
show that every physical pack inherited by THM-4041 belongs to that class.

## 1. Statement

Write `||z||` for distance from `z` to the nearest integer and put

```text
R={0,11,14,22,23,28,33,34,42,45} subset Z/56Z,       (1)
T={1/28,15/28,5/112,61/112}.                          (2)
```

Let `H` be any finite set of integers such that

```text
h mod 56 notin R             for every h in H,        (3)
```

and let `alpha,beta` be any odd integers, with repetition allowed. Then some
`x in T` satisfies

```text
||2h x|| >= 1/14             for every h in H,
||alpha x|| >= 1/14,         ||beta x|| >= 1/14.      (4)
```

Consequently, if `H` is positive and `alpha,beta` are distinct positive odd
integers, the distinct-speed row

```text
2H union {alpha,beta}                                  (5)
```

is `1/14`-lonely. In particular, because

```text
H_0={1,2,...,10,12}                                   (6)
```

obeys `(3)`, every row `2H_0 union {alpha,beta}` with two distinct positive
odd exceptions is lonely. There is no exception-height restriction and no
affine-certificate premise.

## 2. The two pack-safe phase banks

Take

```text
y_1=1/14,       y_2=5/56,
x_(i,j)=(y_i+j)/2,              i=1,2,  j=0,1.       (7)
```

These four lifts are exactly the times in `(2)`. For every integer `h`,

```text
||2h x_(i,j)||=||h(y_i+j)||=||h y_i||.                (8)
```

At `y_1`, the pack inequality holds exactly when `h` is nonzero modulo `14`.
At `y_2`, it holds exactly when the circular residue of `5h` modulo `56` has
magnitude at least `4`. The strict failures of this second test are the
preimages under multiplication by `5` of

```text
0,+/-1,+/-2,+/-3 modulo 56,
```

namely `{0,11,22,23,33,34,45}`. Union with the four multiples of `14` gives
exactly `R`. Thus `(3)` and `(8)` prove all pack inequalities in `(4)` at all
four displayed times.

## 3. The odd-exception mask identity

For fixed `i`, the two lifts in `(7)` differ by `1/2`. Multiplication by an
odd speed preserves that half-turn. Hence one odd exception can be dangerous,
in the strict sense `||wx||<1/14`, at most one label in each bank.

At the first bank, direct reduction modulo `28` gives

```text
w dangerous at 1/28   iff w=+/-1  (mod 28),
w dangerous at 15/28  iff w=+/-15 (mod 28).          (9)
```

If `alpha,beta` do not cover both labels, one first-bank lift proves `(4)`.
Otherwise one exception, call it `w`, lies in the second class in `(9)`.
The second-bank numerators over denominator `112` are `5` and `61`, congruent
modulo `28`, and therefore

```text
5w = 61w = +/-19 (mod 28).                            (10)
```

If `w` were dangerous at either second-bank time, its odd numerator modulo
`112` would have circular magnitude below `8`. Its reduction modulo `28`
would then lie in

```text
{+/-1,+/-3,+/-5,+/-7},                               (11)
```

contradicting `(10)`. Thus `w` is safe at both second-bank labels. The other
odd exception cannot spoil both half-turn labels by itself, so at least one
second-bank lift is safe for both exceptions. Together with the already safe
pack, this proves `(4)`.

## 4. Endpoint and hostile audit

The pack inequality is closed (`>=1/14`) and exception danger is open
(`<1/14`). The named pack genuinely uses the boundary:

```text
||2*12*(5/112)||=1/14.                                (12)
```

Three controls show that the hypotheses are properties of this fixed
certificate, not decoration:

- `H={11}`, exceptions `(1,15)` defeat all four times if the `y_2` pack
  condition is dropped;
- `H={14}`, exceptions `(1,11)` defeat all four if the `y_1` condition is
  dropped;
- the valid pack `H_0` with even exceptions `(22,28)` defeats all four, so
  oddness is load-bearing.

These are hostile rows for the four-time certificate only. None is a Lonely
Runner counterexample.

The discovery companion checks both its exact Fraction masks and an integer
residue implementation on all `56` odd classes modulo `112` and all `1,596`
unordered pairs with repetition. The independent audit imports none of that
code. It evaluates the complete row directly with common denominator `112`,
recovers the `46` allowed pack classes, and finds the safe-time-count
histogram

```text
1 time: 60 pairs;  2 times: 370;  3 times: 760;  4 times: 406.       (13)
```

Both scripts reproduce their frozen outputs byte-for-byte in normal and
optimized Python modes.

## 5. Exact scope and remaining physical gate

THM-4041 shows that the physical `d=2,c_2=9` boundary has two odd exceptions
and a divided eleven-speed pack `H`; this supplies the parity sidecar used
above. It does **not** imply `(3)`. THM-4049 therefore closes every physical
row whose projected divided pack happens to avoid `R`, but it does not close
the entire THM-4041 image, any other `11+2` branch, or LRC(14).

The complementary image is not merely formal. In THM-3818's finite box, the
typed two-component row

```text
u=(1,4,6,8,10,12,14,15,16,18,22),   v=(1,3),
s=1,                                 t=2^45           (14)
```

has component sizes `11+2`, rank eleven, and no bounded crossing row. Its two
odd exceptions are `(1,15)`, while its divided pack has residues

```text
(2,3,4,5,6,7,8,9,11,32,40) modulo 56.               (15)
```

Thus class `11` really occurs. The four times in `T` have full-row
clearances `(1/28,1/28,1/56,1/56)`, but `x=9/19` has clearance `2/19`.
This is a hostile to uniform physical entry into the firewall, not an LRC
counterexample.

The precise remaining connection is:

```text
source:    a physical THM-3818 d=2 producer
target:    (H modulo 56, alpha modulo 2, beta modulo 2)
map:       project to (H,alpha,beta), then reduce the displayed coordinates
preserves: the pack firewall test and odd-exception sidecar
needed:    prove H avoids R, or classify and close the rows that do not
lost data: owner, height, all phases outside T, and the physical tuple origin
test:      execute the projection census in THM-3818's 91^12 box, or sharpen
           that bound enough to make the census tractable.                (16)
```

The existing finite control bank does not replace that unexecuted physical
projection. This residue firewall is a proved sufficient subfamily, not an
extrapolation from exception height `79` and not a proof of LRC(14).

**QED.**

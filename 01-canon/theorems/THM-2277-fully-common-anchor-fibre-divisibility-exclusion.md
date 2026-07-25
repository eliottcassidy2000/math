---
id: THM-2277
title: "Fully common-anchor fibre divisibility exclusion"
status: >
  PROVED + VERIFIED-EXACT. In any of THM-2266's 120 interior first-depth-one
  profiles, the fully common height-20 dependency packet is impossible:
  the guard and all five unit coefficients cannot simultaneously lie on
  THM-2275's bounded multiplier rays over the shallow owner. On an odd
  common-owner fibre, a nonconstant danger comb occupies at most one third
  of the roots, so the two remaining blockers force a common-owner divisor.
  After division, THM-2137 requires boundary budget at least 231, whereas
  the complete bounded packet has budget at most 122. This excludes a
  coefficient locus inside every interior profile; it removes no profile
  and does not prove LRC(14).
source: codex-2026-07-25-common-anchor-fibre
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2137-deep-scalar-tail-boundary-complexity
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2275-mixed-scalar-relation-and-guard-blocker-crossing
related:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2246-depth-one-private-joint-two-step-fibre-cap
  - THM-2274-mixed-scalar-relative-rank-harvest-and-adaptive-pair-crossing
script: 04-computation/lrc14_fully_common_anchor_fibre_exclusion_thm2277.py
output: 05-knowledge/results/lrc14_fully_common_anchor_fibre_exclusion_thm2277.out
script_sha256: 3147a2e6ab5d39f1b6fa88118d6c882c7a47df72863f3f6c1f24ec36641574c7
output_sha256: 79ada529729f323392533d2d3eb08349679c12baa17b9c644fc3192bb01920a0
hash_basis: working-tree bytes (LF)
---

# THM-2277 -- fully common-anchor fibre divisibility exclusion

Retain THM-2266's scalar five-unit/three-blocker cover

```text
D_r={x in R/Z:||rx||<1/14},
C_H={x in R/Z:||Hx||>1/7},

C_H subset union_(i=1)^5 D_(q_i)
             union D_(c_1) union D_(c_2) union D_(c_3)     a.e.       (1)
```

in one of its `120` interior first-depth-one profiles:

```text
c_1=13u,             c_2=13^b u_2,       c_3=13^c u_3,

3<=b<=c-2,           5<=c<=19.                              (2)
```

Here `H,q_1,...,q_5,u,u_2,u_3` are positive thirteen-units, `H` is
odd, and the five `q_i` are pairwise distinct.

THM-2275 isolates the height-20 dependency rays

```text
A={1,3,5,7,9,11,15,17,19},

Q={1,2,3,4,5,6,7,8,9,10,11,12,
   14,15,16,17,18,19,20}.                              (3)
```

This theorem excludes the *fully common* packet

```text
H=au,                       a in A,
q_i=b_i u,                  b_i in Q,                  (4)
```

where the five `b_i` are distinct. Thus a counterexample in an interior
profile cannot place all six guard/unit anchors on the displayed bounded
rays over the same shallow owner `u`.

## 1. Exact occupancy of one common-owner fibre

Let `u,c` be positive integers. On a fibre of

```text
M_u:x |-> ux mod 1,
```

put

```text
g=gcd(u,c),                    n=u/g.                  (5)
```

As the `u` roots are traversed, their `c`-images form a translate of the
uniform `n`-grid, with every grid point repeated exactly `g` times. An
open circular interval of length `1/7` contains at most

```text
ceil(n/7)                                                    (6)
```

points of an `n`-grid. Consequently

```text
#(M_u^(-1)(y) intersection D_c)
 <=g ceil(n/7)
 =u ceil(n/7)/n.                                      (7)
```

The ratio in (7) has the complete all-`n` classification

```text
n=1:       ceil(n/7)/n=1;

n>=2:      ceil(n/7)/n<=1/2,
           with equality if and only if n=2;

n>=3:      ceil(n/7)/n<=1/3,
           with equality if and only if n=3.           (8)
```

Indeed, the values for `1<=n<=6` are

```text
1, 1/2, 1/3, 1/4, 1/5, 1/6.
```

For `n=7q+r`, `q>=1`, `0<=r<=6`, the two strict slacks in (8) are

```text
n-2ceil(n/7)
 =5q                         if r=0,
 =5q+r-2                     if r>0,

n-3ceil(n/7)
 =4q                         if r=0,
 =4q+r-3                     if r>0,                  (9)
```

and all four expressions are positive. This proves (8) for every positive
integer `n`, not merely below a finite cutoff.

## 2. Two blockers force common-owner divisibility

Under (4), oddness of `H=au` and of every `a` in `A` implies

```text
u is odd.                                              (10)
```

Define the divided private set

```text
P_0
 =C_a minus (D_13 union union_(i=1)^5 D_(b_i)).       (11)
```

The original private remainder of THM-2234 is exactly

```text
P
 =C_H minus (D_(13u) union union_(i=1)^5D_(q_i))
 =M_u^(-1)(P_0).                                      (12)
```

In particular Haar invariance and THM-2234 give

```text
measure(P_0)=measure(P)>=2593/90090.                  (13)
```

The cover (1) gives

```text
P subset D_(c_2) union D_(c_3)                        (14)
```

outside a null set. Root-count disintegration shows that for almost every
`y in P_0`, every root in `M_u^(-1)(y)` obeys (14): the integral of the
number of exceptional roots is `u` times the measure of the exceptional
set and is therefore zero.

For `j=2,3`, set

```text
n_j=u/gcd(u,c_j).                                     (15)
```

Every `n_j` divides the odd integer `u`. If neither blocker is divisible
by `u`, then both `n_j` are odd and at least three. Equations (7)--(8)
show that each blocker occupies at most `u/3` roots, so their union
occupies at most `2u/3<u` roots. This contradicts the full-fibre cover.
Therefore

```text
u divides c_2 or u divides c_3.                       (16)
```

If exactly one blocker is divisible, relabel it as

```text
c_2=uk_2.
```

Its indicator is constant on every common-owner fibre:

```text
1_(D_(c_2))(x)=1_(D_(k_2))(M_u x).                   (17)
```

For `y in P_0` outside `D_(k_2)`, the other, nondivisible blocker would
have to cover all `u` roots, contradicting (7). Hence

```text
P_0 subset D_(k_2)                                    a.e.       (18)
```

If both blockers divide, write `c_j=uk_j`. Both indicators are constant
on the fibre, and (14) gives

```text
P_0 subset D_(k_2) union D_(k_3)                     a.e.       (19)
```

Thus the original two-blocker cover becomes either a one-comb or a
two-comb cover after division by the common owner. This is the point where
the full common-anchor hypothesis is used.

## 3. The height-20 boundary budget is too small

The boundary budget of (11) is

```text
B=a+13+sum_(i=1)^5 b_i.                              (20)
```

The five multipliers are distinct. From (3),

```text
B
 <=19+13+(20+19+18+17+16)
 =122.                                               (21)
```

Because `u` is a thirteen-unit, division does not change the thirteen-adic
depth of a divisible blocker. The interior condition `b>=3` in (2) therefore
gives

```text
13^3 divides every k_j that occurs in (18) or (19).  (22)
```

Put `M=13^3=2197`. In the one-comb case, the divided target has measure
`1/7`. Apply THM-2137 to (11), (13), and (18):

```text
B
 >=M (2593/90090)/(1/7)
 =438217/990
 =122+317437/990
 >122.                                               (23)
```

Equivalently the integer budget would have to be at least `443`.

In the two-comb case, THM-1166 gives intersection at least `1/91`, hence

```text
measure(D_(k_2/M) union D_(k_3/M))
 <=2/7-1/91
 =25/91.                                             (24)
```

THM-2137 applied to (19) now gives

```text
B
 >=M (2593/90090)/(25/91)
 =5696821/24750
 =122+2677321/24750
 >122.                                               (25)
```

The integer budget would have to be at least `231`. Both cases contradict
(21), proving that (4) is impossible in a scalar counterexample. QED.

## 4. Exact packet size and scope

There are nine choices for `a` and nineteen choices for each distinct
unit multiplier. Consequently the excluded multiplier atlas contains

```text
9 binom(19,5)=104,652
```

unlabelled five-element packets, or

```text
9(19)_5=12,558,240                                  (26)
```

packets when `q_1,...,q_5` retain their labels. These are multiplier
shapes; the common thirteen-unit scale `u` is unrestricted.

The quantifiers must not be weakened:

1. the theorem applies only to the `120` interior profiles
   `3<=b<=c-2` from THM-2266;
2. all six coefficients `H,q_1,...,q_5` must have the same factor `u` and
   multipliers in (3);
3. a single dependency ray, or any proper subset of the six common rays,
   does not make (12) a full pullback and is not excluded here;
4. the bounded-rank branch of THM-2275 is not classified;
5. every interior valuation profile still contains coefficient choices
   outside (4), so no scalar profile is removed;
6. this theorem does not prove LRC(14).

The source is the fully labelled private remainder. The map is complete
root-fibre disintegration under `M_u`; it preserves simultaneous coverage
and the exact owner divisibility test. The quotient forgets the locations
of individual roots, but the occupancy count (7) is the needed sidecar.
The boundary-complexity transfer of THM-2137 then converts divided deep
coverage into an impossible height invoice.

## 5. Exact verification

The companion checks the all-integer ratio theorem (8) through its seven
symbolic residue classes, rather than through a numerical cutoff. It also
checks the two ray sets, both packet counts, the maximum boundary budget,
and every fraction and strict margin in (23)--(25). Reproduce with

```bash
python3 04-computation/lrc14_fully_common_anchor_fibre_exclusion_thm2277.py
python3 -O 04-computation/lrc14_fully_common_anchor_fibre_exclusion_thm2277.py
```

Normal and optimized transcripts are identical. After platform newlines
are normalized to `LF`, both reproduce the stored output byte for byte.

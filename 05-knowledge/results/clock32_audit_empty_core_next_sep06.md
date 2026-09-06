# Independent audit: three-mask obstruction at clock 32

**Status: PROVED + FINITE-EXACT + INDEPENDENT AUDIT PASS.** This audits
the root agent's clock-32 strengthening. It gives an unbounded-in-height
consequence relative to the already audited hereditary gcd sieve and its
cited lower-runner LRC input. It neither constructs nor classifies actual
counterexamples to LRC(14).

## 1. Exact consequence and inherited universe

Let V be thirteen distinct positive integer speeds with gcd(V)=1 and
M(V)<1/14. Then every eight-speed subset P has gcd(P)<=30.

The [preceding complete sieve](lrc14_recursive_gcd_empty_core_next_sep06.md)
already confines that gcd to

```text
1,2,3,4,5,6,8,9,10,12,15,16,18,24,30,32.
```

Its unique clock-32 five-tail signature is
g=(1,1,2,4,4), where g_i=gcd(32,w_i). Thus its effective orders and
maximum padded labelled sizes are

```text
q=(32,32,16,8,8),       sizes=(5,5,6,8,8).
```

The independent script reads the frozen complete signature JSON at commit
`b77d7fb5b7cf3d5e2ea7831ae8f720e632576f0f`, checks its SHA-256, and confirms
this unique boundary. No primary producer is imported and no new
speed-height census is used.

## 2. A forced three-mask overlap

A padded mask of effective order q is the inverse image in Z/32 of a
cyclic affine unit block of length ceil(q/7) in Z/q. This padding contains
the actual open danger set at every body-safe phase. Units and block
origins may differ between tails.

Take the two order-8 masks A,B and one order-32 mask Z. Each of A,B is a
union of two complete residue classes modulo 8 and has size eight. Z
consists of five terms of an odd-step progression modulo 32, hence its
five terms occupy distinct classes modulo 8.

If A and B meet, their union has size at most fifteen, so
|A union B union Z|<=20. If they are disjoint, they occupy four distinct
classes modulo 8. The five classes of Z must meet these four, since
5+4>8. Again |A union B union Z|<=8+8+5-1=20.

The two remaining masks have sizes five and six. Consequently the full
five-mask union has size at most 20+5+6=31<32, uniformly in all units,
origins, and body-safe phases. Cited lower-runner LRC supplies such a
phase for the divided eight-speed body. Its free label is a time when
all thirteen speeds have clearance at least 1/14, contradicting the
strict-failure assumption. This removes clock 32 and proves the cap 30.

This proof also identifies the information missing from the previous
pair bound: every pair type in this signature admits disjoint masks.
Pairwise minimum overlap zero does not imply a simultaneous disjoint
configuration. A three-mask union bound detects the obstruction.

## 3. The exact clock-30 boundary for this test

The inherited clock-30 signature is g=(2,2,2,3,3), with orders
(15,15,15,10,10). Choose the following maximum blocks:

```text
mod 15: {2,3,4}, {7,8,9}, {12,13,14};
mod 10: {0,1}, {5,6}.
```

Their inverse images in Z/30 partition all thirty labels, each mask
having size six. Thus the unrestricted padded-mask test at the original
clock cannot exclude this signature. This is a statement about a
specified mask relaxation. It is not an unsafe thirteen-speed row, a
claim about prescribed tails at their safe phases, or a proof that every
stronger divisor or body-sensitive test must fail.

## 4. Hereditary transfer, without a new six-tail classification

For a seven-speed body with six complementary tails, adjoining any one
tail gives an eight-speed subset. Its gcd g_i therefore cannot be 32.
In the frozen 1,217-profile six-tail certificate this immediately removes
exactly

```text
(c,g)=(32,(1,1,2,4,4,32)), (64,(1,1,2,4,4,32)).
```

This does not exclude clock 32 or clock 64 when another signature is
present. The old necessary seven-speed gcd bound 96 remains valid.
The two clock-96 signatures have no child gcd 32, so this particular
filter alone does not remove them. Further joint-mask results are a
separate input; this note does not pre-empt their classification.

## 5. Independent finite certificate

The [source](../../04-computation/clock32_audit_empty_core_next_sep06.py)
enumerates every affine unit block at orders 8,10,15,16,32, producing
respectively 16,20,60,64,256 distinct masks. It checks every one of the
34,816 unordered order-8/order-8/order-32 triples, including repeated
order-8 masks. Their exact maximum union is twenty. It also checks all
pair-type minima are zero, the literal clock-30 partition, and the
inherited boundary and child exclusions. These are finite masks, not
actual unsafe rows.

The [output](clock32_audit_empty_core_next_sep06.out) has 54,792 explicit
checks. Normal and optimized Python runs are byte-identical:

```text
python3 -B 04-computation/clock32_audit_empty_core_next_sep06.py
python3 -B -O 04-computation/clock32_audit_empty_core_next_sep06.py
```

SHA-256:

```text
source   9807186898f582d9062e87040f6546d469e466609dc30bb94aae0cc57e794114
output   816b7ea4190fb43204289bfbaac48af35587fae6c39bfeb0075b74b035b9365b
semantic e5bd026c32fd2675df1318577e2aca7b1aec02ddedbfbc753696e0c7101e4b59
```

No earlier certificate or shared navigation file was edited by this audit.

## 6. Final master proof and compiler acceptance

The subsequent [joint-shadow master](lrc14_joint_shadow_empty_core_next_sep06.md)
adds the two clock-96 exclusions and concludes the full necessary gcd
ceilings `(1,2,4,9,30,90)` for subset sizes twelve down to seven.
Its complete proof and source received an independent final audit here:
**PASS**. In particular:

- The clock-96 modulo-eight argument forces either an E-D intersection
  of at least six, an E-A or E-B intersection of at least four, or both
  D-A and D-B intersections of at least two. All six displayed five-edge
  graphs are trees. The CRT credits give union bounds 93 and 95 for the
  two respective profiles, below 96.
- The master compiler's full profile arrays equal the frozen independent
  certificate after removing only c=32 at five tails and, at six tails,
  every row with a child gcd 32 or original clock 96. This independently
  verifies all rows, not merely the counts `1,2,5,19,109,1213` or maxima.
- The common-phase realization lemma has the stated quantifiers. At
  order q the open offset interval for a specified maximum block has
  length `1/7-(ceil(q/7)-1)/q`, at least `1/(7q)`. This includes q divisible
  by seven. With `0<p<L` and L prime, the offsets `Np/L` traverse the
  full L-grid. Adding multiples of L separates the chosen positive tails
  without changing any mask or gcd. The tails are selected; no claim
  about choosing phases for already fixed tails is inferred.
- Both explicit physical boundary rows have a strictly body-safe phase
  whose lifts are all blocked, and are safe at another named phase.
  Their independently bounded subset gcd maxima satisfy every new cap.

The final master source replays with 46 gates, with normal and optimized
outputs byte-identical to its saved output. Its regenerated profile JSON
is also byte-identical. Final SHA-256 pins:

```text
source  3f906146953677c5e1734020e97ef82fb801ee66cf5ae7a6c697ce83b8d21245
output  950acfb7073ec93af9372c6ff41a5a010281a5b9f1cf38df5de62b808478fa21
JSON    935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f
```

The cap-30 and cap-90 statements are necessary conditions for a strict
counterexample, while the actual boundary rows are safe examples showing
failure of a phase-uniform lifting predicate. These types remain distinct.
All three owned audit files are now frozen for root integration.

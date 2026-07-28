---
id: THM-2721
title: "Semantic inner triangle equal following amplitude and current reanchoring no-go"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  On THM-2712's triangle 0->13->26->0, the three frozen following
  pieces are exact translates with one common positive integer weight, so
  their raw current/following arm integrals are equal and positive.  Their
  two nontrivial abstract C3 transforms vanish, and every pure endpoint
  character telescopes because the address gain is a coboundary.  More
  decisively, the entire following base support is disjoint from current
  rail 2: all 304 semantic endpoints have zero fixed-current midpoint and
  whole-cylinder reanchoring.  This is a positive parallel-corolla theorem
  and a chronological-current no-go, not an endpoint current, row exclusion,
  or LRC(14) conclusion.
source: lrc-semantic-inner-cycle-scout-2026-07-28
depends_on:
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
  - THM-2712-semantic-following-congruence-lock-and-address-coboundary-descent
related:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2706-relative-segal-macro-cycle-and-minimal-ghost-midpoint-completion
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
script: 04-computation/lrc14_semantic_inner_triangle_amplitude_thm2721.py
output: 05-knowledge/results/lrc14_semantic_inner_triangle_amplitude_thm2721.out
script_sha256: 83ba5a846805d8bd67eea5c9ed8db43236042a018ac2079ad0dd4efe20111801
output_sha256: b4fdc6ff0f34e9eaa6a500dd80b67cad9c53a354f4713f3ecba31ed95e3e54bb
hash_basis: LF-normalized bytes
---

# THM-2721 -- semantic inner triangle amplitude and reanchoring boundary

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2712 exposes a root-neutral physical support SCC one odometer level below
the nonsplit THM-2657 quotient.  Its smallest displayed triangle has genuine
positive following support, but support does not say whether the old
current can be restarted at a target.  This theorem computes both sides of
that question.  The raw three arms are exactly equal and positive; the
chronological reanchoring fibre is exactly empty.

## 1. Frozen current, following atom, and triangle

Retain THM-2712's exact data

```text
p=13,                         R=p^6=4826809,
T=297836897838480,
x=649039434905733/1304692766858936,
z={13x}=46873542509301/100360982066072,

I=(960117507257/1930018885886,
   324519717452867/652346383429468),
|I|=1/652346383429468.                                  (1)
```

Let `C` be the frozen THM-2680 current atom

```text
(rail,future,j,h,epsilon,kappa)=(2,3,5,2,0,0),           (2)
```

and let `F` be its frozen following atom

```text
(rail,future,j,h,epsilon,kappa)=(0,1,2,6,1,1).           (3)
```

For an address `n`, write

```text
M_n(u)={13u}+7n/R.                                      (4)
```

THM-2712 proves that precisely `304` addresses, all with `n=0 mod13`,
retain `(3)` on the whole cylinder `(1)`.  In particular,

```text
n=(0,13,26),                    k=(91,91,-182)            (5)
```

is a root-neutral support triangle and `sum k=0` in the integers.

## 2. The raw three-arm amplitude is equal and positive

The selected weighted base piece of `(2)` on `I` is

```text
(148163522522400,148163551171020; 27580222516).          (6)
```

The three selected weighted base pieces of `(3)` at `M_n(I)`, in the order
`n=0,13,26`, are

```text
(139104562225860,139104590874480; 27582102210),
(139110177355380,139110206004000; 27582102210),
(139115792484900,139115821133520; 27582102210).          (7)
```

Every interval in `(7)` has length `28648620`, and successive rows are
literal translates by

```text
91T/R=5615129520.                                       (8)
```

The delayed Boolean factors of both atoms equal one on the appropriate
pullbacks of all of `I`; the exact whole-cylinder slacks are independently
rechecked by the companion.  Thus, if `rho_C,rho_F` denote the unintegrated
integer weighted densities of the two frozen atoms, each of the three
parallel arms has the same exact raw overlap

```text
A_n = integral_I rho_C(u) rho_F(M_n(u)) du

    =760720516410855360360/652346383429468
    =2089891528601250990/1792160394037 >0.              (9)
```

This is stronger than equal support.  It is an exact positive equality of
the inherited raw weighted integrands.  It is deliberately not called a
canonical endpoint-current coefficient: the weights in `(6)--(7)` are the
positive atom numerators before the signed/complex endpoint allocation.

There is also a sharp modular reading.  Index the three arms abstractly by
`C_3` and let `omega^2+omega+1=0`.  Since their vector is `(A,A,A)`,

```text
A+omega A+omega^2 A=0,
A+omega^2 A+omega A=0.                                  (10)
```

So the positive scalar supplies only the trivial abstract `C_3` character;
both charged characters vanish.  The word **abstract** is load-bearing:
`(5)` uses gains `(+91,+91,-182)`, not three powers of one physical
order-three translation.

## 3. Pure endpoint phase has trivial holonomy

On every address edge THM-2712 gives the exact potential

```text
phi(n)=7n,                       k_(a,b)=phi(b)-phi(a).    (11)
```

Hence for every integer physical frequency `N`, the pure translation phase
around `(5)` is

```text
product_edges exp(2 pi i N k_edge/R)
 =exp(2 pi i N(91+91-182)/R)=1.                          (12)
```

Equations `(10)` and `(12)` expose two different losses.  The raw positive
arm mass is flat in the comparison `C_3` index, while the available pure
endpoint phase is an exact coboundary.  Neither produces a charged
transition curvature.

## 4. Fixed-current reanchoring is globally empty

Let `Supp_base(F)` be the complete merged base support of `(3)` and let
`Rail_2` be the complete weighted rail support used in `(2)`.  Exact sorted
interval intersection gives

```text
Supp_base(F) intersection Rail_2 = empty.                (13)
```

Since the full current base support is contained in `Rail_2`, also

```text
Supp_base(F) intersection Supp_base(C) = empty.          (14)
```

This is not a midpoint seam.  Every semantic target cylinder `M_n(I)` lies
inside `Supp_base(F)`, so `(13)` implies

```text
#{semantic n : M_n(x) is a current midpoint}=0,
#{semantic n : M_n(I) is a whole current cylinder}=0     (15)
```

over all `304` addresses.  In particular none of the three targets in `(5)`
can be the source of another copy of the fixed current `(2)`.

Categorically, `(9)` is therefore a positive three-arm **corolla**: the maps
`M_0,M_13,M_26` share the source cylinder `I`.  The address graph compares
their targets, but its comparison edges are not composable chronological
current arrows.  A triangle in that comparison graph is not a two-simplex
in the current/following nerve.  This is the same controlled-forgetting
boundary abstracted by THM-2706, now on the exact semantic inner cycle.

## 5. Sharp changed-source control

The failure is specific to target reanchoring, not a blanket absence of
current support near every factorization.  Write a semantic address as
`n=13j`.  Then

```text
M_(13j)(u)=D(u+7j/R),                 D(v)={13v}.         (16)
```

If the source cylinder itself is translated by `7j/R`, the fixed current
atom retains the whole translated cylinder for exactly four semantic arms:

```text
n=0,169,338,507.                                          (17)
```

Equation `(17)` is a positive control and a warning.  It changes the source
gauge before applying `D`; it does not put the target `M_n(I)` back in the
current support, and it does not repair the triangle `(5)`.  A successful
continuation must therefore supply a lawful changed current/rail/future
branch at the target, or an independently typed endpoint-current
intertwiner.  Another address potential, pure endpoint character, or
support-SCC computation cannot supply that missing map.

## 6. Scope

The theorem proves all of the following in the frozen THM-2680/2712 box:

```text
positive equal raw three-arm amplitude;
flat nontrivial abstract-C3 scalar spectrum;
trivial pure endpoint-character holonomy;
zero fixed-current reanchoring on all 304 semantic endpoints;
four sharp changed-source factorizations.                (18)
```

It does not construct or cancel a THM-2334/2625 endpoint coefficient, choose
a target action, transport an owner/source label, build a changed current
branch, exclude a scalar row, or prove LRC(14).  The ledger remains `165`.

## 7. Exact reproduction

Run

```bash
python 04-computation/lrc14_semantic_inner_triangle_amplitude_thm2721.py
python -O 04-computation/lrc14_semantic_inner_triangle_amplitude_thm2721.py
```

Both executions byte-match the stored `15`-line transcript
`05-knowledge/results/lrc14_semantic_inner_triangle_amplitude_thm2721.out`.
The companion contains no optimized-away assertions.  It reconstructs the
`304` semantic addresses without reading THM-2712's stored output; checks
every exact interval, weight, translate, overlap, and `C_3` identity; proves
the global support disjointness `(13)--(14)`; tests both reanchoring counts in
`(15)`; and exhausts the changed-source controls `(16)--(17)`.

Awaiting independent hostile audit before promotion.

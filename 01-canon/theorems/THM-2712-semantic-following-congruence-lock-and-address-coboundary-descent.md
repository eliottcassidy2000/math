---
id: THM-2712
title: "Semantic following congruence lock and address-coboundary descent"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Of the 3,346 THM-2707 whole-cylinder packet addresses, exactly
  304 retain the full frozen THM-2680 following atom.  They are precisely
  n=0 mod13, one in each of 304 following components, and share the whole
  inherited cylinder.  Their required nonzero-root THM-2657 graph is
  edgeless, although root-neutral 13-divisible physical translations form a
  complete graph and its exact-valuation-one descent is a complete directed
  13-partite SCC with 85,274 edges.  The full address gain is an exact
  coboundary upstairs; its failure to descend through n mod13 is the kernel
  class 7.  This is a fixed-current/rail/clock support theorem, not an
  endpoint-current amplitude, row exclusion, or LRC(14) conclusion.
source: codex/thm2707-semantic-cycle-scout-2026-07-28
depends_on:
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
  - THM-2707-full-physical-lift-fibre-common-simplex-and-packet-scc
related:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2710-central-half-phase-literal-word-nilpotence-and-prescribed-clock-invisibility
script: 04-computation/lrc14_semantic_following_congruence_lock_thm2712.py
output: 05-knowledge/results/lrc14_semantic_following_congruence_lock_thm2712.out
script_sha256: 76df65e9a79d6f0bde103b40fa3df3b7ff1467cb3de954e66969620b8e4dd215
output_sha256: 017d5d6d6f5feb08294e53582b30edff8744e034933f6fcbb1a04a2d42daff22
hash_basis: LF-normalized bytes
---

# THM-2712 -- semantic following congruence lock and address-coboundary descent

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2707 proves an enormous fixed-skeleton support SCC, but its frozen
THM-2680 following atom is absent on the displayed outer cycles.  The exact
repair is not merely negative.  Requiring that semantic atom selects one
base-13 congruence class, kills all outer nonzero-root edges, and reveals the
same slope-seven dynamics one odometer level lower.

## 1. Frozen objects

Retain

```text
p=13,                        R=p^6=4826809,
x=649039434905733/1304692766858936,
z={13x}=46873542509301/100360982066072,

I=(960117507257/1930018885886,
   324519717452867/652346383429468),
|I|=1/652346383429468.                                  (1)
```

The THM-2707 packet address is

```text
q_n(x')={13x'}+7n/R.                                    (2)
```

Its six-label skeleton is `(rail,sector,edge,kappa,h,shallow)
=(0,0,0,1,6,1)`, while

```text
carry(n)=2+7n,                   root(n)=6+n             (mod13). (3)
```

Use the exact THM-2680 atoms

```text
current:   (rail,future,j,h,epsilon,kappa)=(2,3,5,2,0,0),
following: (rail,future,j,h,epsilon,kappa)=(0,1,2,6,1,1). (4)
```

The following atom has deep root `r=-h-1=6` and exact positive numerator

```text
867661831383617727737280.                                (5)
```

All "following" statements below refer to this full weighted atom: rail,
present packet, sharp deep root/tooth half, predecessor digit, future half,
and delayed word.  They are stronger than retaining the THM-2707 skeleton.

## 2. The semantic congruence lock

If `(4)` holds at address `n`, its predecessor and root labels must agree
with `(3)`.  Hence

```text
2+7n=2,                    6+n=6                         (mod13),
```

and therefore

```text
n=0 mod13.                                                   (6)
```

Conversely, scan the `13^5=371293` addresses in `(6)` on the exact common
grid.  Precisely `304` retain the packet and the frozen following atom at
the midpoint, and every one retains the **whole** open cylinder `I`.  No
midpoint hit is lost at a boundary.  Relative to the full lift universe,

```text
outer nonzero lifts                         4455516;
label-compatible delta=7 coset               371293;
same-coset support failures                  370989;
semantic whole-I survivors                       304;
other label-incompatible lifts              4084223.      (7)
```

Relative to THM-2707's `3042` terminal packet survivors, exactly `304` are
semantic and `2738` are skeleton-only.

The following atom has `304` weighted pieces and `304` merged support
components.  Each component contains exactly one semantic address.  The
weight census is

```text
27582102210:266,                 27581135604:38.           (8)
```

After subtracting the radius of `I`, the minimum following-base slacks are

```text
left  =368024580/2197,
right =10574833707900/371293.                              (9)
```

The delayed pullback is sharp: its containing interval is exactly the
pullback of `I`, so both residual slacks are zero.  The frozen current atom
also contains all of `I`; its base/prefix slacks are

```text
(344122815960/28561, 80124658752600/4826809,
 28648620,              312931080).                       (10)
```

Thus the `304` restricted current/following supports share a genuine common
base and form a full support nerve

```text
Delta^303.                                                 (11)
```

The positive rail weights in `(8)` certify support only; they are not
promoted to canonical transition-current amplitudes.

## 3. Outer semantic dynamics is discrete

A required nonzero-root THM-2657 edge has translation numerator

```text
k=7(b-a) mod R,                  13 not dividing k.        (12)
```

Every two vertices in `(11)` have `a=b=0 mod13`, so `(12)` is divisible by
thirteen.  Therefore the induced **nonzero-root** graph on the semantic locus
has

```text
304 vertices,             0 edges,             304 SCCs. (13)
```

There is no outer semantic cycle.

The THM-2707 terminal address `n_0=110` is not semantic.  It nevertheless
has one packet edge in each direction to every semantic address, giving a
bidirected `608`-edge star after adjoining it.  Every resulting two-cycle
loses the following atom at `n_0`.  The shortest terminal-to-semantic lifts
are

```text
k=10  -> n=2758288,
k=-3  -> n=2068742.                                      (14)
```

Both have root increment `delta=7`; equation `(14)` does not supply a
semantic loop.

## 4. Root-neutral dynamics descends one level

It would be false to say that the semantic locus has no physical dynamics.
Write a semantic address as

```text
n=13j.                                                    (15)
```

For two such addresses, the physical translation numerator is

```text
k=7(n'-n)=91(j'-j).                                      (16)
```

It is nonzero when the vertices differ, but it is root-neutral and lies in
the kernel of the outer THM-2657 quotient.  Admitting all nonzero kernel
translations makes the semantic graph the complete directed graph on `304`
vertices, with

```text
304*303=92112                                             (17)
```

edges.

Divide `(16)` by the forced factor `13`.  Over denominator `13^5`, the
descended gain is

```text
k/13=7(j'-j).                                            (18)
```

The exact-valuation-one edges are those with `j'!=j mod13`.  The inner
residue parts have sizes

```text
(23,27,22,27,22,23,22,23,23,23,23,23,23).               (19)
```

Hence the valuation-one graph is the complete directed thirteen-partite
graph with

```text
85274 edges,             one SCC,             diameter <=2. (20)
```

The unordered semantic-pair valuation census is

```text
v_13(k):       1:42637,       2:3222,       3:159,       4:38. (21)
```

In particular, the semantic locus has the explicit root-neutral physical
triangle

```text
0 -> 13 -> 26 -> 0,
k=(91,91,-182).                                          (22)
```

All three vertices carry the same current/following cylinder `I`, and the
signed numerator sum is exactly zero.  This is an inner support cycle, not a
canonical reapplication of THM-2657 at depth five.

## 5. Coboundary upstairs, class seven downstairs

On the complete THM-2707 address graph, define the vertex potential

```text
phi(n)=7n mod R.                                          (23)
```

Then every physical edge gain is the exact coboundary

```text
k_(a,b)=phi(b)-phi(a)=7(b-a) mod R.                       (24)
```

Thus every full-address cycle is balanced and every pure circle endpoint
character has trivial translation-phase holonomy.  This does not split the
outer quotient.  Two addresses with the same residue but separated by
thirteen have

```text
phi(n+13)-phi(n)=91,                                     (25)
```

which lies in the outer kernel, while its divided kernel class is

```text
91/13=7 mod13.                                           (26)
```

Therefore `phi` cannot descend through `n mod13`; `(26)` is exactly the
THM-2657 cocycle class.  Equations `(23)`--`(26)` explain both phenomena at
once: the address-resolved groupoid is balanced, but forgetting the address
creates the nonsplit quotient.  On the semantic fibre `(15)`, dividing by
the forced factor recovers the inner coboundary `(18)`.

## 6. No hidden same-current sibling

Fix the current endpoint `h_current=2`, hence following predecessor `j=2`,
and keep rail `0` and future clock `1`.  The sharp graph has exactly four
formal compatible choices

```text
(h,epsilon,kappa,r)
 =(6,1,1,6), (7,0,1,5), (7,1,0,5), (8,0,0,4).           (27)
```

Their whole-`I` address counts are

```text
304, 0, 0, 0.                                            (28)
```

The `h=7,8` delayed prefixes miss `I`; two candidate atoms also have zero
integrated numerator.  Thus `(4)` is the unique whole-cylinder semantic
branch in this fixed-current/rail/future-clock box.  Equation `(28)` says
nothing about other rails, future clocks, current atoms, or configurations.

## 7. Consequence and stopping boundary

The support/semantic distinction is now exact:

```text
fixed skeleton:   Delta^3345 and one outer lift SCC;
frozen following: Delta^303 but no outer nonzero-root edge;
inner kernel:     a 304-vertex physical support SCC one level lower.      (29)
```

This is controlled forgetting in both directions.  Forgetting the semantic
atom creates false outer dynamics; forgetting the full address creates the
class-seven quotient obstruction.  Retaining both does not yet produce a
signed or complex current.  The live target is a common-gauge amplitude or
endpoint coefficient on the inner cycle `(22)`, or a lawful semantic branch
outside the fixed box `(27)`.

No canonical transition amplitude, source/owner identity, target action,
all-row theorem, row exclusion, or proof of LRC(14) follows.

## 8. Exact reproduction

Run

```bash
python3 04-computation/lrc14_semantic_following_congruence_lock_thm2712.py
python3 -O 04-computation/lrc14_semantic_following_congruence_lock_thm2712.py
```

Both executions byte-match the stored `324`-line transcript
`05-knowledge/results/lrc14_semantic_following_congruence_lock_thm2712.out`.
The companion has no optimized-away assertions.  It checks the label gate,
scans all `13^5` compatible addresses, proves whole-cylinder containment,
prints every semantic row, verifies the component/weight/sibling censuses,
and exhausts all `46,056` unordered semantic pairs for `(21)`.

---
id: HYP-2242
status: OPEN method hypothesis with finite order/toy evidence
source: codex-2026-06-05-S667
related:
  - HYP-2241
  - HYP-2240
  - HYP-2239
  - HYP-2237
  - HYP-2171
  - HYP-2167
  - HYP-2165
  - HYP-2164
---

# HYP-2242: Embedded Maximality Order Carrier

## Claim

Maximality should be treated as an embedded predicate:

```text
maximal(object, ambient embedding, allowed extensions).
```

A scalar or quotient maximum is only locally meaningful.  It becomes proof
meaningful only after naming the ambient extension types that can destroy it and
the address coordinate that blocks or detects those extensions.

The order toy is `(Q,<)`.  A finite embedded chain has a last element, but the
ambient dense order can insert a rational above it or inside any open cut.
Finite denominator searches also have best lower approximations, while the
Dedekind cut itself has no rational maximum.  Completion works by adding the
missing cut address.

The LRC14 analogue is:

```text
visible Res_27 maximum
+ C=27 gcd/product shell
+ owner/carry/deletion ambient address
=> embedded maximum or positive tax.
```

So HYP-2241's private-owner deletion flag is not just another invariant.  It is
the first checked embedded-maximality address: it names which visible floor
atoms remain maximal after the owner-deletion ambient is allowed to act.

## S667 Evidence

S667 adds `04-computation/embedded_maximality_order_atlas_s667.py` with stored
output in `05-knowledge/results/embedded_maximality_order_atlas_s667.out`, plus
reflection `07-reflections/embedded-maximality-order-carrier-s667.md`.

The script records three order labs:

- Finite chain `[0,1,2]` has maximum `2`, but one-point dense-order extensions
  realize all four cuts `(-infty,0)`, `(0,1)`, `(1,2)`, and `(2,infty)`.
- The finite-denominator search below `sqrt(2)` up to denominator `12` has best
  lower approximation `7/5`, but denominator `17` supplies the better lower
  approximation `24/17`; the irrational cut has no rational maximum.
- The same rational point `1` is maximal in `{0,1/2,1}` and in `[0,1]_Q`, but
  not in `(0,1)_Q` or in all of `Q`.

This separates three notions that were easy to blur in the LRC work:

1. internal maximum inside a chosen finite image;
2. boundary maximum inside a named compact/closed ambient;
3. embedded maximum after all allowed extension cuts are admitted.

S667 then imports the S666 projection audit:

| Repair key | Groups | Mixed fibers | Max bucket |
|---|---:|---:|---:|
| visible shadow | `3` | `3` | `378` |
| visible + cheap pair | `112` | `3` | `106` |
| visible + owner cover count | `976` | `2` | `2` |
| visible + owner private flag | `1067` | `0` | `2` |
| visible + owner private count | `1134` | `0` | `1` |

The embedded-maximality reading is: the visible shell is the finite chain; local
`+27` carries are ambient extensions; owner-private deletion data is the cut
address that prevents the finite maximum from leaking.

## Transfer Lanes

The S667 lane atlas ranks embedded-maximality carriers by exact evidence,
ambient extendability, address need, derivative leverage, LRC transfer, and
actionability.

Top lanes:

- `LRC14 owner-private stability`: total `30`, vector `(5,5,5,5,5,5)`.
- `Q,< finite-chain maximum`: total `28`, vector `(5,5,5,4,4,5)`.
- `Dedekind cut completion`: total `28`, vector `(5,5,5,4,4,5)`.
- `Endpoint protection core`: total `28`, vector `(4,5,5,5,5,4)`.
- `Tournament endpoint private child`: total `27`, vector `(4,5,5,5,4,4)`.

Other useful lanes are graph reconstruction rooted decks, unit-distance frontier
owners, A000568 marked observer fibers, Rado/Fraisse extension properties,
matroid basis activities, and forcing model maximality.

## Tournament Analysis

Vertices are embedded-maximality lanes, not runners or arcs.  The pairwise
observable is:

```text
(exact evidence,
 ambient extendability,
 address need,
 derivative leverage,
 LRC transfer,
 actionability).
```

The switch is majority comparison of these coordinates, with curated lane order
as the tie Hamiltonian path.

Fingerprints:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}`
- `directed_3cycles=0`
- `scc_sizes=[1,1,1,1,1,1,1,1,1,1,1]`
- `hamiltonian_paths=1`

This transitive result says the current concept is more definition-sharpening
than analogy hunting.  The clean leader is LRC14 owner-private stability.

## Assumption Challenge

S667 did not assume vertices must be runners.  Candidate vertices included
rationals, cuts, finite chains, runners, endpoint cores, owner obligations,
graph cards, point-set frontiers, matroid bases, forcing models, and proof
obligations.  The selected vertices are lanes because the preserved predicate is
"which ambient address makes a local maximum stable, and which extension
destroys it?"  The quotient deliberately destroys raw speed identity, phase
order, and complete geometric position.

## Consequences for LRC14

Embedded maximality suggests a better shape for the next check:

```text
For every allowed extension inside the visible Res_27 fiber,
either the extension changes an owner/deletion/carry address,
or the maximin score rises above 1/14,
or the extension is a globally coherent scalar floor lift.
```

This lets us ignore clocks that only move inside an already separated embedded
fiber.  The clocks that still matter are the clocks that define extension cuts:
carry cocycles, owner-private deletion bits, endpoint protector ownership,
gcd-shell/lift addresses, HYP-2165 route labels, and any coherent scalar lift
that keeps the floor by symmetry rather than by local leakage.

## Next Tests

1. Turn HYP-2241 into an embedded-extension theorem over coherent carry
   subspaces, not just Hamming balls.
2. Build a Dedekind-cut-style completion of the `Res_27` quotient: list all
   extension cuts that remain after visible/gcd/product projection and name the
   smallest address coordinate that separates each.
3. Try the same embedded-maximality audit for tournament deck collisions and
   unit-distance frontier owners: a local maximum should be accepted only if it
   survives the rooted/deleted ambient.

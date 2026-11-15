# Recontact, nine vertices, and exact tail return

This session works on three distinct objectives: a literal recontact that preserves the next reader, placement of consecutive odd counts on nine vertices, and an exact analysis of what a tail can accumulate. It gives new physical windows and scoped obstructions, not a proof of [Chan–Pak](https://arxiv.org/abs/2308.10214).

## A recontact can now keep the continuation intact

- [Single edge and target inverse](SINGLE_EDGE_AND_TARGET.md): reversing only the interior spine edge 2→0 to 0→2, with **every apex contact fixed**, changes 25 to 23 while keeping the right state (h,S,d)=(9,41,7). A common guard gives strong 71 and 73; every common transitive tail preserves the exact gap. A single changed apex arc alone cannot do this, because it changes the positive start weight S.
- [Two-sided recontact](TWO_SIDED_RECONTACT.md): two apex contacts on one regular five-vertex spine change 41 to 43 while keeping **both** ordinal interfaces. Arbitrary common transitive tails on either side preserve +2. Two tail addresses x,y interact by a literal 120xy term. A different matched triple at exactly nine vertices gives 999,1001,1003 and keeps both interfaces.
- [Recontact energy and Bézout bridge](RECONTACT_ENERGY.md): inside a matched mask fiber the count is exactly a quadratic induced-edge energy. Two contact switches interfere through a computable mixed term. Same-spine masks give a three-odd window; two cyclic modules then combine opposite weighted changes 15·4−29·2 to achieve a literal +2 at common continuation state (435,20658,174). A final guard makes the strong seeds 34475 and 34477. This gives a bounded physical Bézout return, not a claim that every gcd-compatible correction can be assembled.
- [Two-tail capacity and decoding](DOUBLE_TAIL_CAPACITY.md): addresses give order X log X presentations, not necessarily distinct values. The target-driven exact inverse costs O(sqrt X log X) candidate tests. The guarded order-nine triple collectively reaches every odd-compatible residue modulo every modulus even at fixed cycle count 45, but exact absolute coverage still does not follow.

## What nine really means

[The exact order-nine census](NINE_ORDER.md) separates the number of odds in the window from the graph order. All nine values on **each of the first two rungs** (1255–1271 and 2275–2291 by twos) first occur on nine vertices; later rungs cannot fit at order nine. The nine counts do not share one order-nine right interface: the first rung has maximum matched width two, or three after neutral-one preparation. Across all marked order-nine tournaments the exact maximum width is four before that preparation and six after it. The six can be completed strongly and survive every common tail. These are exact maxima under the named order and preparation, not a general bound on window width.

## Exact tail behavior: rhythm, folds, and full history

- [Periodic tails](PERIODIC_TAILS.md): a repeating word produces an exponential term **plus an exact linear drift**, with the sharp minimum and maximum words at each length. Matching the periodic exponential rate is not sufficient to make two actual counts return. At the normalized slope level, a finite limit set is exactly eventual periodicity; every binary block recurring gives the whole limiting interval. This is an analogy with rhythm versus long memory, not a biological claim.
- [Fold overlays](TAIL_FOLDS.md): neighboring addresses have an exact, sign-changing difference. On the previous four-count matched seed, different phases of that fold give five, six or eight consecutive odds at infinitely many preparation heights. The first eight-count overlay is 150401–150415 by twos. The two folded addresses end at different right states, so this expands actual coverage but not matched continuation width.

[The fold-closure obstruction](FOLD_CLOSURE.md) makes the next move precise: the full terminal start weight actually reconstructs every chronological mark of the final ordinal tail. The two folded states also cannot merge by appending separately chosen transitive words. The next target is to **close the wider folded window**: recontact after the two tails so that their exact terminal readouts agree without destroying the small gap, or supply a target-dependent two-sided palette whose absolute sums fill an interval. The two existing aims, matching a state and increasing absolute window width, are now both realized but not yet in one unbounded construction.

## Reproduction

The two C++17 census programs generate canonical order-nine representatives and all marked spine presentations. Their receipts include the literal extremal rows. The companion Python scripts check independent path counts, exact finite arithmetic and the named windows. The asymptotic and all-tail statements follow from the proofs, not extrapolation from the replays.

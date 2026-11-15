# Exact carried targets, richer readers, and matched windows

This session starts at what the target needs to read and works backward to actual arcs and cores. There are exact infinite constructions and all-order rejection tests below. [Chan–Pak](https://arxiv.org/abs/2308.10214) remains open.

## Carry and inverse

- [Absolute carry](ABSOLUTE_CARRY.md): in the physical binary-gap gluing, the canonical odd coefficient λ and checksum c satisfy 1≤c≤λ, even when c crosses the gap. A prescribed N and gap g give the exact interval ceil(N/(2^g+1))≤λ≤floor((N−1)/2^g). For each λ the full physical interface is finite at all orders; reconstruct the suffix from its reversed bits, and enumerate only the target-forced core budget. A literal count 3841 exhibits carry 639.
- [Prefix lines](PREFIX_LINES.md) and [equal response](EQUAL_RESPONSE.md): unequal fixed responses meet a prescribed gap at only O(log² X) prefix addresses through X; dyadic slopes give finitely many explicit dyadic rays and finite exceptions. Nonzero start defect gives exact +2 strips at 90·2^s+{7,9} and 150·2^s+{7,9}. Equal responses admit wider scope: true +36 and −758 locks hold for every prefix, including arbitrarily large carries. These are not same-interface +2 locks.
- [Response inverse](RESPONSE_INVERSE.md): a prescribed full prefix response (L,t) leaves at most one physical d candidate at each suffix width. Factoring t and removing singleton SCCs then gives a complete decoder, with no size search growing with L when the t-bank is precomputed. At (252,9), the complete constant list is {207,243}.

## Admit a richer core by changing what the host actually reads

- [Outlet width](OUTLET_WIDTH.md): uniform substitution through at most d outlets reads exactly the first d+1 ordered slot counts, calibrated by literal nested triangle gates. The first new read has a sharp d! excursion multiplier. [Rich core admission](RICH_CORE_ADMISSION.md) gives the stronger all-uniform-host factorial tariff at fixed lower profile. The old order-eight pair gives a sharp uniform gap of 1440, but one selective apex distinguishes the same two cores by **1939 versus 1941**. That switch by itself does not match the whole continuation interface.
- [Mixed-contact windows](MIXED_CONTACT_LOCKS.md): targeting oddpart(S+h) selects actual painted cores. There are **eight** consecutive odd counts 17·2^j+{23,25,…,37} at every j≥4, and **nine** at 255·2^j+{235,237,…,251} for every j≥2. At slope 75 the triple includes a profile absent from the old ordinal one-apex alphabet. A separate physical negative-low construction realizes 21·2^j−17 and 21·2^j−15 with a real borrow at every j≥5. Fixed-k guarded higher-apex cores have slot growth (j+k)^n; the resulting rich isomorphism-type supply is not yet a proved supply of distinct profiles.

## Preserve a window, then try to widen it

- [Matched windows](MATCHED_WINDOWS.md): four **strong** seeds start at 163,165,167,169 with exactly the same h,S,d. This four-odd window survives **every common transitive continuation**, not just a specified pattern. There is also an independent adjustable preparation height, giving the exact two-parameter family 27+17·2^j(2t+1)−30t+15p(t)+{0,2,4,6}, j≥3 and t≥0. A complete tail-alignment criterion explains why differing widths can synchronize and how to choose them backward.
- [Multiplier-three obstruction](MULTIPLIER_THREE.md): trying to widen this window by an ordinal extension of Hamilton multiplier 3 reduces exactly to replacing a two-bit subword by one painted triangle. Once the final start is fixed, the effective word is fixed, and distinct offsets are at least 2h apart. After scaling the four-count window, they do not even produce an adjacent pair at a common final state. This excludes that primitive ternary widening rule at all lengths. A genuinely new extension must change the multiplier, recontact old vertices, or use a different retained interface.

The strongest positive conclusion is a real matched window under arbitrary common growth. The remaining coverage problem is to move that window to an arbitrary absolute target or widen it by an operation that respects the typed endpoint state. An unbounded number of exact heights on a fixed dyadic lane does not by itself do either.

## Replay

Run the companion `replay_*.py` files from this directory. Exact receipts are included; the independent literal controls use subset-endpoint DP. Infinite and impossibility claims come from the stated arguments, not from the bounded experiments.

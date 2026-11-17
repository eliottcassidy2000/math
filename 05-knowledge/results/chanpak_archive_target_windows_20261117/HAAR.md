# A unit repair from the mixed boundary term, not from the margins

The randomly sampled [Haar-product reflection](../../../07-reflections/lrc14-haar-product-discrepancy-tournament-tiling-codex-s165.md) isolates a signed 2×2 term that two marginal reads do not see. The ordinal marked-cover transfer has an exact analogue: switching order of two pieces changes a mixed product even when their marginals are identical.

For a marked block T let (h,P,W,Z,E,S,K) be its usual passport: H, optional two-slot cover count, weighted marked starts and ends, full marked ends and starts, and K=h−U. For X→Y:

    K = h_Y K_X + h_X K_Y − h_X h_Y − E_X S_Y − W_X Z_Y.

The last two terms distinguish the endpoint order: the E/S term is a join across adjacent blocks; the W/Z term is a pair with reversed block roles for the distinguished endpoints. If two tiles A,B share h,P,W,Z, swapping their positions keeps (h,P,W,Z) of the pair fixed. Only −E_A S_B versus −E_B S_A survives in the order difference in K. This is the exact oriented determinant rather than the individual endpoint margins.

Use the same underlying tournament with outgoing masks (20,121,114,245,32,1,49,119) and tags 77 or146 as the two tiles. Both have (h,P,W,Z)=(51,1512,705,534). The endpoint pairs (E,S) are (15,22) and (19,28); their determinant is 15·28−19·22=2. Their differing individual K terms cancel on swapping. Put the two tiles between a transitive two-vertex prefix (first vertex only tagged) and an untagged singleton suffix. The prefix has E=0 and the suffix S=0, so neither adds an order-sensitive adjacency term. The full outer E and S are fixed, the root can reach the leading tag, and the terminal reader returns all vertices to the root.

The two full marked hosts have common W=11348406 and respective K=−2064993,−2064991, certified by independent subset read and ordinal composition. Any common transitive return TT_k gives the exact two-count window 11348406·2^k−2064993 and two more. At k=2 the upper member is the prime 43328633. This all-height pair has only a dyadic baseline; it is not an absolute selector. What transfers to a future selector is the ability to move by two solely through reordering tagged pieces, without a change to their measured marginal data. See `replay_haar.py` for literal masks and strength.

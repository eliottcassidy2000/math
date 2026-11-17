# Prime-compatible transfer of the full-P₂ repair

The [one-arc +12 family](UNIT_SEAM.md) holds P₂ fixed. The reader 7(64P₂−H) transfers it to a count gap but does not give prime quotients on the growing branch: the source has even order, so P₂≡H mod3 and the quotient is divisible by 3. This distinguishes exact response from prime-compatible response.

Instead apply the [physical 7(5P₂−3H) reader](../chanpak_changed_contacts_20261117/HETEROGENEOUS.md), which works for an arbitrary variable source. Its quotient changes by 36 (actual count gap 252). At u=470 the two sources have Hamilton counts 1812827149 and 1812827161 and common P₂=847903347292696. The resulting quotients 4239511297982033 and 4239511297981997 are both prime, with recursive order certificates checked by the replay. Thus the exact invariant can be carried to two literal strong sevenfold prime witnesses, not just to composite examples. There is no infinitude or prescribed-prime claim.

Replay: `python replay_prime_repair.py` calculates the literal word profile and checks full factorizations of p−1 recursively with order witnesses.

# Five deeper targets from last-contact inversion, not a census

Use the [literal last-contact inverse](../depth/PROOF.md): for each specific order-11 or order-12 core, the exact odd-cycle cut law evaluates all its apex contacts. The outer loop proposes only a small number of deterministic seeded cores, stopping at the target. This does not enumerate tournament classes at the output size.

The independently checked literal masks are in `deep_targets.out`. The [order-eleven audit](https://github.com/eliottcassidy2000/math/blob/eea6f39cc271ab196ca99b422ad86fb320d8b903/05-knowledge/results/first_creator_generation_20261108/audit/AUDIT.md) gives global maximum H=95,095 at order 11. Since padding preserves H, every listed target is therefore first possible at order 12. For four of them, the resulting witness itself has order 12, so the lower and upper certificates meet exactly.

| target | interpretation | literal order | certified first-order conclusion |
|---:|---|---:|---|
| 100,003 | prime | 12 | exactly 12 |
| 210,077 | 7·30,011 | 12 | exactly 12 |
| 210,091 | 7·30,013 | 12 | exactly 12 |
| 458,759 | 7·65,537 | 12 | exactly 12 |
| 700,021 | 7·100,003 | 13 | 12 or 13 from these certificates |

The first three successes occur at zero-based attempts 430, 188, 25 respectively on order-11 deletions. For 458,759 use the separate near-circular/Paley proposal followed by exact contact inversion; it hits at attempt 229 on an order-11 deletion. For 700,021 an order-12 deletion gives success at attempt 199; no order-12 absence is claimed. With at most 16,000 near-circular/Paley samples, the order-12 outputs 681,163 and 700,021 did not land, which is only a sampled miss. “Exactly 12” concerns first order, not uniqueness or a complete list of first creators. The arithmetic confirms a prime and three independently chosen 7p values at the new frontier.

Reproduce from this directory:

    g++ -O3 ../depth/random_apex.cpp -o /tmp/pair_apex
    /tmp/pair_apex 11 431 100003
    /tmp/pair_apex 11 189 210077
    /tmp/pair_apex 11 26 210091
    g++ -O3 near_apex.cpp -o /tmp/near_apex
    /tmp/near_apex 11 230 458759
    /tmp/pair_apex 12 200 700021
    python deep_replay.py

The replay is separate ordinary subset DP over the literal completed graphs, and uses trial division for the indicated primes.

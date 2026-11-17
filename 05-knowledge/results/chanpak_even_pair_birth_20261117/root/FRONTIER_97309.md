# An exact hard birth with two certified odd midpoints

The [order-eleven class audit](https://github.com/eliottcassidy2000/math/blob/eea6f39cc271ab196ca99b422ad86fb320d8b903/05-knowledge/results/first_creator_generation_20261108/audit/AUDIT.md) certifies the maximum H=95,095 at order 11 and H=15,745 at order 10. Adding a transitive singleton does not change H, so each of these is also the maximum up to its order. The first-birth [order-nine atlas](../chanpak_birth_atlas_metagraph_20261117/metagraph/witnesses.tsv) gives maximum H=3,357 at order nine.

[Session 244's certified depth-three witness](../session244_even_gap_anchoring_20261110/logic/forest_census/release/DEPTH_EXCEPTION_CHECK.json) is

    (3910,2044,3576,3057,4001,4033,3728,3841,3648,3076,2056,2).

It has 97,309 paths at order twelve. Combining these two inherited certificates gives the exact global first order **12**. This does not require a census of all order-twelve tournaments. The factorization 31·43·73 gives a simpler depth-one construction at order twenty; it is not birth-optimal.

Now delete the pair {5,8} from this one difficult structure (indices are zero-based). The result has H=5,803 and order ten. Removing just 5 or just 8 gives 23,395 and 24,041 at order eleven. Both are above the certified order-ten maximum; the order-ten child is above the certified order-nine maximum. Thus the paired deletion diagram contains **four genuinely first-order vertices**:

                order 12:   97309
                           /       \
                order 11: 23395   24041
                           \       /
                order 10:    5803

The two odd-order midpoints are not arithmetically averaged from the even levels: their marked correction is Γ=97,309+5,803−23,395−24,041=55,676. About half the upper count sits in the pair interaction. On this literal witness, seven of twelve one-deletions clear the previous-order maximum and twenty-seven of sixty-six two-deletions clear theirs. The other deletions are not classified as births by this test alone.

This is a targeted deep first-order certificate and an application of the marked even-pair identity, not a claim that every difficult birth has paired first-birth ancestors. The new signed-cut template in this packet independently reaches 97,309 at order 15; that is a construction control and does not improve the global minimum.

The [depth-first witness](../depth/targets7.jsonl) for 681,163=7·97,309 together with the [literal order-13 search](../depth/PROOF.md) gives order ≤13, so its known-first-order interval from these sources is [12,13], not an exact claim.

# Local audit repair: split source charges are not split node currents

23 October2026. **CORRECTED BEFORE COMMIT.** The root agent independently
identified a missing diagonal term in this packet's first draft debt ledger.
All state counts, complete-fibre tests, paths, connectivity and *total*
current preservation were already correct. The names and values of the
individual old16/helper contributions were not.

The first draft called `2(old_X_sum-1016)` and `2(helper_X_sum)` physical
target-owned currents. These are source-charge contributions. In an
intermediate state the old16 group has seven V-slots instead of eight,
while the helper group has one instead of zero. The individual groups
therefore change the target degree diagonals, even though their combined
margin remains fixed.

For a partial source change w, actual paired edge transport U to V gives

    DeltaL_ss=0 at a source s,
    DeltaL_UU=-2 sum(w), DeltaL_VV=2 sum(w),
    DeltaL_Us=2w_s, DeltaL_Vs=-2w_s.

Consequently

    (DeltaL X)_U = 2 sum_s w_s (X_s-X_U),
    (DeltaL X)_V = -2 sum_s w_s (X_s-X_V).

Only when the partial margin itself is zero may the X_U/X_V terms be
dropped. Here X_U=X_V=146, so a helper of X=139 contributes **-14 at U**,
and the old group contributes **+14**. The former -278/+278 entries are
retained under explicitly named source-charge fields. A helper of X=173
instead contributes **+54 at U**, balanced by **-54** from the old group.

The same correction applies to the four binary effective character rows:
subtract the target code times the partial margin. In characteristic2 this
is an XOR when that margin is odd. The physical old/helper binary residuals
are respectively13/13 for helper1923 and8/8 for helper1532; the separate
source-code contributions are0/0 and5/5.

The repaired producers retain source-charge sums, partial V-counts, target
degree changes, physical node currents and binary residuals separately.
The checker recomputes every one and verifies their sum at the actual U/V
nodes. All receipts were regenerated. No altered connectivity, route,
resource capacity or principal theorem was needed.

**Mechanism of failure:** using a balanced-block formula on an unbalanced
piece of that block. **Strongest survivor:** full-block cancellation and all
proved paths. **Missing sidecar:** each piece's actual target-degree margin.

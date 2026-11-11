# Independent narrow audit: fixed-q unit-offset seed decision

**Verdict: PASS, conditional only on the explicitly cited strong-order floor input.** The floor itself was not re-audited. The statement decides exact U=q within the selected-set/transitive-tail constructor; it does not decide all ray constructors or every seed with odd part c=q.

## Proof and implementation review

1. **Bound.** W counts the ordered covers with both the first end and second start selected; U includes all such covers and also other first-end states and an empty first path. Thus W≤U. Acceptance H−W=1,U=q forces H≤q+1. SCC decomposition is an ordinal join and H is the product of the positive SCC counts, so each SCC count is ≤q+1. The supplied Busch floor therefore bounds its order. In the script this is `bound(q+1)`; distinguish this from a convention R(q) that already includes the +1.

2. **Nine-state meaning.** Index 3i+j stores the first path's end category i and the second path's start category j, each with 0=empty,1=unselected,2=selected. Every ordered cover of an ordinal join restricts uniquely to an ordered cover on each block, allowing empty paths. Conversely these restrictions concatenate uniquely. On concatenation the first end becomes the local first end when nonempty, while the second start retains the earlier start when nonempty. `DEST` implements precisely these operations. The transfer does not omit paths crossing an SCC boundary.

3. **Observables.** H=v01+v02; W=v22; U=v02+v12+v22. These are exactly the intended cover families. The empty word is not accidentally accepted for positive q.

4. **Pruning.** H multiplies by H(new block)≥1. U cannot decrease: append any chosen Hamilton path of the new block to the first path of each old U-cover, using it as the first path when that path was empty. The second path and its selected start remain fixed. This is an injection, and in fact gives at least H(new block) images per old cover. Therefore H>q+1 or U>q can never recover.

5. **Saturation.** Coordinatewise clipping at K=q+2 commutes with addition and multiplication of nonnegative integers after final clipping. Every transfer is a nonnegative integer polynomial. Hence clipped transitions are exact quotient transitions, including zero coordinates. At an accepted state H≤q+1 and U=q force each coordinate used by H or U below K; W is one of the U coordinates. Consequently acceptance cannot be created by saturation. Discarded information in other coordinates is harmless because every later transition still respects the quotient.

6. **Completeness and length.** All marked strong SCCs up to the derived order bound, together with both marks on a singleton, form a finite complete alphabet. Transfer-identical letters can be merged. BFS over saturated states decides existence and finds a shortest **block word**. A shortest accepted path has no repeated state, giving at most (q+3)^9−1 blocks and at most `bound(q+1)*((q+3)^9−1)` vertices. This is not asserted to be the minimum vertex order.

## Independent executable path

`audit_independent.py` checks all 146 distinct alphabet transfer vectors against a separate literal algorithm: enumerate full permutations and every split into two directed paths, rather than subset endpoint DP. It independently implements transfer concatenation and exhausts the entire pruned saturated closure for each q=1,...,13 (not stopping at the first winner). Every production witness is recounted by the literal permutation/split method. Explicit RuntimeError checks remain active under `-O`.

Reproduce:

    python3 audit_independent.py > audit_receipt.json
    python3 -O audit_independent.py > audit_receipt_optimized.json
    cmp audit_receipt.json audit_receipt_optimized.json

## Scope caveats

The positive q values are {1,2,4,6,8,11,12}; the negatives are {3,5,7,9,10,13}. These negatives exclude exact U=q for every base order within this constructor. They do not exclude U=q·2^j with j>0, nor a doubling-vertex ray of odd slope q, nor individual Hamilton counts 1+2^t q. In particular no q=7 tail seed says the required simultaneous 15,29,... ray cannot arise from one such seed with U=7, even though 15 and 29 separately occur. A c=9 ray elsewhere is consistent with no exact-q9 tail seed.

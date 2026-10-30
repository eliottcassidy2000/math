# Independent audit: global phase completion (Session 116)

## Verdict

**PASS within the stated scope.** The finite all-entry optimization theorem, primitive equality classification, primitive-root reduction, and rotation corollary follow from the inherited Session 114 ordinary-mean theorem. No universal controller bound, finite-word bound, or Collatz convergence follows. This audit does not establish literature novelty.

Read in full: `../collatz/PHASE_COMPLETION_THEOREM.md`; also checked the membership-depth and ordinary-limit arguments in the Session 114 `collatz/PERIODIC_CONTENT_RATES.md` source.

## Proof-critical checks

1. **Composition and cosets.** With chronological concatenation, `F_xy=F_y F_x`. Writing `U=F_u` and `V_i=F_v_i`, the actual prefix is `V_i M^n U`, so its identity condition is `M^n=(U V_i)^(-1)`. Thus the relevant phase condition is `U V_i in C_e`, not `V_i U in C_e`. If `U V_j in C_e`, then `U=(U V_j)V_j^(-1)` and membership is equivalent to `V_j^(-1)V_i in C_e`. These are left cosets `V_j C_e`. No normality is assumed or needed. The second-window conjugated subgroup `D C_e D^(-1)` has the correct orientation.

2. **Maximum depth, infinity, strictness.** The finite phase set has an attained maximum depth, even if this is infinity. For a maximizing phase j, completion `S_j=M V_j^(-1)` preserves every lower-layer hit set. Above a finite maximal depth, the old entry has no hits while completion keeps phase j forever. The positive convergent tail gives a strictly larger ordinary rate. If maximal depth is infinite, all layers agree. This proves an attained finite-candidate maximum, not merely a supremum approximation.

3. **Infinite depth forces word completion.** Membership at every depth implies a rational commutator vanishes (its numerator is divisible by every power of 5). Affine faithfulness and the free-monoid commuting-word theorem then apply. For primitive w, `v_j u=w^b`; because j is in `1..d`, the left word is nonempty and b is at least one. Cancellation yields `u=s_j w^(b-1)`. Conversely, multiplying a completion entry by powers of w preserves every subgroup test. There is no need to assert that every rational centralizer element belongs to every cyclic subgroup.

4. **Primitive phases and duplicates.** The word `v_i s_j` has length `d+i-j`, strictly between zero and `2d`. If it commutes with primitive w, its length must equal d, forcing i=j. Thus all off-diagonal depths are finite. The entry languages for distinct primitive phases are disjoint: their lengths are distinct modulo d. Empty entry is correctly included by j=d and exponent zero. For a nonprimitive controller, passing to its primitive root avoids erroneous phase uniqueness claims; `S(N;u,z^k)=S(kN;u,z)` is exact and the rate scales by k.

5. **Symmetry and ultrametric.** Inversion in the subgroup gives symmetry. Multiplying `V_i^(-1)V_j` and `V_j^(-1)V_k` gives transitivity at each fixed layer. Therefore nested partitions give the stated depth ultrametric inequality. This is a depth similarity, not the usual distance convention; a distance can be formed as `5^(-Delta)`, with infinity mapped to zero.

6. **Ordinary mean and rotation.** Finite depths give bounded discrepancies. For infinite depths, positive commuting words are powers of a common primitive word, and LTE bounds terminal heights by `O(log N)` with a nonzero positive exponent. This validates the inherited ordinary mean and avoids an unjustified infinite-sum/limit exchange. The rotated and completed words share the same infinite letter stream, with only fixed-length end boundaries; their terminal charges are `O(log N)`. Consequently rotation equality is valid for ordinary means, but does not imply equal finite-N sums or pointwise entry domination.

## Independent computational evidence

`check_independent.py` imports no production checker. It constructs rational affine maps directly, enumerates finite affine cyclic groups, and separately traverses chronological prefixes over full periods. The arithmetic depth formula is implemented independently in rational coordinates, and checked against literal membership; the direct chronological test provides an additional check on composition orientation.

Exhaustive domain: all 63 nonempty positive words of clock at most 6 and all 32 entries of clock at most 5. Membership checks use moduli 5, 25, 125, and 625. Direct chronological period counts use moduli 5, 25, and 125, with a smaller deterministic set of entries plus every completion.

**48,736 explicit exception-gated checks pass**, including:

- 24,576 literal group-membership checks;
- 1,716 direct chronological full-period charge counts;
- 8,064 simultaneous-completion layer comparisons;
- 2,016 optimizer-language equivalences;
- 2,016 exact primitive-root scaling identities;
- 1,777 strict finite-depth improvement checks;
- 192 rotation-score identities;
- symmetry, infinite diagonal, primitive separation, and ultrametric tests;
- both exact two-reset counterexamples stated in the theorem.

Reproduction:

```sh
python session116_research/audit_collatz/check_independent.py > /tmp/audit-normal.json
python -O session116_research/audit_collatz/check_independent.py > /tmp/audit-optimized.json
cmp /tmp/audit-normal.json /tmp/audit-optimized.json
cmp /tmp/audit-normal.json session116_research/audit_collatz/receipt.json
```

Normal and optimized runs matched byte for byte. Finite tests are evidence against implementation/sign/indexing mistakes, not a replacement for the all-depth proof. Infinity is represented only for an exact zero rational quantity, never inferred from a finite modular cutoff.

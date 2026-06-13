# Positivity Hierarchy for the Redei-Berge Function of Tournaments

**Instance:** opus-2026-03-06-S10
**Status:** Partial results proved, conjectures verified computationally

## Summary of Results

For a tournament T on [n], the Redei-Berge symmetric function U_T admits
expansions in various bases of the ring of symmetric functions. The positivity
properties of these expansions form a strict hierarchy:

```
e-positive => s-positive => h-positive => p-positive
    FAIL         FAIL          FAIL         HOLDS (=OCF)
```

Additionally, we identify a NEW intermediate property:

```
p-positive => hook-s-positive (n>=4) => ... => general properties
   HOLDS         HOLDS (conj)
```

## 1. Even Cycle Vanishing Theorem (PROVED)

**Theorem.** For any tournament T on [n], the power sum coefficient [p_mu]U_T = 0
whenever the partition mu has any even part.

**Proof.** For a permutation sigma with an even k-cycle c, define sigma' by
reversing c. Both have cycle type mu. If c is a directed cycle in T, then the
reversed cycle is directed in T^op. The change in the sign exponent phi is k-1,
which is odd for even k. So sigma and sigma' contribute with opposite signs
and cancel. This involution pairs all contributing permutations with cycle type
mu (when mu has an even part), giving zero total.

**Consequences:**
1. U_T lives in the span of {p_mu : all parts of mu are odd}
2. At n=4: only p_{(1^4)} and p_{(3,1)} are nonzero
3. At n=5: only p_{(1^5)}, p_{(3,1,1)}, p_{(5)}
4. At n=6: only p_{(1^6)}, p_{(3,1,1,1)}, p_{(3,3)}, p_{(5,1)}
5. At n=7: only p_{(1^7)}, p_{(3,1,1,1,1)}, p_{(3,3,1)}, p_{(5,1,1)}, p_{(7)}

**Geometric interpretation:** The involution sigma <-> sigma' is the same
T <-> T^op involution that generates the perpendicular grid reflection in the
tiling model. Even-length directed cycles are "parity-sensitive" to this
reflection and cancel, while odd-length cycles are "parity-insensitive."

## 2. Hook Schur Positivity (CONJECTURED for n >= 4)

**Conjecture.** For n >= 4 and any tournament T on [n], [s_{(k,1^{n-k})}]U_T >= 0
for all k = 1, ..., n.

**Proved at n=4.** The only nonzero p-types are (1^4) and (3,1) (by the Even
Cycle Vanishing Theorem). The hook characters at these types are:

| Hook | chi at (1^4) | chi at (3,1) |
|------|-------------|-------------|
| (4)  | 1           | 1           |
| (3,1)| 3           | 0           |
| (2,1,1) | 3        | 0           |
| (1^4)| 1           | 1           |

ALL entries are non-negative. Combined with p-positivity (OCF), each hook
coefficient is a sum of non-negative terms. QED.

**Verified at n=5** (all 11 isomorphism classes), **n=6** (all 40 classes).
Testing n=7 (242 classes, in progress).

**Why n=3 fails:** All partitions of 3 are hooks. The standard representation
(2,1) has chi^{(2,1)}((3)) = -1. For the 3-cycle tournament, p_{(3)} = 2,
giving [s_{(2,1)}] = 2/6 - 2/3 = -1/3 < 0.

**Why n >= 5 is subtler:** chi^{(n-1,1)}((n)) = -1 (general fact). So
hook positivity requires quantitative bounds: the (3,1,...,1) contribution
must outweigh the (n) contribution. This is a tournament cycle count
inequality.

**Pattern of negativity:** The hook characters at the (n-1,1,...,1) cycle type
(the "single 3-cycle" type) satisfy:
- chi^{(n-j,1^j)}((n-1,1,...,1)) >= 0 for all j (hooks are non-negative here)
- chi^{non-hook}((n-1,1,...,1)) < 0 (non-hooks are negative here)

This is why non-hook Schur coefficients are ALWAYS negative for non-transitive
tournaments: the dominant p-type is (3,1,...,1), and non-hook characters are
negative there.

## 3. Non-hook Negativity (PROVED at n=4)

**Theorem (n=4).** [s_{(2,2)}]U_T < 0 for every non-transitive tournament on 4 vertices.

**Proof.** [s_{(2,2)}]U_T = 2/24 * 1 + (-1)/3 * p_{(3,1)} = 1/12 - 2c3/3.
This is negative whenever c3 >= 1 (which holds for all non-transitive tournaments).

## 4. Palindromic Property (VERIFIED)

For all tested tournaments (n=3,...,6):
[s_{(n-j, 1^j)}]U_T = [s_{(j+1, 1^{n-j-1})}]U_T

This follows from U_T = U_{T^op} (omega-invariance), since conjugating the
hook (n-j, 1^j) gives (j+1, 1^{n-j-1}).

## 5. Connection to Tiling Grid Geometry

The triangular pin grid has a perpendicular axis of symmetry corresponding to
T <-> T^op. The key connections:

1. **Transfer matrix symmetry M[a,b] = M[b,a]** is equivalent to OCF and
   corresponds to the grid's perpendicular symmetry at the level of the
   generating function.

2. **Even cycle vanishing** (p_mu = 0 for even-part mu) is the symmetric
   function version of the "even r-powers vanish" property of the transfer
   matrix (kind-pasteur-S23's even-r-powers conjecture).

3. **Hook Schur positivity** means that the exterior power decomposition of
   the tiling count is non-negative. Geometrically, this says that certain
   "descent-structure" weighted tiling counts are always non-negative.

4. **Non-hook negativity** means that more complex representation-theoretic
   decompositions of the tiling count CAN be negative. The perpendicular
   grid symmetry is "not strong enough" to force positivity for all Schur
   coefficients, only for the hook-shaped ones.

## 6. n=7: Hook Positivity FAILS

**At n=7, hook positivity fails for 11/242 tournament classes.** The ONLY
failing hook is (4,1,1,1) — the "middle" hook (j = 3 = (n-1)/2).

The worst case is the regular tournament (Paley T_7) with scores (3,3,3,3,3,3,3):
- p_{(7)} = 48 (24 directed Hamiltonian cycles, each contributing twice)
- p_{(3,1,1,1,1)} = 28, p_{(3,3,1)} = 28, p_{(5,1,1)} = 84
- [s_{(4,1,1,1)}]U_T = 1/252 + 56/72 + 56/18 + 0 - 48/7 = -83/28 ≈ -2.96

**Structural analysis of failure:**
- chi^{(4,1,1,1)}((7)) = -1 (negative at n-cycles)
- chi^{(4,1,1,1)}((5,1,1)) = 0 (5-cycles don't help!)
- chi^{(4,1,1,1)}((3,1,1,1,1)) = 2, chi^{(4,1,1,1)}((3,3,1)) = 2
- The 48 Hamiltonian cycles contribute -48/7 ≈ -6.86, overwhelming the
  positive 3-cycle terms

**Which hooks are always positive at n=7:**
- (7) and (1^7): chi = [1,1,1,1,1] at odd types → always ≥ 0
- (5,1,1) and (3,1,1,1,1): chi = [1,0,0,3,15] → always ≥ 0
- (6,1) and (2,1,1,1,1,1): chi = [-1,1,0,3,6] → quantitative bounds hold

Only (4,1,1,1) (middle hook, self-conjugate) fails: its positive terms
(dim/n! and 3-cycle contributions) are too weak relative to its negative
n-cycle contribution.

**Refined conjecture:** Hook positivity holds for "outer" hooks (j close to 0
or n) but can fail for "middle" hooks (j ≈ n/2) when n ≥ 7.

## 7. Open Questions

1. **Characterize which hooks are always positive.** Is it true that hooks
   with chi^{hook}(mu) ≥ 0 at ALL odd-part mu are always positive? (Pure sign
   argument.) What about hooks where the quantitative bounds suffice?

2. **Does the failure persist at n ≥ 8?** Check if the middle hook fails at
   n=8 (where the middle hook has even j=3 or 4, so chi at (n) is +1; but
   n=8 is even, so (8)-cycles don't contribute by even cycle vanishing).

3. **Geometric interpretation of the middle hook.** The middle exterior power
   Λ^{⌊n/2⌋} of the standard representation has the largest dimension. Its
   failure at n=7 may relate to the grid geometry in a deep way.

4. **Is hook positivity equivalent to a tournament cycle inequality?** At n=7,
   positivity of (4,1,1,1) requires: 2*c3/72 + 2*c33/18 ≥ 2*h7/7 - 1/252,
   where h7 = #Hamiltonian cycles. This is false for regular tournaments.

## 8. Computational Data Summary

| n | #classes | hook-positive | failing hook | min hook coeff |
|---|----------|--------------|-------------|----------------|
| 3 | 2        | 1/2 (50%)    | (2,1) (all are hooks) | -0.333 |
| 4 | 4        | 4/4 (100%)   | none        | 0.042          |
| 5 | 11       | 11/11 (100%) | none        | 0.008          |
| 6 | 40       | 40/40 (100%) | none        | 0.001          |
| 7 | 242      | 231/242 (95.5%) | (4,1,1,1) only | -2.964    |

**Key pattern:** The transitive tournament always has the SMALLEST positive hook
coefficients (1/n! for s_{(n)}), while the regular tournament (when it exists)
has the most negative middle hook coefficient.

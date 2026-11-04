# The same return tree for every suffix with SCCs of order at most four

2026-11-04. **PROVED complete normalization and finite search bound; FINITE-EXACT complete tree certificate.** This extends the separate Session184 `spectrum/PROOF.md` without changing it. Allowing the unique strong four-vertex tournament, with all its actual port subsets and arbitrarily many copies, adds **no return at all** to the seven-state tree.

## 1. Statement and scope

Keep the root attachment and clock-uniform digit relation of the companion packet. A suffix B now may be **any tournament whose strongly connected components have order at most four**, with any actual port subset meeting its terminal SCC. The outer prefix is TT_t. Its signature still gives

    H_t=c+2^t F,       c=h−J,

and a return with fixed integer clock shift q and digit d in {2,4} requires

    2^q F'=3F,        c'=3c+d.                         (1)

**Theorem.** Starting from marked TT3 with every vertex a port, the complete address language in this enlarged grammar remains

    {empty, 2, 4, 22, 42, 422, 4222}.                 (2)

Each address has the same unique normalized representative as in the companion table. Neither first-generation sibling has a +4 successor, and no lineage has more than four steps, regardless of arbitrary fixed integer clock shifts and arbitrary replacement by any suffix in this grammar.

Consequently a clock-uniform +4 return from either original sibling in this same root-attachment architecture requires **some suffix SCC of order at least five**. This does not exclude a one-clock return or another attachment architecture, and does not claim that introducing a larger SCC suffices. No native guarded-operation construction or old-name retention is asserted.

The component sum law is inherited from Session182; all figures concerning the smaller grammar are credited to the companion proof. The extension is the complete larger-bank normalization and its exhaustive return certificate.

## 2. Complete finite component bank

No two-vertex tournament is strong. The only strong three-vertex tournament is C3. Up to isomorphism, the only strong four-vertex tournament is

    Q=(6,12,8,1),

in outgoing-neighbor-mask notation on names0,1,2,3. For completeness, the checker enumerates every one of the64 labelled four-vertex tournaments, checks strongness by reachability, and exhibits an isomorphism to Q for each of the24 strong ones. This finite classification, together with the linearly ordered SCC decomposition, makes the following22-symbol bank exhaustive: S_0,S_1; C_0,...,C_3; and Q_W for all16 actual subsets W of Q.

Here Q's complete seven-coordinate table is displayed in mask order. F and G denote the two single-boundary cover counts, J their directed joint, and e,s the Hamilton endpoint/start counts.

| W mask | h | p | F | G | J | e | s |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 5 | 32 | 0 | 0 | 0 | 0 | 0 |
| 1 | 5 | 32 | 5 | 7 | 0 | 1 | 1 |
| 2 | 5 | 32 | 6 | 9 | 0 | 1 | 2 |
| 3 | 5 | 32 | 11 | 16 | 4 | 2 | 3 |
| 4 | 5 | 32 | 9 | 6 | 0 | 2 | 1 |
| 5 | 5 | 32 | 14 | 13 | 2 | 3 | 2 |
| 6 | 5 | 32 | 15 | 15 | 6 | 3 | 3 |
| 7 | 5 | 32 | 20 | 22 | 12 | 4 | 4 |
| 8 | 5 | 32 | 7 | 5 | 0 | 1 | 1 |
| 9 | 5 | 32 | 12 | 12 | 4 | 2 | 2 |
| 10 | 5 | 32 | 13 | 14 | 2 | 2 | 3 |
| 11 | 5 | 32 | 18 | 21 | 10 | 3 | 4 |
| 12 | 5 | 32 | 16 | 11 | 4 | 3 | 2 |
| 13 | 5 | 32 | 21 | 18 | 10 | 4 | 3 |
| 14 | 5 | 32 | 22 | 20 | 12 | 4 | 4 |
| 15 | 5 | 32 | 27 | 27 | 22 | 5 | 5 |

All entries are independently obtained by directed subset DP and explicit two-slot convolution. The existing S and C signatures remain unchanged.

## 3. Normalization of arbitrary shifts is still finite

Let a target constant be C=3c+d. Then C is nonzero modulo3. Remove its maximal initial zero-port component prefix P, leaving a suffix D whose first component is marked. The inherited assembly law gives

    c_(P→D)=h_P c_D,          F_(P→D)=p_P F_D.           (3)

An unmarked triangle in P would make C divisible by3, so none is possible. The prefix P therefore contains only S_0 and Q_empty. Let r be the number of Q_empty components and s the number of S_0 components. Then

    h_P=5^r,       p_P=2^(5r+s),
    c_D=C/5^r,     0<=r<=v_5(|C|).                    (4)

Their order within this zero-port prefix does not affect these coordinates. The target may have arbitrarily many S_0 blocks, but they alter only its dyadic F scale, already allowed by the clock shift.

Every **marked first component** in the22-symbol bank satisfies

    0<=v_2(F)<v_2(p),          v_2(F)<=4.               (5)

This follows immediately from the displayed table and the S/C tables. It persists when an arbitrary further component is appended: the update is

    F_new=F_old h_block+p_old F_block,

where h_block is odd and the second summand is divisible by a strictly higher power of2 than the first. Hence v_2(F_new)=v_2(F_old), and p_new has still larger valuation. Therefore every normalized D has v_2(F_D) in {0,1,2,3,4}, independent of its length and later markings.

Let F_0 be the odd part of3F_source. Combining (1), (3), and (4), every normalized target D must occur in one of the finite parameter pairs

    (F_D,c_D)=(2^a F_0, C/5^r),
    0<=a<=4,       0<=r<=v_5(|C|).                    (6)

Conversely each accepted D from (6), preceded by r unmarked Q components and any chosen number of S_0 blocks, gives a target satisfying (1) with the necessary integer clock shift. There are legal sufficiently large clocks. Thus no arbitrary-shift target is lost.

## 4. Complete search at fixed (F,c), without an order cutoff

For each pair (F,c) in (6), recursively append any component in the22-symbol bank, requiring the first to be marked. Evaluate its actual seven-coordinate signature by the inherited ordered-sum law.

- Discard a prefix if its current F exceeds the target F. The update F_new=F_old h_block+p_old F_block is nondecreasing because h_block>=1 and F_block>=0.
- Test a word whenever its F equals the target and its final component is marked; retain it only if h−J equals the target c.
- Recurse to a longer word only when the current unrestricted two-slot count p is at most the target F.

The last pruning rule is complete. If a word is extended to a longer word with a terminal marked component Z, then that final component has F_Z>=1, so the total F is at least p_before_Z F_Z>=p_current. Thus a prefix with p_current>F cannot have any successful longer completion. This remains valid if one or several unmarked blocks come before the eventual terminal marked block.

Every component has p>=2; consequently each recursed prefix has at most floor(log2 F) blocks. A tested terminal word may have one additional block. The procedure therefore terminates and covers **all word lengths**, not just a selected vertex-order window. It is elementary finite decision for the exact pair within this bank.

Sections3–4 are the all-order proof underlying the finite return certificate. Positive recurrences, not extrapolation from sampled orders, justify every omitted word.

## 5. Complete result and reproducibility

The self-contained `verify.py` embeds the seven source words, the full component bank, the finite normalization, and the complete recursion. It tests both digits from all seven states. The exact output has exactly the six edges already in the companion tree, each with one normalized marked word and r=0; every other edge has no candidate. Thus (2) follows by exhaustive successor closure starting at the root.

For each digit test, `receipt.json` records every normalization pair from (6), its visited-prefix count, terminal candidate count, and all accepted literal component words. This includes the otherwise easy-to-miss unmarked-Q prefix case for target C=−35: r=0 and r=1 are both tested. No high-order unmarked prefix or integer clock shift is silently discarded.

The script independently verifies the complete Q port table, every labelled strong four-vertex tournament's membership in the Q isomorphism class, and all22-bank words through total vertex order7 by literal graph subset DP. The normalized 2-adic bound is separately checked on every such word beginning with a marked component. These are controls on the implementations; Sections3–4 prove completeness beyond them.

Run `python verify.py` and `python -O verify.py`; output must match. The proof imports no code or receipts from the smaller packet, so the extension is reproducible separately. No software installation is needed.

**What remains:** a broader physical carrier with larger strong suffix components, or a genuinely different return mechanism. Even such a construction still needs recursive closure, source/entry coverage, and any required native-operation certificates before it can prove the unrestricted odd spectrum.

# An exact one-arc repair on a branching set, with full inverse

This follows the [rooted migration](../chanpak_bank_fusion_20261117/migration/PROOF.md) and [uniform-reader](../chanpak_passive_pubsub_20261117/reader/PROOF.md) constructions. The background pass found the three invariant-preserving blocks and the physical source. This follow-up removes the restriction on allowed words, proves an exact decoder for the entire language, and determines its counting exponent.

## Preserve only the advertised signed coordinate

Start from a transitive pair rooted at its last vertex. Write (h,R,S) for the Hamilton count, the marked terminal cover count with an optional second run, and the full Hamilton count starting at the root. The initial state is (1,2,0), with exactly one full Hamilton path ending at the root. A migration of width k with padded binary address u, leading bit zero, keeps S=0. On (h,R) it has matrix

    [[F(u), u+1], [2^k−rev_k(u), 2^k]],

where F(0)=1, F(2v+1)=F(v) and F(2v)=F(v)+v+[v odd] on positive inputs. Matrices act on column vectors. Here are three blocks preserving 3h−R=1:

| Letter | Width and address in chronological order | h on that line |
|---|---|---|
| A | (6,20) | 82h−21 |
| B | (4,2), (8,84) | 5976h−1579 |
| C | (6,8), (10,340) | 96264h−24407 |

All finite words are legal. Finish by (2,2), with binary address 10, giving (h_f,R_f)=(12h−3,15h−4). The final rooted graph I is strong: the new root reaches the spine, the spine dominates the old graph, and the unique old path ending at its root allows all old vertices to reach the return arc to the new root. In the uniform reader, use I's converse before the source, then append an exterior singleton which loses to all of I's converse except that it beats the distinguished vertex. The unique return edge gives the exact read R_f P₂+h_f H from the [single-return formula](../chanpak_passive_pubsub_20261117/reader/PROOF.md#the-sharp-cone).

The literal source masks (174,124,120,240,65,80,1,118) give (H,P₂)=(187,3736). Reversing just arc (5,6) gives masks (174,124,120,240,65,16,33,118) and (197,3728). Thus the same arc changes every output by exactly 10h_f−8R_f=2. The lower endpoint is N=58284h−15505. Literal DP verifies the initial pair 42779 and 42781. At chronological word CAB the pair is 2052302173293611 and 2052302173293613, certified by Lucas's criterion in the replay. This is one twin-prime witness, not an infinitude claim.

## Complete inverse: two moduli carry different information

A modulo-72 decoder alone identifies the last letter only under the no-AAA restriction. The restriction is unnecessary. Use the shifted coordinate x=27h−7. Then

    A: x → 82x,
    B: x → 5976x−808,
    C: x → 96264x+14852.

Initially x=20. Throughout the grammar x is divisible by 4; A increases v₂(x) by one, whereas B resets it to 3 and C resets it to 2. Separately, h mod 3 is unchanged by A, set to 2 by B, and set to 1 by C. The starting residue is 1, the same as after C. Therefore a target endpoint fixes h, and h fixes the length of its terminal A-run:

    if h≡2 mod3: number of trailing A's = v₂(27h−7)−3;
    if h≡1 mod3: number of trailing A's = v₂(27h−7)−2.

Undo that many A's by exact affine division. If h=1, stop. Otherwise the old mod-3 residue says whether to undo B or C. Reject any failed division, nonpositive state, invalid residue, or negative number of A's. Every accepted inverse can be run forward again, and every legal word is decoded. The maps strictly increase on h≥1, so the inverse terminates and there are no word collisions.

This illustrates a nonuniform way to split information: the 2-adic valuation carries the length of the suffix, the mod-3 residue carries its preceding letter. Neither channel has to preserve the other channel's data.

## Exact size of this bank, and the obstruction that survives

Let σ be the unique positive solution of

    82^(−σ) + 5976^(−σ) + 96264^(−σ) = 1.

It is approximately 0.1450940571. For any word, h lies between fixed positive multiples of the product of its letter multipliers. Indeed for maps m h−b, h divided by the product of multipliers equals 1 minus a sum bounded above by max(b/m) times the geometric series with ratio 1/82; that is strictly less than 1. The corresponding normalized count of endpoints is comparable as well. Assign each letter probability m^(−σ), and cut the full ternary prefix tree at the first product over X. Every leaf then has probability between (96264X)^(−σ) and X^(−σ), with the leaf probabilities summing to 1. A full ternary tree has a number of internal nodes proportional to its leaves. Injectivity of the inverse proves exactly Θ(X^σ) different lower endpoints below X. This improves the restricted grammar exponent but is still sparse.

In fact the fixed source and invariant impose a single progression: h always stays odd, so every lower endpoint is 42779 mod 116568. Removing the word restriction cannot fill other phases. A possible route forward must change the source profile or the signed line, not merely admit more words which preserve the same line.

## Why widening a single block does not produce new exact letters

Among one-step padded leading-zero migrations with width k≥1, preserving the whole line 3h−R=1 requires the second component of (3,−1)M to equal −1, hence 2^k=3u+4. So k≥2 is even and u=(2^k−4)/3. For these addresses rev_k(u)=u/2. Writing u_(k+2)=4u_k+4 and using that u_k is divisible by 4 gives the exact recurrence F(u_(k+2))−F(u_k)=5u_k/2+4. With F(u₂)=1 this yields F(u_k)=(5·2^k+6k−14)/18. Thus (3,−1)M=(k−3,−1), and only k=6 preserves the line. Wider individual letters shear it. The successful second blocks must do something genuinely different from widening a single letter.

Run `python replay.py` for the literal first reader pair and two rooted states, all three transfers, 9841 decoded words (including unrestricted repetitions), Lucas checks, and the arithmetic controls on single-step widths.

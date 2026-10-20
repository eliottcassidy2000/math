# One adjacent swap destroys the repetition factor

Date: 20 October 2026. Status: **PROVED elementary lemma**, independently audited. This statement applies to every repeated nonconstant binary word, not just to mechanical words.

Let v be a binary word of length q_0 with p_0 ones, 0<p_0<q_0. Let d>=2 and c=v^d, and obtain w from c by one cyclic adjacent swap of unequal bits. Define

    D_0=2^(q_0)-3^(p_0),
    D=2^(d q_0)-3^(d p_0)=H D_0,
    H=sum_(j=0)^(d-1) 2^[q_0(d-1-j)] 3^(p_0 j).

Then H>1 and gcd(H,6)=1. For the chronological carry B,

    gcd(H,B(w))=1,
    gcd(|D|,B(w))=gcd(|D_0|,B(w))<=|D_0|<|D|.          (1)

In particular w cannot encode an integer 3x+1 cycle of either sign.

## Proof

Choose a rotation cut not between the two edited positions. Then the swap is internal to the displayed word. The unedited rotation is still a d-fold repetition of a rotation v' of v. The exact repetition law gives

    B((v')^d)=H B(v').

An internal adjacent swap 10->01 changes the carry by +2^a3^b, and the reverse swap changes it by -2^a3^b, for nonnegative integers a,b determined by its position and suffix one-count. Thus

    B(w')=H B(v') +/- 2^a3^b.

The last term is a unit modulo H, giving gcd(H,B(w'))=1. Consequently

    gcd(|H D_0|,B(w'))=gcd(|D_0|,B(w')).

For a one-place rotation with first bit epsilon, the carry identity is

    2 B(rot(w))=3^epsilon B(w)+epsilon D.

Both 2 and 3 are units modulo D and all its divisors. Therefore gcd with H, with |D_0|, and with |D| is invariant under rotation, transferring both conclusions back to the original cut. Finally H>1 gives the strict inequality in (1). The proof does not assume D>0. ∎

## Mechanical consequence and sharp boundary

If d=gcd(p,q)>1, the zero-intercept mechanical word of shape (q,p) is exactly the d-fold repetition of the reduced-shape mechanical word with (q_0,p_0)=(q/d,p/d). Thus this lemma excludes every one-swap noncoprime mechanical shape immediately. Combining it with the separately audited primitive positive-clock exclusion removes the coprimality restriction from that positive mechanical-family result.

The small-base content need not be one. For v=001 and d=3, the edited word w=010001001 has D_0=5, H=97, D=485, and g=5. Thus the exact reduction can retain all of the base clock while destroying the entire repetition multiplier.

Adjacency and a single actual edit matter. Two swaps can undo one another and return to a repeated integer-cycle word. An arbitrary edit does not necessarily have a 2,3-unit carry difference. The conclusion is not asserted for an arbitrary word merely sharing the repeated word's shape.

## Exact checks

`verify_repeated_word_obstruction.py` checks every nonconstant binary base of length 2 through 8, every repetition count 2 through 5, and every unequal cyclic adjacent swap occurrence. It verifies both identities in (1), not merely absence of integer cycles. The receipt records 25,088 swaps and 1,016 crossing-cut cases; normal and optimized Python outputs agree.

Reproduce:

    python3 verify_repeated_word_obstruction.py > repetition-receipt.replay.json
    python3 -O verify_repeated_word_obstruction.py > repetition-receipt.optimized.replay.json
    cmp repetition-receipt.json repetition-receipt.replay.json
    cmp repetition-receipt.json repetition-receipt.optimized.replay.json

The inherited repetition and local-swap identities are in [local carry and prime-power gates](../collatz_local_prime_power_gates_20261020/LOCAL_CARRY_AND_PRIME_POWER_GATES.md). The present lemma is their direct joint application, with the crossing cut and full gcd retained. No historical-priority claim is made.

# Independent proof: the alphabet-three safe-prefix shield

2026-10-19. **PROVED.** Independent argument using maximum rigidity, with exhaustive finite-state continuation controls kept separate. This is an alternative to the producer's simultaneous extension induction.

## Conventions and theorem

For a word a=(a_0,...,a_(F-1)), its left apices are the first entries of all rows of its finite absolute-difference triangle, including a_0 itself. Call a good when every such entry is at most one.

**Shield theorem.** Let a be a good word over {0,1,2,3}, and suppose its terminal entry belongs to {1,2}. Appending any finite or infinite binary word preserves goodness at every depth.

**Maximum rigidity used.** If a triangle with seeds in [0,3] has apex three, then every seed in its cone is zero or three. Therefore a cone containing a seed one or two cannot attain three.

## Backward-zero lemma

Suppose a first bad left apex appears at depth N and has value two. Suppose also that every cell on the backward diagonal

    (N-j,j), 1<=j<=J,

has value at most two, where J<=N. Then at depth N-J the row begins with J zeros followed by a two.

Proof by induction in j. At j=0 the row at depth N starts with a two, which is the assertion with zero leading zeros. Suppose the row at depth N-j+1 starts with j-1 zeros followed by a two. The row above must start with j equal entries b, since its first j-1 consecutive differences are zero. Its first entry is a previous left apex, hence b is zero or one. Write c for its next entry. Then |b-c|=2, and by the diagonal hypothesis c<=2. The only possibility is b=0,c=2. This is the next induction step.

## Shield proof

Assume a binary continuation has a first bad left apex at time N. Goodness of the prefix gives N>=F. All values are between zero and three. The apex cone contains the terminal prefix marker a_(F-1) in {1,2}, so maximum rigidity excludes value three: the first bad apex must be two.

For every j with 1<=j<=F-1, the cone of the cell (N-j,j) consists of seeds from j through N and contains that marker. Hence this cell is not three. For j=F, its cone lies entirely in the binary tail, so its value is zero or one. Thus the backward-zero lemma applies with J=F and forces a two at (N-F,F). This contradicts binary evolution of the tail. No first bad apex exists.

The proof uses no probabilistic assumption. It is pointwise for every binary continuation.

## The stronger one-symbol extension rules

Let w be any good word over {0,1,2,3}.

A. The word w1 is always good.

B. The word w2 is good if and only if w is not all zero.

C. The word w0 is good whenever w is all zero or its rightmost nonzero entry is one or two.

For A and B, only the new final apex needs checking. Let c=1 or2 be the appended seed and let N be the new depth. The marker c lies in every backward-diagonal cone, so three is impossible throughout that diagonal. If the new apex were two, the backward-zero lemma with J=N would force w to be all zero and c=2. Conversely a zero word followed by two does have apex two. This proves both directions of B and all of A. For C, the all-zero case is immediate. Otherwise truncate at the rightmost nonzero marker, which is one or two; that truncated prefix is good. Its remaining zeros, including the newly appended zero, are a binary continuation covered by the shield.

This proves the producer's A/B/C rules without any simultaneous induction or exceptional constant-word analysis.

## Sharp alphabet boundary and scope

The cap three is essential. The word (1,2,3,4) is good, but appending one gives left apices

    (1,1,0,0,2).

Thus the append-one rule already fails at alphabet maximum four. The backward-diagonal proof fails exactly because a cone containing one can now attain three, allowing the alternative |1-3|=2.

The shield does not say that a safe finite prefix with arbitrary terminal value is safe under all binary continuations. A terminal one or two, or an equivalent verified extension condition, is needed. No assertion about deterministic prime gaps follows without the required bounded alphabet and safe-prefix hypotheses.

## Independent exhaustive check

Run:

    python3 verify_prefix_shield.py > prefix_shield_receipt.json
    python3 -O verify_prefix_shield.py > prefix_shield_receipt_optimized.json

For every width from one through eight, the checker builds the complete state graph on {0,1,2,3}^F with both binary boundary transitions. Reverse reachability from states with first coordinate at least two identifies exactly the states that can ever fail under any finite binary input word. Every good prefix ending in one or two is checked against this complete set. This checks all future continuation lengths, not merely a fixed finite suffix length, inside those bounded widths. The checker also tests A/B/C on every good prefix in the same universe and includes the alphabet-four hostile witness. Universal width validity comes from the proof above.

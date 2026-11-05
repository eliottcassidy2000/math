# Independent sharper insertion / compression audit

2026-11-05. **PASS for the stronger inequalities, exact profile compression,
sharper factorwise order bounds, and odd prime unit amplitudes through 29.**

Frozen proof hash: `0d594bd035b99dd3c7b1dfc29cd8688da337c56d983cb72fa4c8a962cef5e826`.
Verifier hash: `b3f3b6ea5248da7c1a9cc955ea19d46cc88f92bd0b0934ff635a69e1e86c6489`.
Receipt hash: `d002434ca56c9d9a8e49e54391e900d00717bf9e2afe0b703ed3f9dedecded57`.

## Proof audit

For a root with an outneighbor, every deletion Hamiltonian path has a
first such outneighbor. Inserting the root immediately before it is legal
and never puts the root last. Deleting the root recovers the input path,
so `h-e>=d` and `U>=d+f`. The outneighbor/positive-U hypothesis is essential;
universal sinks are correctly excluded. Combining this inequality with the
previously audited unit determinant dichotomy proves `2d<=U+beta` for
general two-piece unit synthesis, not merely singleton partners.

Before the first deletion SCC meeting the root's outneighbor set, every
component beats both the root and every later component. These components
are an actual ordered prefix of the full tournament. Removing its singleton
components preserves h,e,R,d separately: each singleton contributes one
to the unique permissible path slot. The retained graph is induced and
no arc is changed or component moved.

The SCC starting-cover formula gives `U>=2^s` for the total number s of
vertices after the first occupied component. This uses
`P_C(2)>=2^|C|`, including both empty slots via H(empty)=1. The prefix's
retained nontrivial vertices plus any nontrivial first occupied component
number at most d. A singleton first occupied component adds at most one,
the suffix adds at most floor(log2 U), and the root adds one. Thus the
displayed bound `|Q|<=d+floor(log2 U)+2` is valid for every positive U,
without bounded slack.

The sharper source bound follows immediately. On the partner converse,
positive v=L-t gives `|Y|<=t+floor(log2 v)+2<=L+1`. When v=0, the unit
determinant forces L=g=t=1 and a singleton replacement is valid. All
quantities needed by the physical star identity remain unchanged.

The complete prime-amplitude statement through 29 is correctly an
existence classification. Every absent pair uses only known deletion rows
through 13. The one unclassified candidate d=15 occurs only for the already
attained pair (29,+1), so its possible additional representations are
irrelevant to whether that pair exists. The explicit source has profile
(23,4,40,11); its positive K multiplicities, not the source at m=0, yield
59,117,233,465,929. No formal zero-multiplicity interpretation is used.

The shorter H15 unit-return obstruction is also correct: at d=3,U=8,
the new bound forces f<=5, hence only f=4 remains. This supplement does
not itself supply the separate complete H15 one-step classification.

## Independent verification

The audit uses the independently complete tournament census and checks
every root of every representative through order six, with exact
automorphism orbit weights. The 400 positive-U representative roots cover
all 195,508 positive-U labelled roots. All insertion and profile-compression
claims hold. The remaining 6,505 zero-U labelled roots are excluded as
required. Literal insertion is additionally verified on 194 deletion paths
in the complete smaller representative universe.

There are 678 larger random root controls and twelve adversarial
multi-component, singleton-chain, and shuffled-label compression checks.
All 42 small unit-profile pairs satisfy the sharper actual-piece bounds,
with 126 independent physical star counts. Prime-budget checks explicitly
record the only unneeded unknown row, d=15 at (29,+1). All five actual
p29 K extensions agree.

Normal and optimized source runs are byte-identical and equal to the frozen
receipt. No mathematical or implementation issue was found. Results are
in `AUDIT.json`; reproduce with `python audit.py`, using the sibling
independent primary, factor, and tournament-census audit helpers.
Author source files were not edited.

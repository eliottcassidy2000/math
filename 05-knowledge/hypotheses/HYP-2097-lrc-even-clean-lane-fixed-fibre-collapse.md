# HYP-2097: On the even clean lane, the n=14 target collapses from 190 merged round nodes to 64 fixed fibres

**Session:** codex-2026-06-03-S577  
**Status:** SUPPORTED by exact funnel bookkeeping; per-class certificates open.

## Claim

HYP-2092's clean-composite n=14 route should not treat the whole 190-node converse-merged round seam as the primary proof target.  After merging Opus S568/HYP-2093 and Oracle S576o/HYP-2094, the target sharpens:

1. Label the 64 self-converse fixed round classes first.
2. Use the remaining 126 non-fixed converse-paired nodes as controls for the open/generic round body.
3. Attach D/U/N private-pivot labels and THM-397 pair-sum endpoint-owner labels only after this fixed-fibre split.

In short: the n=14 proof table should be a 64-row fixed-round certificate table, not a 190-row seam labelling problem with equal priority.

## Exact S577 Funnel

`04-computation/lrc_even_clean_lane_fibre_collapse_s577.py` merges:

- S576/HYP-2092: HYP-2091 gives clean/wall and unit/nonunit proof lanes.
- S568/HYP-2093: even-n measure-zero rows are expected to be floor-tight.
- S576o/HYP-2094: even-ladder worry is contained in self-converse round classes.
- S569/HYP-2095: the cheap certificate is an unblocked small pair; blocked small pairs split into paired shields or endpoint anchors.
- HYP-2096: n=14 has a nine-runner unit spine plus four composite slack runners after the unit-shell gate is forced.
- S574/HYP-2090 and THM-397: residual certificates need D/U/N private pivots and endpoint owners.

For n=14:

- all runner tournaments `A000568(13) = 48,542,114,686,912`
- open round body `A000016(13) = 316`
- converse-merged round seam `= 190`
- self-converse fixed worry nodes `= 64`
- non-fixed converse-paired controls `= 126`
- D/U/N labels per speed set: `D=12`, `U=9`, `N=13`
- pair-sum endpoint cells: `C(13,2)=78`
- nonunit summand shells mod `27`: `4`
- fixed/all compression: `64 / 48,542,114,686,912 ~= 10^-11.88`
- fixed/merged compression: `64 / 190 ~= 0.3368`

Thus the open body is already an extreme compression, but HYP-2094 still removes two thirds of the converse-merged seam from the primary worry target.

## Merged Lemma Queue

**E1 containment.** A counterexample's optimal class on the even clean lane is self-converse round.  This is the HYP-2086/HYP-1998 open-round direction sharpened by HYP-2094.

**E2 certificate table.** Each of the 64 n=14 fixed round classes gets a concrete witness family: n-clock, pair-sum pinch, antipodal/unit-shell clock, or endpoint-owner certificate.

**E3 realization independence.** The class certificate survives all speed realizations of that fixed round class, not just bounded representatives.

**E4 nonunit descent.** The gcd-3/gcd-9 shell defects mod `27` either lift to a second clock or force a THM-397 endpoint owner.

**E5 owner compatibility.** HYP-2090 D/U/N private pivots and THM-397 endpoint owners must agree on the same labelled speed or pair-cell obstruction.

## Why The 126 Paired Nodes Still Matter

The non-fixed converse pairs are not irrelevant.  They are the controls that explain why the generic/open round body is already lonely:

- time reversal sends a round class to its converse;
- if a class is not fixed by converse, the binary dihedral/cosine obstruction is already off the fixed boundary;
- labels on these 126 controls can show which D/U/N or endpoint features disappear before the fixed seam is reached.

But they should be used as comparison cases after the 64 fixed classes are isolated.

## Next Computation

Build `04-computation/lrc_n14_fixed_round_certificate_s578.py` with one row per fixed self-converse round class.  Suggested columns:

- canonical self-converse round node
- gap necklace / round d-vector
- anti-automorphism witness
- unblocked small-reduced-sum pair candidate from HYP-2095
- n-clock candidate
- best pair-sum candidate and value
- unit-spine normalization status from HYP-2096
- D/U/N owner multiplicities
- private obligations
- THM-397 endpoint-owner pair cells
- unit/nonunit shell state mod `27`
- certificate status: floor-tight, pinch-certified, endpoint-certified, open

The script should also include a control table for the 126 non-fixed converse-paired nodes, but only after the 64-row target is explicit.

## Tournament Analysis

For S577 the tournament vertices are even LRC rungs, not runners.  The observable is:

`(self-converse fixed node count, nonunit shell count, D/U/N total, converse-merged seam size)`.

The switch orients toward larger residual proof burden.  The fingerprint is transitive, so S577 is a routing table, not the final cyclic object.  The next nontrivial tournament should use the 64 fixed n=14 nodes as vertices, with switches comparing certificate difficulty, owner-label conflict, or nonunit-shell defect.

## Assumption Challenge

Do not assume the 190 converse-merged nodes are the natural vertices for the next proof object.  That was the right S576 target before HYP-2094 landed.  After the even-ladder self-converse collapse, the LRC predicate preserved by the quotient is sharper on the 64 fixed nodes.  The 190-node view preserves the open round body and time-reversal quotient, but it destroys the fixed-boundary distinction that now carries the putative counterexample.

## Honest Status

This is a reduction and proof-workflow claim, not a proof of n=14.  The open gaps are exactly the ones named by HYP-2094:

- bounded-speed evidence must become all-realization evidence;
- class-level containment must be proved, not only observed;
- every one of the 64 fixed classes needs a uniform witness or descent certificate.


# Independent audit: mobile-triangle physical scalar fibers

Verdict: **PASS** on the frozen producer `inputs/PROOF.md`, sections 1–5. No mathematical correction required. The all-order fiber theorem is proved by the exact combinatorial transfer and its checked annihilating polynomial, not by extrapolation of a finite census. The computational evidence below is FINITE-EXACT. No Chan–Pak coverage conclusion follows.

## Frozen scope and independence

Source checkout supplied: `session228-work`, commit `1a178a44f13f9602a3be9e84714be91cbfc7d7a6`. Audited producer files are copied into `inputs/`; `SHA256SUMS` identifies the exact bytes. The independent `verify.py` imports no producer code and uses only Python's standard library. It separately builds adjacency matrices, computes subset-DP Hamilton paths and rooted start paths, and enumerates local path assignments and permutations to derive the six-state matrices. No assertions are used; failure checks remain active under `-O`.

Read repository AGENTS, RESEARCH-PROTOCOL, ACTIVE-GUARDRAILS, targeted Chan–Pak correction entries MISTAKE-715–717, and Session227 adaptive-arithmetic PROOF. This audit does not edit the repository or producer. Final integration/checkpointing belongs to the parent session.

## Mathematical reconstruction

For word b of length m, replace adjacent bits p,p+1 by C3, direct earlier deletion blocks into later blocks, and let the root beat k=2b_p+b_(p+1) vertices of C3. Singleton neighborhoods encode remaining bits. There are m+2 vertices. The C3 rotation action is transitive on subsets of each fixed size 0,1,2,3, so no hidden neighborhood-type choice remains. Deletion has exactly three Hamilton paths.

An independent concise derivation of U tracks x=H(processed deletion graph) and U for the root plus that graph. Initially (U,x)=(1,1). Processing a singleton b yields

    (U,x) -> (2U+(b−1)x, x).

Processing a triangle with k permitted starting vertices yields

    (U,x) -> (12U+(3k−9)x, 3x).

To see the second formula, the two-path ordered cover count of C3 is 12, its Hamilton count is 3, and each allowed first vertex contributes three nonempty first-path choices. The formula equals three times the effect of processing the two replaced singleton bits. Consequently U=3(1+binary_value(b)) for every word and position. The same argument proves the stated multiple-disjoint-replacement extension with factor 3^t; it says nothing about distinct scalar outputs.

For E, a Hamilton path ending at the root first traverses the entire deletion graph; it must exit its last block. If that block is an unreplaced terminal singleton with bit1, its edge points away from the root, giving E=0. This argument and the hypothesis are sound.

For D, deleting the root from a Hamilton path gives two ordered paths P,Q, with P empty or ending at an inneighbor and Q empty or starting at an outneighbor. Within an ordered sum, each path traverses each block at most once, so local path pieces concatenate uniquely. The producer's six-state automaton retains exactly the necessary endpoint data. Independent enumeration agrees with all four displayed S0,S1,T1,T2 matrices. A=S1*S0 satisfies

    A^4−6A^3+9A^2−4A=0.

The even and odd placement expressions have the stated ordering of factors. Their four-by-four grids exactly match the independently generated ones in normal.json. Both proposed formulas satisfy the induced recurrence in each variable; the boundary delta at j=0 is killed by the shift factor. Induction first in j then k proves equality for every j,k≥0. No invertibility is needed and the zero eigenvalue does not invalidate forward induction.

Thus for r≥2 and p=0,...,2r−3 the physical state is

    (B_r,0,4^r+2,3) for p=0;
    (B_r+(-2)^p−2,0,4^r+2,3) otherwise,
    B_r=(4^r+6r−1)/3.

The scalar values are pairwise distinct: distinct nonzero positions have distinct signed powers, and equality with p0 would require (-2)^p=2, impossible for integer p≥1. This proves at least 2r−2 values in each stated fiber, not exactly 2r−2 in the unrestricted fiber. The order is 2r+2 and deletion count is fixed at 3; the common U varies with r. It is therefore unbounded fiber cardinalities across physical sidecars, not a claim that a single fixed numerical sidecar has infinitely many values. The minimum r2 witnesses are D=9,5 with sidecar (0,18,3).

## Computation and hostile boundaries

Commands, run from this audit directory:

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

Both runs PASS and their outputs are byte-identical.

- All binary words of lengths 2–7 and all adjacent replacement positions: 1,284 literal tournament profiles. D independently matches the six-state count, U matches the binary identity, and d=3.
- 516 retained-terminal-singleton bit1 cases all have E=0.
- 1,892 alternating-family transfer checks, r=2,...,44, verify every formula and distinctness claim at every allowed placement.
- Exact matrix annihilator and all 32 initial-grid values independently checked.
- Terminal replacement control: word01 replaced entirely by C3 gives (D,E,U,d)=(5,2,6,3), although the word ends in1. Thus removing “unreplaced singleton” from the E statement would be false. For length2 words 00,01,10,11 the profiles are respectively (3,3,3,3), (5,2,6,3), (5,1,9,3), (3,0,12,3). The U identity correctly survives all these boundaries.

## Novelty and scope audit

The closest inherited proof is `05-knowledge/results/session227_chan_pak_constructions_20261108/adaptive_arithmetic/PROOF.md`: it supplies two-point collisions and proves a two-value bound for a specific four-word menu, even with arbitrary physical seeds and positive-depth iteration. Its common-graft sensitivity result is also inherited, not overturned. Mobile placement is outside that menu and is not a common suffix applied to two existing states. The new theorem therefore escapes a narrow construction barrier without contradicting it.

Targeted searches for physical-sidecar, unbounded fiber, and mobile-triangle statements in the session-results tree and nearby tournament canon found no prior identical construction. This supports a repository-local distinction, not an external literature-priority claim or an exhaustive corpus novelty claim. THM-4904's root-change **word multiplicity** fibers are a different quotient; many encodings of the same D do not establish many D values at fixed (E,U,d).

All graphs here have a global terminal sink, so no strong-tournament theorem is proved. Their distinct D values are geometrically separated and do not furnish a filled interval, an effective smaller predecessor for a prescribed odd target, or a recursively closed digit grammar. The result is a rigorous new physical-family mechanism and a many-choice fiber theorem. Chan–Pak completeness remains open.

# Arbitrary response does not remove the absolute-target obstruction

**PROVED exact target exclusion for the two fixed sources, with physical positive controls elsewhere; PROVED all-order source-response criterion.** This is a boundary theorem for the proposed decoder, not a counterexample to Chan–Pak.

## Every response and mask still lies inside two fixed languages

Let F_p(W)=h+tW-dK(W)-e(W mod2). Regardless of the coefficient A, the mask, or the amount of shared continuation, a paired output at target N requires F_L(WL)=F_R(WR)=N. Therefore it belongs to the intersection of the two ordinary rooted languages. Choosing A adaptively does not enlarge those languages.

There is a finite COMPLETE test for each fixed source with t>d. Remove trailing zero letters, which do not change H. If a nonempty word of length k ends in1, its final weighted start is at most2N-d, whereas the binary recurrence gives at least2^k(t-d)+d. Thus

    2^k <= (2N-2d)/(t-d).

Test the seed and all addresses of length at most the resulting integer bound. This excludes the entire infinite rooted language, not just words up to a chosen experimental length.

For L=(15,3,19,5) and R=(15,3,21,5), this test rejects BOTH 32,561 and32,563 from BOTH languages. Consequently no coordinated edit or response coefficient between these sources can place either absolute center.

The targets are physically attainable elsewhere. Start with outgoing masks (126,108,56,112,98,64,4), root3, and its version reversing arc2→5. Append the chronological word1000101101. Independent literal subset path counting gives respectively32,561 and32,563 on17 vertices. These inherited neutral-switch constructions are credited to `ans_inverse/ANS_INVERSE.md`; the contribution here is their use as explicit hostile controls for the new proposed selector. They show precisely why additional source supply is necessary.

## A response-orbit criterion for selecting new sources

Take any two physical sources sharing an odd deletion count d>=3, and let D=tR-tL with0<|D|<d. An edit of width n has distant-bit coefficient

    A=2^n D+dS,

where S is an odd integer with |S|<=2^n-1, namely the reversed signed edit address. Conversely every such S has a legal signed-bit realization, hence arises from physical binary words on the given sources (this alone does not ensure a positive collision tail).

**Theorem.** A unit response A=+1 or-1 occurs at arbitrarily large widths if and only if

    D mod d belongs to ±<2> in the unit group modulo d.

Equivalently, 2^nD=±1 modulo d for some n. Necessity is immediate. For sufficiency take infinitely many n in the appropriate residue class modulo ord_d(2), and put S=(±1-2^nD)/d. It is an odd integer. Since |D|<d, the bound |S|<=2^n-1 holds for all sufficiently large n. A signed binary representation provides the required edit.

More generally, for any FIXED odd A, under the same strict interior assumption |D|<d, arbitrarily large physical masks with response A exist exactly when A modulo d lies in the doubling orbit of D. This does not require A or D to be units. If |D|>=d, no unit response exists because |A|>=2^n(|D|-d)+d>=d. The boundary |D|=d can support other fixed responses and is not included in the general interior theorem.

This gives a finite arithmetic classifier for which new source pairs can supply each response. In particular gcd(D,d)=1 is necessary but not sufficient for unit gain: e.g. d=17 has ±<2>={1,2,4,8,9,13,15,16}, excluding D=3. Neither a favorable response orbit nor local residue coverage forces the absolute equation F_p(W)=N; the preceding physical hostile controls demonstrate that distinction.

Replay: `python3 finite_edit_recontact/absolute_target_boundary_replay.py`. It exhausts the derived finite language bounds, independently recounts the17-vertex positive controls, and checks response-orbit equivalence for every odd d from3 through39 and0<|D|<d.

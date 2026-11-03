# FRACTRAN priority scanner on a marked tournament path

This is a construction and finite control based on the exact S159 transfer/FRACTRAN note at commit `7589463e54e20585aff6a2e68019ee4a4af582cd`. It is not a Hamilton-path counting identity, a universality claim, or a simulation by the unmarked metagraph. It combines the ordered FRACTRAN program from that note with the separate fixed-path chord cube. The status tournament below is a new 15-vertex observation gadget, not the growing tournament family whose Hamilton-path count is H_m.

## The priority theorem

Let a finite ordered FRACTRAN program consist of reduced positive fractions f_i=a_i/b_i, i=1,...,m. For an input integer N>0 put e_i(N)=1 if b_i divides N and 0 otherwise. On any fixed directed Hamiltonian path with m designated ordered stations, place unit tiles with boundary colors U (undecided) and S (selected):

* if e_i=0: U→U, weight 1;
* if e_i=1: U→S, weight f_i;
* for either input bit: S→S, weight 1.

Impose the left boundary U and permit either right boundary. There is exactly one tiling: it ends in U iff N halts, and otherwise the unique U→S station is the first applicable fraction. Its weight is f_i, so multiplying N by the weight gives the next FRACTRAN integer. Later enabled stations cannot select another rule because S never returns to U. This is a color-matching tiling of the *marked directed path*, not an assertion that tournament arcs alone impose the guards.

Equivalently, with matrices acting on row vectors indexed U,S,

 M_i(N) = [[1-e_i(N), e_i(N) f_i], [0,1]].

Then [1,0] (∏_i M_i) has first coordinate 1 precisely at a halt, and second coordinate the selected fraction otherwise. For blocks A then B, the two-state product M_A M_B gives an associative compositional scanner: B acts only if A did not select. This explicitly retains ordinary first-applicable priority, rather than treating fractions as an unordered set.

For the example (3/2, 5/3, 1/5), N=6 enables the first two fractions and the scanner returns 3/2, hence 9. Swapping the first two stations returns 5/3, hence 10. The same underlying path can support either program; its marked station order and attached labels are indispensable.

## Actual S159 machine as an observation gadget

Take the S159 ordered list

 (17/143, 133/85, 17/19, 23/17, 1015/69, 23/29, 31/23,
  38073/62, 31/37, 41/31, 86/287, 41/43, 13/41).

Use the marked path 0→1→⋯→14. It has 91 nonconsecutive chords. Assign e_i to chord (i−1,i+1), i=1,...,13, interpreted as forward iff e_i=1. Assign every remaining chord the value e_1. Thus every e∈{0,1}^13 gives a labeled tournament C(e) with that same forward Hamiltonian path, and C is injective. The stations 1,...,13 scan those 13 distinguished chord bits, in order. Their fractions are decorations; encoding alone does not make the graph choose a rule.

The native partner reverses all 91 chords and fixes the path, so K_P(C(e))=C(1−e). Scanning it chooses the first fraction *inapplicable to the original N*, if there is one. This native complement does not reverse program order and is not the last transition or a previous computation step. It may not correspond to any input integer at all. Its line colour is a separate predicate: under THM-796 it is blue exactly when the chord tiling is fixed by path reflection, and black otherwise.

An exact criterion settles this: for any bit pattern e, put L=lcm{b_i:e_i=1}, with lcm(empty)=1. Then e is realized by a positive integer iff no b_j with e_j=0 divides L. Necessity holds because L divides any realizing N; sufficiency holds by choosing N=L. This criterion also handles repeated/shared prime factors, without assuming independence of divisibility tests.

The 13 actual denominators are

 143,85,19,17,69,29,23,62,37,31,287,43,41.

Their only constraints are four independent Boolean implications

 e_2≤e_4, e_5≤e_7, e_8≤e_10, e_11≤e_13,

with five other unconstrained bits. Each implication pair admits 00,01,11, hence the marked cube contains 3^4·2^5=2,592 arithmetic images out of 2^13 patterns. A pattern and its bitwise complement are both realizable iff each implication pair is 00 or 11, giving 2^4·2^5=512 such patterns. For example N=41 realizes 0000000000001. Its native complement enables denominator 287=7·41 but disables denominator 41, and so is unrealizable. This example is black: reflection exchanges chords (0,2) and (12,14), whose bits are respectively0 and1.

The duplicated unused chords are important: choosing their orientations as constants would generally stop the native mate from being literally in the encoded subcube. These duplicates carry no additional arithmetic information.

## State and quotient warnings

Even the full 13-bit decorated tournament does not determine the *next* status. For the actual program N=23 and N=46 both have pattern 0000001000000 and both choose 31/23, but their successors 31 and 62 have patterns 0000000001000 and 0000000101000 respectively. Thus a projected deterministic transition on the observation tournament is not well defined. Retain the integer, or exact relevant prime-exponent registers, for a genuine dynamical state. The observation can decorate a state transition system, not replace that system.

Forgetting the path, station labels, and edge markings loses still more. In the order-four fixed-path cube take chord order (0,2),(0,3),(1,3) with 1=forward. The patterns (0,0,0) and (0,0,1) define isomorphic unmarked tournaments via (0,1,2,3)↦(3,0,1,2), but scanning the chords sees no enabled bit in the former and the third in the latter. This finite isomorphism example is for an arbitrary three-guard system, not claimed to occur inside the S159 arithmetic image. A marked path has a unique order-preserving vertex identification, whereas an unmarked isomorphism can exchange which path is viewed as the reference; black moves or quotienting are not automatically semantics-preserving.

## Checks

`python3 check.py` independently enumerates all 2^13 bit patterns using the lcm criterion, checks the native-complement image counts, compares the two-state scan with explicit FRACTRAN execution for N=1,...,10,000, and finds the unmarked four-vertex collision via permutation search. `python3 collisions.py` searches for the actual 23/46 status collision. These checks corroborate the finite claims; the scanner and feasibility proofs above hold for all inputs. The colour terminology was repaired during the [coverage-machine integration audit](../../coverage_machine_induction_20261103/INCOMING.md); the arithmetic counts and scanner construction are unchanged.

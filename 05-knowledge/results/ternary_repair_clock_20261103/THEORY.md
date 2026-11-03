# Repairing the two offsets: exact clocks, residue clocks and integer carries

Status: **PROVED scoped mechanisms; FINITE-EXACT source certificates.**
The full odd-spectrum conjecture and an all-size closed source bank are OPEN.
Words use source-component order, as in the inherited golden seam packet.
The actual path complement keeps the supplied spine and reverses every other
arc. Counts and ports refer to the resulting literal tournament.

## R1. Three q2 probes force the old exact affine states

Retain v=(h,x,y), tau=(2,-1,0), and P=P_q2 from
[the inherited transfer](../golden_seam_recurrence_20261103/THEORY.md).
The observation matrix, with columns e_h,P e_h,P^2 e_h, is

    O=[[1,5,73],[0,8,144],[0,6,120]],  det O=96.        (R1)

For any two full states v,u and j=1,2, the following are equivalent:

1. h(q2^k U)=3h(q2^k V)+2j for k=0,1,2;
2. u=3v+j tau;
3. the same count identity holds for every k>=0.

The first three observations determine u-3v uniquely, and tau P=tau
supplies the solution. This is the exact meaning of commuting a replacement
with the *same* q2 insertion clock. It does not cover changing the clock
origin or merely preserving its reduction modulo 7.

The complete source-word saturations give:

| Maximum supplied strong-atom order / alphabet | H95 states | H97 states | States supporting both exact offsets |
|---|---:|---:|---:|
|Inherited selected 65|11|8|0 in either fibre|
|All charts through 6:960|31|30|0 in either fibre|
|All charts through 7:31733|74|84|0 in either fibre|

These are all-word-length results at the stated H values. Thus changing
the 95/97 representative anywhere inside these alphabets cannot restore
both offsets under the exact q2-equivariance requirement. In particular
the original four failed states

    (287,122,225), (289,121,225),
    (293,119,225), (295,118,225)                       (R2)

remain absent after adding every strong seven-vertex chart. This is an
operation-relative obstruction, not a global absence theorem for their H.

## R2. Complete compilation, and a further eight-vertex boundary test

There are 30773 strong seven-vertex charts with the fixed legal spine.
They extend the inherited 960 charts to 31733, with 10553 distinct counting
matrices. Repeated matrices still have distinct labelled source charts.

For any supplied atom c, the incoming
[universal three-coordinate theorem G1--G4](../golden_hamilton_transfer_20261103/RECURRENCE.md)
proves the existence of its matrix. Three literal calibrations determine it:

    row0=v(c), row2=v(c a)-v(c), row1=v(c aa)-v(c a).  (R3)

The primary compiler computes Hamilton/end/two-path-cover counts by subset
DP. The independent audit expands literal permutations and ordered cuts
over the entire orientation cube. It verifies every entry of all 30773
matrices, not a sample. Strongness is independently checked by absence of
an all-forward interval cut. The recurrence

    A_n=2^((n-1)(n-2)/2),
    S_n=A_n-sum_(1<=k<n) S_k A_(n-k)                 (R4)

follows from the unique ordered strong-component decomposition of a chart.
It gives S_7=30773 and S_8=2032504, independently confirming the universes.

Every supplied matrix is entrywise at least P_a. Consequently a word of
l letters has H>=h(a^l), which grows tribonaccially without bound. At
H<=2671 there are at most 14 letters. A generated set is certified exact
by an actual predecessor for each state and closure under every matrix
transition that remains below the bound. No vertex-order cutoff or
old C16 cone is imposed. Indeed the new strong source chart t0 has
v(t0)=(1,32,32), outside C16; using that old bank-specific guard would be
an incorrect pruning rule.

There is one further exhaustive boundary test. Remove any sequence of
leading 65 letters arithmetically from each state in (R2), retaining integral
nonnegative physical defects but imposing no C16 bound. The complete
reverse fringe has 14 states, listed in the receipt. None is the initial
state of any of the 2032504 strong eight-vertex charts. Both subset DP and
an independent full-cube permutation/cut audit agree. Therefore no word
of the form w C8, with w in the 65 alphabet and C8 a single strong
eight-vertex source chart, realizes(R2). This does not exclude an internal
eight-vertex block, several such blocks, or arbitrary larger atoms.

A useful all-order necessary guard comes from the inherited canonical
insertion mechanism. If p is not a sink, the number E of paths ending
at p satisfies 2E<H. Given such a path, remove p and reinsert it immediately
before the first vertex it beats. This vertex exists, is not the old last
vertex, and the result is a Hamilton path not ending at p. Deletion and
appending p invert this injection. Hence H>=2E, and the known oddness of
tournament H makes the inequality strict. In a nontrivial supplied-spine
chart, port 0 is not a sink. Thus the fringe state(43,19,41), with E=22,
is impossible at every size, independently of the atom bank. Other fringe
states are not excluded by this argument.

## R3. Actual repairs preserving the residue clock

Let t_m denote the strong seven-vertex chart with nonadjacent lexicographic
source bits m; bit 1 means smaller vertex points to larger. Keep all 65 old
symbols and add only t18413 and t12083. The following four source replacements
are literal:

| Input word / v | Output word / v | j | kappa in u=3v+j tau+7kappa |
|---|---|---:|---|
|a a q1 a a / (95,41,75)|t18413 a / (287,206,239)|1|(0,12,2)|
|a a q1 a a / (95,41,75)|r21 a a b / (289,191,239)|2|(0,10,2)|
|a a a a a q0 / (97,40,75)|r30 r2 / (293,210,281)|1|(0,13,8)|
|a a a a a q0 / (97,40,75)|t12083 a / (295,216,260)|2|(0,14,5)|

All have H'=3H+2j exactly, and all are strong: the retained spine and the
reversed last-to-first cross-component arc form a Hamilton cycle. The
integer carries are nonzero; the failed exact states have not been declared
realizable. Two of these residue repairs already use the old bank, while
the 287 and 295 repairs require an atom of order at least 7 in this source
language, by the complete 960-bank saturation.

This 67-letter bank is not closed. Its next residue-preserving targets fail
at 863,871,887,889, with the required (x,y) residues given in the receipt.
These failures cover all word lengths in that bank. Supplying two useful
atoms repairs one stage, not the induction.

Using the full 31733 alphabet instead gives a **complete four-level finite
construction from H31**:2+4+8+16=30 offset branches, each with its correct
parent state, literal source word and integer carry. Endpoints use at most
12 vertices. The receipt certifies every full state by independent subset
path and cover counts. This is finite depth; no fifth or arbitrary level
is claimed, and the separate 3H branch is not supplied by this tree.

## R4. The finite repairs are legal native programs

Apply the incoming
[one-temporary-vertex compiler, METAGRAPH Section 5](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md).
After padding the smaller endpoint with named sinks, flip each differing
pair using a prescribed birth, two legal path complements and named
deletion; finally delete padding belonging only to the target extension.
All names, incident cuts and both path words are retained. A subset search
on units, with the required arc contracted into a two-ended unit, finds
the traced path. For an untraced arc, use the inherited inverse gadget.

The four 67-bank repairs above use respectively 10,38,16,24 native moves,
with peak orders 9,11,11,10. The full four-level tree is also compiled and
replayed:1112 native moves in total over its 30 separate programs, with
maximum peak order 13. These are supplied certificates, not minimum-length
claims. Birth and deletion clocks are reported separately. Intermediate
H and strongness are not preserved, and these are not paths in a
fixed-order complement-only metagraph. The endpoint maps are exact;
no equality-of-H class reset is assumed.

## R5. An explicit residue/carry decomposition retaining all Q4 phases

Reduce v modulo 7 and set

    Y=y, X2=2x+h-y, g=h-y,
    alpha=Y theta+4X2 in F49,
    theta^2=4theta+2.                                (R5)

This is an invertible F7-linear encoding of v mod 7 as (alpha,g).
q2 multiplies alpha by theta and fixes g. The inherited seed has
alpha0=theta^38, and Frobenius on its nonzero clock is k->7k+36 mod 48.

Every repair in R3 and every edge of the finite tree satisfies

    alpha'=3alpha,  g'=3g+2j mod 7.                    (R6)

Since 3=theta^40, a branch shifts the field clock by 40. That shift commutes
with Frobenius. All branches at depth d have the same alpha=3^d alpha0;
the charge set fills F7 by depth 3. The phase-aware residue dynamics thus
has six field phases and seven charge states, or 42 states after saturation.
This is a finite control clock, not an all-size realization theorem. If
only H mod 7 is observed under the two digit operations, seven residues
suffice; retaining the field phase is an additional requested output.

For an integer representative rho in {0,...,6}^3 write v=7z+rho. A source
matrix P acts by the exact Euclidean carry rule

    rho'=(rho P) mod 7,
    z'=z P+floor((rho P)/7).                          (R7)

An actual ternary repair has

    rho'=(3rho+j tau) mod 7,
    z'=3z+kappa+floor((3rho+j tau)/7).                 (R8)

Floors are coordinatewise Euclidean floors, including negative numerators.
These identities retain every Q4 matrix without identifying phases. The
inverse/actual source certificate still decides which repairs are legal;
the residue equations alone never grant realization.

The loss is visible on actual equal-order sources:

    a t11547: v=(293,98,218),       order 8,
    t16510 a: v=(293,147,204),      order 8.

They have equal H, equal F49 data and equal charge modulo 7. Prefixing q2
gives H3557 and H3865. Their difference is 7(0,7,-2); the next observation
sees 308. Thus even H plus vertex clock plus the entire residue clock
cannot replace the integer carry.

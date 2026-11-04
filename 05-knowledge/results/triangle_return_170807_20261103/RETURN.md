# The missing56935 branch and its exact two-module reuse boundary

2026-11-03. **FINITE-EXACT literal source interface; PROVED all-module
response law.** This reaches170807 but does not close a ternary source
grammar or prove the full {7,21} spectrum theorem.

## R1. An actual strong interface retaining every old name

Start with the inherited16-vertex source of H=56935, in the exact order

    [1,2,4,6,7,8,14,21,22,31,32,33,39,40,50,51].

Reverse the two old edges on pairs {39,51} and {31,33}, in the legal order
saved in interface.json. The completed intermediate has H=67599. Add a
vertex54 with in-neighbor mask64641 in that old order. Equivalently its
old in-neighbors are {1,21,32,33,39,40,50,51}. The result is strong and has

    H=170807=3*56935+2.

Each edge reversal is compiled using one temporary helper, two actual
Hamilton-path complements, and deletion of only that helper. Thus the
entire program uses four legal native complements, two temporary births
and deletions, and one permanent birth. Peak order is17. All sixteen old
names remain. The paths have a common frame at both physical root50 and
spare6; every completed repair is strong. The star of port51 changes and
its new joint responses are retained, rather than identified with its old
responses. The whole operation is a count-changing rewrite followed by a
birth, not a fictitious1-ready source at H56935.

The exact named program, full matrices, intermediate counts, joint-port
responses and phase codecs are in interface.json and audit.json. The prior
complete no1-ready theorem at56935 is unchanged.

The current F49/charge coordinates pass through
(1,4;0),(0,3;0),(3,5;2),(3,0;0). The final F49 pair is Frobenius-fixed,
while the F343 vector codec has the nontrivial three-step orbit
(0,3,0)->(5,6,2)->(6,5,5)->(0,3,0). Both are retained. They are different
observer field actions; this is not a computation of tournament automorphisms
or a license to forget the final F343 phase.

## R2. Complete finite interface exclusions

For the single pinned source above, exact search gives:

* No target170807 using at most three arbitrary births alone. Monotonicity
  of H under birth prunes65536 first cuts to46, then all second cuts to192;
  every third cut of each retained prefix is tested. There are zero weak
  hits as well as zero strong hits. Thus an edge-preserving extension at
  this count requires at least four added vertices, if it exists.
* No strong target from one-old-edge reversal followed by at most two arbitrary births.
  All120 reversals and6583 retained first-birth prefixes are tested.
* No zero/one/two old-edge reversals followed by an arbitrary re-marking
  and the specific two-birth ready reset. The two-edge menu has7140 pairs.
* Two old-edge reversals followed by one arbitrary birth succeed at the
  displayed witness, found at pair6771 of the fixed search order. The rest
  of this discovery menu is not claimed exhausted.

The first three bullets even allow arbitrary edge flips without demanding
native guards, so their exclusions remain valid in the corresponding
guarded sublanguages. The successful flips are separately guard-compiled.
The witness minimizes the number of old-edge changes **among programs
with one arbitrary birth and otherwise only old-edge reversals**. No
minimal number of actual complements, helpers, or arbitrary operations
is asserted.

The pruning theorem is the inherited positive cut formula: every birth
has H(new)>=H(old). Therefore a prefix already above170807 can never
return to it using births alone. It does not justify pruning a later
complement or deletion. The underlying named first/second cuts are saved.

## R3. Arbitrary replacements in two modules have an exact bilinear law

There are two untouched transitive modules in both the source and target:

    A={6,8}, ordered6->8;  B={14,7}, ordered14->7.

Every external out-neighbor of A lies in B. The sole external out-neighbor
of B is31. These assertions hold after the rewrite and the permanent
birth:54 dominates both modules. Replace A and B simultaneously by any
nonempty tournaments, holding every exterior incidence fixed. Write
h_A=H(A), c_A=C2(A), where C2 counts ordered nonempty two-path covers;
define h_B,c_B similarly. Both full substituted tournaments remain strong.

No Hamilton word can visit B more than twice, since each nonfinal visit
must exit through the single physical vertex31. The same holds for A:
between two A visits the path must leave A through B, then leave B
through31. Consequently the maximal-run expansion is bilinear in
(h_A,c_A) and(h_B,c_B), with nonnegative integer coefficients. Four exact
TT1/TT2 substitutions determine those coefficients, giving

    H_source =12631 h_A h_B +4449 c_A h_B
              +8805 h_A c_B +4449 c_A c_B,
    H_target =37943 h_A h_B +13297 c_A h_B
              +26541 h_A c_B +13297 c_A c_B.             (R3)

This is an all-size theorem: the run bound proves the functional form;
the four finite exact coefficient computations are supplied by two
independent Hamilton-count implementations. Additional C3 substitutions
test a nontransitive input, not another fitted degree of freedom.

Subtracting gives the exact offset:

    (H_target-3H_source)/2
       =25(h_A-c_A)h_B+(63h_A-25c_A)c_B.                 (R4)

## R4. The successful+2 base is isolated in this entire substitution bank

For any n-vertex tournament X, cutting each Hamilton path at one of its
n-1 internal positions injects (n-1)H(X) objects into C2(X). Hence
C2(X)>=(n-1)H(X). This positive bound completely classifies offsets2,4
in(R4), even when A and B have arbitrary internal structure:

    H_target=3H_source+2 iff A=TT2 and B=TT2;
    H_target=3H_source+4 never occurs.                   (R5)

Proof. If |A|>=4, then c_A>=3h_A, and both coefficients in(R4) are
negative. If A=TT1 its value is25h_B+63c_B>=25. If A=TT3 it is negative.
If A=C3, the value is-75h_B+39c_B, a multiple of3, hence neither1 nor2.
The remaining case A=TT2 gives13c_B-25h_B. For |B|>=4 this is at least
14h_B>=14. For B=TT1,TT2,TT3,C3 its values are respectively-25,1,53,3.
These are all cases. The two orders of a transitive pair describe the
same module isomorphism class, with physical port incidences retained.

In particular this is stronger than failure of a few clone tests. If
vertex6 and vertex7 are inflated to transitive blocks of sizes a,b>=1,
put x=2^(a-1)-1 and y=2^(b-1)-1. The exact offset is

    H_target-3H_source=2-600x+104y-800xy.                 (R6)

It equals2 only at x=y=0, and is always2 modulo8, so it never equals4.
The mixed term is the missing joint incidence in a separate-port summary.
All twelve untouched singleton clone directions also fail already at size2;
these are finite controls, while(R5) handles arbitrary replacement of A,B.

## Next return target

The170807 interface is now concrete, but its output must be returned to a
new grammar before either ternary branch can iterate. R5 says that keeping
this quotient and merely enriching its two sink-side modules cannot do it.
The next useful change is an additional exit from B, or a context that can
observe three module visits, with its full joint response retained.
Incoming [gate-flow work](../circuit_flow_tournaments_20261103/README.md)
independently finds a gate-order commutator measured by a three-run term.
TR1 and that commutator have the same order-of-observation boundary; their
operations and coefficients are different and have not been identified.

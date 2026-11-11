# Odd line seeds and the symmetry of extension counts

**Status: FINITE-EXACT constructions and PROVED IN NOTE restrictions;
same-agent review with independent counting algorithms. The universal seed
construction and Chan-Pak remain OPEN.**

The first checkpoint gives explicit unit-offset roots of order eleven for
18 of the 33 lines left without a seed in the precisely bounded bank of
[HYP-9658](../../hypotheses/HYP-9658-unit-offset-seeds-of-the-two-adic-lines.md).
They include the least missing line, c=245. These are additional witnesses,
not a census or an assertion of minimum order. Every first and second tail
point is independently recounted by subset dynamic programming and path
backtracking in `verify.py`; literal adjacency lists are in `search.jsonl`.

For c=245 the deletion base has rows
`(552,517,49,534,67,346,655,63,223,436)` (bit j means an arc to j), selected
set O={1,4,5}, H(B)=6965, W(O)=6964 and U_B(O)=31360=245*2^7.
Thus the inherited tail construction gives H=1+245*2^s for every s>=8.

The inherited mechanism is
[THM-5126](../../../01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md):
for a tournament B and a selected set O, add B -> TT_t -> v with v -> O
and B minus O -> v. Its Hamiltonian count is

    H = H(B) - W_B(O) + 2^t U_B(O), t >= 1.

Here W sums ordered two-path covers with first end and second start in O;
U sums covers whose second path starts in O and whose first path may be
empty. Offset one is exactly W=H(B)-1.

## Symmetry restriction

Let G=Aut(B,O), the automorphisms of B preserving O as a set. G acts freely
on Hamiltonian paths and on each of the above ordered covering-path
families: an element fixing an ordered cover fixes every vertex in its
specified list position. Consequently |G| divides H(B), W_B(O), U_B(O),
and their difference H(B)-W_B(O). **Every unit-offset seed has trivial
selected-set stabilizer.** The base B itself need not be asymmetric.

This is a necessary condition, not a construction or converse. In
particular a transitive two-vertex base with both vertices selected has
trivial stabilizer but offset -1. The same proof applies to every integer
offset: a preserved symmetry of order g forces g to divide that offset.

## Reproduction and search scope

`bases.txt` contains 1,347 distinct nine-vertex deletions of tail seeds of
order ten in the inherited `line_seeds.txt` bank. Search seed 20261110
alternates a random labelled ten-vertex tournament and a random one-vertex
extension of one of those bases. It tests **every** selected subset of each
sampled base, with no selected-cardinality filter. The 200,000 sampled
presentations need not be distinct and are not uniformly distributed over
isomorphism classes. The stopping bound was chosen before this run;
negative output has no all-order implication.

    g++ -O3 -std=c++17 search.cpp -o search.exe
    ./search.exe bases.txt 200000 > search.jsonl 2> search-summary.txt
    python3 -B verify.py > verification.out

The executable is temporary and is not a research artifact. `verify.py`
also checks symmetry divisibility on every labelled tournament through
order five and every selected subset, the exact qualifying-path cylinder
formula for small extensions, the indicated inherited number seeds, and
the seven-vertex Paley controls. All computations use integer arithmetic
or exact rational means. No Lean formalization or independent agent audit
is claimed.

Further synthesis of the extension compiler and the two supplied papers
will be recorded in this packet. Reserved THM-5130 and its exploratory
sign-flip recursion are not used as proved dependencies.

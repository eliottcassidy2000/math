# Odd line seeds and the symmetry of extension counts

**Status: FINITE-EXACT constructions and PROVED IN NOTE restrictions;
same-agent review with independent counting algorithms. The universal seed
construction and Chan-Pak remain OPEN.**

This packet gives explicit unit-offset roots of order eleven for
all 33 lines left without a seed in the precisely bounded bank of
[HYP-9658](../../hypotheses/HYP-9658-unit-offset-seeds-of-the-two-adic-lines.md).

**2026-11-11 continuation:** the [aspect-sequence census](../metagraph_aspect_sequences_20261111/PROOFS.md)
independently replays the complete nine-vertex-base exclusion and proves source
padding preserves a seed. Hence these 33 roots have minimum order eleven
within the selected-set unit-offset tail construction. This does not assert
minimum order for arbitrary realizations of individual ray values.

Thus HYP-9658(B) is settled, with order eleven rather than twelve.
The universal statement HYP-9658(A) remains OPEN. These are additional witnesses,
not a census; the minimum-order refinement above is subsequent work. Every first and second tail
point is independently recounted by subset dynamic programming and path
backtracking in `verify.py`; all literal adjacency lists are in `verification.json`.

For c=245 the deletion base has rows
`(552,517,49,534,67,346,655,63,223,436)` (bit j means an arc to j), selected
set O={1,4,5}, H(B)=6965, W(O)=6964 and U_B(O)=31360=245*2^7.
Thus the inherited tail construction gives H=1+245*2^s for every s>=8.

**Infinite consequence.** Every positive odd N>1 whose N-1 has odd part
at most 4097 is a tournament Hamiltonian-path count, except 7 and 21.
This extends the inherited 65 certified lines to 2,049 lines; it is still
a density-zero subset of the odd integers, not an eventual threshold.

Proof: the inherited bank supplies 2,016 of these lines; the 33 new roots
complete it. The largest point below the chosen ray starts is 63,857,
on c=3991. Every such point except 7 and 21 is already covered by the
independently audited finite range through 2^32 in
[THM-4875 D](../../../01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md).
The later, single-run 2^35 extension is unnecessary. `line_thresholds.json`
and `verification.json` check every line and its below-ray bound; the
inherited 2,016 certificates and 2^32 computation are used as proved inputs,
not claimed as newly rerun. The general exclusions are inherited from
[THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).

The chosen starts s0 of the new rays are:

| c | s0 | c | s0 | c | s0 |
|---:|---:|---:|---:|---:|---:|
|245|8|2023|4|3427|1|
|853|6|2037|4|3545|3|
|909|4|2085|3|3713|4|
|1049|4|2343|2|3733|2|
|1407|5|2393|4|3753|4|
|1467|4|2703|5|3773|3|
|1515|6|3001|1|3907|3|
|1577|2|3033|3|3909|1|
|1629|3|3123|3|3985|4|
|1829|3|3313|1|3991|5|
|1887|4|3411|5|4023|4|

These starts are certified choices, not minimum possible exponents.

As a corollary, every binary-weight-three number
1+2^b+2^(b+d), with b>=1 and 1<=d<=12, is a Hamiltonian count except
7 and 21. Indeed its line parameter is c=1+2^d<=4097. The inherited
corollary covered d<=7; unbounded d is still open.

## Construct a certified tournament

Run `python3 -B construct.py --c 245 --exponent 8` in this directory.
It emits the literal twelve-vertex adjacency rows and expected H=62721.
Every odd c from 1 through 4097 is accepted; omitting `--exponent` uses
its certified first ray point. Larger exponents append the prescribed
tail or inflate a certified doubling vertex. The constructor explicitly
rejects exponents below that chosen start; those points use the inherited
finite construction in the proof above. It performs no new search.

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
sampled base, with no selected-cardinality filter. The 1,000,000 sampled
presentations need not be distinct and are not uniformly distributed over
isomorphism classes. Sampling was expanded after a successful 200,000
checkpoint; no statistical prediction or preregistration claim is made.
Negative output has no all-order implication. A second run enumerates all
1,347*512=689,664 one-vertex extensions of the frozen base bank, testing
every selected subset. That restricted exhaustive run finds 32 target lines;
the sampled run supplies the remaining c=3753. It is exhaustive only over
the declared base presentations, not over all ten-vertex tournament classes.

    g++ -O3 -std=c++17 search.cpp -o search.exe
    ./search.exe bases.txt 1000000 > search.jsonl 2> search-summary.txt
    ./search.exe bases.txt 0 exhaustive > exhaustive.jsonl 2> exhaustive-summary.txt
    python3 -B verify.py > verification.out
    python3 -B check_constructor.py

The executable is temporary and is not a research artifact. `verify.py`
also checks symmetry divisibility on every labelled tournament through
order five and every selected subset, the exact qualifying-path cylinder
formula for small extensions, the indicated inherited number seeds, and
the seven-vertex Paley controls. All computations use integer arithmetic
or exact rational means. No Lean formalization or independent agent audit
is claimed.

Read [the qualifying-path compiler](EXTENSIONS.md) for the degree 2k
extension polynomial, the exact target-support obstruction, and the
inherited sufficient threshold 95 induction. Read
[symmetry, prime geometry and the papers](PRIMES_AND_PAPERS.md) for the
7/1039 doubling action, the fixed-valuation Collatz character, the
3/13/27 field obstruction, and the two explicit paper transfers.
[Research board](BOARD.md) records inheritance, scope and stopping points.
Reserved THM-5130 and its exploratory sign-flip recursion are not used as
proved dependencies.

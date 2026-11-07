# Asynchronous physical escape from the old paired-carrier synchronization wall

Status: PROVED physical two-residue infinite families and synchronization criterion; FINITE-EXACT literal witnesses.
No all-odd threshold, F(8,137), or novel Hamilton-count value is claimed. In particular, old coverage already extends far
beyond these counts. Novelty is the J-containing asynchronous construction
route, relative to the four fixed Session 203 carriers and Session 214's common
word obstruction, not the scalar Hamilton-count spectrum.

## Inheritance and concept board

Pinned main: 0f8750677dd20b594a0e032a3b56aa23f58a950c.
Closest mechanism: Session 203 physical nine-state O/Z/J transfer and its
four paired carriers, frozen in Session 214's paired-coverage/dependencies.
Hostile: a shared J irreversibly separates each original pair (Session 214).
Corrected near miss: matching X is not matching the continuation state.
Least-used sidecar: independent word histories with independent source choice.
Concept board: asynchronous scalar collision; persistent four-state
synchronization; pure-Z rays; O^kJ rays; cyclic module substitution.
The last route has a complete inherited walk-content law in
THM-3121-path-cover-walk-content-substitution-kernel.md. It is not a new
composition theorem here and cannot be replaced by scalar H substitution.

## New persistent physical pairs: reset the missing endpoint coordinate

Theorem. The following two pairs are actual bounded-start paired carriers at
new bases X=4225 and X=4357, and every common O/Z word preserves their equal
adjusted bases. Thus each supplies an infinite family of simultaneous
certificates for the two residues in F(8,137).

For X=4225 start with these rooted eight-vertex graphs and apply JO:
- A: root 6; rows (2,20,41,67,77,27,167,63).
- B: root 4; rows (6,12,24,49,195,23,47,111).
The output states are respectively
(703,25,0,3531,1422,195,31816,7488,24501),
(537,11,0,3693,1086,195,28056,7824,21699).
Take j_A=8,j_B=3. Both have (X,e,V,d)=(4225,0,8910,195).
The compiler counts at the empty continuation are 12677 and 12679.

For X=4357 start with these rooted seven-vertex graphs and apply ZJO:
- A: root 3; rows (10,20,65,22,37,15,59).
- B: root 2; rows (2,48,11,67,45,77,23).
The output states are respectively
(511,13,0,3851,1084,251,26344,8752,21333),
(533,14,0,3830,1132,251,27056,8704,21858).
Take j_A=4,j_B=4. Both have (X,e,V,d)=(4357,0,9836,251).
The compiler counts at the empty continuation are 13073 and 13075.

Proof. Direct exact counting from the displayed adjacency rows, followed by
the all-input inherited transfer, gives these states. The starts have forms
1+3j_A and 2+3j_B and all j lie in [0,8]. The synchronization theorem below
then applies to every common O/Z word. Counts from the fixed four-vertex
context are exactly 3X+2 and 3X+4. QED.

In particular pure Z continuations give, for every k>=1,
X=8715*2^k-4295 and X=9585*2^k-4977, respectively.
These yield the two paired Hamilton-count families 3X+2,3X+4.
The full common binary O/Z language is certified, not only these rays.

Mechanism: re-pair entirely different sources before J. After their common
J-containing prefix the X,V,d coordinates agree but e need not. A final O
resets e to zero and doubles the already equal V, making the pair stable
under the complete binary continuation language. This is a positive use of
the information erased by O; no signed coordinate is treated as realizable.
It escapes the four old sources' common-J wall by changing the source pair.
Before the reset, the 4225 pair has (X,V,d)=(4225,4455,195) but e=187
versus 189; the 4357 pair has (X,V,d)=(4357,4918,251) but e=220 versus 218.
Thus immediately using Z would break equality, whereas O converts the
imperfect matches into permanent binary carriers.

These two initial bases are absent from the four old O/Z scalar orbits.
Indeed for M=4357 the inherited order cap is floor(log2 M)+2=14; all old
sources have order 8. After removing terminal O moves, any represented base
<=M therefore occurs within six letters. Exhaustive exact evaluation of all
four old binary trees to depth 6 excludes 4225 and 4357. Moreover the new
four-states and every common O/Z descendant lie outside the old four-state
orbits, since their d values 195 and 251 differ from all old values 49,29,109,33.
This does not prove every new ray base avoids the old scalar image.
Nor does it make the displayed Hamilton counts new to the inherited spectrum.

Discovery used a supplied bank of 7412 tournament rows through order 8, without
requiring its class-completeness. Filter: nontrivial sources with 1<=s<=26 and
s nonzero modulo 3. All 12893 retained rooted sources were extended by all
words of length <=3. Hashing (X,e,V,d) yielded 341727 keys and 20 keys with both
residues and at least one retained J-containing representative. This count
is about stored representatives, not all word presentations or isomorphism
classes. NEW_BANK.txt gives the exact matches; NEW_CARRIERS.json gives the
four selected full literal graphs, compiler graphs, and independent counts.
The witnesses and the theorem rely only on the four displayed seeds, not on
any completeness claim about the bank or search.

## Explicit positive constructions

Use coordinates p=(h,s,e,L,R,d,C,b,T), X=h+L-s+2j, and literal physical
O/Z/J definitions in ROOTED_CLOSURE.md. Words act from left to right.
All starts and j labels remain fixed under these nontrivial source extensions.

1. Start from carrier 1547 A, root 0, outgoing rows
(2,148,97,7,45,75,27,125). Apply JJ. The resulting state is
(793,13,267,8169,951,425,30400,11476,18507), j=4, r=1.
Thus X=8957.

2. Start from carrier 759 B, root 1, outgoing rows
(2,52,25,35,105,197,15,95). Apply ZOOZ. The resulting state is
(2691,11,29,6271,4331,29,246624,10240,204157), j=3, r=2.
Again X=8957.

Both j values are in [0,8]. The inherited four-vertex rooted context with
rows (2,8,3,5) therefore gives Hamilton counts 26873 and 26875.
This replaces one source and uses asynchronous histories; Session 214's
common-word exclusion is respected, not contradicted.

There is also an escape with no source replacement:
677 A under ZZJJJ has state
(4607,10,1552,149996,5408,2353,561640,208792,339616),
and 677 B under ZOOJJ has state
(3989,11,1796,150615,4386,2193,494812,190944,287735).
Their fixed j values are both 3 and their adjusted bases are both 154599.
The two literal compiled counts are 463799 and 463801.

These are rigorous finite constructions: inherited exact transfer plus the
stated source rows proves them, and verify_witnesses.py independently recounts
all nine output coordinates and both compiled graphs. WITNESSES.json includes
full literal rooted and compiled adjacency rows; no implicit realizability
assumption is used.

## Exact universal O/Z synchronization criterion

Theorem. Two physical carriers with fixed start labels and adjusted bases
have equal X after every common word over O,Z (including the empty word)
if and only if their initial (X,e,V,d) agree, where V=R+b.

Proof. Sufficiency follows from the inherited exact maps
O:(X,e,V,d)->(X,0,2V,d),
Z:(X,e,V,d)->(X+V-e,d,2V-d,d).
For necessity, equality at the empty word gives delta X=0. Equality after
Z gives delta V=delta e. Equality after OZ gives 2 delta V=0, hence both
are zero. With these zero, equality after ZZ gives -2 delta d=0.
Thus every required coordinate agrees. No positivity or unproved state
realization condition is used. QED.

Corollary. If two carriers have different d, no common O/Z prefix can make
them universally synchronized: d is invariant. If d and e agree, X agrees,
and V differs, then every common word containing Z separates X, since
at each stage delta V=2^k delta V_initial and delta e=0 after every move.
The sum of increments has its fixed nonzero sign. Words containing only O
preserve the scalar X but yield no new base.

This distinguishes the two finite asynchronous collisions above from the
new persistent pairs: each finite collision has different d on its two sides.

## Finite searches and attempted persistent positive extensions

search.py enumerates every O/Z/J word of length at most 9 from each of the
eight sources. It merges only identical full nine-states and stores a shortest
word; each state has identical continuations. It discards terminal-O recipes
from the scalar inventory because deleting that suffix preserves X. For each
source and X it retains one first-discovered representative. This is not an
enumeration of every collision recipe, and the J-containing-filter count must
not be interpreted as the count of all possible J-containing intersections.
It yields 122 selected A/B intersections whose stored pair has a J somewhere,
including cross-source pairs; collisions.json records each exact recipe.
Among those 122 selected collisions, only X=881319 has equal d. There,
e_A=e_B=d=253 but V_A=1232967 and V_B=1267879, so the corollary excludes any
common word containing Z from maintaining or restoring equality. No collision
in this selected list can become a universally synchronized carrier by any
common O/Z continuation. This is a statement about the list, not all words.

Two separate exact, bounded searches sought positive infinite families.
For any state, Z^k (k>=1) has base
X_k=(X-V+2d-e)+2^k(V-d).
Here V-d>0. Two such rays have an infinite intersection exactly when their
constant offsets and odd parts of positive coefficients agree. Indeed unequal
offsets allow at most finitely many equalities, by 2-adic valuation after
factoring the smaller power of 2; equal offsets reduce to a power-of-two ratio.
The search families.py checks all states represented by words of length <=10
from each source and retains representatives by (offset,oddpart(coefficient)).
No matched A/B signatures involving a stored J-containing history occur.
The finite search does not classify arbitrary asynchronous suffixes.

Similarly, O^kJ for k>=1 gives
X_k=(X-2L)+2^k(R+C).
The analogous signature search jrays.py again finds zero selected
J-containing prefix matches through length 10. Appended J is inherent in this
ray even when the prefix has no J; the script's reported count excludes those
pure O/Z prefix matches, so it must not be read as a blanket no-ray theorem.
These are failed positive-family probes, not promoted global obstructions.

## Reproduction

g++ -O2 -std=c++17 new_bank.cpp -o new_bank
./new_bank BANK.txt > NEW_BANK.txt
python3 verify_new_carriers.py > NEW_CARRIERS.json
python3 search.py > SEARCH.txt
python3 families.py > FAMILIES.txt
python3 jrays.py > JRAYS.txt
python3 verify_witnesses.py > WITNESSES.json
python3 -O verify_witnesses.py > WITNESSES_OPTIMIZED.json
cmp WITNESSES.json WITNESSES_OPTIMIZED.json

inherited.py is an exact copy of Session 214's paired-coverage/verify.py.
Only its definitions are imported; its main program is not run by the new
scripts. Provenance and hashes are in PROVENANCE.json.

Remaining frontier: widen the successful source-pair/reset search to derive a target-dependent family covering a tail, rather than merely adding finitely many binary carrier orbits. The new infinite families do not imply any coverage theorem.

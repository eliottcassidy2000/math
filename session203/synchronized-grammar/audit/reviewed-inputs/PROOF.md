# Eight short tests characterize paired synchronization under all O/Z/J words

Status: proof candidate. This is an exact reduction for a possible constructive bridge, not a claim that the required physical seed pair exists or that its orbit has a tail.

Use the physical nine-coordinate states and transfers in the independently reviewed `dependencies/ROOTED_COVER_CLOSURE.md`, in coordinate order `(h,s,e,L,R,d,C,b,T)`. Inputs are nontrivial rooted tournaments with positive start count, so the singleton Z exception is excluded. Each operation preserves s. Define `G=h+L`, so the staggered adjusted base for residue r is

    X= B+2j = G-(s+2r)/3,  where s=r+3j.

## Exact synchronization criterion

Fix two such physical sources, assigned to residues r=1 and r=2, with fixed admissible j values. Assume their adjusted bases agree initially. The following conditions are equivalent:

1. Their adjusted bases agree after every common finite word in O,Z,J.
2. Their adjusted bases agree after each of the eight words

       Z, J, OZ, OJ, ZZ, ZJ, JJ, JOZ.

3. Their eight coordinates `(h,e,L,R,d,C,b,T)` are identical, and their starts satisfy `s_1=s_2+2`.

In condition3, the j relation is `j_1=j_2+1`. Thus for width8 the larger start can be any of4,7,...,25 and the smaller start is two less. The original small source pair still must be physically realized; arithmetic equality of a freely chosen state is not sufficient.

**Proof.** Let delta be the difference of the two nine-coordinate column states. Since starts are preserved, equality of their adjusted bases after a word w, together with equality at the empty word, says exactly

    (l M_w-l) delta=0,  l=(1,0,0,1,0,0,0,0,0).

Here M_w is chronological, so its rightmost matrix is the first physical letter. For the eight words in the listed order, the row matrix of these conditions is

    [0 0 -1  0 1  0 0  1 0]
    [0 0 -1 -2 1  0 1  0 0]
    [0 0  0  0 2  0 0  2 0]
    [0 0  0 -2 2  0 2  0 0]
    [0 0 -1  0 3 -2 0  3 0]
    [0 0 -1 -4 3 -2 2 -1 2]
    [0 0  0 -2 2  0 1  0 2]
    [2 0 -1 -2 3  0 3  0 0].

Deleting the zero s column gives determinant32. Consequently all eight other differences vanish. Initial adjusted-base equality now says `(s_1+2)/3=(s_2+4)/3`, giving the stated start difference2. This proves2 implies3. Condition1 immediately implies2.

Conversely the displayed physical transfers of all eight non-s coordinates do not involve s. Thus equality of those coordinates is preserved by every letter. Starts remain constant, and the relation `s_1+2=s_2+4` proves adjusted-base equality for every continuation. This proves3 implies1. QED.

## Constructive consequence and remaining obligations

A physical seed pair satisfying3 would lift the paired-carrier strategy from the O/Z language to the full O/Z/J language without losing either staggered residue. Every common continuation would retain the same two allowed j values, and the inherited four-vertex context would construct actual counts3X+2 and3X+4. To prove an eventual threshold, one would still need an all-target occupancy theorem for the resulting common X spectrum. Seed-pair existence alone would not establish it.

The eight tests provide a finite decisive test of an all-word pairing assertion. If initially synchronized sources fail to match any non-s coordinate, some listed word of length at most3 separates their adjusted bases. Thus searching longer common words cannot rescue a claim of synchronization under *every* word.

## Reproduction and scope

The theorem is checked independently of any discovery search by exact symbolic multiplication in `verify.py`, including determinant32 and the complete nullspace. Normal and optimized Python runs agree. It requires no graph census or scalar-coverage computation. Neither existence nor nonexistence of a physical seed pair is asserted.

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

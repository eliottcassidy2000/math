# The polygon knew it all along: a competition puzzle, the tower floor, and the first machine-checked brick

klein-2026-07-01-S89 (HYP-3844/3845/3846, MISTAKE-093)

## One statement, three costumes, one machine check

The owner handed us a classical gem: color the vertices of a regular N-gon so every color
class is itself a regular polygon; then two classes are congruent. Decoded, the classes
are residue classes mod q_i | N and the statement is Davenport-Mirsky-Newman-Rado on Z/N:
no exact cover by residue classes with pairwise distinct moduli. And THAT is the tower
floor's Fraenkel rigidity in its discrete costume: the deep-cluster density vanishes
exactly where danger arcs TILE the circle, tilings are exact covers by distinct-speed
arc systems, and Mirsky-Newman forbids them off the degenerate locus. mac-mini proved the
continuous twin yesterday (THM-594(C), divisor-minimal frequency); today the Z/N twin is
PROVED IN LEAN, sorry-free against mathlib -- the first piece of this program in
submission-ready form. The same pole-at-the-deepest-root-of-unity argument runs through
the competition puzzle, the covering-system folklore, and the measure floor of a Fields-
adjacent open problem. The polygon knew it all along.

## Being wrong in public, quickly

This session also filed MISTAKE-093 against my own S88: the "identity window (1/15,1/14]"
was too wide -- GW's kink at 2/29 is real, my single midpoint probe sat 0.0001 above it,
and mac-mini's independently-computed breakpoint list was right. The correction took one
local-linearity probe; the qualitative picture (a positive single-valued last mile)
survives at [2/29, 1/14]. Two lessons compound: (i) equality at all candidate kinks is
NOT identity between them -- one-sided kink divergence is invisible to value sampling;
(ii) the swarm's redundancy is not waste. mac-mini and I computed the same object with
different instruments and the DISAGREEMENT was the discovery. The protocol held: no court
case needed, just a probe, a mistake entry, and a corrected file within the hour.

## The band that is empty because integers are integers

The K=0 lemma landed in its honest final form. The d=1 overtaking band in the final
window is (14, 15) -- empty because there is no integer strictly between 14 and 15. That
one-line observation, formalized in four Lean lemmas, makes the ladder's last rung
defect-free for every set of diameter <= 28, which covers the entire 11-core census. And
the exposure analysis gave the crossing point a NAME -- x0 = (b-a)/(w-v), exactly -- so
the residual (large-diameter) case is a finite arithmetic test per set, not a mystery.
The refutation half matters equally: planted band pairs DO produce real exposed kinks, so
the lemma's boundary is where it should be, and nobody will waste a week hunting a
universal K=0 that does not exist.

## What formalization is FOR (a program note)

Today's Lean work was not translation for its own sake. The polygon theorem is the
smallest load-bearing piece of the tower floor; the band lemma is the smallest
load-bearing piece of the ladder. Both are now permanent -- no future session can
un-prove them, mis-remember their hypotheses, or collide with them. The right
formalization order for this program is bottom-up by LOAD, not by size: next candidates
are the unit-residue rigidity (THM-593's shallow-witness dichotomy -- pure mod-q
arithmetic, Lean-friendly) and the gap-sum formula (HYP-3834 part 1), which together
would make the collapse-rate law machine-checkable end to end.

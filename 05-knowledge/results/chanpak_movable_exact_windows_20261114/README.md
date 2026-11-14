# Moving exact Chan–Pak windows with physical core changes

**New exact families, not universal window placement.** The
[previous complementary tails](../chanpak_tight_mod_six_20261114/README.md)
gave windows with a shared binary high prefix. This session
allows different sources and independently chosen readers,
and proves exactly what arithmetic must match to put a window
at the target.

1. **[Two independently movable three-point windows](MOVABLE_TRIPLES.md).**
   Two banks of three literal cores have identical six-coordinate
   boundary data and ordinary counts spaced by two. Thus both
   wing addresses can move independently, with no symmetry or
   high-prefix restriction. The six-vertex bank switches even
   between a prime and a decomposable puncture without changing
   the exposed data. At addresses (1,12) it realizes **7605,
   7607, 7609**, where 7607 and 1087 are prime and 7609=7·1087.
   The entire earlier complementary chart rejects all three.
   The exact inverse for these banks is implemented, and the
   full literal counts are checked.
2. **[Four- and five-point windows by a tilted reader](DIAGONAL_WINDOWS.md).**
   A linear combination of the two boundary responses suffices
   on the tilted line. The remaining offset is an actual count
   of incompatible path pairs. One slope gives five consecutive
   odd counts, with mod-10 selection and error at most 8.
   At its first named height this gives **5953 through 5961**,
   odd values only, again outside the entire earlier
   complementary chart. Replay includes literal seven- and
   eight-vertex puncture banks and an independently checked
   finite scan. The finite maxima do not imply any all-order
   maximum window length.
3. **[Synchronize source size and reader width](ANCHOR_LAW.md).**
   A common quadratic coefficient and common linear coefficient
   can align a window even across different two-path-cover
   counts and different wing widths. The new phase formula
   gives exact signed offsets without changing that baseline.
   An exact worksheet for x=2,y=1 reduces prescribed placement
   to a submodular quadratic subset sum with coefficients from
   literal endpoint pairs. These coefficients are not freely
   selectable: four contacts cannot have all six pair penalties
   equal to one.

On any constructed window the mod-6 or mod-10 choice is a real
absolute proof: the error is a multiple of the modulus and
strictly smaller than it. The missing step is a physical choice
of anchor at every prescribed odd target. Stacking identical
full-profile windows grows the number of options but multiplies
their spacing; freely varying dyadic tilted scales of a fixed
bank is still sparse. Neither is the missing all-target argument.

Run `python3 replay_independent.py 7607` for a literal witness
using independent wing addresses, or
`python3 replay_relocation.py 5953` for a tilted witness.
The companion replays test the named targets and exact formulas.

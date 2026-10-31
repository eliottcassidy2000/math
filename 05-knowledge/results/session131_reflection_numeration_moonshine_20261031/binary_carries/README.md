# Session131: carry-correct 1–2–4 / Hamiltonian bridge

Main result: for a tournament with an involutory anti-automorphism,

    (b0,b1,b2)=(1, q mod2, (floor(q/2)+i) mod2),
    H=1+2q+4i+8(j+p+t) mod16.

Here q counts fixed literal odd cycles; i/j count intersecting/disjoint
nonfixed mirror-cycle pairs; p counts nonfixed orbits of disjoint
TWO-CYCLE PACKINGS; t counts fixed THREE-CYCLE PACKINGS. The actual fourth
bit is `(floor((q+2i)/4)+j+p+t) mod2`. See `PROOFS.md` for complete proofs,
the exact positive lift, the central-cut version, inherited sources, and
why neither arithmetic 7=1+2+4 nor Zeckendorf indices identify these counts
with binary digits.

`normal.json` and `optimized.json` each pass **106,736 explicit gates**:
all **4,695** reversal-compatible labelled tournaments at n=1..7; 128
seeded n8 and 48 seeded n9 samples; four triangle stacks (n3,6,9,12);
and two inherited OCF-collision fixtures. This is not a full labelled
all-tournament census. `hostile.normal.json` and `hostile.optimized.json`
independently check all eight distinct hostile tournaments by literal
Hamiltonian permutations and a separately implemented odd-cycle census,
with **87 explicit gates**. The normal/optimized files are byte-identical
within each pair.

Key hostile: the transitive n6 tournament has H1 and the two-triangle
ordered stack has H9, yet both have (q,i)=(0,0). Thus q and i do not
recover H modulo16. At n5, H9 with(q,i)=(2,1) requires carrying, while
H15 with(q,i)=(3,2) realizes residue7 despite i even. At n9 the three
triangle stack has H27 and(q,i,j,p,t)=(1,0,1,1,1), requiring the t term.

Reproduce using Python's standard library only:

```sh
python3 -B check_binary_carries.py > normal.replay.json
python3 -B -O check_binary_carries.py > optimized.replay.json
cmp normal.replay.json optimized.replay.json
cmp normal.json normal.replay.json
python3 -B verify_hostiles.py > hostile.normal.replay.json
python3 -B -O verify_hostiles.py > hostile.optimized.replay.json
cmp hostile.normal.replay.json hostile.optimized.replay.json
cmp hostile.normal.json hostile.normal.replay.json
```

Sources are pinned by commit, Git blob, and SHA256 in `SOURCE_PINS.json`.
The OCF attribution is Grinberg–Stanley; Irving–Omar v3 Corollary20 is the
checked restatement. The existing THM-466 correction, rather than the
stale early Zeckendorf prose, governs coefficient nonuniqueness. Nothing
here newly proves forbidden7/21, full H-spectrum completeness, an
all-order presentation-fibre gap, or a Collatz correspondence.

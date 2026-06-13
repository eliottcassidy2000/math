# Images Index

Images produced during the research process. Primarily TikZ diagrams from the LaTeX paper, preserved here for future use in the paper rebuild.

| ID | Filename | Description | Source | Used in |
|----|----------|-------------|--------|---------|
| IMG-001 | grid-n4.png (to be extracted) | Pin grid Grid(4), m=3 tiles, with σ-fixed diagonal r=c highlighted | parity_tournaments_fixed.tex §setup | Paper §setup |
| IMG-002 | grid-n5.png (to be extracted) | Pin grid Grid(5), m=6 tiles, with σ-fixed diagonal highlighted | parity_tournaments_fixed.tex §setup | Paper §setup |

---

## Notes

Images in the LaTeX source are TikZ diagrams. To extract them:
1. Compile `03-artifacts/drafts/parity_tournaments_fixed.tex` with pdflatex
2. Use pdfcrop + pdftoppm to extract individual figures
3. Save as PNG here and update this index

Alternatively, each TikZ figure can be compiled standalone for quicker extraction.

## Adding Images

Add a row to the table above with:
- A sequential IMG-NNN ID
- The filename (place the file in this directory)
- A description clear enough that someone who hasn't seen it knows what it shows
- Which source it came from
- Where it is (or will be) used in the paper

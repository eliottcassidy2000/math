# Independent review: quotient cycle parity and two-cover correction

PASS. No correction was needed. The concise all-order proof and every retained cycle/coefficient control were independently checked against the frozen source pins.

Rooting at the unique actual port x linearizes a Hamilton cycle without dividing by its length. It also prevents a child run from straddling the cut. Filling each ordered run list uses exactly one complete child-cover factor. Every repeated-child term is even, while all one-run Hamilton factors are odd. Thus the bare quotient-extension cycle parity is correct for arbitrary nonempty children, including singleton children and even total order.

The augmented physical order is at least three, so 2^N vanishes modulo eight. Division by the odd Hamilton product preserves the terms with coefficient four. This proves the displayed normalized q2 formula with the actual product retained. The coefficient relations B0 modulo eight and B_i,C_i modulo four follow with the stated factors; in particular C_i is 2Z modulo four, not an unscaled equality of integer cycle counts. Matching bare cycle parities removes the odd-order correction, but does not itself supply a nonzero local contrast. The fixed ambient must retain the prescribed physical exterior port incidences.

The independent checker imports no producer code. Its cycle counter uses predecessor subset recursion rooted at physical vertex zero, rather than the producer's forward recursion rooted at x. It reproduces all 2,060 exact retained cycle counts, with the closing arc included. Every complete coefficient row is checked against the two independently audited affine source banks, with no omitted rows: all 4,360 B/C relations and every base relation agree.

All 663 actual child-collapse records have their named physical graphs and block masks rebuilt independently. The normalized two-cover values and actual odd units are checked. A separate Hamilton-subset/ordered-convolution implementation reproduces the entire negative core profiles, not only their residues, and the actual ambient counts. The hostile PPO boundary has bare cycle counts 1 versus 2, core first-three gaps -68,-396,-1116 and ambient Hamilton gap -600, so the modular-loss warning is retained exactly.

Normal and optimized independent receipts agree byte for byte. The complete values remain in the bound source packet. This is a structural compression and sufficient ambient-correction test, not arbitrary-quotient non-annihilation or reconstruction.

Reproduce with python audit.py and python -O audit.py. Exact six-file audit allowlist: REVIEW.md, audit.py, cycles.cpp, normal.json, optimized.json, SHA256SUMS.

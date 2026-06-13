# Spectral Analysis: The Personality of a Tournament

## What is Spectral Analysis?

"Spectral analysis" sounds technical, but the idea comes from something physical: the way light breaks into a rainbow when passed through a prism.

When white light hits a prism, it splits into its component colors — red, orange, yellow, green, blue, violet. Each color corresponds to a specific wavelength. The spectrum of light is the fingerprint of its composition.

Mathematicians do something analogous with graphs and matrices. You encode a mathematical object (like a tournament) as a matrix of numbers, then look at the "eigenvalues" of that matrix. These eigenvalues are like the characteristic frequencies — the "color spectrum" — of the object.

---

## What is a Matrix?

A matrix is just a rectangular grid of numbers. For a tournament on 5 players, you might build a 5×5 grid where the entry in row i and column j records something about the relationship between player i and player j.

The entries we care about: position-weighted path counts. For each pair of players (i, j), count the number of Hamiltonian paths where player i starts at a specific position and player j starts at another specific position. This builds a matrix M(T).

---

## What are Eigenvalues?

A matrix "acts" on lists of numbers — it transforms them. An eigenvalue is a special number λ such that when the matrix acts on a particular list of numbers v (called an eigenvector), the result is just λ times v: the list gets scaled, but not rotated or mixed.

Think of it this way: most transformations are like rotations — they mix everything up. But some special directions (the eigenvectors) just get stretched or shrunk — and the stretch factor is the eigenvalue.

The eigenvalues of a matrix tell you its fundamental behavior. A matrix with eigenvalues all close to 0 "damps" everything — multiply by it repeatedly and everything shrinks. A matrix with large eigenvalues amplifies things.

---

## Why Apply This to Tournaments?

We built the matrix M(T) — where entry M[i,j] counts paths weighted by the starting positions of players i and j — and studied its eigenvalues.

The main discoveries:

### 1. The Matrix is Always Symmetric

M[i,j] = M[j,i] for every pair of players i and j.

This means the number of paths "starting with player i in position k and player j in position l" equals the number of paths "starting with player j in position k and player i in position l." This is not obvious — it's an algebraic miracle that comes from a hidden pairing of paths (technically, from an involution on permutations).

A symmetric matrix has a beautiful property: all its eigenvalues are real numbers (not complex).

### 2. The Eigenvalues Follow a Cosine Formula

For rotationally symmetric tournaments (like Paley and cyclic interval tournaments), the eigenvalues of M(T) follow a perfect cosine pattern:

**λⱼ = 3 + 2cos(jπ/K)**

where K = ⌊(n+1)/2⌋ and j ranges over 1, 2, ..., K.

This is extraordinary. The "spectrum" of the tournament matrix — its characteristic frequencies — forms a perfect cosine wave. The tournament's algebraic structure is as regular as a pure musical tone.

### 3. The Spectral Radius Approaches √5

As the tournament size grows, the largest eigenvalue approaches **√5 ≈ 2.236...**

This constant √5 appears naturally in the golden ratio (φ = (1+√5)/2 ≈ 1.618). The connection isn't coincidental — the path-counting recursion in tournaments has a structure similar to Fibonacci-type recurrences, which have the golden ratio built in.

### 4. Transitive Tournaments Give a Pure Signal

The transitive tournament (complete ranking, exactly one Hamiltonian path) gives the simplest possible matrix M: it's a scalar multiple of the identity matrix. Every eigenvalue is the same number (H/n). The spectrum is a single point — just one frequency, perfectly pure.

This is the tournament analog of a pure monochromatic light beam.

More complex tournaments have richer spectra with multiple eigenvalues — the more symmetry is broken, the more "colors" appear in the spectrum.

---

## What the Spectrum Reveals

The eigenvalue spectrum is a compressed description of a tournament. Different tournaments have different spectra, and the spectrum encodes:

1. **How many Hamiltonian paths the tournament has**: H(T) is related to the sum of all eigenvalues (the "trace").
2. **How symmetric the tournament is**: high symmetry = fewer distinct eigenvalues.
3. **How "spread out" the path structure is**: the ratio of largest to smallest eigenvalue (condition number) measures how evenly distributed the paths are across starting positions.

---

## The Band-Limited Structure

One particularly deep discovery: the path-count matrix M(T) is **band-limited** in a specific sense.

In signal processing, a band-limited signal is one that contains only frequencies up to some maximum. Radio stations broadcast at specific frequencies; a band-limited signal doesn't bleed into other stations.

Tournament matrices turn out to be band-limited in the space of "Walsh functions" — a particular set of mathematical wave patterns. The degree (highest frequency present) is at most 2⌊(n-1)/2⌋.

This means: out of the enormous number of possible matrix entries, the tournament matrices are constrained to a much smaller subspace. This compactness is why eigenvalue formulas like the cosine pattern above are possible — the structure is forced to be simpler than it looks.

---

## Connection to the Krawtchouk Framework

The band-limited structure connects tournaments to a coordinate system used in coding theory — the **Krawtchouk functions** — which measure how evenly spread a pattern is across all sizes of subsets.

In this coordinate system:
- The first two coordinates strongly predict H(T) (correlation ≈ -0.94)
- Self-complementary tournaments have a special symmetry: their odd-indexed coordinates are exactly zero

This connection says: tournament path-counting is, in a precise sense, a coding theory problem. The "frequency content" of a tournament (its Krawtchouk spectrum) determines how many Hamiltonian paths it has. This opened connections to the theory of error-correcting codes.

---

## Key Words

- **Matrix**: a rectangular grid of numbers that can "act" on lists of numbers
- **Eigenvalue**: a stretch factor — the matrix acts on a special vector by just scaling it
- **Spectral analysis**: studying the eigenvalues of a matrix associated to a mathematical object
- **Cosine formula**: the eigenvalues of rotationally symmetric tournament matrices follow a cosine pattern
- **Band-limited**: contains only "frequencies" up to some maximum, like a radio signal constrained to one station's range
- **Krawtchouk functions**: a coordinate system for signals on binary strings, borrowed from coding theory

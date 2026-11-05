# Independent review: constructive ratio-lattice cancellation

Status: **PASS**, no correction. The result is an exact decision and construction theorem for cancelling the binary mode in the stated ordered-join module language. It does not prescribe the ternary amplitude, constant offset, primality or general spectrum coverage.

## Frozen source

Under `../../session197-odd-cycle-path-polynomials/`:

- `RATIO_LATTICE.md`: SHA256 `a6f9c0b09d715865d01bbd62d95f2261c59e1f4a7eea2ed05cc09ac90aff60a3`.
- `verify_ratio_lattice.py`: SHA256 `3b40fd4e066f3b9d9eb3b6be44645a140035fd6df1d7adf9d9c29e757814a2e8`.
- `ratio_lattice_normal.stdout` and `ratio_lattice_optimized.stdout`: SHA256 `a5931452d19b5fd6df234d277fda12d9856c18a0ddd60276bd9454129d788814`.

The inherited arbitrary-module coefficient formula is `UNIT_OFFSET_RIGIDITY.md`, SHA256 `3c116086bfd697eaeaf5fcbf46308f1e31cdf653eb2a651053f01bce46ce374b`, whose verifier is `verify_unit_offset_rigidity.py`, SHA256 `756ea2c09f32d4830a05ebf92bee079812b3b013d9526b173a04ca2f70bbab0f`. That physical identity was already independently reviewed. The present review checks its specialization and the new inverse arithmetic theorem.

## 1. Physical reduction and nonempty realization

Labelled-bin paths in an ordered join split uniquely into the first block's path followed by the second block's path, allowing either part to be empty. Conversely their concatenation is valid. This proves full-bin multiplicativity and hence multiplicativity of rho=B(2)/H. The graph operation itself is correctly not called commutative.

Dividing the binary coefficient only by positive h_U h_V gives exactly the stated two-term equation. The both-zero, one-zero and same-sign branches are complete and valid. In the remaining opposite-sign branch kappa is positive and the ratio direction is rho_U/rho_V=kappa, with no inversion error.

Valuations of reduced positive rationals, including denominator valuations and primes appearing only in the target, give the exact integer system Mz=t. A solution z is equivalent to nonnegative multiplicities on the two sides by taking its positive and negative parts. If either side is empty, adding the same bank graph to both leaves the ratio unchanged and makes both modules nonempty. This correctly changes amplitude and offset; their preservation is not asserted.

The seven-vertex quotient has the literal Hamilton cycle 0,3,4,6,2,1,5,0. Thus substituting any nonempty modules preserves strongness. Positive leading coefficient and independence of bases three, two and one establish the eventual-to-exact statement.

## 2. Smith form and certificate directions

The source maintains D=E M F by recording row operations in E and column operations in F. Therefore t'=E t and z=F y are the correct conventions. In particular F, not its inverse, maps a diagonal-coordinate solution back to an exponent vector.

The Euclidean row/column steps preserve integrality and unimodularity. Whenever a nonzero remainder is swapped into the pivot position, the pivot's absolute value strictly decreases. Once its row and column are clear, a lower-block entry not divisible by the pivot is exposed by adding its row to the pivot row; the following Euclidean column step then decreases the pivot. This proves termination of each stage. Divisibility of the entire remaining block is preserved by later integer operations, giving the claimed divisibility chain.

The obstruction rows are exact certificates: a failed diagonal congruence supplies a functional divisible by the modulus on every generator but not on the target; a failed zero row annihilates all generators but not the target. Neither is merely a rational-rank heuristic.

Physical bank ratios are greater than one, so a nonempty physical bank has a nonempty prime-support union. The implementation's nonempty matrix shapes therefore cover this physical application, including target kappa=1 and rank-deficient matrices. Its generic rectangular tests include zero matrices.

## 3. Independent exact algorithm checks

`check_smith.py` compares the source routine against a separate exact Smith normal form implementation in SymPy over ZZ. It also independently recomputes E M F, both determinants, and every returned witness or obstruction directly on the original matrix and target.

The audit covers 2,250 systems:

- Every one of the 625 two-by-two integer matrices with entries from -2 through 2, each with two targets.
- 500 additional rectangular matrices of sizes one through six, each with a known integer target and an arbitrary target.

All 1,479 solvable and 771 obstructed results carry valid exact certificates, and all Smith invariant factors agree. No floating-point calculation is used. This external algebra package is an audit dependency only; the source algorithm itself uses the standard library.

## 4. Independent physical checks

`check_physical.cpp` independently constructs the bank, ordered joins, quotient and positive graph witnesses. Hamilton counts use inclusion-exclusion of fixed-length directed walks, rather than the source's endpoint subset DP.

It confirms the profiles of K1, Tr3 and the four-vertex H5 graph, checks multiplicativity on all nine ordered pairs from that bank, verifies the explicit quotient Hamilton cycle, and counts the actual sixteen- and seventeen-vertex graphs as 46,015 and 138,175. These agree with 15,360*3^m-65 for m=1,2. Signed 128-bit integers are ample: even the crude bound 2^17*17^17 is below 2^88.

The missing-prime example and full-support index obstruction are both correct. For the latter the functional v2+5v5 vanishes modulo three on both generator columns but equals four on the target. The positive solution (4,-1) realizes U=Tr4 and V=the H5 graph. The stated zero-sign physical examples and common-padding case are valid.

## 5. Replays and artifacts

Source normal and Python `-O` runs reproduce the frozen receipt. Independent Smith tests also agree normally and under `-O`. The physical checker agrees at `-O3` and `-O0 -D_GLIBCXX_ASSERTIONS`. No source file was changed.

- `check_smith.py`: `64a2b1c0831266d980d5d1d9f76a7347b935de04af4685989df10751d73e29da`.
- `smith-normal.json`: `55739057771b80b0808bedbfb68eeb9aea30a8a9c21a7da68b2e55f61555ef4c`.
- `check_physical.cpp`: `8f8a43dacc57822aefb58a59919aaaae643d89a120b2f0c1c41d88cbb1fd6a8d`.
- `physical.txt`: `24f3b79634e6508d3816cbcedd00bb4cf7094827e002ec3831b9761f662a0b00`.

The separate unbounded transitive-leaf corollary discussed during review is not part of this frozen source or signoff.

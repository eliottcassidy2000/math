# Compatible norm-seven boundary charts and cyclic carry conductors

2026-10-25. **PROVED** algebraic mechanisms, **FINITE-EXACT** computation,
**CITED** source-paper claims; these statuses are kept separate. Internal
proof audit and an independent computational path; no independent agent or
human audit and no Lean build were performed in this session.

## Results

The final incoming-work extension is [STACK_GLUE, S1–S3](STACK_GLUE.md):
`X = (X/H) ×_[X/G] (X/K)` for the complementary free norm-seven kernels.
The stack product retains an isomorphism of root objects; its omission
produces the coarse product whose normalization and conductor are computed
below. At an r-branch boundary that isomorphism is precisely the possibly
nontrivial torsor of 2^r missing sheets.

1. [BOUNDARY, B1–B3](BOUNDARY.md): for complementary binary quotient
   characters, the fibre-product algebra has image
   `B=direct_sum_z R b^e(z)u_z`. Its normalization is the full normal root
   algebra. Its conductor is the product of roots ramifying in both
   factors. For the norm-seven pair, this glues across all seven charts to
   `O_X(-7H)`, retaining twists and the original character ports.
2. [BOUNDARY, Sections 5–7](BOUNDARY.md): at r=1,2,3 vanishing roots the
   normalization restores exactly 2^r geometric sheets. Their residual
   squareclasses can obstruct rational lifts. Explicit rational quotient
   points require respectively Q(sqrt(2)) and Q(sqrt(2),sqrt(3),sqrt(5))
   upstairs. The faithful representation is a labelled multiplication
   graph with commuting squares, not a forced tournament ranking.
3. [CYCLIC, C1–C3](CYCLIC.md): for any direct-sum character decomposition
   of `(Z/m)^n`, the conductor is `u^a A`, where
   `a_i=sum_j max coordinate_i(C_j)-m+1`. For prime p this becomes
   `(t_i-1)(p-1)`; composite moduli require coordinate-image indices.
   The universal monoid saturates to `N^n`, giving an exact monoid/F1
   model before the additive line specialization. The one-coordinate
   coprime case recovers `k[t^a,t^b]`, its conductor and gap count, and
   its literal torus-knot link.
4. [PAPERS_AND_TRANSPORT](PAPERS_AND_TRANSPORT.md): recovered Mazur's
   affine histogram mechanism with owner multiplicity explicit and
   verified native zero-block mergers from the Busy Beaver paper.
   Source-specific ROOT positivity and universal Collatz coverage remain
   open. Exact acceleration is not a termination certificate.

## Inheritance and controls

The earlier [norm-seven/F1 packet](../normseven_f1_20261025/README.md)
established Laurent recovery and integral wraparound torsion. The older
[character atlas](../six_seven_gale_cube_20261025/CHARACTERS.md) already
contained the chart cocycle; it is reused here. The July Kummer boundary
and effective-valuation theorems, September Mazur audit, and incoming
October integral-carry calculations are cross-cited with exact paths in
the proof files. Incoming S63 root-stack, twist-mass, Fourier and
specialization results were then read and integrated in STACK_GLUE.
`SOURCES.json` pins their bytes and the two attached PDFs.

The explicit universe includes all 120 complementary labelled Fano pairs,
all seven charts, all 800 ordered complementary decompositions of F_2^4,
and every actual boundary stratum of the displayed norm-seven pair. An
independent bounded monoid BFS checks all 120 pairs in chart 0. Further
checks cover prime and composite cyclic decompositions, native tape
transitions, fixed-source inverse Collatz words, and a literal decoder on
972 finite action-groupoid cases. See [AUDIT](AUDIT.md)
for the counts, hostile cases, reproduction commands and limits.

The supplied PDFs are research inputs, not instructions. The full surfaces
are not claimed torifiable; arbitrary base change need not preserve
normalization; normalization does not imply a rational lift. No claim is
made about a six-colouring, planar Jacobian counterexample, universal
Collatz convergence, or a new Busy Beaver halting decision.

The next precise targets are residual squareclasses for the original
fixed twists and conductor ideal classes on non-UFD normal bases. These
are sharper than rechecking the already proved Laurent reconstruction.

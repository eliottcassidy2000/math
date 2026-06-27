#!/usr/bin/env python3
"""S212: expanded tournament-matrix atlas.

This script expands the S210 tournament matrix dictionary.  Each row is a
classic matrix result or object from some area of mathematics, followed by a
useful tournament-facing statement.  The point is not that every row is already
a theorem in the repo.  The point is to give future agents a large menu of
matrix proof carriers and quotient guardrails.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations


BASE_S210_ROWS = 135


DATA = """
linear spectral|Courant-Fischer min-max|Eigenvalues of a Hermitian matrix are variational extrema over subspaces.|Apply to iS for a tournament sign matrix: extremal modes are the strongest rank or cyclic directions, but the eigenvector support must be retained before using the mode as a proof carrier.
linear spectral|Rayleigh quotient|x*Mx / x*x gives a one-vector spectral certificate.|Use signed vectors on vertices, edges, or proof obligations to certify a cut, cycle, or residual direction in a tournament carrier.
linear spectral|Schur-Horn theorem|The diagonal of a Hermitian matrix is majorized by its spectrum.|Score or diagonal sidecars are only shadows of the full spectrum; any score-only tournament quotient needs a cyclic residual sidecar.
linear spectral|Lidskii inequalities|Eigenvalues of sums obey majorization constraints.|Split a tournament operator into ranking, cyclic, and packet-error blocks; impossible spectral majorizations become certificate failures.
linear spectral|Ky Fan extremal principle|Sums of top eigenvalues optimize traces over projections.|Choose a low-dimensional proof-carrier subspace by maximizing retained cyclic or owner payload, not by raw variance alone.
linear spectral|Davis-Kahan sin-theta theorem|Spectral subspaces are stable under perturbations with a gap.|An edge-flip or tail-add update changes the tournament spectral carrier predictably when a gap exists; no gap means the quotient is fragile.
linear spectral|Bauer-Fike theorem|Eigenvalues of a diagonalizable matrix move within condition-scaled perturbation disks.|A poorly conditioned tournament operator is a bad proof carrier: tiny hidden coordinates can move route labels.
linear spectral|Pseudospectrum|Non-normal matrices can have large spectral instability.|Adjacency A is non-normal; use pseudospectra to detect tournaments where eigenvalues understate transient cyclic pressure.
linear spectral|Spectral projectors|Contour integrals isolate invariant spectral subspaces.|Project tournament data onto ranking, cyclic, and exceptional modes while keeping the projector as the retained sidecar.
linear spectral|von Neumann trace inequality|Trace pairings are bounded by singular-value alignments.|Pair a tournament certificate matrix with a residual matrix; equality cases identify aligned hidden obstruction directions.
linear spectral|Loewner order|PSD order compares quadratic forms.|A matrix tournament certificate should state which quadratic tests it dominates; scalar domination is not enough.
linear spectral|Matrix monotone functions|Functions preserving Loewner order have strong structure.|If a tournament invariant is obtained by applying f(A), prove f preserves the order relevant to the LRC predicate.
linear spectral|Hadamard determinant inequality|The determinant is bounded by row norm product.|Extremal tournament sign matrices near conference/Hadamard objects may maximize volume while not maximizing Hamilton paths or LRC usefulness.
linear spectral|Fischer determinant inequality|Principal minors of PSD matrices satisfy log-submodularity.|Subtournament determinant data may have submodular structure; failures mark non-PSD or skew residual debt.
linear spectral|Sylvester rank inequality|rank(AB) is bounded below by rank(A)+rank(B)-n.|Composition of two tournament quotient maps loses at most a quantifiable dimension only when the matrix ranks are known.
linear spectral|Frobenius inner product|tr(A*B) is the Euclidean pairing on matrices.|Turn pairwise tournament comparisons into an inner-product geometry of proof carriers; orthogonal residuals are invisible to the chosen invariant.
linear spectral|Matrix sign function|The sign of a matrix separates positive and negative spectral subspaces.|For Hermitian iS, sign(iS) is a coarse cyclic phase carrier; test whether it preserves route/status labels.
linear spectral|Polarization identity|Quadratic forms determine bilinear forms.|If a tournament certificate is quadratic, recover the cross terms before quotienting, because those are often the hidden coupling payload.
matrix scaling positivity|Birkhoff-von Neumann theorem|A doubly stochastic matrix is a convex combination of permutation matrices.|A softened tournament adjacency can be decomposed into orderings; residual mass outside near-permutation orderings measures cyclic ambiguity.
matrix scaling positivity|Sinkhorn-Knopp scaling|Positive matrices can often be diagonally scaled to doubly stochastic form.|Scale weighted tournament comparison matrices to remove degree bias, then study the remaining cyclic residual.
matrix scaling positivity|Perron complement|Eliminating states in a nonnegative matrix gives a corrected Perron block.|Deleting tournament vertices or proof states should use Perron/Schur complements, not raw deletion, if ranking pressure matters.
matrix scaling positivity|M-matrix theory|Matrices with nonpositive off-diagonal and positive inverse encode monotone systems.|A Laplacian-like tournament carrier with M-matrix structure gives monotone repair certificates; failure means feedback/cycle debt.
matrix scaling positivity|P-matrix theorem|All principal minors positive has implications for complementarity uniqueness.|If a tournament packet LCP has a P-matrix, route selection is unique; otherwise multiple residual routes may coexist.
matrix scaling positivity|Total positivity|All minors positive force variation-diminishing behavior.|A totally positive tournament feature matrix would forbid oscillating route changes; non-total-positive minors locate cycle witnesses.
matrix scaling positivity|Oscillatory matrices|Totally nonnegative matrices with positivity in powers have ordered real spectra.|If a tournament quotient becomes oscillatory, it is essentially order-like; cyclic information has probably been destroyed.
matrix scaling positivity|Stieltjes matrices|Symmetric positive definite M-matrices have inverse positivity.|Use for energy forms on tournament cut or proof-obligation graphs where inverse entries should be influence weights.
matrix scaling positivity|Hilbert projective metric contraction|Positive maps contract projective distance.|Repeated stochastic tournament updates may converge to a ranking pressure vector; contraction says nothing about exact cycle sidecars.
factorizations canonical forms|Rational canonical form|A matrix over a field decomposes into companion blocks.|Tournament recurrences over finite fields, especially mod p sidecars, should report invariant factors instead of only eigenvalues.
factorizations canonical forms|Smith-McCoy for polynomial matrices|Polynomial matrix equivalence exposes invariant factors.|Dynamic tournament transfer matrices can have hidden periodic factors that eigenvalue samples miss.
factorizations canonical forms|Kronecker canonical form of pencils|Matrix pencils classify generalized eigenvalue and singular structure.|Use pencils lambda B - C for two competing tournament carriers; singular blocks are hidden coordinates neither carrier sees alone.
factorizations canonical forms|QZ decomposition|A pair of matrices is triangularized for generalized eigenvalues.|Compare adjacency and Laplacian, or score and cycle operators, without forcing them to commute.
factorizations canonical forms|Autonne-Takagi factorization|Complex symmetric matrices diagonalize by unitary congruence.|Symmetric lifts of tournament signs can expose phase-aligned cycle packets distinct from Hermitian iS spectra.
factorizations canonical forms|Williamson normal form|Positive symplectic matrices have canonical symplectic eigenvalues.|If a tournament carrier is built from paired edge/cycle coordinates, symplectic eigenvalues separate canonical coupled modes.
factorizations canonical forms|LDL transpose factorization|Symmetric indefinite matrices factor into pivots and signs.|Boundary or conflict forms for tournaments may be indefinite; pivot signs identify obstruction directions.
factorizations canonical forms|Rank-revealing QR|QR with pivoting exposes numerical rank and important columns.|Pick sidecar columns, such as owner, sector, primitive deck, by pivoted rank rather than intuition.
factorizations canonical forms|CUR decomposition|A matrix can be approximated by selected rows and columns.|Choose representative tournament vertices, edges, or proof fibers that explain most matrix signal while preserving exact residual checks.
factorizations canonical forms|Nonnegative matrix factorization|A nonnegative matrix decomposes into additive latent parts.|Decompose walk or danger-count matrices into latent route families; nonunique factors must not be mistaken for proof labels.
factorizations canonical forms|Independent component analysis|A matrix factorization seeks statistically independent sources.|Feature banks over tournaments may split into ranking, cycle, endpoint-owner, and arithmetic-clock components.
factorizations canonical forms|Butterfly factorization|Fourier-like matrices admit fast multiscale factors.|Circulant and residue tournaments should exploit butterfly/FFT structure before brute-force enumeration.
graph combinatorics|Lindstrom-Gessel-Viennot lemma|Path families are determinants of path-count matrices.|Nonintersecting path packets in a tournament carrier can become determinant certificates with visible path addresses.
graph combinatorics|All-minors matrix-tree theorem|Laplacian minors count directed forests connecting boundary sets.|Use minors to count redundancy of proof routes between source, obstruction, and certificate vertices.
graph combinatorics|Matrix forest theorem|The inverse of I plus Laplacian counts rooted forests.|Forest weights in carrier tournaments measure how strongly residual obligations attach to certificates.
graph combinatorics|Markov chain tree theorem|Stationary probabilities equal weighted arborescence sums.|Random-walk ranking in a tournament has a tree expansion; high rank pressure should be backed by arborescence sidecars.
graph combinatorics|Critical group and sandpile matrix|Laplacian cokernels define finite abelian groups.|Smith normal form of tournament carrier Laplacians exposes p-adic obstruction clocks.
graph combinatorics|BEST theorem via matrices|Euler tours in directed graphs are counted by arborescences and degrees.|For edge-state tournament carriers, Eulerian route counts factor into local degree data plus tree redundancy.
graph combinatorics|De Bruijn transfer matrices|Words and overlaps are counted by adjacency matrices of state graphs.|Automatic tournament shadows should keep state-transition matrices, not just accepted word labels.
graph combinatorics|Goulden-Jackson cluster matrix|Forbidden-word counts are rational functions from overlap matrices.|Forbidden tournament patterns can be counted by overlap matrices; overlaps are the hidden address coordinates.
graph combinatorics|Line digraph adjacency|Edges become vertices and adjacency is incidence-compatible.|When node perspectives lose information, pass to directed-edge matrices; this is exactly the edge-perspective lift.
graph combinatorics|Godsil-McKay switching matrices|Certain switching operations preserve spectra.|Cospectral tournaments or conflict graphs need switching sidecars; spectrum alone is not a complete proof state.
graph combinatorics|Spectral excess theorem|Distance-regularity is characterized by spectral excess.|A tournament carrier that looks distance-regular spectrally should be tested for exact distance or proof-depth regularity.
graph combinatorics|Walk-regularity|Closed walk counts are independent of base vertex.|If a tournament is walk-regular but not vertex-transitive, trace data has forgotten stabilizer or orbit information.
graph combinatorics|Incidence matrix of a hypergraph|B B^T records co-incidence of blocks.|Cycle, edge, or owner-block incidence matrices can replace raw tournament vertices in LRC proof carriers.
graph combinatorics|Krawtchouk matrices|Hamming-scheme eigenmatrices diagonalize subset convolution.|Boolean edge-flip layers of tournaments can be analyzed in Krawtchouk coordinates; low layers are not automatically safe quotients.
finite fields coding|Generator matrices of codes|Rows span a linear code.|Treat selected tournament invariants as codewords; low distance means hidden coordinates are vulnerable to collisions.
finite fields coding|Parity-check matrices|Syndromes detect errors.|A tournament quotient needs a parity-check sidecar that detects route/status-changing hidden errors.
finite fields coding|MacWilliams transform|Weight enumerators of a code and its dual are related by Krawtchouk matrices.|Edge-flip weight distributions and dual obstruction checks should be paired, not used separately.
finite fields coding|Delsarte linear programming|Association-scheme eigenmatrices give bounds on codes.|Tournaments in symmetric schemes, especially Paley/cyclotomic families, may admit LP bounds on forbidden cycle packets.
finite fields coding|Reed-Solomon Vandermonde generator|Evaluations of low-degree polynomials form MDS codes.|Clock-sampling sidecars can be MDS-like: a few exact phases may reconstruct the hidden polynomial packet.
finite fields coding|BCH parity matrices|Consecutive roots enforce distance.|For cyclic tournament residue words, root constraints can forbid small hidden defects.
finite fields coding|LDPC Tanner matrix|Sparse checks define iterative decoding.|Large tournament proof ledgers should use sparse sidecar checks so local residual errors can be decoded.
finite fields coding|Expander code matrices|Graph expansion boosts distance.|Carrier tournaments with expansion in their check matrix can tolerate local forgotten coordinates.
finite fields coding|Hadamard code matrix|Rows of a Hadamard matrix encode all linear characters.|Conference and Paley tournament signs live near Hadamard codes; character rows should be exact sidecars.
finite fields coding|Rank-metric codes|Distance is matrix rank.|Compare tournaments by rank distance of sign matrices to detect structural jumps not visible in Hamming edge distance.
number theory arithmetic|SL2 continued-fraction matrices|Continued fractions are products of 2 by 2 integer matrices.|Farey and denominator clocks in LRC are matrix products; tournament state should keep the product, not only the final denominator.
number theory arithmetic|Stern-Brocot matrices|Left/right products enumerate rationals.|Tournament threshold walks over Farey cells can be represented by L/R matrix words with exact ancestry.
number theory arithmetic|Modular group generators S and T|PSL2Z acts by integer matrices on the upper half-plane.|Structured tournament transformations on rational phases may be modular actions with hidden cusp data.
number theory arithmetic|Brandt matrices|Quaternion ideal classes are connected by Hecke correspondences.|Tournament isomorphism classes can be moved by prime correspondences; Brandt-like matrices record class-to-class transfer.
number theory arithmetic|Hecke operator matrices|Arithmetic correspondences act on modular forms.|Prime-denominator tournament packets should be tested as Hecke-like transitions on class or sidecar spaces.
number theory arithmetic|GCD matrices|Entries gcd(i,j)^alpha have determinant factorizations.|Divisor-sensitive tournament clocks should use gcd matrices to retain exact-period intersections.
number theory arithmetic|LCM matrices|Entries lcm(i,j) encode divisibility joins.|A quotient that stores only lcm tails may lose prime-power coordinates; Smith form can expose them.
number theory arithmetic|Divisor zeta matrix|The poset matrix 1_{d divides n} has Mobius inverse.|Exact-period packet ledgers are divisor-poset matrices; Mobius inversion should be explicit before scalarizing.
number theory arithmetic|Redheffer matrix|A determinant identity is equivalent to a Mertens sum.|Divisibility/tournament incidence matrices can hide number-theoretic sums; determinant signs may flag hard arithmetic clocks.
number theory arithmetic|Dedekind sum matrices|Cotangent and modular reciprocity produce structured finite matrices.|Two-large residue tails in LRC can be organized as Dedekind/cotangent channel matrices instead of absolute sums.
number theory arithmetic|Kloosterman sum matrices|Exponential sums over inverses form trace kernels.|Reciprocal tournament clocks should retain inverse-residue matrices, especially when raw residue counts collide.
number theory arithmetic|Character circulants|Multiplicative characters define diagonalizable circulants.|Quadratic/cyclotomic tournament orientations are exactly character-circulant sign matrices.
number theory arithmetic|Period matrices|Integrals over cycles produce matrices of arithmetic geometry.|When tournament packets are boundary cycles, period-like matrices can record which cycles are genuinely independent.
algebra representation|Macaulay matrix|Polynomial systems are solved by large coefficient matrices.|Tournament feasibility constraints can be translated into polynomial systems; Macaulay rank detects hidden impossibility.
algebra representation|Bezoutian|Common roots and stability are encoded by a bilinear matrix.|Two tournament clocks collide if their Bezoutian degenerates; useful for scalar quotient failures.
algebra representation|Multiplication matrices in quotient rings|Algebraic solutions are eigenvalues of multiplication operators.|Finite tournament packet equations can be solved by multiplication matrices, with eigenvectors retaining labels.
algebra representation|Jacobian matrix|Local rank gives dimension and singularity.|A family of tournament deformations is locally rigid when its Jacobian over sidecar constraints has full rank.
algebra representation|Hessian matrix|Second derivatives classify critical points.|Extremal tournament invariants, such as Hamilton paths or energy, need Hessian directions to detect flat hidden coordinates.
algebra representation|McKay adjacency matrix|Tensoring with a representation gives a graph matrix.|Automorphism groups of symmetric tournaments act on edge/cycle modules; McKay matrices classify transfer between irreducible packets.
algebra representation|Fusion matrices|Products in a fusion ring are matrices with nonnegative integer entries.|Combining tournament proof carriers may obey fusion rules; noncommuting fusion means route order matters.
algebra representation|Modular S and T matrices|A modular tensor category is controlled by two matrices.|If tournament packet sectors have anyonic/modular behavior, S/T matrices separate charge sectors better than scalar counts.
algebra representation|Auslander-Reiten matrix|Representation categories have translation quivers.|Tournaments as complete quivers can be studied by module mutation; AR matrices may classify persistent cycle modules.
algebra representation|Cartan determinant theorem|Finite type algebras have constrained Cartan matrices.|A tournament mutation category with forbidden determinant values could certify impossible obstruction packets.
algebra representation|Exchange matrix mutation|Cluster algebras mutate skew-symmetrizable matrices.|Tournament sign matrices are skew inputs; mutation gives a controlled deformation grammar for edge and proof-coordinate flips.
algebra representation|Incidence algebra zeta matrix|Poset convolution is matrix multiplication.|Transitive subtournaments are total-order incidence algebras; nontransitive debt is exactly failure of a poset zeta model.
topology geometry|Cellular boundary matrix|Homology is kernel of boundary modulo image.|Tournament path, arc, or conflict complexes should expose obstruction cycles as nullspace classes.
topology geometry|Coboundary matrix|Cohomology classes are quotient cocycles.|A hidden coordinate can be a cocycle; if it is not a coboundary, no scalar quotient can kill it safely.
topology geometry|Cup product matrix|Cohomology multiplication gives bilinear structure.|Cycle interactions in tournament complexes need product data, not only Betti numbers.
topology geometry|Intersection form|Manifold middle homology has a bilinear matrix.|Pairing tournament cycles by intersection or overlap may identify owner-essential cycles.
topology geometry|Seifert matrix|A knot's Seifert surface gives an asymmetric form.|Directed cycle carriers are naturally asymmetric; Seifert-like matrices may distinguish orientation-reversal debt.
topology geometry|Alexander matrix|Fox derivatives present knot modules.|Tournament mutation relations can be differentiated into presentation matrices whose minors are invariants.
topology geometry|Goeritz matrix|Checkerboard surfaces produce link signature matrices.|Planar conflict diagrams from tournaments can use Goeritz forms to certify parity or signature obstructions.
topology geometry|Linking matrix|Pairwise linking numbers assemble into a symmetric matrix.|Cycle families in a tournament may have linking/overlap matrices; off-diagonal entries are interaction debt.
topology geometry|Dirac matrix|A first-order operator squares to a Laplacian.|Path complexes of tournaments may be better studied with Dirac operators than with Laplacians alone.
topology geometry|Hodge star matrix|Metric duality maps k-forms to complementary forms.|LRC endpoint arcs need a metric sidecar; the Hodge star records weights lost by pure incidence.
topology geometry|Persistent reduction matrix|Column operations pair births and deaths.|As thresholds move, tournament or danger-arc complexes emit barcodes; unpaired bars are residual proof debt.
topology geometry|Cayley-Menger determinant|Distances determine simplex volume.|Geometric realization of a tournament sign pattern is constrained by distance determinant feasibility.
topology geometry|Rigidity matrix|Edge constraints have infinitesimal motion kernels.|If a geometric tournament is rigid, only named wall crossings can change its orientation.
topology geometry|Stress matrix|Self-stresses certify rigidity or tensegrity.|Dual stress matrices can prove no continuous deformation removes a tournament obstruction.
analysis operators|Galerkin stiffness matrix|Weak PDE forms become sparse PSD systems.|LRC safe-set energy can be discretized on arcs; zero modes are boundary equality atoms.
analysis operators|Mass matrix|Finite elements encode inner products by a mass matrix.|A tournament carrier with weighted cells must keep the mass matrix or risk false orthogonality.
analysis operators|Finite-difference differentiation matrix|Derivatives on grids become banded matrices.|Wall-crossing sequences of tournaments can be differentiated to find first active hidden coordinates.
analysis operators|Chebyshev differentiation matrix|Spectral collocation gives dense high-accuracy derivative matrices.|Smooth parameterized tournament families can use spectral derivatives to locate exact threshold changes.
analysis operators|Fourier matrix|DFT diagonalizes cyclic convolution.|Residue tournaments should be Fourier-diagonalized before using any scalar residue statistic.
analysis operators|Wavelet transform matrix|Multiscale basis matrices localize time and frequency.|Tournament evolution over thresholds may have local jumps and global trends; wavelets separate them.
analysis operators|Convolution matrix|Toeplitz/circulant matrices implement convolution.|Pairwise danger patterns can be convolved exactly when translation symmetry exists; otherwise convolution is a lossy approximation.
analysis operators|Nyström matrix|Integral operators are approximated by sampled kernels.|Infinite or graphon tournaments can be sampled, but final residuals need exact finite sidecars.
analysis operators|Ulam matrix|A dynamical system transfer operator is discretized by boxes.|Threshold dynamics of tournament states can be approximated by a Markov matrix over quotient boxes.
analysis operators|Floquet monodromy matrix|Periodic linear systems are classified by one-period transfer.|Periodic tournament or automaton clocks should report monodromy; fixed points and multipliers are hidden clock data.
analysis operators|Carleman linearization matrix|Nonlinear dynamics embed into infinite linear systems.|Nonlinear tournament updates can be lifted to linear but infinite moment coordinates; finite truncations need residual bounds.
analysis operators|Resolvent identity|Inverse shifts satisfy a precise difference formula.|Compare tournament operators before/after an edge flip using resolvent updates; poles mark fragile proof carriers.
optimization algorithms|KKT matrix|Optimality conditions combine Hessian and constraints.|Tournament ranking, sidecar fitting, or LRC packet certification should emit the KKT residual, not just the optimum.
optimization algorithms|Barrier Hessian|Interior-point methods use curvature of log barriers.|Polyhedral tournament quotients have barrier curvature; flat directions are hidden coordinates.
optimization algorithms|Linear complementarity matrix|LCP uniqueness and algorithms depend on matrix classes.|Route selection among tournament proof obligations can be phrased as complementarity; P-matrix sidecars give uniqueness.
optimization algorithms|Monge array|Quadrangle inequalities enable fast optimization.|If tournament cost matrices are Monge after an ordering, the ordering is a real structural certificate.
optimization algorithms|Totally balanced matrices|Certain set-cover matrices have integral polyhedra.|Obstruction-cover matrices with total balance can be solved without fractional residual ambiguity.
optimization algorithms|Spectral sparsification|A sparse PSD matrix approximates a Laplacian.|Sparse proof-carrier matrices can approximate bulk behavior, but exceptions must list omitted edge/cycle sidecars.
optimization algorithms|Leverage scores|Diagonal of projection matrices measures row influence.|Rank tournament sidecar rows by leverage to find which hidden coordinates control route purity.
optimization algorithms|Matrix multiplicative weights|PSD weights are updated to satisfy constraints.|Iteratively search for tournament dual certificates while maintaining a PSD witness matrix.
optimization algorithms|Sum-of-squares moment matrix|Polynomial nonnegativity is certified by PSD moment matrices.|Forbidden tournament/LRC packets may have SOS certificates over edge or endpoint variables.
optimization algorithms|Lasserre hierarchy|Moment matrices give increasingly tight polynomial relaxations.|Use low levels as scouts, but record which tournament labels are lost at each relaxation.
optimization algorithms|Copositive programming|Hard combinatorial problems can be exact over copositive cones.|Feedback arc or clique-style tournament problems have copositive exact formulations; SDP is only a relaxation.
optimization algorithms|Interior-point normal equations|Newton steps solve structured linear systems.|The numerical conditioning of a tournament proof LP tells which sidecars are nearly dependent.
probability random|Absorbing Markov fundamental matrix|N equals (I-Q)^{-1} and counts expected visits.|In proof-state tournaments, expected visits before certificate absorption measure residual difficulty.
probability random|Detailed balance matrix equation|pi_i P_ij = pi_j P_ji characterizes reversibility.|A tournament random walk is generally irreversible; reversibilization destroys cyclic current and must be tracked.
probability random|Additive reversiblization|P plus time reversal creates a reversible chain.|Use reversible spectra as a prefilter, but keep antisymmetric current for Condorcet and LRC debt.
probability random|Fisher information matrix|Curvature of log likelihood controls estimability.|In noisy tournament models, parameters with low Fisher information are hidden coordinates; do not quotient them.
probability random|Cramer-Rao bound|Estimator variance is bounded by inverse Fisher information.|No tournament statistic can reliably recover a sidecar below this information threshold.
probability random|Marchenko-Pastur law|Sample covariance eigenvalues have a universal bulk.|Feature matrices of random tournaments should match MP baselines; outliers signal structure.
probability random|Circular law|Eigenvalues of random non-Hermitian matrices fill a disk.|Adjacency A of random tournaments has non-normal random-matrix behavior; compare before reading eigenvalue clouds as structure.
probability random|Girko hermitization|Non-Hermitian spectra are studied through Hermitian block matrices.|Tournament adjacency spectra should be hermitized when using rigorous random-matrix tools.
probability random|Free convolution|Large random matrix sums have additive free laws.|Tournament core plus random tail may be modeled by free convolution, with deviations as structured sidecars.
probability random|Dyson Brownian motion|Eigenvalues evolve under matrix noise.|Random edge flips induce spectral diffusion; barriers or stuck modes indicate rigid tournament structure.
probability random|BBP phase transition|Spiked random matrices reveal signals past a threshold.|A latent ranking spike in a random tournament is detectable only past a signal threshold; below it, cycle sidecars dominate.
probability random|Nonbacktracking spectral detection|Community signals can emerge in nonbacktracking spectra.|Edge-state tournament matrices may reveal block or proof-route structure invisible in adjacency spectra.
cs machine learning|Spectral clustering matrix|Eigenvectors of Laplacians cluster data.|Cluster tournament proof states by Laplacian modes, then verify that each cluster is route-pure.
cs machine learning|HITS authority-hub matrices|A^T A and A A^T rank hubs and authorities.|Tournament winners and validators can differ; compare hub/authority modes with cyclic residuals.
cs machine learning|Graph neural message matrix|Message passing repeatedly applies normalized adjacency.|A GNN on tournaments is a learned WL/matrix quotient; inspect failures for hidden sidecars.
cs machine learning|Attention matrix|Rows of softmax scores mix tokens.|A tournament can seed attention by pairwise dominance, but softmax loses antisymmetric exactness unless S is retained.
cs machine learning|Neural tangent kernel matrix|Training dynamics linearize through a kernel.|Learned tournament invariants can be audited by NTK rank: low rank means likely hidden-coordinate collisions.
cs machine learning|Kernel ridge inverse|Predictions use (K+lambda I)^{-1}.|A tournament kernel predictor should expose influence weights; high influence residuals deserve exact proof checks.
cs machine learning|Graph transformer bias matrix|Structural biases are added to attention logits.|Encode tournament edge, cycle, and owner sidecars as separate bias matrices rather than one scalar feature.
cs machine learning|Matrix completion with one-bit observations|Signs of low-rank matrices can be recovered under assumptions.|Near-transitive tournaments are one-bit low-rank data; cyclic tournaments violate the model in the residual.
cs machine learning|Permutation equivariant matrices|Equivariant linear maps have constrained block forms.|Tournament neural features must respect isomorphism by simultaneous permutation conjugacy.
games social choice|Payoff matrix minimax theorem|Zero-sum games have equal max-min and min-max values.|A tournament sign matrix is an antisymmetric game; mixed equilibria identify balanced cyclic cores.
games social choice|Replicator dynamics matrix form|Population shares evolve by payoff differences.|Dominance dynamics on a tournament can cycle; equilibria expose nontransitive support.
games social choice|Condorcet majority matrix|Pairwise majority margins form a skew matrix.|Social choice tournaments should keep margin magnitude, not only win orientation, when proving stability.
games social choice|Borda score vector|Row sums of the majority matrix rank alternatives.|Borda is the score vector of A or S; it is a first-order shadow and needs cycle sidecars.
games social choice|McGarvey realization matrix idea|Any tournament can be realized as a majority relation.|Matrix proof: majority tournaments are universal, so social-choice realizability alone gives no structural proof unless voter profile matrices are retained.
games social choice|Maximal lottery equilibrium|A skew payoff matrix defines a mixed equilibrium.|Balanced support of a tournament game can be a canonical cyclic perspective, complementary to rooted node views.
games social choice|Pairwise comparison consistency index|Reciprocal matrices measure inconsistency by eigenvalues.|Tournament inconsistency is not noise by default; high inconsistency may be the central cyclic object.
games social choice|Elo update matrix|Pairwise outcomes update latent ratings.|Elo extracts a ranking pressure channel and discards exact cycle witnesses unless match history is retained.
physics control|Hamiltonian matrix|Symplectic systems have paired eigenvalues.|Edge-cycle coordinates may form Hamiltonian pairs; paired eigenvalues suggest conserved cyclic quantities.
physics control|Lindblad superoperator matrix|Quantum Markov evolution acts on density matrices.|Uncertainty over tournament classes can evolve by measurement; decoherence corresponds to forgetting labels.
physics control|Scattering S-matrix unitarity|Incoming and outgoing amplitudes preserve probability.|A proof quotient should conserve obligation mass: every incoming residual exits as certificate, descent, annihilation, or named debt.
physics control|Green function resolvent|Poles of (zI-M)^{-1} identify modes.|Tournament resolvent poles reveal persistent cycles or bottlenecks in proof-state dynamics.
physics control|Controllability Gramian|Reachability energy is a PSD matrix.|Edge flips or sidecar additions with low controllability energy are efficient proof moves.
physics control|Observability Gramian|Output energy detects hidden states.|A scalar invariant set is proof-safe only if its observability Gramian detects route/status-changing coordinates.
physics control|Algebraic Riccati equation|Optimal feedback solves a matrix quadratic equation.|Adaptive tournament querying should balance comparison cost against reduction in cyclic uncertainty.
physics control|Transfer matrix method for Ising models|Partition functions are traces/products of local matrices.|Layered tournament families can be counted by transfer matrices with weights for cycles, owners, or endpoints.
physics control|Kasteleyn orientation matrix|Planar dimer counts are Pfaffians.|Planar subcarriers inside tournaments may have Pfaffian certificates even when the whole tournament is nonplanar.
physics control|Dirac gamma matrices|Clifford generators anticommute.|Small parity-rich tournament signs can be represented by Clifford generators; anticommutation tracks cycle parity.
"""


DOMAIN_SCORES = {
    "linear spectral": (7, 5, 8, 7, 4),
    "matrix scaling positivity": (8, 5, 7, 7, 5),
    "factorizations canonical forms": (7, 7, 7, 6, 5),
    "graph combinatorics": (9, 9, 8, 9, 6),
    "finite fields coding": (8, 6, 9, 8, 6),
    "number theory arithmetic": (8, 6, 10, 9, 6),
    "algebra representation": (7, 7, 8, 7, 5),
    "topology geometry": (9, 8, 8, 10, 6),
    "analysis operators": (7, 6, 7, 7, 5),
    "optimization algorithms": (9, 7, 8, 9, 7),
    "probability random": (6, 5, 6, 6, 5),
    "cs machine learning": (6, 5, 6, 5, 7),
    "games social choice": (8, 8, 6, 7, 6),
    "physics control": (8, 7, 7, 8, 6),
}


def parse_rows() -> list[dict[str, str]]:
    rows = []
    for line in DATA.strip().splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) != 4:
            raise ValueError(f"bad row: {line!r}")
        domain, name, matrix, tournament = parts
        rows.append(
            {
                "domain": domain,
                "name": name,
                "matrix": matrix,
                "tournament": tournament,
            }
        )
    return rows


def directed_three_cycles(adj: list[list[int]]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def ham_paths_adj(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if seen & (1 << u):
                    continue
                if adj[v][u]:
                    dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def domain_tournament(rows: list[dict[str, str]]) -> tuple[list[str], list[list[int]]]:
    domains = sorted({row["domain"] for row in rows})
    counts = Counter(row["domain"] for row in rows)

    def score(domain: str) -> tuple[int, int, int, int, int, int, str]:
        exact, incident, arithmetic, lrc, computable = DOMAIN_SCORES[domain]
        return (
            lrc + exact + incident + arithmetic,
            lrc,
            exact,
            incident,
            arithmetic,
            computable + counts[domain],
            domain,
        )

    n = len(domains)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if score(domains[i]) > score(domains[j]):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return domains, adj


def main() -> None:
    rows = parse_rows()
    counts = Counter(row["domain"] for row in rows)
    domains, adj = domain_tournament(rows)
    out_scores = [sum(row) for row in adj]
    order = sorted(range(len(domains)), key=lambda idx: -out_scores[idx])

    print("# Expanded Tournament-Matrix Atlas - codex-S212")
    print()
    print("This expands the S210 tournament matrix translation layer.")
    print(f"New rows in this expansion: {len(rows)}")
    print(f"Combined named hooks with S210: {BASE_S210_ROWS + len(rows)}")
    print(f"Domains covered: {len(counts)}")
    print()
    print("## Core Readout")
    print()
    print(
        "A tournament is not just an adjacency matrix.  It can be read as a "
        "skew sign matrix, Hermitian operator iS, Laplacian, Markov kernel, "
        "boundary matrix, incidence matrix, transfer matrix, comparison game, "
        "or sidecar observability system.  Each matrix theorem below suggests "
        "a tournament proof carrier, but every carrier must name what it "
        "forgets."
    )
    print()
    print("## Domain Counts")
    print()
    for domain, count in sorted(counts.items()):
        print(f"- {domain}: {count}")
    print()
    print("## Expanded Atlas")
    print()
    current = None
    number = 1
    for row in rows:
        if row["domain"] != current:
            current = row["domain"]
            print(f"### {current.title()}")
            print()
        print(f"{number}. **{row['name']}.**")
        print(f"   Matrix statement: {row['matrix']}")
        print(f"   Tournament statement: {row['tournament']}")
        print()
        number += 1

    print("## Tournament Analysis Over Domains")
    print()
    print("Vertices are matrix-result domains, not tournament runners.")
    print(
        "Pairwise observable: retained exactness, incident/cycle payload, "
        "arithmetic hidden-clock payload, LRC sidecar usefulness, and "
        "computability."
    )
    print(
        "Switch: orient toward the domain that keeps more proof-relevant "
        "sidecar information before scalarizing."
    )
    print(f"vertices={domains}")
    print(f"score_hist={dict(sorted(Counter(out_scores).items()))}")
    print(f"directed_3_cycles={directed_three_cycles(adj)}")
    print(f"hamiltonian_paths={ham_paths_adj(adj)}")
    print("one_hamiltonian_path=" + " -> ".join(domains[i] for i in order))
    print()
    print("## Next Formal Pulls")
    print()
    pulls = [
        "Build a sidecar observability matrix whose rows are residual packet pairs and whose columns are candidate hidden coordinates.",
        "For HYP-3047 edge perspectives, emit four-sector block matrices, skew-cycle traces, and Schur-complement deletion corrections.",
        "Run Smith normal form on integer boundary/incidence sidecars when real or complex spectra merge packets.",
        "Treat each low-rank edge flip as a rank-2 skew update and compute determinant/resolvent/Pfaffian sensitivity.",
        "Use matrix-tree and forest counts on proof-carrier tournaments to measure redundancy of certificate routes.",
        "Use KKT/Farkas/SOS matrices to turn impossible residual packets into dual certificates.",
    ]
    for i, pull in enumerate(pulls, 1):
        print(f"{i}. {pull}")
    print()
    print("## Assumption Challenge")
    print()
    print(
        "The rows and columns need not be runners.  They may be directed edges, "
        "cycles, gaps, wall crossings, endpoint owners, denominator clocks, "
        "proof obligations, quotient fibers, cohomology classes, sidecar fields, "
        "or low-rank update directions.  The preserved predicate is whether "
        "the matrix carrier keeps enough information to decide the tournament "
        "or LRC proof obligation.  The destroyed information must be retained, "
        "reconstructed, dual-annihilated, descended, or routed to named debt."
    )


if __name__ == "__main__":
    main()

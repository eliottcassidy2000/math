# Session 28: recovered prime-density, Hessian, and representation threads

Status: PROVED scoped deductions / FINITE-EXACT computations / historical recovery. No RH proof, no general lattice existence result, no claim that three fractions uniquely specify a sequence.

## 1. The three fractions have a precise natural interpretation

Let pi(x) count primes <=x. The proportion of k-digit base-b integers that are prime is

    d_b(k) = [pi(b^k-1)-pi(b^(k-1)-1)] / [(b-1)b^(k-1)].

For b=10 this yields exactly 4/9,21/90,143/900,1061/9000,8363/90000,68906/900000,586081/9000000. A complete independent Eratosthenes sieve below 10^7 verifies these numerators. An exact-string web search identifies the same interpretation in rxiv:2311.0050v1, but no other claim from that manuscript is used.

Primary table: Andrew Granville, *A good new millenium for the primes*, Table 1, printed p.2:
https://dms.umontreal.ca/~andrew/PDF/PrimesICM1.pdf
It gives pi(10^8)=5761455, pi(10^9)=50847534, pi(10^10)=455052511, pi(10^11)=4118054813. All decimal endpoints here are composite, so replacing b^k-1 by b^k changes no counts.

The classical prime number theorem implies, for fixed integer b>=2,

    d_b(k) ~ 1/(k log b).

This follows by subtracting the two PNT asymptotics at X=b^k and X/b. Since b is fixed, their leading difference is a fixed positive proportion of X/log X. The PNT source can also be checked in the formalized proof: https://isa-afp.org/entries/Prime_Number_Theorem.html .

Define D_b(z)=sum_(k>=1) d_b(k) z^k. An epsilon comparison with sum z^k/k gives

    D_b(r) ~ [1/log b] log[1/(1-r)], r ->1-.

Consequences (PROVED): radius of convergence is 1; D_b is not rational (a rational function near 1 is bounded or has a pole of integral order, neither logarithmic); and D_b cannot be a positive-integral-weight holomorphic level-one modular form in its native q-coordinate. Indeed its constant term at infinity is zero, so modular inversion would force exponential vanishing, times a power, along q->1, contrary to the positive logarithmic divergence. This does not exclude transforms, completions, vector-valued or nonholomorphic constructions.

A logarithmic singularity is not an isolated eigenvalue pole of a scalar self-adjoint resolvent. It does not by itself exclude a continuous spectral measure: multiplication by x on L2([0,1]) has scalar response integral_0^1 dx/(1-zx)=-log(1-z)/z.

## 2. Exact prime densities FAIL positive self-adjoint moment realization

Let m_n=d_10(n+1). If m_n=<v,A^n v> for a self-adjoint operator, with v in the domains of the needed powers, then for every real polynomial p,

    L(p^2)=||p(A)v||^2 >=0, where L(x^n)=m_n.

The exact eleven digit-shell densities violate this. Independent 6x6 Hankel determinant:

    det(m_(i+j))_(0<=i,j<=5)
    = -5505704589894950389 / 132860250000000000000000000000000000.

A small exact certificate, independently found in this lane, is

    p(x)=x-3x^2+5x^4-3x^5 = -x(x-1)^2(3x^2+x-1),
    L(p^2)=-5833441/45000000000 <0.

Thus these native moments admit no positive scalar self-adjoint realization, regardless of spectral support. This is not an obstruction to signed measures, cross-matrix-elements <u,A^n v>, a different sequence, or RH spectral programs in general.

An even shorter shifted-tail obstruction uses L1(x^n)=d_10(n+2):

    q(x)=x(x-1)^2(2x-1)=-x+4x^2-5x^3+2x^4,
    L1(q^2)=-1152073/2250000000 <0.

This needs only d_10(4),...,d_10(10). The shifted tail itself cannot be a Hamburger moment sequence. It should not be confused with the statement that an arbitrary odd shift of a Hamburger sequence remains Hamburger; that is false if its measure has negative support.

### The smooth Li comparator actually DOES have a positive measure

For k>=2 define the continuous prime-density model

    ell_b(k)=1/(b-1) integral_1^b du / [(k-1)log b+log u].

This equals the logarithmic-integral mass of the digit interval divided by its length. Direct positive integral interchange gives

    ell_b(n+2)=integral_0^1 x^n w_b(x) dx,
    w_b(x)=(b x-1)/[(b-1) log(b x)] >0,

with removable value w_b(1/b)=1/(b-1). To prove it, use

    1/[(n+1)log b+log u]
      =(1/log b) integral_0^1 x^(n+log(u)/log(b)) dx

and integrate u from 1 to b. Hence the Li shell model is genuinely a Hausdorff moment sequence, while the exact prime-density tail is not. This is a precise way that prime fluctuations destroy a positive spectral model; it is not an RH implication.

## 3. Three different Hessians must not be merged

### 3.1 Canonical unbiased full-arc Hessian (PROVED)

Fix edge coordinates y_ij=-y_ji, i<j, and define the independent-edge expected Hamiltonian count

    F_n(y)=2^(-(n-1)) sum_(permutations v) product_(a=1..n-1)(1+y_(v_a,v_(a+1))).

This is the unique multiaffine extension of H on the signed Boolean cube. Let B be an oriented incidence matrix of K_n, one column e_i-e_j per edge. For n>=2, set c=(n-2)!/2^(n-2). Then

    Hess F_n(0)=c(2I-B^T B).

Proof: diagonal entries vanish by multiaffinity. Disjoint edge pairs have signed path contributions cancelling under reversing one two-vertex block. For two edges sharing a vertex, the relevant directed paths consist of a contiguous three-vertex block. There are 2(n-2)! contributing orders, with sign minus the incidence inner product; division by 2^(n-1) gives the formula. Since BB^T=nI-J, the score/cut space has eigenvalue -(n-2)c and dimension n-1, while the circulation kernel has eigenvalue 2c and dimension (n-1)(n-2)/2. At n=2 the sole eigenvalue is zero.

Independent verification: enumerate all n! paths for each n=2..7, accumulate exact signed quadratic coefficients, compare every matrix entry and exact SymPy eigenvalue multiplicity. The output is `check_recovery.json`.

This does not assert a global maximum, nor describe the Hessian at a transitive fixed-path tiling corner, nor the Hessian of an arbitrary interpolating function of an adjacency spectrum.

### 3.2 Signed versus unsigned representation

The tangent space of oriented edge signs is the signed permutation module Lambda^2(R^n), decomposing as standard plus Lambda^2(standard), dimensions n-1 and (n-1)(n-2)/2. Divergence B is the score projection. The unsigned permutation module on unordered pairs instead decomposes as trivial plus standard plus the two-row module, dimensions 1,n-1,n(n-3)/2 (n>=4). Both modules are legitimate, but cannot be substituted for one another. In particular vertex transpositions flip some orientation coordinates.

Closest prior provenance: `07-reflections/representation-lens-fundamentals-2026-05-30.md` uses the unsigned pair module, while speaking about score coordinates; `01-canon/theorems/THM-584-complement-is-antipodal-map-level-parity-spectrum.md` has the explicit dated signed-action repair. Empirical n=5 signed Hessian entries +/-0.75 are already in `05-knowledge/results/ramanujan_expander_tournaments.out`; it does not establish all the information-geometric conclusions surrounding them. No priority claim is made for the elementary all-order calculation above.

### 3.3 Fixed-path transitive Hessian: recovered old refutation

`07-reflections/n-2-negative-eigenvalues.md` extrapolates exactly n-2 negative eigenvalues from n<=7. Current `01-canon/theorems/THM-300-quadratic-signature.md` ALREADY REFUTES this at n=9, with n<=12 exact verification. This lane independently reproduced n=3..10 by Hamiltonian DP on every one- and two-tile flip, followed by exact rational congruence, obtaining:

    n=3 (0+,1zero,0-); n=4 (1+,1zero,1-); n=5 (2+,1zero,3-);
    n=6 (6+,0,4-); n=7 (10+,0,5-); n=8 (15+,0,6-);
    n=9 (20+,0,8-); n=10 (27+,0,9-).

This is replication of current canon, not a new refutation. The old reflections need visible historical warnings. The inherited all-order replacement remains OPEN.

### 3.4 Continuous spectral interpolation does not fix a Hessian

THM-134 (`01-canon/theorems/THM-134-paley-local-max-parseval.md`) states a continuous Parseval-simplex local-maximum conclusion using fitted coefficients and says symmetry forces every critical point uniform. The latter inference is false: on x1+x2=1, t=x1-1/2, the symmetric polynomial t^4-t^2/8 has three interior critical points t=0,+/-1/4.

More generally, let a finite spectral sample contain the uniform point u, and let R be the finite set of distinct nonzero squared distances ||x-u||^2 of sample points. The symmetric polynomial

    G(x)=||x-u||^2 product_(r in R)(||x-u||^2-r)^2

vanishes on every sampled point, while its tangent Hessian at u is 2 product_(r in R)r^2 times the identity. Adding any real multiple of G preserves all sampled values but changes the sign of the Hessian at u arbitrarily. Thus discrete fitting plus symmetry cannot certify an intrinsic continuous Hessian. A prescribed extension or a proved functional identity on an open set is necessary. This diagnoses a proof gap and disproves the quoted general symmetry implication; it does not refute a specified exact polynomial's directly computed local maximum.

## 4. Broad inheritance search and surviving connections

Initial physically present tracked-file scan: 72,887 UTF-8-readable files; 40,709 tracked paths absent; 12 existing files over5MB excluded. Query match-file counts: fraction substring768, Hessian485, representation theory230, literal {2,3,11}130, repunit173,1729 904. Fraction substring search has false positives such as 321/9016; a second digit-boundary check found no exact21/90 or143/900 occurrence in the recorded hits. This is not an exhaustive all-history or all-remote-file claim. `inventory.json` preserves paths, line numbers, queries, and snippets.

Selected underused prior sources:
- `THM-134-paley-local-max-parseval.md`: continuous-interpolation gap above.
- `THM-135-paley-not-h-max-p19.md`: discrete Paley-vs-interval counterexample remains separate from continuous interpolation.
- `THM-138-ising-decomposition.md`: OCF disjoint-cycle interactions are a better retained sidecar than scalar spectral flatness; its informal causal extrapolations are not inherited as proofs.
- `THM-300-quadratic-signature.md`: current refutation overrides old Hessian reflections.
- `THM-584-complement-is-antipodal-map-level-parity-spectrum.md`: signed S_n action preserves Fourier degree, not Hamming weight; its metagraph is the full single-arc-flip quotient, not the fixed-path complement wiring.
- `THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels.md`: exact classical Vandermonde cancellation links scores with alternating projection. It loses all intransitive signed mass, so cannot alone recover H; do not inherit its broader SL_n/bracket wording without type checks.
- `THM-3367-berggren-spinor-pencil-hessian-gauge-and-affine-line-keller-closure.md`: an orthogonal Hessian-gauge/representation lane; no direct lattice identification claimed.
- `THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md`: the actual prior {2,3,11}/341 mechanism is modular strong-tournament counts and Camion pseudoprimes. This is a structural arithmetic mechanism, unlike merely observing88=8*11.
- `07-reflections/information-geometry-of-tournaments.md`: valuable Fourier/cube provenance but includes unsupported topology and Fisher deductions; not a source for those statements.

Reproduction:

    python3 check_recovery.py > check_recovery.json
    python3 check_transitive_hessian.py
    python3 certify_inertia.py > transitive_hessian_certified.json
    python3 check_prime_hankel.py > prime_hankel.json

The independent first script needs SymPy; the transitive probe uses NumPy only for a non-authoritative preliminary signature. The rational congruence certificate is authoritative. No finite positivity test is promoted to an all-order realization.

### Reproducible inventory replay

A later replay while the shared checkout continued changing reports114,041 tracked paths:72,889 readable UTF-8 files,40,718 missing,12 oversized,422 unreadable. Match-file counts are fraction769, Hessian485, representation230, triple130, repunit175,1729 905. The original inventory is retained; `inventory-replay.json` supplies the new tracked-manifest SHA256. `recover_inventory.py` never runs git grep or reads unavailable git blobs. `REPUNIT-INVENTORY.md` expands all175 repunit file matches and their recorded line snippets. This is extensive present-checkout coverage, not all historical versions.

`check_li_comparator.py` also verifies the measure identity for n=0..8 at50decimal precision and evaluates the same squared polynomial positively under the Li model, negatively under actual prime densities. Quadrature is a sanity check only; positivity and the density formula follow from the exact integral proof.

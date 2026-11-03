# Session 143 — independent mathematical audit

**Disposition: PASS for the explicitly scoped results below.** This is a
proof/finite-certificate audit, not a claim that the proposed universal chain
exists or that an open problem has been solved. Repository source baseline:
`8286d710498a3b4e9b9f1ac683c6eb77716f3f4b`.

## 1. Tournament and sporadic lane

Independent ordered-arc enumeration recovers 12 classes, eight self-converse,
ten modulo converse and nine spectra. Newton identities independently recover
characteristic polynomials and agree with directed-cycle enumeration. The
regular class, its C5 automorphism group and D10 extended action, and the
order-20 affine stabilizer of its four labeled circulant orientations all pass.
The elementary Burnside/uniqueness/normalizer proofs are sound. Two imaginary
magnitudes within one spectrum are not two regular tournament spectra.

The generic complex-representation theorem for `rho(g)+rho(g)^2` is sound.
Exact Fourier arithmetic separately gives 5A multiplicities 891 and 870, and
5B multiplicities 871 and 875, **conditional on the cited degree-4371 character**.
The nonzero traces exclude loopless adjacency operators. Rational conjugacy
classes, ordinary character values, actual eigenvalues and group-algebra
operators are correctly distinguished. The Q(zeta5)/Q(zeta17) field-intersection
obstruction, transcendence-degree obstruction for Q(X0(13)), and weight-two
cuspidal restriction are valid at their stated scopes.

[McKay's catalogue](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html) and
[ATLAS B](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B) were independently
opened. The catalogue confirms the graph counts; ATLAS supplies the subgroup
occurrence and class data, not a representation intertwiner.
[Breuer–Magaard–Wilson](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/htm/ctblbm.htm)
sections 6–7 support the rational classes and character transcription.
[Steel's thesis](https://magma.maths.usyd.edu.au/~allan/reps/AllanSteelPhD.pdf),
section 5.11, supports minimal faithful ordinary degree 4371.
[Daniels et al.](https://math.mit.edu/~drew/mcom3213.pdf), printed page 446,
supports genus zero and the j-map.
[Conway–Elkies–Martin](https://arxiv.org/html/math/0508630v1), section 2.2,
supports the 13-point puzzle and non-group M13 semantics. These citations
were checked, not independently reconstructed. No local GAP execution,
Baby Monster matrices, subgroup straight-line execution or intertwiner is
claimed. A characteristic-zero 24-dimensional nontrivial Baby Monster
representation is excluded; modular and nonlinear actions are outside that claim.

## 2. Magic-cycle lane, including complete cubic classification

The elementary central symmetry and third-centered-moment arguments are
correct over the reals. Independent coefficient arithmetic gives
`(44/3) A(65 A^2+63 B^2)` for the nine quadratic samples. Therefore a genuine
quadratic cannot be rearranged into a real 3-by-3 magic square. The all-real
AP-root-square corollary and separate center/maximum proof are sound, without
assuming distinct or integral values. The general odd-packet formula is
algebraically correct and its m=1 boundary is correctly separated.

The complete cubic result is independently verified, not merely inherited:
for `f(k)=A k^3+B k^2+C k+D`, `k=-4,...,4`, a magic arrangement exists exactly
when `B=0` and either `A=0` or `C=-15A`.

- A separate exact Gaussian checker enumerates all 945 center-index/perfect-
  matching systems. Its nullity census is 914 zero-dimensional, 30 lines,
  one plane. Rational row operations describe the entire real kernels.
- The only nonodd rays are `(1,±3,-18)` and `(1,±3/2,-39/2)`. Their centered
  absolute quartets are `{20,20,34,36}` and `{10,27,28,35}`. Each violates
  the necessary largest-versus-two-smallest inequality, and independently
  fails all signed quartet assignments.
- On the odd plane, all 384 signed-permutation linear systems are checked:
  eight survive, all at `C/A=-15`; 88 have a nonzero constant obstruction
  and 288 incompatible parameter constraints. Zero/repeated values do not
  evade the indexed-pairing argument. Linear and genuine cubic witnesses
  verify sufficiency. This is a classification of numerical assignments,
  not squarehood of those assignments.

The nine-cycle is genuinely +4 in Z/9, not an eight-vertex cube. Its digit
carry formula is correct. The homogeneous-matrix proof excludes every affine
order-nine permutation of F3^2 under any relabeling; independent enumeration
of all 432 affine maps confirms the order census. The odd-cycle bipartite
obstruction is correct for loopless graph homomorphisms to a hypercube.

The inherited nonmagic square-total example must stay distinct from an
MSOS. The extra demand that the common total itself be an integer square is
impossible for nonzero integer-square entries; the integer/rational domain is
essential. The real all-ones grid has total 3, itself a real square. This audit
requested explicit arithmetic qualification, which the producer adopted.
[Auel–Singer v1](https://arxiv.org/html/2609.09351v1), dated September 2026,
was independently opened and retains the distinct-integer-square question as
open. This is a primary status check, not exhaustive literature certification.

## 3. Octonion and conditional cube lane

An independently implemented Hamilton-quaternion doubling (not the producer's
recursive signed-basis code) verifies all 64 multiplication entries and 512
associators in the stated convention. The associator **ratio** is the
binary determinant sign; the additive associator is 0 or ±2 times a basis
unit. Its nonzero support is 168 ordered independent triples. A basis gauge
can flip a raw ±2 coefficient while preserving the ratio; the report now
keeps these notions separate.

The 16 binary-supported subalgebras form the F2^3 subspace lattice. Each
chosen basis picks an eight-vertex Boolean sublattice, with 28 choices.
Dimensions 1,2,4,8 and R/C/H/O types are justified by closure and associativity,
not vertex-count analogy. The diagonal character action is an actual algebra
-automorphism action. The scalar-field correction requested by this audit is
adopted: `r_x -> e_x` is Q-linear from the multiquadratic field to the rational
span O_Q, followed by real extension as needed; it is not an algebra map or a
Q-vector-space isomorphism to all real O.

The perimeter left-unit action has square -I and gives an order-four signed
lift of the square-grid half-turn after coordinate squaring. It preserves the
actual eight-line magic predicate because it permutes grid lines. It is a
norm-isometry, not an octonion automorphism, and changes no squared multiset.
Only the specified two XOR translations preserve the magic equations in this
labeling; general binary translations fail. No existence proof follows.

[Albuquerque–Majid](https://arxiv.org/html/math/9802116v1) gives the classical
twisted-group-algebra framework. [Baez §2.1](https://math.ucr.edu/home/baez/octonions/node4.html)
explains the binary/Fano subalgebras. [Barendregt](https://homepages.inf.ed.ac.uk/wadler/papers/barendregt/pure-type-systems.pdf)
provides the type-system inclusion cube. All were independently opened. The
last comparison stays **conditional**: the user's intended lambda-cube referent
is unresolved. No term, typing judgement or normalization transport is claimed.

## 4. Verification and limitations

`check_independent.py` executes 2,188 explicit checks. `check_full_cubic.py`
separately exhausts the 945 full-cubic systems. Neither imports producer code;
both use only the standard library and survive Python `-O`. Frozen ordinary
and optimized outputs match byte-for-byte. Wrong spectrum, moment, octonion
sign and full-cubic ray mutations must fail in both modes. The portable
manifest verifier reruns these positive and negative controls from any cwd.

The sporadic lane's original bare assertions were an audit finding, not a
mathematical refutation. They were replaced with explicit gates; the repaired
positive/source/transcription-mutation regression passes in both Python modes.

Named repository blobs are checked against source OID, SHA256 and size, with
no broad hydration. A source pin is provenance, not endorsement of every
historical theorem or experiment in that file. In particular, inherited LRC
experiments and malformed associator normalization are not promoted to new
transport theorems. No general MSOS, LRC, Baby Monster correspondence,
Hecke-equivariant construction, or priority claim is established.

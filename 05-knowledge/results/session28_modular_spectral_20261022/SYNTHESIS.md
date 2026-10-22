# Session28 synthesis: which mathematical structures survive the proposed bridges?

Research date2026-10-22. Status: audited proofs, exact computations and explicitly cited external inputs. No RH, Collatz, global packing, or complete tournament-spectrum claim. This report is a research result packet, not a claim of priority over classical methods.

## 1. Prime-density coefficients have a decisive positivity obstruction

The sequence4/9,21/90,143/900 is the proportion of primes in successive decimal digit blocks. Its generating function has logarithmic boundary growth, by PNT, and cannot be rational or a positive-weight holomorphic level-one modular form in its native q-coordinate.

More decisively, if L(x^n)=delta_(n+1), then

    L([x(x-1)^2(x+1)(2x-1)]²)=-2447633/22500000000.

So these coefficients are not moments of any positive measure on R or of a self-adjoint operator observed through one vector. This is not a refutation of RH. A logarithm alone would not suffice: continuous positive spectrum produces logarithmic resolvents. The Li approximation actually does have an explicit positive Hausdorff measure, while the exact arithmetic sequence fails the test.

[Proof and source boundary](PRIME_DENSITY_SPECTRAL_OBSTRUCTION.md), [independent recovery and alternate certificates](session28-recovery/REPORT.md), [cross-audit](session28-lorentzian/AUDIT_CROSS.md).

## 2. Rank88 supplies a genuine four-dimensional Hecke response

The cited2026 rank88 constructions use Q(sqrt(-23)), with the structural dimension factorization88=2*4*11. We independently verified a catalogue Gram matrix is even, integral, positive definite and determinant1; its minimum8 remains cited from the construction paper, not independently re-enumerated here.

The extremal theta is

    E4^11−2640E4^8Δ+1813680E4^5Δ²−244992000E4²Δ³.

It is modular but not a Hecke eigenform. Its T2 orbit spans all four dimensions of M44. The three cusp eigenvalues have an explicit positive trace Gram, giving an honest self-adjoint representation. In contrast, the recovered level22 bad-prime U2 has eigenvalues−1±i. 'Hecke' is not, by itself, a promise of a real spectrum.

[Exact independent checkpoint](HECKE_CHECKPOINT.md), [full rank88 report](session28-rank88/REPORT.md).

## 3. Rank88 Gaussian geometry has one fixed angular operator

For any actual extremal rank88 lattice, all determinant-one Gaussian energy Hessians have the form

    Hess_t=a(t)I+b(t)A, b(t)>0,

where A is one fixed quartic minimum-shell operator on the3915-dimensional trace-free symmetric-matrix space. Scalar theta fixes a,b but loses A. The Hessians commute and keep the same eigenvalue ordering as temperature varies. For the cited perfect shells, sufficiently hot and cold temperatures are stable; any transitions form a finite reciprocal set. This does not prove stability at all temperatures.

The dimension residue explains part of the8/24 pattern: for extremal dimensions n=24m+r, r=0,8,16, the quartic harmonic quotient has weight r/2−8. It vanishes for r=0,8 and can be a constant for r=16. Thus ranks24m and24m+8 have scalar Gaussian Hessians in this setting, whereas24m+16 can retain a quartic angular operator. This is a modular-weight mechanism, not numerical resemblance.

A late hostile probe changed the scalar story: the self-dual temperature t=1 is a rigorously certified strict LOCAL MAXIMUM of r(t)=a(t)/b(t), not the minimum suggested by a coarse grid. Exact modular differentiation and rational bounds prove negative curvature. Therefore at least two distinct reciprocal off-center global minimizers exist. Their uniqueness and locations, and the actual angular eigenvalues of the88-dimensional lattices, remain open here.

[Full theorem](session28-rank88/REPORT.md), [independent audit](session28-cm-audit/AUDIT.md), [certified maximum and minimizing-pair theorem](session28-cm-audit/scalar-followon/ADDENDUM.md).

## 4. D=-163 really can be glued into E8, but is not exceptional for this map

Four copies of the binary Gram matrix[[2,1],[1,82]], with determinant163 each, acquire an explicit isotropic discriminant-plane overlattice which is E8. The CM endomorphism survives and satisfies alpha²−alpha+41I=0 and alpha*alpha=41I. The associated J=2alpha−I is skew-adjoint, with J²=−163I; it is not a lattice isometry.

The construction extends to every positive d congruent3 modulo4, including63. Thus it gives a genuine CM/lattice bridge while disproving the inference that the bridge itself detects Heegner class number one or explains tournament repunit rarity. Gluing changes the theta series and modular level. A level163 CM newform is not thereby turned into a level-one CM eigenform.

[Construction, matrices and generalization](session28-rank88/REPORT.md), [independent240-root enumeration and composite controls](session28-cm-audit/AUDIT.md).

## 5. Lorentzian score geometry has a sharp symmetry boundary

For25-vertex tournaments, the score lift v=(s;70) satisfies v²=−2c3. Yet the full origin-fixing integral Lorentzian isometry group preserving ALL such score lifts is only S25. Consequently only coordinate-swap root reflections preserve the entire score set. A specific Leech-root reflection exits the score slice while preserving Lorentzian norm.

Converse is a valid affine lattice symmetry but has no linear integral Lorentzian extension in this fixed embedding. These qualifications prevent an overbroad no-go claim.

For cyclic tournaments the hyperbolic radial identity is

    cosh d(center,T)=sqrt(650/c3(T)).

A cusp is at infinite hyperbolic distance. The ranking pairing is a horospherical statistic, not a literal finite distance to that cusp. The geometry still forgets Hamiltonian information: Session27's equal-score657/643 example remains a hostile.

[Proof, metric and exact controls](session28-lorentzian/GEOMETRY.md).

## 6. A precise tournament Hessian connects representation theory and lost information

For the unique independent-edge multiaffine extension F_n of the Hamiltonian count at the unbiased point,

    Hess F_n(0)=((n−2)!/2^(n−2))(2I−BᵀB),

with B the signed incidence matrix. Its negative eigenspace is the(n−1)-dimensional score module, and its positive eigenspace is the(n−1)(n−2)/2-dimensional circulation module. The latter is absent from the score norm. The correct signed edge representation differs from the unsigned pair permutation module.

This Hessian is proportional to the natural expected-cyclic-triangle Hessian. But extending the Boolean score identity as a quadratic norm gives a different Hessian: equality at all Boolean tournaments does not determine off-cube derivatives.

Disjoint triangle circulations give an all-order exact response polynomial and Gamma integral. At order25 one fully biased triangle raises expected H by exactly1% without changing expected scores. This is an ensemble statement, not a deterministic score-fiber move.

[Representation proof](session28-weil/REPORT.md), [independent verification](session28-recovery/REPORT.md), [exact circulation polynomial](session28-lorentzian/TRIANGLE_FIBRES.md).

## 7. Fixed-H padding counts have Jordan, not positive-moment, behavior

A nonconstant eventually polynomial positive sequence violates eventual even-subsequence log-convexity, so cannot be positive scalar moments of a self-adjoint operator. Session27's fixed-H polynomial padding counts inherit this obstruction. The63 sequence2(m+1) has an exact nonnormal Jordan realization instead. A determinant's repeated poles and a scalar vector resolvent's poles obey different rules; eigenvalue multiplicity must not be confused with pole order in every observable.

[Proof and survivors](session28-weil/REPORT.md).

## 8. Explicit formulas yield a rigorous sampling map, not an RH proof

For normalized Chebyshev-psi blocks at bases such as6 or10, finite zero truncations give poles at b^(1−rho). The real part of rho controls radius; height controls phase modulo2pi/log b. Aliased zeros cannot cancel in this particular finite model. Two multiplicatively independent bases determine a single matched zero uniquely, but do not solve matching or stable recovery.

A standard sharp explicit-formula remainder gives locally uniform convergence of the truncated generating functions inside the unit disk, with error O_(b,r)((1+log²T)/T) on |z|<=r<1. No interchange is taken on faith. Pole persistence or continuation beyond that disk remains a separate problem; the finite circle-of-poles statement is not a global RH equivalence.

[Proof, source pin and remaining boundary](session28-weil/Log_Scale_Explicit_Formula.md).

## 9. Guinand-Weil positivity needs the omitted directions and their coupling

We recovered a concrete Cayley normalization error in a Conway-cusp audit: (z−i)/(z+i) tests the upper half-plane, not the right half-plane. The surviving Blaschke example shows contractivity alone does not exclude interior zeros.

For block forms, a rigorous sufficient positivity budget is a*gamma>=beta² when the two diagonal blocks are bounded below by a and gamma and the coupling has norm at most beta. Positive finite compressions alone do not control the complement. Infinite-dimensional transfer needs the domain, dense form core and tail estimates. None of these elementary criteria supplies the currently missing global Weil positivity proof.

[Source audit](session28-weil/REPORT.md), [tail and coupling criteria](session28-weil/Tail_Criteria.md).

## 10. Historical recovery and repairs

The broad recovery replay inspected72,889 readable present files among114,041 tracked paths, with absent/oversized/unreadable counts preserved. It is not an exhaustive all-history claim. The repunit subset records175 matching present files and their snippets.

Three narrow repairs preserve the old calculations while correcting current status: THM134's unsupported continuous Parseval-simplex claim becomes OPEN; the old n−2 reflection points to its existing THM300 refutation; a historical finite reciprocal 'L-function' script gains an AST-preserving warning that its claimed zeros are poles.

[Recovery and inventory](session28-recovery/REPORT.md), [correction scope](session28-recovery/CORRECTIONS.md), [exploration record](EXPLORATION.md).

## Next substantive obligations

1. Compute or sharply bound the fixed quartic operator's least eigenvalue using the rank88 automorphism representation, retaining multiplicity spaces.
2. Certify the scalar Gaussian ratio's global minimum or stationary-point count; symmetry alone has now been decisively ruled out as a minimum proof.
3. Establish outer-domain convergence/pole persistence for logarithmically sampled explicit formulas, with matching and separation if two bases are used.
4. Retain marked tournament response data when transporting score, lattice or spectral information toward metagraph wiring. Scalar projections are demonstrably insufficient.

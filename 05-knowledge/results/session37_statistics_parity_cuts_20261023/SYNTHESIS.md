# Session37: exclusion, exchange signs, and arithmetic completion

Status: proved identities and obstructions, independently checked finite certificates, and explicitly cited structural inputs. The session treats 'fermionic' and 'bosonic' as precise algebraic constructions where justified, not as a claim that tournaments are physical quantum systems.

## 1. The main answer: several statistics, not one interchangeable duality

The useful distinction is between three operations:

- **Excluding repeated vertices:** a hard support constraint, naturally encoded by commuting square-zero variables.
- **Changing signs under exchange:** exterior/fermionic antisymmetry, which is not the same constraint.
- **Allowing repeated components:** symmetric/bosonic occupation, whose interaction with exchange signs can be classified exactly.

The arithmetic distinctions are another axis. Rational versus irrational partitions the real numbers; algebraic versus transcendental is a different partition. A transcendental real is still real, and an algebraic irrational is neither rational nor transcendental. Dedekind cuts describe all real numbers, not one side of these partitions.

The positive results below connect these axes through explicit maps. The counterexamples say which tempting identifications do not survive those maps.

## 2. A matrix dictionary prevents false mergers

Throughout the linked papers:

- A is the 0/1 tournament adjacency matrix.
- S=A−A^T is its integral skew matrix.
- K=−S² is a positive semidefinite energy/weight matrix.
- X contains commuting nilpotent vertex labels in path/cycle formulas.
- P_g is a permutation matrix acting on the vertex representation.

The adjacency determinant/permanent pair, the skew/energy spectrum, and det(I+Λ²P_g) in a character computation are DIFFERENT invariants. In particular a five-vertex skew-cospectral H=3/15 pair does not contradict recovery from the different adjacency determinant/permanent pair through eight vertices.

## 3. Hard exclusion has an optimal paired-exterior realization

[Exclusion and mixed statistics](exclusion/RESULT.md) recovers the exact path algebra Q[x_1,...,x_n]/(x_i²). Its variables commute; repeated use of a vertex annihilates a term without adding exchange signs.

The substitution x_i=θ_iη_i embeds this algebra faithfully into the even part of an exterior algebra. For n>=1, two odd generators per vertex are not an arbitrary convenience: 2n is the exact minimum for a faithful even-sector realization. If parity preservation is dropped, the exact minimum is2n−1. A four-vertex hostile shows that simply replacing x_i by one odd variable does not compute the unsigned path count.

The determinant/inverse-determinant cycle filter has an exact formal two-species interpretation. Setting both multiplicities to1 cancels even-length cycles and gives odd-length cycles weight2; square-zero labels retain disjoint support. The determinant/permanent path identities themselves are inherited from the cited Irving–Omar and Grinberg–Stanley work, not claimed as new.

## 4. A genuine component Fock model for tournament enumeration

[Exterior tournament character](parity/EXTERIOR_TOURNAMENT_FOCK.md) proves that the ungraded S_n-character of Λ(Λ²Q^n) is exactly the permutation character of labelled tournaments. Its invariant dimension is therefore A000568.

In the invariant or coinvariant space, a graph-monomial orbit survives precisely when every graph automorphism has determinant+1 on oriented edge space. This sign is edge-reversal sign TIMES edge-permutation sign. It is not the older Royle automorphism-even condition.

Connected surviving components with even numbers of edges repeat freely; those with odd numbers of edges occur at most once. Their graph algebra is an actual Sym(even)⊗Λ(odd) construction. If b_n and f_n count those two connected species, then

A(z)=product_n (1−z^n)^(-b_n)(1+z^n)^(f_n),
b_n−f_n=δ_(n,1)−δ_(n,3).

An exact dyadic recurrence relates b_n,f_n to the inherited connected Royle-even counts. Direct graph and character checks support the proof through the stated finite ranges.

This is a linear representation/statistics bridge. It does not identify tournament converse with exterior parity: they already disagree at n=3. Nor does it supply a canonical multiplication-compatible set bijection between tournaments and graphs.

## 5. Vertex parity produces a different, valid superalgebra

[Super-PBW and signed curl](parity/SUPER_PBW_AND_SIGNED_CURL.md) separately grades strong tournament generators by vertex number. Odd Lie generators have nonzero squares in the tensor algebra; super-PBW moves those squares into additional even generators. Omitting that correction would already erase the actual two-vertex transitive tournament. Here Sym(even)⊗Λ(odd) describes the PBW associated-graded algebra, not an associative-algebra identification of the original SCC word algebra.

In the super-cocenter, exactly even repetitions of odd-weight primitive SCC necklaces vanish. The note gives the signed-curl dimension and its alternating-group orbit interpretation. Ordinary H is not a supertrace: the odd singleton a has supercommutator[a,a]=2a², but H(a²)=1.

The two statistics models must remain distinct. Their connected-space total dimensions agree at orders3–7, but already differ at2 (0 versus1) and differ again at8 (6414 versus6413). Their parity assignments already disagree at1. This is an explicit warning against turning a short numerical match into an identification of structures.

## 6. A transparent interacting cycle gas and exact generating function

[The recovered real-root boundary](root/CYCLE_GAS_REAL_ROOT_BOUNDARY.md) credits the old THM-020/025 threshold: odd-cycle gas polynomials are real-rooted through order8 but need not be at9. A simple new witness joins three directed triangles and reverses one first-to-last arc, giving

I(z)=1+49z+7z²+z³,   H=I(2)=135.

Its nonreal roots prevent a positive free-fermion realization det(I+zK) for finite Hermitian K>=0. This does not forbid interacting, signed, nilpotent or non-Hermitian constructions.

[The one-back-arc family](root/ONE_BACK_ARC_CYCLE_GAS.md) extends the construction to m triangles. All new cycles share the reversed arc, yielding the exact formula

I_m(z)=(1+z)^m + z/2[9(z+10)^(m−2)−(z−2)^(m−2)].

Every member with m>=3 is strong and has nonreal cycle-gas roots. At the Hamiltonian fugacity2, one parity channel vanishes and

H(U_m)=3^m[1+4^(m−2)]  for m>=3,

while m2 is the exceptional value17. A rational generating function and exact recurrence are proved, including the correct starting index. The formula comes from a retained common-arc interaction, not independent-particle factorization.

## 7. Two older claims are corrected, not quietly repeated

[Independent inheritance audit](../session37_inheritance_audit_20261023/REPORT.md) repairs current canon and records the error mechanisms.

1. Adjacency characteristic plus permanental polynomials determine H through order8, not merely through7. Two independent complete computations cover all2,097,152 fixed-path presentations at8, with3217 invariant fibres and no H split. Failure first occurs at9: an explicit equal-polynomial pair has H2843 and2841. Source extension propagates failure to all larger orders. The order8 assertion is finite-exact/computer-assisted; the order9 witness and propagation have explicit proofs.
2. The old rank obstruction at8 was real, but obstructed AFFINE reconstruction, not arbitrary functional recovery. Three actual tournaments exhibit that distinction. At9 the realized missing-information polynomial is z^9(u³−u)+z^8(v²−u²), invisible to both determinant/permanent evaluations but not to H.
3. The bosonic THM-1810 unsigned Vandermonde identity was false: the permanent counts transitive tournaments, while the unsigned edge product counts all tournaments. They differ already by6 versus8 at three vertices. Positive Wick matching signs also do not prevent polynomial coefficient cancellation. The valid earlier detection-depth argument and attribution are retained.

A concurrent incoming repair of THM-1970 independently withdrew 'counting hardness implies invariant non-determination'. This is the same logical discipline: information sufficiency, linear decoding and computational difficulty are distinct questions.

## 8. Exact arithmetic maps and certified real cuts

[Arithmetic cut completion](arithmetic/ARITHMETIC_CUT_COMPLETION.md) starts from nonzero eigenvalues μ of K=−S². They are totally positive algebraic; Q(μ) is totally real and Q(μ,√−μ) is CM. The Cayley phase and inherited Coulomb charge satisfy

t=(1+i√μ)/(1−i√μ),
c=13−6(t+t^(-1))=(1+25μ)/(1+μ).

This preserves specified fields, not merely approximate numerical values. The full rational Cayley matrix retains S by an inverse transform. Squaring to K already identifies S and −S; taking spectra forgets further data.

The arithmetic types do not pair naively. At order4, irrational energies3±2√2 give eighth-root phases and rational normalized angles. At order5, the rational energy7 gives charge22 and a transcendental normalized angle θ/π. Five is the sharp first order for this normalized-angle phenomenon. The raw nonzero angle θ is already transcendental even in the root-of-unity cases, a different statement.

At positive rational fugacity within the stated domain, det(I+zK) and det(I−zK)^(-1) are rational, while their nonzero real logarithms are transcendental by the cited Hermite–Lindemann theorem. Explicit rational series bounds define convergent lower and upper cuts. The logarithms’ half-sum and half-difference isolate the odd and even trace-power sectors, respectively. The transcendence input is classical, not a new transcendence theorem.

The P163 example reaches Q(√−163), phase(−81+i√163)/82 and charge1019/41 through this exact map. The construction does not use the Heegner class-number-one property and implies no RH or modular/VOA identification.

## 9. Arithmetic richness does not restore the missing path information

For any regular odd base T, add two singleton blocks and allow all eight uniform quotient orientations. They share the full skew and energy spectra and every phase/charge/determinant/log observable derived from those spectra. Six transitive choices have H(T); two cyclic choices have3H(T)+2pc_T(2), where pc2 counts unordered spanning covers by two nonempty directed paths, with singleton paths allowed.

Cutting Hamiltonian paths gives2pc_T(2)>=(n−1)H(T), so the ratio is at least the extended order N=n+2 and is unbounded. A skew/energy-spectrum-only positive estimator therefore has symmetric multiplicative error at least√N on one member. This is an explicit all-order information-loss theorem, not an inference from #P-hardness.

The cyclic extension is strong and the transitive one reducible. No unbounded-ratio claim within a fixed strong class, score class, or regular target class is made. Adjacency characteristic/permanental data are not asserted equal here.

## 10. What was recovered, and what remains open

The session recovered March cycle-gas obstructions, old determinant/permanent and Gaussian-moment notes, the even-graph and SCC/Lyndon constructions, Dedekind/order-observer blindness (THM4088/4095), and the more recent CM/Coulomb and curl work. Standard representation theory, super-PBW and transcendence inputs are cited in their packets. Repository-relative advances are not asserted to be first discoveries in the literature.

Incoming work during the session also supplies exact metagraph reconstruction from authenticated actual ports, while leaving classification from small scalar inventories, black-core connectivity and faster general enumeration open. These inputs are preserved, not overwritten or claimed as this session's results.

The common direction is precise: first specify occupancy, exchange sign, retained labels and allowed operations; then identify which quotient or analytic completion is being taken. A powerful new representation can expose a real structure while still destroying the invariant we wanted. The counterexample and inverse-map obligations are as important as the analogy.

No complete metagraph classification, full Hamiltonian-count spectrum, uniform153N coverage, LRC, Collatz, Yang–Mills or RH proof is claimed.

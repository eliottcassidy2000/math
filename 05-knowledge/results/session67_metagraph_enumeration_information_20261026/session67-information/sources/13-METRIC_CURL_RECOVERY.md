# The augmentation metric repairs the Heisenberg quotient

Session36. Status: PROVED elementary identities and fibre classification; FINITE-EXACT certificates separately identified. This is a new connection in this session, not a literature-priority assertion.

## Inheritance and the exact missing coordinate

- Session35 `session35-recursions/AUGMENTATION-LATTICE.md`: the alternating form on the integral augmentation lattice has odd determinant in odd tournament order; its integral alternating/Heisenberg type already loses Hamiltonian count at order5. The two explicit Smith-(1,1) witnesses have H=11 and H=5.
- THM-4739, `01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md`, §1, records score/gradient and triangle/curl decomposition. The file is scoped PROVED (pre-audit) and has other, separate finite claims; the needed elementary Hodge algebra is reproved below. Its per-vertex curl-profile compression is not the full labelled triangle curl.
- Session35 `session35-zeta/CONFERENCE_COLLAPSE.md`: all128 specified Paley8 switches have identical Q=8I and marked spectral pair (p,F), yet H=189 or333.
- THM-868, `THM-868-e8-bridge-score-lattice.md`, puts score half-vectors in an E8 slice; Session27's score-root collision shows scores do not determine H. We do not infer a general exceptional-lattice theorem from the numerical rank.
- Historical Hessian caution: THM-134, `THM-134-paley-local-max-parseval.md`, is explicitly OPEN after the2026-10-22 correction. Finite spectral interpolation does not specify a continuous H extension or its Hessian. No claim here uses that historical Hessian.

The new repair keeps the standard Euclidean metric AND the integral lattice together with the alternating form. Keeping its Smith type alone, or its real metric-operator spectrum alone, is insufficient. The simultaneous integral isometry class is a different, richer invariant.

## 1. Objects

Let T be a tournament on [n], n>=1, and M its skew matrix: M_ij=1 for i->j, M_ij=-1 for j->i, M_ii=0. Set

L=A_(n-1)={x in Z^n: sum x_i=0},
g(x,y)=x^T y,
omega_T(x,y)=x^T M y.

With root basis b_i=e_i-e_(n-1), i=0,...,n-2, write B=(b_i), G=B^T B=I+J and K_T=B^T M B. Thus

K_ij=M_ij-M_i,r+M_j,r, r=n-1.

This is the labelled triangle curl on triangles through r. The general curl is recovered by K_ij+K_jk+K_ki. All results below hold in even as well as odd n. For even n the augmentation lattice has odd dimension, so the restricted alternating form is necessarily degenerate; no symplectic/nondegeneracy assertion is needed.

## 2. Fixed-label fibre theorem

THEOREM. Two tournament matrices M,N have the same restricted alternating form on L if and only if either M=N or there is a nonempty proper subset X of vertices such that all cross-arcs of M point X->X^c, N agrees with M internally on X and X^c, and N reverses exactly those cross-arcs.

Proof. Put D=M-N. If B^T D B=0, then D_ij-D_i,r+D_j,r=0. Put a_i=D_i,r and a_r=0; hence D_ij=a_i-a_j. Every entry of D belongs to {0,2,-2}. Thus every pairwise difference a_i-a_j belongs to that set. If a is not constant, it has exactly two levels separated by2. Let X be its upper level. Then D_ij=2 for i in X,j outside X, forcing M_ij=1,N_ij=-1; within either level D=0. This gives the stated cut. Conversely its difference matrix is a1^T-1a^T, which vanishes on L x L. QED.

COROLLARY (complete fibre). Let the ordered strong components of T be C_1->C_2->...->C_k. Its fixed-label K-fibre consists of exactly the k tournaments obtained by cyclically rotating this ordered component list, retaining the labelled subtournaments C_i unchanged. In particular, a strongly connected tournament has a singleton labelled fibre.

Proof. A uniform cut X->X^c cannot split a strong component: a return path from the second side to the first would be impossible. Hence X is exactly an initial segment of the ordered SCCs. Its cut reversal rotates that segment to the end. The k rotations are distinct as labelled tournaments when k>1, since a nontrivial rotation reverses some cross-block arc. Include the empty rotation. QED.

COROLLARY (Hamiltonian recovery). K determines H(T), even though it need not determine T. Indeed H(C_1->...->C_k)=product_i H(C_i): every directed Hamiltonian path traverses each component contiguously and in this order, and arbitrary component Hamiltonian paths concatenate. Rotating factors preserves the product. The same argument preserves every multiplicative invariant of the ordered SCC decomposition.

This is not a claim that the SCC list can be reordered arbitrarily: the full K-fibre permits only cyclic rotations.

## 3. Unlabelled integral-metric classification

LEMMA. Every integral isometry of A_(n-1) is the restriction of a signed coordinate permutation, with a single global sign. For n>=3 the group is {+/-P: P in S_n}; at n2 it is the two-element isometry group of A1 (the two descriptions overlap); n1 is the zero lattice.

Elementary proof. The minimal nonzero vectors are exactly roots e_i-e_j of squared norm2. An isometry maps the basis b_i to n-1 roots with pairwise inner product1. For n>=3 such a family has a common positive endpoint or a common negative endpoint: take two roots with inner product1, which share one same-sign endpoint; any further root pairing1 with both must share that endpoint, since the two other endpoints differ. The resulting family is precisely a signed coordinate permutation of the standard star basis. The rank-one and rank-zero cases are immediate. QED.

THEOREM. The metric alternating lattices (L,g,omega_T) and (L,g,omega_U) are isometric if and only if the ordered SCC lists of T and U agree up to a cyclic rotation and componentwise tournament isomorphisms. In particular, this invariant classifies strong tournaments and always determines H.

Proof. The isometry is +/-P. The global sign appears twice in the bilinear form and cancels, so this is exactly equality of restricted forms after relabelling by P. Apply §2. Conversely a componentwise isomorphism followed by one cut rotation gives such an isometry. QED.

Global -I does NOT reverse omega. Tournament reversal sends omega to -omega, and is not automatically identified. Allowing anti-isometries of the alternating form would additionally identify reversal of the SCC cycle with reversal of each component; H is still unchanged, but this is a distinct convention.

## 4. Precisely why the metric must be integral

Let p_M(t)=det(tI-M), q_M(t)=1^T adj(tI-M)1. Then

det(tG-K_T)=q_M(t).

Proof. Choose an orthonormal basis U for 1^perp and append u=1/sqrt(n). The cofactor in the u position of tI-M is det(tI_(n-1)-U^T M U), hence q_M(t)=n det(tI_(n-1)-U^T M U). Write B=UR. Since det(R)^2=det(G)=n,
det(tG-K)=det(R)^2 det(tI-U^T M U)=q_M(t). QED.

Thus the singular/eigenvalue data of the alternating operator relative to g is already in the marked spectral pair (p,F=q/p). On the order8 conference switching family q_M(t)=8t(t^2+7)^3, so every restricted real skew operator is orthogonally equivalent: one zero eigenvalue and three equal 2-planes of frequency sqrt(7). Yet H is189 or333. The real metric pair loses the orientation's placement relative to the integral root lattice. That placement is exactly what §3 retains.

Similarly, dropping g and retaining only integral alternating congruence returns Session35's Smith-(1,1) order5 H=11/5 collision. Neither form separately repairs the loss; the simultaneous integral pair does.

## 5. Finite decorated-theta recovery, without scalar-theta sleight of hand

For g>=1 and x=(x_1,...,x_g) in L^g, retain BOTH matrices

Gamma(x)_ij=<x_i,x_j>, Omega_T(x)_ij=omega_T(x_i,x_j).

Let N_T^(g)(Gamma,Omega) count tuples with those two matrices. The count is finite when Gamma is fixed because its diagonal bounds every vector. These are the joint coefficients of a decorated theta series, for example

Theta_T^(g)(Z,W)=sum_(x in L^g) exp(pi i tr(Gamma(x) Z)+2pi i sum_(i<j) W_ij Omega_T(x)_ij),

with Im Z positive definite and W real. There is no assertion that this is a classical scalar Siegel modular form. The independent alternating variables are the orientation-sensitive sidecar. The metric Gram variable alone is blind here because L is the same A lattice for every tournament.

THEOREM (finite recovery at genus n-1). For n>=2, equality of the joint coefficients at the single Gram matrix G=I+J in genus n-1 is equivalent to isometry of the metric alternating lattices, and hence determines H and classifies strong tournaments.

Proof. Every tuple with Gram G is an integral basis of L: its Gram determinant equals detL=n, so its index is1. If T and U have the same coefficients, the coefficient at (G,K_T) is positive for T (use B), hence also for U. A U tuple with that Gram and alternating matrix gives an integral metric isometry carrying omega_U to omega_T. Conversely any isometry bijects all tuples. QED.

This uses only squared-norm2 roots. Equivalently retain the finite multiset of alternating matrices evaluated on all ordered star-root bases. For n>=3 there are2n! such bases, with the global sign duplicating each alternating matrix; at n2 there are2 bases, and n1 is trivial. This is a finite tensor/tuple depth upper bound, not a minimal-depth theorem for all n.

## 6. The complete two-root response is exactly (c3,c4)

THEOREM. For any tournament of order n>=2, put A=n(n-1)(n-2), and let c3,c4 count directed simple3- and4-cycles modulo cyclic rotation. For ordered norm2 roots x,y, the joint distribution of g=<x,y> and w=omega_T(x,y) is:

- g=+2 or-2: N(g,0)=n(n-1), separately for each g.
- g=+1 or-1: N(g,+3)=N(g,-3)=6c3, and N(g,+1)=N(g,-1)=A-6c3.
- g=0: N(0,+4)=N(0,-4)=4c4;
  N(0,+2)=N(0,-2)=A(n-3)/6+8c3(n-3)-16c4;
  N(0,0)=2A(n-3)/3-16c3(n-3)+24c4.
- All other counts are0.

Thus at fixed n the entire genus2 root-shell response is equivalent to the pair(c3,c4), not a stronger orientation invariant.

Proof. Adjacent roots give a signed triangle curl. Each cyclic triangle contributes six ordered pairs to each of the fixed g=+/-1 and w=+/-3 bins; the remaining adjacent-pair counts follow by root sign symmetry. A disjoint-root pair with w=4 is a directed4-cycle, and each cycle gives four such pairs. For the second moment, put r=M1 and P=I-J/n. Direct orthogonal decomposition gives

||PMP||_F^2=n(n-1)-2||r||^2/n=(n-1)(n-2)/3+16c3/n.

Here ||r||^2=(n^3-n)/3-8c3 follows from sum_i binom(outdegree_i,2)=binom(n,3)-c3. The root tight-frame identity sum_x xx^T=2nP gives total alternating second moment4nA/3+64nc3. The adjacent pairs contribute4A+192c3, leaving4A(n-3)/3+64c3(n-3) on disjoint supports. Subtract the w=+/-4 terms, divide the remainder by8 to obtain each w=+/-2 count, then subtract those from the total A(n-3) disjoint pairs to obtain w=0. The collinear cases g=+/-2 are immediate. Conversely the bins w=3 and w=4 directly recover c3 and c4. QED.

### Conference corollary: universal two-root blindness

THEOREM. Let M be a skew conference tournament matrix of even order n>=4, so M^2=-(n-1)I, and put A=n(n-1)(n-2). For ordered roots x,y of A_(n-1), the complete joint distribution of g=<x,y> and w=omega_T(x,y) is independent of M:

- g=+2 or-2: w=0 has n(n-1) occurrences for each g.
- g=+1 or-1: w=+3 and-3 each have A/4 occurrences; w=+1 and-1 each have3A/4, for each g.
- g=0: w=+4 and-4 each have A(n-4)/16 occurrences; w=+2 and-2 each have A(n-2)/4; w=0 has3A(n-4)/8.
- All other pairs have count0.

Proof. The score row-sum vector r=M1 satisfies ||r||^2=n(n-1). The elementary triangle identity gives c3=(n^3-n)/24-||r||^2/8=A/24. Root pairs of inner product+/-1 share exactly one endpoint; their alternating pairing is a signed triangle curl. Each cyclic triangle contributes six ordered pairs of each given g and each sign of3, yielding6c3=A/4. By sign symmetry, the remaining counts at each sign of1 are3A/4.

Disjointly supported roots are exactly the g=0 pairs. Pairing w=4 is a directed4-cycle on their four endpoints, with unrestricted diagonals; each such cycle gives four ordered root pairs, so N(0,4)=4c4. The conference adjacency characteristic polynomial is

lambda(lambda-(n-2)/2)(lambda^2+lambda+n/4)^(n/2-1).

Since there are no directed2-cycles, tr(adjacency^4)=4c4. Its displayed spectrum gives c4=A(n-4)/64, hence N(0,+/-4)=A(n-4)/16.

For the remaining counts use the tight root frame identity sum_(x root)xx^T=2nP, P=I-J/n. Consequently sum_(x,y roots)omega(x,y)^2=4n^2||PMP||_F^2. The conference compression has one zero eigenvalue and n-2 nonzero eigenvalues of squared absolute value n-1, so this is4n^2(n-1)(n-2)=4nA. The g=+/-1 cases contribute12A. Thus the g=0 second moment is4A(n-3). If B=N(0,4)=N(0,-4) and C=N(0,2)=N(0,-2), then32B+8C=4A(n-3). Solve to get C=A(n-2)/4. Subtract2B+2C from the total A(n-3) disjoint root pairs to get N(0,0)=3A(n-4)/8. The cases g=+/-2 are y=+/-x and have zero alternating value. QED.

This is an all-order theorem about the norm2 root shell, not about every coefficient of the full genus2 decorated theta series. It explains why one must pass beyond two-root incidence to recover the Hamiltonian split in the next example.

## 7. Exact root-shell depth in the Paley8 switching class

Take the128 switches defined in Session35 CONFERENCE_COLLAPSE. Genus1 decorated data on norm2 vectors cannot see omega at all. Genus2 norm2 joint data N^(2) is IDENTICAL for all128 switches. Its entire distribution is recorded in `verify_metric_curl.json`; H is still189 on16 members and333 on112.

Genus3 norm2 data separates those two H-values. One single coefficient is enough:

Gamma=[[2,-1,0],[-1,2,0],[0,0,2]],
Omega_12=-3, Omega_13=-4, Omega_23=2.

Its value is84 for all16 universal-vertex switches and36 for all112 nonuniversal switches. Within precisely this family, H=441-3N^(3)(Gamma,Omega). Therefore the first root-shell genus that can distinguish the two H-values is exactly3. This is a FINITE-EXACT statement, not a claim that full genus2 decorated theta series coincide, and not a general formula for arbitrary conference matrices or tournaments.

A combinatorial reading: the Gram condition uses a pair of roots sharing an endpoint and a third root supported on two other vertices, so it is a five-vertex incidence statistic. It detects how directed triangles are wired to external vertices, not merely how many triangles or simple5-cycles exist.

## 8. Recovering the lost cut, and its sharp boundary

No root marking that depends only on the fixed standard lattice can recover the cut: all labelled roots AND K are already identical throughout a nontrivial fixed-label fibre. In particular, adding just a marked root is insufficient even for transitive tournaments.

A marked vertex together with that vertex's score (outdegree) does recover the cut. If v lies in component C_j, then s_v equals its fixed internal score plus the total sizes of the later components. Across the k cyclic rotations these later-size sums are k different partial sums of positive component sizes. Thus only one rotation has a specified score. Intrinsically one may mark the dual simplex weight e_v-(1/n)1 and attach s_v. This statement is about the residual fibre; it is not an efficiency claim for constructing the candidates.

The same selected genus3 coefficient is not a universal H-reconstruction. On the Paley11-derived order12 conference switches with masks0 and1 it equals660 in both cases, while H is95095 and196815. This is an exact hostile to extrapolating the particular coefficient, not to the whole genus3 response.

## 9. Verification and scope

`python3 session36-lattice/verify_metric_curl.py` is the consolidated certificate. It checks every labelled tournament through n6, including complete K-fibres, exact fibre size=#SCC, H invariance, every SCC rotation and two-level cut differences. There are33867 tournaments, with distinct form counts1,1,4,38,728,26704. Hamiltonian subset DP is independently checked against full permutation enumeration through n5. The general two-root distribution is directly checked on every labelled tournament through n5.

It also checks all128 Paley8 switches with independent root-coordinate and five-vertex-pattern computations of the selected genus3 coefficient; root-star basis counts for n2..6;16 conference controls at n4,8,12,20;132 integer matrix-pencil evaluations; Session35's Smith-(1,1) hostile; and the order12 failure of the selected single coefficient. Exact JSON/output and normal/optimized replay receipts accompany the source. The general theorems are proved algebraically, not extrapolated from this census.

Connection contract: tournament -> integral metric alternating augmentation lattice. Preserved: exact triangle curl and all multiplicative SCC invariants, in particular H. Lost: a cyclic choice of cut between ordered strong components. Restoring sidecar: marked vertex score. Coarser quotient losses: Smith congruence discards root geometry; real orthogonal equivalence discards integral placement; undecorated theta discards omega entirely. The finite decorated coefficient family restores the simultaneous pair at genus n-1.


## 10. The horizon beyond this recovery

`FINITE_GENUS_OBSTRUCTION.md` proves that no fixed root-shell genus cutoff recovers H uniformly: for every cutoff, every sufficiently large order has strong tournaments with identical labelled scores and all those responses but different H. The proof recovers Moon's interval-chain product law and compares independent arithmetic H factors against the polynomial information budget of a fixed number of response bins. At the concrete genus2 level, IC(2,2,4) and IC(2,3,3) already have order9, identical scores and(c3,c4)=(7,10), but H=81 versus75.

Thus the complete integral metric pair is faithful on strong tournaments, while every uniformly bounded root-shell truncation eventually ceases to be faithful even to H. No statement about the full unbounded-norm theta series follows.

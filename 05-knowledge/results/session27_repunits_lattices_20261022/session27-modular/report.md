# Session27: typed modular, lattice and tournament connections

Status: independent synthesis; no new repository theorem IDs. Exact elementary propositions below are PROVED; imported structure is CITED; script checks are FINITE-EXACT. Repo read at dca0424aed4c, then relevant THM-4749 refreshed at 54bf9664293. No repository files changed.

## Executive results

1. There is a genuine chain: positive even unimodular lattice → theta modular form; bosonic oscillator Fock space → eta reciprocal; their product → lattice character. This gives E8 rank8 → 240 roots → Lie dimension248 via the explicit coefficient 8+240, and rank24 Leech → Theta/eta^24=j−720. No numerology is needed.
2. The native fixed-H tournament generating function is rational and cannot itself be a nonzero positive-weight holomorphic level-one modular form. This is a proved, precisely scoped obstruction, not a ban on every modular transform.
3. Prior work already supplies a concrete tournament-to-E8 map, T→score-deviation vector. At n=24 the same coordinate map lands in D24+, not automatically in the Leech lattice. It loses virtually all orientation data.
4. The historical 1729 reflection has a valid hostile against a split-prime pattern, but its blanket claim of no modular structure and its proposed e-limit are not consequences of that hostile.
5. Ordinary Construction A of the extended binary Golay code is NOT the Leech lattice: its theta series has 48 roots. Additional construction/gluing is essential.

## Inheritance and live concept board

Closest proved mechanisms: THM-868-e8-bridge-score-lattice.md (not the colliding geometric-ladder THM-868), THM-872-milgram-discriminant-residue.md, THM-489-eta-discriminant-bridge.md; fixed-H SCC polynomial law from Session24 and THM-4749 (repunit-degree-law file). Canonical hostile: E8⊕E8 versus D16+ have the same genus-one theta and different root systems. Corrected near miss: 1729 split-prime ratio fails at Paley p=19; historical claims about universal nonmodularity do not follow. Least-used sidecar: marked score hyperplane/half-integral coset, which the full theta series forgets.

Board: (i) fixed-H SCC arithmetic, (ii) score lattice, (iii) theta/eta, (iv) Lorentzian null quotient/deep holes, (v) CM j-values. Anchor: classify lawful invariant transfers, without implying LRC(14). Niche: rational-modular obstruction. Wildcard: 1729's distinct cubic and Carmichael mechanisms. No conclusion settles LRC.

## 1. Normalizations and exact dictionary

Set q=exp(2πiτ), Imτ>0. Define theta_L=Σ_{v∈L}q^{(v,v)/2}. Then eta(τ)=q^(1/24)Π_{m≥1}(1−q^m), and Delta=eta^24=qΠ(1−q^m)^24. The q^(1/24) is not a spatial dimension. It is the normalization that produces eta's multiplier eta(τ+1)=exp(2πi/24)eta(τ); for a single chiral boson its reciprocal supplies the vacuum energy shift −1/24.

Normalized Eisenstein series have constant term one:
E4=1+240Σσ3(n)q^n; E6=1−504Σσ5(n)q^n.
Since their weights are4 and6, E4³−E6² is a weight12 cusp form. Its first coefficient is 720−(−1008)=1728; normalization gives Delta=(E4³−E6²)/1728. The factor1728 is fixed by normalization, not derived from the nearby integer1729. The one-dimensional spaces of weights8,10,14 give E8=E4², E10=E4E6, E14=E4²E6. Here E8 is an Eisenstein series, not the E8 lattice or Lie group. The lattice theta_E8 is E4, not E8.

Classical modular forms can be viewed algebraically as sections of powers of the Hodge line bundle over elliptic-curve moduli; Delta detects degeneration at the cusp. Fractional eta is a multiplier/line-bundle object, not simply a weight-zero algebraic function of q. Its 24th power removes the multiplier. j=E4³/Delta has j(i)=1728, and j(exp(2πi/3))=0. Avoid mixing normalized Delta with the analytic Weierstrass discriminant, which carries (2π)^12 factors.

Sources: Milne, Modular Functions and Modular Forms, §§4,8,12, https://www.jmilne.org/math/CourseNotes/MF.pdf (definitions, structure and algebraic interpretation); Tong, String Theory, §6.4, https://www.damtp.cam.ac.uk/user/tong/string/string.pdf (eta and oscillator normalization). These are classical imported identities, not newly proved from finite prefixes.

## 2. E8, 248, 24 and the exact character bridge

PROVED from the standard coordinate root model: E8 has 112 integer roots (two coordinates ±1) and128 half-integral roots (eight coordinates ±1/2 with even parity), hence240 roots. The Cartan decomposition of its complex simple Lie algebra adds rank8 to one dimension per root: dim e8=8+240=248. Root-datum reconstruction supplies the Lie algebra; the scalar rank alone does not.

A stronger graded identity explains the same248:
 theta_E8/eta^8=q^(−1/3)(1+248q+4124q²+…).
The q coefficient is exactly240 lattice roots plus8 oscillator states. This is a genuine representation-theoretic join, not a coincidence between an eight-dimensional lattice and a248-dimensional Lie algebra. The denominator exponent is rank/central charge8, not248.

For a positive-definite even unimodular rank24 lattice N, theta_N is weight12 and
 theta_N=E4³+(R_N−720)Delta,
where R_N counts norm2 roots. Therefore
 theta_N/eta^24=j+R_N−720.
For the rootless Leech lattice, theta_Λ=E4³−720Delta=1+196560q²+16773120q³+…, and theta_Λ/eta^24=j−720=q^(−1)+24+196884q+… . In particular j−744 (the normalized monstrous J) is NOT the Leech lattice character: it differs by24. A further orbifold construction is not the identity map.

For D24+, R=2·24·23=1104, so theta=E4³+384Delta. For E8³ and D16+⊕E8, R=720 in both cases, so both theta series equal E4³ although their root decompositions differ. Thus even complete genus-one theta coefficients cannot classify Niemeier lattices. This is a statement about ambient genus-one theta only: do not infer an assertion about higher-genus theta, its Jacobian restriction, or hyperelliptic restriction.

The eta/code bridge THM-489 has W_Golay=W_Hamming³−42P24 and theta(P24)=16Delta in its stated convention, hence ordinary Construction A gives theta=E4³−672Delta, with48 roots. Those are the coordinate ±sqrt(2)e_i roots, present in ordinary binary Construction A. This is the A1^24 Niemeier case, not the rootless Leech case.

The E8 dimension root decomposition is independently documented in Kostant's lecture summary https://math.ucr.edu/home/baez/kostant/summary.html. General lattice/vertex-algebra construction: Borcherds, Vertex algebras, Kac-Moody algebras, and the monster, linked at https://math.berkeley.edu/~reb/papers/index.html . The displayed coefficient identity is also derived directly by multiplication in check.py, without requiring a VOA theorem to check it.

## 3. 24 versus26, II25,1 and deep holes

CITED: critical bosonic string theory has26 spacetime coordinates and24 transverse oscillator polarizations in light-cone gauge. Its chiral oscillator factor is eta^−24=Delta^−1, whose coefficients count24-coloured partitions after the q^−1 shift. This does not say physical spacetime is24-dimensional; superstring critical dimension10 is a separate theory. Tong §6.4 gives the full closed-string integrand, including conjugate/zero-mode/measure factors; 1/Delta alone is not that full partition function.

A Euclidean even unimodular lattice has rank divisible by8. For an indefinite even unimodular lattice, signature difference p−q must be divisible by8. II25,1 has rank26 and signature difference24: it is Lorentzian, not a26-dimensional positive-definite packing lattice. With U the hyperbolic plane, II25,1≅Λ⊕U. For primitive isotropic v, v-perp/Zv is a positive even unimodular rank24 lattice. Isometry orbits yield the24 Niemeier classes. The Leech's deep holes form23 affine-isometry types, corresponding to the23 rootful Niemeier classes, not24 deep-hole types. Its covering radius is sqrt(2) when minimum squared norm is4.

Primary sources: Borcherds, The Leech lattice, https://math.berkeley.edu/~reb/papers/leech/leech.pdf ; Borcherds–Conway–Queen, The cellular structure of the Leech lattice, Theorem1, https://math.berkeley.edu/~reb/papers/shallow/shallow.pdf . The first also provides explicit simple roots (λ,1,λ²/2−1) in the Lorentzian model. This is a precise geometric construction connecting26 and24.

Sphere packing optimality is an additional theorem, not a consequence solely of rootlessness or Niemeier classification. Cohn–Kumar–Miller–Radchenko–Viazovska prove optimality in dimension24 and uniqueness among periodic packings: https://annals.math.princeton.edu/2017/185-3/p08 . Do not silently strengthen that uniqueness claim to all arbitrary nonperiodic configurations.

Langlands duality exchanges character/cocharacter lattices and roots/coroots in a root datum. E8 is simply laced with unimodular root lattice and is self-dual at this level. That is not the same operation as Euclidean lattice duality L→L*, despite agreement in this case. The dual group alone is not a Langlands correspondence, nor does it identify all automorphic forms. Source for precise dual-root-datum definition: https://www.math.columbia.edu/~chaoli/docs/AutomorphicForm.html .

## 4. The actual tournament score map

For an n-vertex tournament, scores s_i satisfy Σs_i=n(n−1)/2. Let v_i=s_i−(n−1)/2. For n divisible by8, v lies in the half-integral coset of D_n+: every coordinate is half-integral, total sum is0, and the standard even-unimodular D_n+ gluing criterion is met. For n=8 this is E8; for n=24 it is D24+.

Preserved: score norm, degree variance, the half-integral coset and Σ=0 hyperplane. Lost: arcs, Hamiltonian-path count, automorphism group, most cycle information. Needed sidecar: actual tournament/realization multiplicity and Landau inequalities. Cheap test: at n=8 the score floor consists of70 balanced half-roots, whereas E8 has240 roots in total; the full theta is not the score enumerator. At n=24 a score half-vector has squared norm≥24/4=6; the Leech's first nonzero shell has squared norm4. This forbids treating the same marked score vectors as the complete Leech shell system, but does not prohibit some unrelated embedding or transform.

Prior files: 01-canon/theorems/THM-868-e8-bridge-score-lattice.md; THM-869-milgram-corona-sedenion.md; THM-872-milgram-discriminant-residue.md. Their speculative Cayley–Dickson rhetoric is not needed for these elementary membership statements.

## 5. New scoped proposition: rational native generators are not modular forms

PROVED. Let R(q) be nonzero rational, holomorphic on |q|<1. If f(τ)=R(exp(2πiτ)) is a holomorphic modular form of positive integral weight k for SL2(Z), a contradiction follows.

Proof. Write a0=R(0). If a0≠0, for every reduced a/c choose γ∈SL2(Z) with γ∞=a/c. Along τ=γ(iY), Y→∞, modularity gives f(τ)=(ciY+d)^k(a0+O(exp(−2πY))). Thus f is unbounded near exp(2πia/c). A rational function bounded near any point that is not a pole must have a pole at each such root of unity; infinitely many poles are impossible.
If a0=0, let m≥1 be the first nonzero Taylor exponent. Along τ=it, t→0+, the S-law gives f(it)=(it)^(-k)f(i/t), of order t^(−k)exp(−2πm/t), up to a nonzero constant. But R(exp(−2πt)) has a finite Laurent order at t=0 if R is nonzero rational. Finite-order behavior cannot equal exponential vanishing. Contradiction. ∎

This proof uses only modular transformation and rational local behavior. It does not assume eigenform Hecke relations or positivity. Weight0 constants are excluded by the positive-weight hypothesis; zero generators such as h=7,21 are excluded by nonzero. It makes no statement about weakly holomorphic, mock, quantum, vector-valued, transformed, or two-variable generating functions.

Application: fix h attainable as H(T), and let N(n,h) count isomorphism classes. SCC decomposition and the inherited finite strong-component size bound give
 F_h(q)=Σ_n N(n,h)q^n=Σ_{ordered strong sequences A1,…,Ar; product H=h} q^(Σ|Ai|)/(1−q)^(r+1).
The sum is finite; singleton runs occupy r+1 slots. For the empty sequence h=1 this is1/(1−q). Hence every attainable fixed-h native generator is nonzero rational and is not a positive-weight holomorphic level-one modular form. Examples: F3=q³/(1−q)² and F63=2q⁸/(1−q)² (the latter inherits the complete two-atom inventory). The repunit lane's symmetry-weighted W_h inserts factors1/Π|Aut Ai| and remains rational, so the same obstruction applies. It also applies to the exponential generating function of labelled tournaments at fixed h, since that is W_h.

Source for finite SCC law and bounds: current THM-4749 (audited MISTAKE-657), and 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md. This is a new local consequence of inherited mechanisms, with no claim of literature novelty. It does not require conjecturing that only7,21 are excluded for every h.

## 6. 1729 and CM: distinguish mechanisms

FINITE-EXACT plus elementary proof:1729=7·13·19=1³+12³=9³+10³. Exhaustion of all positive cube pairs up to1729 verifies its minimal two-representation property. It is Carmichael: squarefree and each p−1 (6,12,18) divides1728; Fermat modulo each prime and CRT give a^1728≡1 mod1729 for all1296 units. These are real arithmetic structures. Neither implies a string-theory relation.

Prior tournament ratio: H(Paley11)/|Aut|=95095/55=1729. This has a native quotient interpretation because an automorphism fixing an ordered Hamiltonian path fixes every vertex, so the action on paths is free. The p=19 ratio in the prior computation contains inert primes5,11,23 relative to Q(sqrt(−3)), refuting the proposed universal split-prime property. It does NOT show every modular interpretation impossible, prove convergence to e, or establish Paley as universal Hamiltonian-path maximizer. Historical prose asserting these should be scoped down.

The Heegner radicands are exactly the user's {1,2,3,7,11,19,43,67,163};1 belongs. Their fundamental discriminants are−4,−8,−3,−7,−11,−19,−43,−67,−163. This is a classification of imaginary quadratic fields with class number one, not a list of prime dimensions. CM attaches j-values to elliptic curves with endomorphisms by their integer rings. No canonical map from tournament path orbits to those CM elliptic curves has been supplied. For each proposed map one must specify the elliptic curve or period, prove independence of choices, and check more than one scalar coincidence. In particular j(i)=1728 concerns discriminant−4, not the factorization1729 in Q(sqrt(−3)). The CM mechanism is sourced by Milne §12; completeness of the nine-value classification is a classical imported theorem, not recomputed here.

## 7. Proposed {7,21}→{8,24} bridge: cheap hostile first

The arithmetic mapping7→8,21→24 is multiplication by8/7 on the two chosen values, but no tournament operation, lattice functor or invariant-preserving transformation is thereby defined. Adding1 gives8,22, not8,24. If the intended map is prime-field size→extended quadratic-residue-code length q+1, then7→8 and23→24, not21→24;21 is not a prime power. Conversely n=8,24 score-lattice membership above is lawful but does not characterize forbidden H values7,21, because it uses vertex order n, not H. This is a type obstruction, not proof that every future relationship must fail.

## Reproducibility, results and stopping conditions

Run python3 check.py in this directory. check.json records all exact integer outputs. It tests Delta via two independent definitions through q^60; the E8/E10/E14 identities; reciprocal oscillator coefficients; theta candidates;248; all positive cube representations≤1729; and all1296 Carmichael units. Higher coefficients confirm consistency, not standalone all-order proofs or lattice existence. It does not simulate string theory, enumerate all Niemeier lattices, or calculate actual Leech vectors.

Coordination status: startup packet was read-only; message processor could not run because .machine-id was absent. No identity was invented. No source/Git changes were made. Retain the no-go, normalization controls, and score map; stop numerical coincidence hunting unless a proposed source-target map and preserved predicate are supplied.

### Primary source pin for Heegner completeness

H. M. Stark, *There is no tenth complex quadratic field with class-number one*, PNAS57 (1967),216–221, DOI10.1073/pnas.57.2.216; original paper archived at https://pmc.ncbi.nlm.nih.gov/articles/PMC335491/ . This pins the nine-field completeness input separately from Milne's CM dictionary.
Access qualification: the search index verified Stark's bibliographic record; opening the PMC full text encountered a browser/CAPTCHA check. No CAPTCHA was completed and the proof itself was not reread. Accordingly this source pin is bibliographic, not a fresh full-proof audit.

## 8. Fresh positive arithmetic bridge: Luca2026 → Camion criterion

CITED (Theorem1, statement and full HTML proof inspected): Florian Luca, *Yet another curious property of1729*, Ramanujan Journal70, article74, published19August2026, https://link.springer.com/article/10.1007/s11139-026-01450-1 . For prime p and positive integer n, if2^n p³+1 is Carmichael, then p=3. The example1729=2^6·3³+1 establishes occurrence. This does not say every exponent at p=3 works or that1729 is the only such number. The theorem restricts the prime-cube odd part; it is not a classification of general p^q r forms.

PROVED composition with inherited THM-4746: write s(m) for strong labelled tournaments, c(m)=2^((m−1)/2) for odd m, D=s−c. The inherited quasi-period law D(m+p−1)≡(2/p)D(m) modulo odd prime p and D(1)=0 show that every Carmichael N passes D(N)≡0 modulo N: each p|N has N≡1 mod(p−1), and squarefreeness permits CRT. Thus1729 is a Camion pseudoprime with s(1729)≡c(1729)≡1 modulo1729. Direct independent Wright-recurrence computation in check_camion.py confirms both residues. The map is Carmichael integers→the Camion congruence locus; it preserves the passed congruence, forgets universal-base Fermat information, and is not reversible:21 passes with s(21)≡c(21)≡16 mod21, yet2^20≡4 mod21, so21 is not Carmichael. This is an exact bridge to tournament-count arithmetic, not to the1728 normalization of Delta. The modular no-go in §5 concerns fixed-H class inventories, a different sequence from s(m).

Inherited source: 01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md, §§1,4,5. Reproduce with python3 check_camion.py; output check_camion.json. All computation uses integer residues, not floating approximations or a giant exact strong count.

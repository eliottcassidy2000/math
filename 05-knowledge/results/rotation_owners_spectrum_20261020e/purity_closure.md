# Zero black defect closes the finite pure-blue atom window

Research continuation, October 20, 2026.

**PROVED composition of imported results, with their CITED inputs and
FINITE-EXACT certificates:** the strong pure-blue tournament classes are
exactly the singleton, the directed triangle C3, and the regular cyclic
pentagon C5. The two imported packets separately leave an apparently
unresolved finite window; their bounds and certificates together close it.
This is an independent proof route, not a priority claim. It neither
retracts nor claims to supersede the separate
[uniform-prefix classification](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md).

The consequence is the exact all-order ordinary class generating function

    B(z) = z+z^3+z^5,
    P(z) = (1+B(z))/(1-z^2-z^6).

B counts strong pure-blue classes, including the singleton. P counts all
pure-blue classes, including the empty tournament. Neither series counts
labelled tournaments, Hamiltonian paths, path orbits, or all tournaments.
This does not establish black-core connectivity or the ordinary H spectrum.

## 1. Definitions and dependency boundary

For a strong tournament T of order n>=3, let H(T) be the number of directed
Hamiltonian vertex words. The automorphism group A=Aut(T) acts freely on
these words: fixing a complete ordered word fixes every vertex. Write

    h = H(T)/|A|,
    b = number of blue A-orbits of Hamiltonian words,
    Delta = h-b.

A word P=(p_0,...,p_(n-1)) is blue precisely when its positional reversal
rho_P(p_i)=p_(n-1-i) is an anti-automorphism of T. This property is invariant
under A. Consequently Delta is exactly the number of black ordinary path
orbits, and T is pure-blue if and only if Delta=0. Delta is not the feedback
deficiency used in owners.md, and no paired-black graph degree is silently
substituted for it.

The inheritance chain read and audited for this composition is:

* **Imported PROVED bound:** commit 946e63daee,
  [effective black-defect order bound, theorem and Sections1--5](../effective_black_defect_order_bound_20261020/PROOF.txt).
  For integer d>=0 with Delta<=d, put s=d+2 and let r be the least power of
  two at least 4s(d+1). Then n<=2sr. Its only external theorem is CITED
  Moon vertex-pancyclicity, which also supplies Hamiltonicity.
* **Imported reduction, with CITED input:** at commit 8313f3db06,
  [strong cycle reduction, Sections1--2](../merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md).
  A nontrivial pure-blue strong tournament has odd order and is a regular
  circulant when labelled along any Hamiltonian cycle. The even-order
  exclusion uses CITED Lu2020, Theorem7; its exact access and convention
  are recorded in [Lu's source note](../merged_metagraph_purity_and_strong_core_20261020/LU_CROSSING_SOURCE.md).
* **Imported FINITE-EXACT certificates:** commit 8313f3db06,
  [finite atoms, Section5](../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md),
  its main cycle certificate JSON and all eight extended gzip JSON files.
  Every skew circulant presentation at each odd order7 through31 has a
  verified nonblue Hamiltonian word.
* **Imported PROVED grammar:** [THM-2454, pure-blue alphabet and center law,
  Section2](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md),
  together with the all-order rigid-strong classification in
  [finite atoms, Section3](../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md).
  The latter uses the explicitly CITED Benhocine--Wojda bypass theorem.

No Cuckler asymptotic, Dixon group bound, Feit--Thompson theorem, or
Thomassen cycle-count estimate is needed for this new composition.
The cited published inputs have not been reproved here. This pass read
their inherited precise statements and source-access records.

## 2. Audit of the zero-defect cutoff

Setting d=0 in the imported bound gives s=2, r=8, hence n<=32.
The bound does not assume d>0 or nonpurity. Its zero case has the following
direct specialization of the proof.

Choose a Hamiltonian cycle and let q be the number of ordinary path orbits
met by its cuts. If the automorphisms among its cyclic rotations have
index q, the blue cuts occupy at most two of these orbits: composing cut
reversals gives rotation by twice their displacement. Thus Delta=0 forces
q<=2. Each rotation orbit has size m=n/q and induces a regular tournament.
The imported Section2 proves, with the two-orbit boundary included,

    minimum outdegree >= n/4-1/2,
    kappa(T) >= ceil(n/3) if q=1,
    kappa(T) >= ceil(m/2) if q=2.

Assume n>32. Both alternatives give kappa(T)>8. Moon's pancyclicity gives
a directed eight-cycle on a set S, and deleting S leaves a strong
tournament. If every outside vertex were beaten by at most one vertex
of S, summing outdegrees over S would give

    2n-4 <= 28+(n-8) = n+20,

and therefore n<=24, a contradiction. Some outside vertex v is beaten by
two distinct vertices of S. Start a Hamiltonian word W of T-S at v,
using a cut of a Hamiltonian cycle of that strong remainder. Cut the
eight-cycle immediately after either of those two vertices and append W.
This gives two Hamiltonian words P,Q.

The positional map u with Q=uP is a nonidentity power of an eight-cycle
on S, fixing everything else. It has even order and hence cannot be a
tournament automorphism. The words therefore lie in different ordinary
path orbits. If both were blue, then

    rho_Q rho_P = u (rho_P u^(-1) rho_P)

would be an automorphism. The two factors have disjoint supports in the
first and last eight positions of P, since n>32. Both have the same even
order, so their product has that even order. An even-order tournament
automorphism would have a nonidentity involutory power, which reverses
the arc of a transposed pair. This is impossible. At least one word is
black, contradicting Delta=0.

Thus n<=32, including every composite-order and two-orbit case. This
argument is a proof, not an extrapolation from the certificate census.

## 3. Odd circulants, complete coverage, and small orders

For a pure strong tournament, index a Hamiltonian cycle by Z/nZ.
Every cut reversal r_i(j)=2i-1-j is an anti-automorphism, so
r_(i+1)r_i is rotation by2 and is an automorphism. If n is odd, rotation2
generates rotation1: the cycle labels identify T with a circulant.

If 4 divides n, rotation2 has even order, which is impossible. For the
remaining even case n=2k, k odd and n>=6, rotation2 has two regular vertex
orbits U,W. Their constant scores a,b satisfy a+b=n-1 and a!=b, because
a cut reversal exchanges the orbits. Every Hamiltonian path has endpoint
score sum n-1, so its endpoints must lie in different orbits. Lu's cited
crossing theorem excludes such a strong tournament at order>=6. There is
no strong order2 tournament. This is the exact imported parity argument.

Together with n<=32, this leaves only odd orders3 through31.
For each such n, every skew connection set S on Z/nZ is represented by
one binary choice from each pair {d,n-d}, 1<=d<=(n-1)/2. There are exactly
2^((n-1)/2) masks. No primality, isomorphism filter, randomness, or condition
1 in S is imposed. In particular every possible cyclic generator and
every presentation obtained from a pure tournament's Hamiltonian cycle
is included.

The independent checker reads certificate data only, authenticates each
data file against the imported SHA256 manifest (LF-normalized text, exact
binary gzip bytes), and requires every mask
exactly once at every specified order. At orders7 through31 it checks a
directed Hamiltonian cycle, a genuine cyclic cut P, and an explicit pair
whose orientation is not reversed by rho_P. These 65,528 certificates
exclude purity directly; no converse of the cycle reduction is used.
The successor's failed automorphism pair is independently checked too.

At orders3 and5 the checker does not merely accept absent negative
certificates: it enumerates every Hamiltonian word in all six skew
presentations and verifies every reversal. Their H values are3 and15.
The two order3 presentations give C3; all four order5 presentations lie
in the multiplier orbit of S={1,2}, giving the regular cyclic C5.
The singleton is pure separately. Therefore the complete strong alphabet
is exactly {1,C3,C5}.

The imported verifier used Python assertions. This new checker uses
explicit exceptions throughout, so its complete acceptance conditions
remain active under python -O.

## 4. The symmetric SCC grammar fixes the full numerator

Every tournament has a unique ordered strong-component word. For a
pure-blue tournament, converse reversal pairs its outer components; the
middle is either absent or one pure-blue strong component. If the ordinary
path-orbit counts on the left half are t_i and the center has counts h_0,b_0,
the inherited center law gives

    total path orbits = (product t_i^2) h_0,
    blue path orbits  = (product t_i) b_0.

Use h_0=b_0=1 for an absent center. All t_i>=1 and 0<=b_0<=h_0.
Equality of these products holds precisely when every t_i=1 and b_0=h_0.
Thus the halves consist of rigid strong atoms and the center is pure-blue.
The imported all-order rigid-strong classification says these atoms are
exactly singleton and C3. Both are self-converse, so the right half is
the reversed left word without an additional choice.

An arbitrary left word in sizes1 and3, counted with both halves, therefore
has generating function 1/(1-z^2-z^6). The optional center contributes
1+z+z^3+z^5. Unique ordered SCC decomposition prevents multiple counting.
It follows that

    B(z)=z+z^3+z^5,
    P(z)=(1+z+z^3+z^5)/(1-z^2-z^6).

If c_0=1, c_j=0 for j<0 and c_j=c_(j-1)+c_(j-3) for j>=1, then

    p_(2m)=c_m,
    p_(2m+1)=c_m+c_(m-1)+c_(m-2).

Equivalently p_n=p_(n-2)+p_(n-6) for n>=6, with
(p_0,...,p_5)=(1,1,1,2,1,3). All pure-blue classes are self-converse,
so ordinary versus converse-merged class counts agree on this subset.
The empty and orders1,2 here are class conventions, not claims about
the exceptional small-order complement-line conventions.

## 5. Scope, controls, and reproduction

The scalar all-order count now has a completely known numerator.
This closes the residual of the finite-atom packet by composition and
gives an alternative to the separate uniform-prefix classification.
It does not turn within-fibre relocation into a black metagraph edge,
join different black terminal basins, or prove all-order black-core
connectivity. The regular-endcap bypass in regular.md remains a separate
statement with a different target.

Run:

    python 05-knowledge/results/rotation_owners_spectrum_20261020e/purity_check.py
    python -O 05-knowledge/results/rotation_owners_spectrum_20261020e/purity_check.py

Both modes must produce the same receipt. The finite universes are all
65,534 skew presentations at odd orders3 through31, all Hamiltonian words
of the six positive order3/5 presentations, and the explicit palindromic
SCC-word grammar through total order32. The last check verifies uniqueness,
the rational-series coefficients, and the prior order3..13 counts; the
all-order identity follows from the preceding grammar proof.

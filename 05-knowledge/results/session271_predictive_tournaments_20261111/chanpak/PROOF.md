# Occupied-slot phases: a predictive classification of fixed-core growth

**Status: PROVED IN NOTE** for the stated carrier and iid/periodic regimes;
**FINITE STOCHASTIC PROBES** for numerical tables; **OPEN** for arbitrary-target
inversion and infinite Chan–Pak. No external novelty claim or global tournament
classification is made. This packet inherits, rather than rechecks, Session270's
36-state positive transfer and minimal 17-dimensional rational representation.

## Inheritance, portfolio and concept board

Closest mechanisms: Session270 `chanpak/PROOF.md` (slot transfer and averaged
cover transform), Session270 `carrier/PROOF.md` (marked class count), Session268
`spectrum/PROOF.md` (full q-vector and sparse product obstruction), Session268
`automata/PROOF.md` (rooted seams and genuine interactions). The old rank17
lineage is THM-4878, `01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md`.

Anchor: large H and witness order. Niche: saturation of occupied slots rather
than generic random-matrix asymptotics. Wildcard: a tropical support invariant
coupled to exact residue automata. Six live concepts are occupied-slot cone,
quenched versus annealed growth, scalar range, q tomography, modular inversion,
and symmetry/SCC sidecars. Hostiles are constant-color words, treating a mean as
coverage, cancelling signed 17-state coordinates, and treating marked words as
ordinary classes. The corrected near miss is a generic lognormal model:
log H has bounded, not linearly growing, variance in every fixed iid carrier
below. This makes endpoint geometry more important than bulk noise.

## 1. The support-rank invariant

Fix a tournament C on k>=1 core vertices. A word w over colors {0,...,2^k-1}
adds an ordered transitive spine; bit j of a letter is one exactly when core j
beats that spine vertex. Let H_C(w) count directed Hamiltonian paths. The
inherited slot representation considers each permutation p of the core and
its k+1 intervening spine slots. Vertices assigned to a slot occur in increasing
spine order. Only the first and last vertex of each occupied slot constrain
cross orientations.

Let S be a nonempty set of allowed letters. For a permutation p, declare:

* slot 0 allowed if some c in S has bit p[0]=0;
* slot k allowed if some c in S has bit p[k-1]=1;
* internal slot j allowed if some c has bit p[j-1]=1 **and** some (possibly
  different) c has bit p[j]=0.

A subset A of allowed slots is feasible when every empty internal slot j has
the core edge p[j-1] -> p[j]. Define

    r(C,S) = max |A| over feasible (p,A).

Every nonempty S has r>=1, since a tournament with any nonempty spine has a
Hamiltonian path and its occupied slots give a feasible subset. A maximizing
subset can include **all** allowed slots: adding an allowed slot only removes
an empty-core-edge requirement. Computing r costs O(k! k |S|) using this fact;
the supplied small implementation deliberately enumerates all subsets too.
This invariant is combinatorial, not an empirical matrix eigenvalue.

### Theorem 1: bounded-endpoint growth

There is a constant B<=k! 2^(k+1), depending only on k, such that

    1 <= H_C(w) <= B r(C,S)^m                 (w in S^m, m>=1).

There also exist words u,v in S*, each of length at most k+1, for which

    H_C(a u z v b) >= r(C,S)^|z|             (a,z,b arbitrary words over S).

**Upper proof.** Each path specifies a core permutation and its occupied slot
subset. The subset is feasible over S. For that pair there are at most |A|^m
assignments of m spine vertices to the occupied slots. Summing over at most
k! 2^(k+1) pairs proves the bound. The sharper universal bound
H_C(w)<=k!(k+1)^m follows without summing over occupied subsets.

**Lower proof, including arbitrary outside words.** Choose a maximizing
permutation and let A be all its allowed slots. Put the letters of a into the
earliest active slot. This is always a legal first entrance: if slot0 is
active it has no previous core; otherwise its failure forces bit p[0]=1 for
all S. Inductively, each successive inactive internal slot forces the next
core bit to be identically one until the earliest active slot is reached.
Its previous-core bit is therefore always one.

Construct u with one letter per active slot, in any fixed order, choosing
its previous-core bit to be one when needed. This fills every active slot;
slot0 can use any letter. All vertices of z now have r independent choices
of an active slot. Construct v with one letter for each active slot except
slot k, choosing its next-core bit to be zero. This fixes every required
last-edge condition. Put b in the latest active slot. Symmetrically, if that
slot is not k, inactivity of later slots forces its next-core bit identically
zero, so b cannot spoil acceptance. Empty internal slots have the required
core edges by feasibility. Distinct assignments of z give distinct paths.
The argument does not require u,v to be contiguous in a Hamiltonian path;
they are contiguous in the scanned spine word. QED.

For full support S all slots are allowed for every core permutation, so
r=k+1, independently of the orientations of C. For k=2 one may use u=333
and v=00; the sharper full-support sandwich is

    3^|z| <= H(a 333 z 00 b) <= 2*3^(|a333z00b|).

## 2. Exact random and periodic phase classification

### Theorem 2: iid colors

Let W_1,W_2,... be iid with probabilities p_c>0 exactly on S. Set r=r(C,S).
Then, as m grows:

    log H(W_1...W_m)/m -> log r almost surely;
    E log H(W_1...W_m) = m log r + O(1);
    Var(log H(W_1...W_m)) = O(1);
    lim (1/m) log E[H(W_1...W_m)^s] = s log r, for every real s>0.

Constants can depend badly on k and the probabilities. These are fixed-core,
fixed-law limits, not uniform statements as k grows or probabilities vanish.

**Proof.** Waiting for u near the left end and v near the right end has
geometrically decaying tails: test disjoint blocks of the appropriate fixed
lengths. The two waiting times need not be independent for the union bound.
If the occurrences are disjoint, Theorem1 gives
log H >= (m-T_left-T_right)log r. If they overlap, H>=1 suffices. Thus the
deficit D_m=m log r-log H is bounded below by -log B and satisfies a uniform
exponential upper-tail bound P(D_m>x)<=K exp(-c x), for 0<=x<=m log r,
with adjusted constants covering short words. For r=1 the upper bound itself
makes everything immediate. Integrating this tail proves bounded first and
second moments of D_m; Borel–Cantelli gives D_m=O(log m) almost surely.
For positive moments, the upper bound is B^s r^(sm). The event that u is the
literal prefix and v the literal suffix has fixed positive probability,
independent of m once disjoint, and contributes at least a constant times
r^(sm). This proves the pressure formula. QED.

Consequences: quenched and annealed *exponents coincide*. There is no nonzero
sqrt(m)-scale central-limit variance for log H in this ensemble. The theorem
does not assert convergence of H/r^m to a constant or even almost surely;
changing suffixes retain O(1) random fluctuations. Finite-sample means can
still be unstable when probabilities are very skewed.

### Theorem 3: periodic words

If a finite block b has letter support S and W=b^infinity, then
H(W[1:m])=Theta(r(C,S)^m). To see this, initialize each active slot using a
suitable letter in a bounded number of periods, allow arbitrary assignments
in the middle, then fix each final constraint in a bounded number of terminal
periods. Theorem1 gives the upper bound. Consequently the output-visible
spectral growth per letter is the integer r. This does not say every eigenvalue
of a period matrix is an integer or that every transient block has this radius.

`TABLES.md` gives all 15 nonempty supports for the two-core language. There
are exactly three growth phases, r=1,2,3. In particular S={2} already has base2,
S={0,1} has base2, and S={0,3} or {1,2} has base3. Entropy of the letter source
alone cannot predict H growth: a deterministic color2 word grows exponentially,
while other deterministic colors give H=1.

## 3. Spectral, projective, tropical, and rare-event interpretations

The positive 36-state machine is stratified by the subset of occupied slots.
Occupancy only grows. In a saturated feasible stratum with r occupied slots,
each letter has exactly r outgoing choices per state, so after division by r
the internal transfer preserves column mass. Start and finish constraints are
boundary layers, not a source of bulk Lyapunov noise. The integer r is a
max-over-feasible-strata, or tropical, growth invariant. This explains the
spectral phase classification without diagonalizing products of matrices.

Projectivize the positive cone on each stratum. Its mass functional is exact,
and initialization/acceptance words connect it to the observable. The cone
need not be strictly positive or uniformly contracting in Hilbert's metric;
zero entries and disjoint permutation sectors persist. Generic strict-positive
matrix theorems are unnecessary. The rational17 quotient is still excellent
for algebraic equality and future-context compression, but its signed
coordinates conceal this cone. Use the 36-state realization for inequalities;
transport the cone to the quotient only with its reachable-lattice sidecar.

For full-support iid colors, let P(s) be a pressure limit when it exists.
Theorem2 gives P(s)=s log(k+1) for all s>0. This is an unusually flat positive
moment spectrum, not evidence that all H values are common. In the two-core
uniform law the all0 and all3 words have H=1 and probability4^-m each. Thus
for s<0,

    liminf m^-1 log E H^s >= max(s log3, -log4),
    limsup m^-1 log E H^s <= 0.

For nonuniform full support replace -log4 by max(log p0,log p3); additional
low-growth supports can strengthen the bound. More generally, conditioning
all letters into A subset S costs P(A)^m and yields the rigorous lower bound

    liminf m^-1 log E H^s >= log P(A) + s log r(C,A),    s<0.

Here the upper bound H<=B r(C,A)^m gives the lower bound after raising to a
negative power. Maximize over A. This is a finite variational lower envelope
for rare-event pressure, not an equality theorem: temporal ordering patterns
can beat every support-only mechanism. Negative moments and lower large
-deviation tails are the natural place to seek genuinely nonlinear pressure.

## 4. Predicting H scale and preserving the q signal

For uniformly random colors, the inherited averaged-cover transform gives

    E H_C(W_m) = (k!/4^k)(k+1)^m + O(k^m),

because q_k(C)=k!. The leading coefficient and exponential base are universal
across all k-core tournaments. If two k-cores have largest differing cover
coordinate q_j, their mean difference has leading term

    (Delta q_j/4^j)(j+1)^m.

This follows from the coefficient of y^(j+1) in
(y+1)^2/4 * sum_r q_r((y-1)/4)^(r-1); larger equal coordinates cancel.
For same-order cores q_k and q_(k-1) are universal, so j<=k-2 whenever there
is a difference. The relative mean signal is exponentially suppressed by
at most ((k-1)/(k+1))^m. More random spine vertices improve the bulk growth
prediction while hiding the core's distinguishing q coordinates. Short
triangular tomography and long growth prediction are complementary observers.

### Theorem 4: biased mean coefficient depends only on bit marginals

For any iid color law let alpha_j=P(bit j=1). If every alpha_j is strictly
between zero and one, then r=k+1 even if many colors have probability zero, and

    E H_C(W_m) = k! product_j[alpha_j(1-alpha_j)] (k+1)^m + O(m k^m).

**Proof.** For each core permutation, choose an assignment of the m spine
vertices to its k+1 slots uniformly. Except for assignments having a slot of
size zero or one, all constrained first/last vertices are distinct. Independence
of colors at those vertices gives exactly product_j alpha_j(1-alpha_j),
regardless of correlations between bits of the *same* letter. The probability
of a slot of size zero or one is at most (k+1) times
(k/(k+1))^m + m/(k+1)*(k/(k+1))^(m-1). Multiply this error by (k+1)^m and sum
over k! permutations. Inactive empty-slot core-edge conditions occur only
inside this error event. QED.

Thus the skewed probe with both alpha_j=.2 has leading coefficient .0512,
while perfectly correlated binary03 colors with alpha_j=.5 retain coefficient
.125. Their fluctuation distributions are very different, despite the same
leading mean. The sharp uniform-color remainder above comes from the inherited
exact transform; this general elementary bound need not be sharp.

For k=2, E H=(3^m+3*2^m+3)/8 for m>=1. Hence the mean-scale witness-order proxy
is n approximately 2+log_3(8H). The universal language-specific lower bound is
n>=2+ceil(log_3(H/2)). Neither proves an upper bound for a prescribed target.
`TABLES.md` extends the scale to a million decimal digits without constructing
any such graph. A scalar H evaluation at that size still needs linear word
length and large integers; graph output is quadratic in n.

The probes use five distributions, lengths16..512, 1,200 samples per law.
Uniform sd(log H) stays roughly1.28..1.38 while m increases32-fold; the skewed
law stabilizes near3.47. Means H/3^m for the uniform law remain close to1/8.
These are predictions consistent with the theorem, not an independent audit
of Session270 or estimates of all-target success probability.

## 5. Coverage, entropy, symmetry, and SCC distinctions

At k=2 there are 4^m words but at most O(3^m) possible H values. The gap
log4-log3 is average multiplicity entropy, not a coverage certificate. A
local-limit hypothesis with order3^m support could suggest about(4/3)^m
representations per bulk value and typical target-order log_3 H+O(1). This is
**HEURISTIC ONLY**: arithmetic holes and strongly nonuniform multiplicities
can persist, and 7,21 remain the inherited excluded values. Exact small residues
and broad random support do not repair this logical gap.

For fixed core C, let G=Aut(C). Under uniform independent cross orientations,
a nonidentity g with c(g) core cycles preserves the word with probability
2^(-m(k-c(g))). Tournament automorphisms have odd order, so k-c(g)>=2 for
nonidentity g. A union bound gives

    P(nontrivial marked automorphism) <= (|G|-1) 4^-m.

This concerns automorphisms preserving the core. Every automorphism stabilizing
the core fixes its transitive complement pointwise and embeds into G. For an
unmarked graph T in this carrier, orbit-stabilizer therefore gives
|Aut(T)|<=binom(k+m,k)|G|. Unmarked automorphisms moving the core are not ruled
out by the exponential estimate. Similarly, forgetting a designated core has
at most binom(k+m,k) possible core subsets (and at most k! identifications
with the fixed labelled core). Thus word-class entropy survives unmarking up
to polynomial factors for fixed k, even though exact ordinary counts change.
This still says nothing about scalar H coverage.

SCCs and substitution require different sidecars. An ordinal sum has H equal
to the product of its component H counts, so log H is additive across SCCs.
The present random law is over cross-edge words, not over SCC factorizations.
Keep SCC order and rooted endpoint/2-cover data when allowing seam reversal;
the Session268 seam law contains subtraction and is not recoverable from
log H or the support-rank phase. The q-vector is sufficient for uniform
substitution; nonuniform seams demand rooted data. None of the scalar growth
phases classifies these interactions up to graph isomorphism.

## 6. Target-directed mechanisms: what to try instead of exhaustive coverage

1. **Scale first.** Start near m=ceil(log_3(8N)) for full support; compare nearby
   lengths because the bounded multiplicative endpoint fluctuations are broad.
   A graph with H close to N is not a witness for N.
2. **Cone-constrained suffix search.** Keep v(prefix) in the positive36-space
   and a^T M_suffix as a dual covector. H(prefix suffix) is their exact pairing.
   Precompute bounded suffix trees and use min/max pairings to prune only when
   their certified interval excludes N. Merging prefixes solely by H is invalid.
3. **Meet in the middle with congruences.** For modulus q, store the full36-state
   vector modulo q and suffix covector modulo q. Their pairing must equal N
   modulo q. This is an exact finite residue sieve. Moduli can be combined by
   CRT only for the same compatible prefix/suffix states, not independently
   projected scalar residue sets. The rational17 matrices may have denominators
   noninvertible modulo q; using them blindly loses information. The integer36
   realization avoids that problem.
4. **Arithmetic obstruction certificates.** A finite reachable residue-state
   closure disjoint from N modulo q proves impossibility in that grammar. A
   bounded search missing a residue does not. Conversely hitting all odd
   residues for several q does not imply hitting every odd integer.
5. **Digit algorithm on a returning cone.** Search finite macros with affine
   action on an integral slice (h,sidecars), and inverse branches reducing a
   certified height while covering residue classes. The required certificate
   has four parts: invariant physical cone, integral inverse choices, all
   residue classes covered, and strictly decreasing height. Spectral radius3
   alone supplies none of the inverse integrality or coverage clauses.
6. **Interaction enlargement.** When residue states stagnate, add the inherited
   rooted seam and track its end/start/two-cover coordinates. Recompute the
   physical image cone rather than assuming arbitrary integer sidecars exist.
   Product-only grammars are already sparse; a useful enlargement must create
   new residues or change the returning cone, not just add more products.

`target_probe.py` compares equal budgets of iid guesses, local integer-distance
mutations, and coarse distance plus residue tie-breaking, at targets10^12+1,
10^30+1,10^100+1. Every scored H is computed by exact transfer. This is a cheap
hostile test of search heuristics, not an exhaustive solver or a calibrated
probability study. It saves the best word and exact H whether or not a hit
occurs. Close misses are reported as misses. The important distinction is
finding the right exponential scale versus solving exact arithmetic.

## 7. Scope, reproduction, and remaining questions

Run from the repository root:

    OPENBLAS_NUM_THREADS=1 python3 05-knowledge/results/session271_predictive_tournaments_20261111/chanpak/predict.py
    python3 05-knowledge/results/session271_predictive_tournaments_20261111/chanpak/target_probe.py

`predict.py` imports only the inherited transfer constructor from Session270.
It performs 10,000 modest random coverage samples, not the previous million
-target exhaustive census. Floating normalized transfer is used for growth
statistics; target scores and stored H values are integers. Random sampling
has no confidence guarantee for unobserved tails. The high-index table is a
scale extrapolation supported by an exact growth theorem, not a claim of
witness construction. No repository-global files or proof IDs are changed.

Open next objects: limiting boundary distribution of H/r^m; exact negative
pressure and optimal rare words; finite integral cones for inverse macros;
residue-state compression respecting the physical lattice; uniformly growing
k regimes; and rooted-seam interaction phases. The complete classification
here is of exponential growth by finite support for fixed cores and iid or
periodic words. Infinite Chan–Pak and complete recursive tournament structure
classification remain open.

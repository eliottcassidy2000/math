# Predict enormous tournaments; classify what recursion can actually see

**Session271, 11 November 2026. Research packet complete. Chan–Pak remains OPEN.**

This session changes the objective from recovering exact A000568 integers to
predicting much larger scales, then explains why numerical census accuracy,
recursive graph classification and exact Hamiltonian-target coverage are
three different problems. No exact A000568 cache was queried or extended.
Inherited results were used in their stated scope rather than re-audited.

The strongest new junction is an explicit theory of **which observer closes
under which operation**. The full graph has a canonical modular tree with an
infinite prime alphabet; H needs a path-cover vector; adjacency spectrum needs
only its characteristic polynomial under uniform substitution; seams expose
rooted states; fixed-core random spines have an integer growth classification.
These are actual formulas and scoped completeness theorems, not just labels.

## 1. Predict A000568 at indices up to 10^1000

Let N=binom(n,2), b_n=2^N/n!, and a_n=b_n(1+delta_n). The identity predictor
b_n already becomes extraordinarily accurate in relative terms as n grows.
The complete toolkit combines logGamma, bounded Stirling expansions, finite
odd-cycle support, inherited defect polynomials, adaptive truncation and local
ratio forecasts. All calculations retain logarithms and scientific notation.

| n | A000568 prediction, rounded for display | log10 of identity-model relative-error bound |
|---:|---:|---:|
| 300 | 5.122742023739228 × 10^12886 | −172.46398… |
| 1,000 | 7.554225356196997 × 10^147796 | −592.33429… |
| 10,000 | 1.509893164724477 × 10^15014335 | −6007.87304… |
| 1,000,000 | 1.190087511215613 × 10^150509281608 | −602041.26433… |
| 10^12 | 2.039399725886988 × 10^150514997820274377090934 | −602059991291.23539… |

The error column applies to the **unrounded mathematical model**, not the
16-digit mantissa displayed here. Numerical evaluation is high-precision
mpmath with a precision-doubling diagnostic, not directed interval arithmetic.
Rounding the displayed mantissa introduces its own much larger relative error.
The committed JSONL also covers 10^30,10^100 and10^1000.

The reusable envelope is simple. Every odd-cycle permutation moving m vertices
loses at least m(2n−m−1)/3 pair-orbit bits. Convexity then bounds every omitted
support m>=M by

    U_M(n)=(n−M+1) max{n^M 2^[−M(2n−M−1)/3], n^n 2^[−n(n−1)/3]}.

Summing a finite head and adding this tail gives an all-index symbolic
relative-error bound. Its decimal rendering is approximate. At fixed desired
accuracy, huge n actually makes the model easier, although input/output bit
costs still grow and exact integer printing has quadratic-size output.

The leading correction is (n)_3 2^(−2n+4)/3. Disjoint 3-cycles are not an
independent Poisson gas: the k-cycle term has a pair interaction factor
64^binom(k,2), apart from finite-population depletion. A single5-cycle does
not even dominate its defect layer: the two3-cycle contribution is larger by
10(n−5)/9. The defect transseries itself is inherited from Session98, credited
through THM-4845; the practical envelope, adaptive implementation and gas
comparison are the contributions here.

**Absolute error remains gigantic.** At n=1000 the identity model's absolute
error bound has log10 about147204.544, despite relative error below roughly
10^−592.334. This does not license exact rounding, exact low bits, or a claim
that the integer was computed. The Session199 warning is respected.

[Forecast methods, proofs, full tables and executable](forecast/README.md).

## 2. A complete recursion backbone, and a rigorous limit to finite catalogues

Accepting the inherited modular decomposition theorem, every finite tournament
has a unique canonical tree of ordered SCC nodes and prime-quotient nodes,
with singleton leaves. Isomorphism is exactly equality of the decorated tree
modulo the actual quotient automorphisms. The packet specifies the complete
representation, automorphism extension, converse action and species interface.
The prime decorations form an **infinite alphabet**.

For s fixed prime types, unrestricted recursion and arbitrary ordered sums
produce at most [16(s+1)]^n ordinary n-vertex classes. By contrast the global
count is at least 2^N/n!. Thus any fixed finite prime catalogue is exponentially
thin on the n² entropy scale, even with unbounded recursive depth.

The stronger interaction barrier uses coding theory. For every alpha<1/2,

    Pr[d_arc(T,F)<=alpha N]
       <= n! [16(s+1)]^n 2^[−N(1−h_2(alpha))],

where T is uniform labelled, F is that finite-alphabet family, and the distance
minimizes over all labellings. Consequently d_arc(T,F)/N tends to1/2 in
probability. A transitive ordering gives the matching deterministic upper
bound. Even o(n²) arbitrary seam edits cannot make the finite catalogue
represent a typical tournament.

This is not a reason to abandon recursive families. They are excellent
constructive languages, but the global census is dominated by prime atoms.
A direct module calculation gives

    Pr(nonprime) <= [4 binom(n,2)+4n+64]2^(-n),  n>=5.

Twin-module rarity is of order n²2^−n in expectation; automorphism corrections
are of order n³2^−2n. The two vertices of a twin module cannot be swapped because
their internal edge is directed. Module recursion and symmetry corrections
therefore describe different rare populations.

[Canonical atlas and interaction boundaries](structures/ATLAS.md).

## 3. The spectral observer closes more tightly than expected

For any tournament child T_i, regular or irregular, write f_i(z)=chi(T_i;z)
and g_i(z)=(-1)^n_i f_i(-z-1)-f_i(z). For quotient adjacency B,

    chi(Q[T_i];z)=det(diag(f_i)-B diag(g_i)).

The reason is the tournament identity A+A^t=J−I: it determines the all-ones
resolvent response from chi itself. For C3 this becomes f1 f2 f3−g1 g2 g3.
Thus adjacency-cospectral children are interchangeable in every uniform
substitution context for this observer. No regularity assumption is needed.
This is not graph isomorphism and does not recover finite-field Jordan data.

A seam changes the necessary state. If E selects named ports and an edit is
A'=A+EKE^t, the port response C=E^t(zI−A)^−1E obeys

    chi(A')=chi(A)det(I−KC),    C'=(I−CK)^−1C.

Uniform block connections use the all-ones column; seam edits use port
columns. The resulting determinant/Woodbury network is exact on the actual
tournament adjacency, not a loose physical analogy. The same unrooted TT2,
with source versus sink chosen as its port, produces z³−1 versus z³ after
one seam. Rooted information is essential.

The Hamiltonian counterpart is a path-cover tensor retaining ordered
start/end port labels. It is sufficient for arbitrary named cross-edge
interactions, but its run index grows with graph size. Fixed port number alone
does not prove a universal finite-state H machine.

[Proofs and executable reference compiler](structures/README.md).

## 4. Chan–Pak growth has finitely classified phases in a fixed carrier

Fix a k-vertex core and an ordered transitive spine. A color letter records
its k cross orientations. For a nonempty letter support S, define r(C,S) as
the maximum number of feasible occupied spine slots over all core orders.
The feasibility condition is explicit and finite.

The new theorem gives, for iid letters positive exactly on S,

    log H/m -> log r almost surely,
    E log H=m log r+O(1),       Var(log H)=O(1),
    lim m^(-1)log E[H^s]=s log r   for every s>0.

Periodic words with support S have H=Theta(r^m). Full support gives r=k+1,
regardless of the core's orientations. All15 nonempty supports for the two-core
language fall into exactly three phases, with bases1,2,3.

The mechanism is important: occupied slots only accumulate. Once the active
slots are filled, every new letter has r legal bulk assignments. Fixed
initialization and acceptance words control the endpoints. Bulk growth is
therefore deterministic at the exponential scale; bounded endpoint noise
replaces a generic sqrt(m)-fluctuation story. In probes from m16 to512, the
uniform law's sd(log H) stays roughly1.28–1.38.

For an iid color law with bit marginals alpha_j strictly between0 and1,

    E H = k! product_j alpha_j(1-alpha_j) (k+1)^m + O(m k^m).

Within-letter bit correlations and core orientation do not affect this leading
coefficient. They do affect the boundary distribution. For uniform colors,
the inherited exact transform sharpens this to k!/4^k times(k+1)^m+O(k^m).
Meanwhile any distinguishing lower q-coordinate is exponentially suppressed
in the averaged bulk. Long spines predict growth well but hide the core.

For a target of scale10^100, the two-core language has a rigorous order lower
bound211 and an ensemble-mean order proxy214. At scale10^1000000 these become
2095905 and2095908. **The proxy is not a witness-size upper bound or an
existence guarantee.** Three search heuristics were compared at equal finite
budgets. A local search reached relative error2.05×10^−6 at10^30+1, but none
of the nine runs hit its target exactly. Residue tie-breaking helped one case
and badly hurt another. Near-target prediction is not Chan–Pak.

[Phase classification, tables, target experiments and next algorithms](chanpak/README.md).

## 5. Census prediction now controls other metric distributions

Let U be uniform ordinary classes and L the distribution induced by uniform
labelled tournaments. There is an exact positive mixture

    U=L/(1+delta)+delta D/(1+delta),

where D weights a class by1−1/|Aut|. Hence TV(U,L)<=delta/(1+delta).
Every event, including primeness and the edit-distance barrier, transports
between these ensembles with this absolute error. Rare conditional events
can amplify it; the weights have not been silently discarded.

The full labelled path-cover mean is

    E_L Q_T(t)=n! 2^(1-n)(1+2t)^(n-1).

For its r-th coordinate, the ordinary-class mean has downward relative error
at most delta/(1+delta) and upward error at most

    delta(2^(n-r)-1)/(1+delta).

For example, ordinary-class mean H at n=1000 is predicted by
7.51066780758×10^2266, with upper relative model error below roughly10^−291.605.
The n=10^6 prediction is1.66937049362×10^5264679, with bound roughly10^−301011.57.
These refer to means, not concentration or individual tournament counts.
Positive multilinear substitution propagates the coordinate intervals through
independently sampled marked modules. Uniform unmarked parent sampling still
requires its quotient-cycle-index weights.

A niche connection to the repository's Gaussian gcd field explains which
stochastic predictor to avoid. Naive sampling of the identity sector has
relative variance2^(n²)−1. An exact Cameron–Martin identity tilt makes that
sector constant. A retained3-cycle correction then has relative variance63,
and two such log contributions have correlation factor64, exactly matching
the cycle-gas interaction above. The full growing-support variance remains
open; deterministic summation is still best for this unweighted sequence.

[Ensemble transport](ensemble/PROOF.md) · [Gaussian prediction bridge](ensemble/GAUSSIAN_PREDICTOR.md).

## 6. What is completed, and the precise open frontier

Completed: a giant-index predictive toolkit; complete finite-tournament normal
form with an infinite prime alphabet; exact graph/symmetry/spectral/port
interfaces; a finite support-phase classification for fixed-core iid/periodic
languages; and quantified links between census, distributions and recursion.

Not completed: a finite classification of every prime or every dynamical
behavior, an all-context finite H state, a limiting boundary distribution,
or the infinite Chan–Pak realization theorem. The most concrete next task is
an integral returning cone with inverse macros that cover residues and strictly
decrease a height. The17-state quotient helps equality; the positive36-state
realization helps inequalities and lattice-safe residues. Those jobs should
not be forced onto the same reduced state without its missing sidecars.

[Ranked open problems](OPEN_QUESTIONS.md) · [Reproduction and validation scope](REPRODUCE.md).

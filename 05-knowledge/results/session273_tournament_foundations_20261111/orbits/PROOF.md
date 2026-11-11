# Orbit conditioning, coherent growth, and an infinite symmetry reversal

**PROVED IN NOTE**, except explicitly inherited results and open targets.
Session273, 2026-11-11. No external novelty claim. This extends Session271's
positive orbit-mixture to an exact growth process and infinite-volume bridge.
It does not supply a new proof of Chan–Pak or Yang–Mills uniqueness.

## 1. Three measures and a boundary that cannot be omitted

Put `m_n=binom(n,2)`, `a_n=A000568(n)`, `b_n=2^m_n/n!`, and
`delta_n=a_n/b_n-1`. A tournament's automorphism order is `g(T)`.
Let `P_n` be fair independent arcs on the named vertices `[n]`. Let `U_n` be
uniform on ordinary isomorphism classes. Lift `U_n` to a named law `Q_n` by
choosing a uniform labelling of the selected class. Orbit–stabilizer gives

    Q_n(T) = g(T)/(a_n n!),
    dQ_n/dP_n(T) = g(T)/(1+delta_n),
    E_P g = 1+delta_n.                                      (1)

Thus `Q_n`, not `P_n`, is the literal labelled realization of the uniform
ordinary-class census. S271's positive decomposition also works before
forgetting labels: since `g>=1`,

    Q_n = [P_n + delta_n D_n]/(1+delta_n),
    TV(Q_n,P_n) <= delta_n/(1+delta_n),                       (2)

where `D_n` is the normalized excess density `g-1` if delta is positive.

**The laws Q_n are not projectively consistent.** There are four ordinary
classes at order four: transitive, a cyclic triangle with a source, a cyclic
triangle with a sink, and the strong four-vertex class. Their numbers of
cyclic induced triples are respectively `0,1,1,2`. A uniform three-subset
under Q_4 is therefore cyclic with probability `(0+1+1+2)/(4*4)=1/4`.
Under Q_3 the cyclic and transitive classes each have probability `1/2`.
The discrepancy is not a choice of gauge: uniformly label first and delete
vertex four, or delete a uniform vertex of an unlabelled class, gives the
same calculation. Consequently there is **no infinite law whose [n] marginal
is Q_n for every n**. This is an obstruction to a proposed exact model, not an
obstruction to convergence of Q_n as n tends to infinity.

## 2. Exact finite-horizon growth by a conditional automorphism martingale

Fix a final size N. Under the fair process reveal vertices one at a time,
with `F_n` containing all arcs within `[n]`. Define the positive rational
function on named n-tournaments

    h_n^(N)(t) = E_P[g(T_N) | T_n=t],  0<=n<=N.             (3)

Then `h_0^(N)=1+delta_N`, `h_N^(N)=g`, and `h_n^(N)>=1`.
For the `2^n` possible incident rows x of vertex n+1,

    h_n^(N)(t) = 2^-n sum_x h_(n+1)^(N)(t+x).              (4)

The row transition

    K_n^(N)(t,t+x) = 2^-n h_(n+1)^(N)(t+x)/h_n^(N)(t)      (5)

is therefore a probability kernel. Multiplying its probabilities telescopes:

    Pr(T_N=t) = 2^-m_N g(t)/(1+delta_N) = Q_N(t).           (6)

More generally its n-prefix density relative to P_n is
`h_n^(N)/(1+delta_N)`. This is an exact positive Doob transform, not a signed
sampling prescription. The missing state is the **terminal horizon N and
its conditional automorphism response**, not just the present class or its
current automorphism order. Formula (5) is invariant under relabelling of the
current vertices and so can be quotiented to a class-valued branching law
with the actual multiplicities of incident rows retained.

### A finite, exact completion formula

For a permutation sigma of `[N]`, act on unordered pairs, retaining whether
the chosen orientation of each pair is reversed. If an orbit returns with
opposite orientation, it supports no fixed tournament. Otherwise let
`r_sigma` be its number of unordered-pair orbits, and let `q_sigma(t)` count
the orbits touched by prescribed arcs of t. If prescribed arcs on an orbit
force conflicting orientations, the contribution is zero. In every other
case it is

    Pr(sigma in Aut(T_N) | t)
       = 2^[r_sigma-q_sigma(t)-(m_N-m_n)].                  (7)

Summing (7) over permutations gives (3). The proof is simply that each
untouched consistent orbit has one independent orientation bit. This is a
conditioned Burnside sum; it does not assert a fast algorithm. It reveals
exactly which cross-boundary symmetry information a growing state needs.

### Why the boundary disappears locally

By conditional expectation and `g-1>=0`,

    E_P[h_n^(N)-1] = delta_N.

Consequently, simultaneously for every `0<=n<=N`,

    TV(Q_N restricted to [n], P_n) <= delta_N/(1+delta_N).   (8)

This follows also by data processing from (2). It is uniform even when n
increases with N. It is an averaged statement, not a uniform bound on
`h_n^(N)(t)` for rare prescribed t. Conditioning on a rare cylinder of
probability p can amplify errors by a factor on the order of `1/p`.
The full N-vertex histories can be coupled to the fair history with failure
probability at most the bound in (8); separate per-step error sums are not
needed.

## 3. A complete local limit and a reversal of symmetry

The proven Burnside envelope of
[S271 forecast](../../session271_predictive_tournaments_20261111/forecast/PROOF.md)
implies `delta_N -> 0`. Thus all fixed-prefix marginals of Q_N converge,
without choosing a subsequence, to P_n. The P_n are mutually consistent and
define the unique probability law P_infinity on the countable named arc
space `{0,1}^{binom(Natural,2)}` with independent fair coordinates. One may
extend each Q_N arbitrarily outside `[N]`; every such extension has the same
weak limit in this compact product space, because cylinders eventually lie
inside `[N]` and continuous functions are uniformly approximable by cylinder
functions. Arbitrary extensions need not be exchangeable at finite N.

There are two different uniqueness statements:

1. On the **named coordinate space**, the cylinder probabilities determine a
   unique limiting probability measure.
2. After forgetting names, almost every infinite sample belongs to the same
   isomorphism class, the countable homogeneous random tournament.

For the second statement, for any finite disjoint A,B each fresh vertex has
probability `2^(-|A|-|B|)` of satisfying `A -> v -> B`; disjoint fresh choices
are independent. Infinitely many satisfy it almost surely. Intersect over
countably many pairs A,B. The usual explicit back-and-forth construction
then extends every finite partial isomorphism, proving uniqueness of the
countable extension-property tournament and its homogeneity. This elementary
argument is included to state the precise meaning of the limit, not as a new
Fraïssé theorem.

**Symmetry reversal.** Finite random tournaments are rigid with probability
approaching one, yet the infinite limit has nontrivial automorphisms: any
one-vertex map extends by back-and-forth. For finite labelled tournaments,
`Pr_P(g>1)<=delta_N/2`, because tournament automorphism groups have odd order,
so nontrivial ones have at least three elements. For Q_N use this and (2).
An infinite automorphism need not preserve any named finite prefix. Hence no
continuity theorem transfers finite rigidity to infinite rigidity. This is
a concrete answer to how symmetry can emerge at an infinite scale even when
A000568 is asymptotically governed by rigid finite objects.

## 4. Chan–Pak: age coverage versus bulk target sampling

Let h be a positive integer. If a k-vertex tournament S has exactly h
Hamiltonian paths, then the infinite random tournament contains S induced,
indeed infinitely many copies. Conversely an induced finite witness in that
random tournament is an ordinary tournament witness. Therefore, almost surely,

    {H(S): S finite induced subtournament of T_infinity}
      = {H(S): S any finite tournament}.                   (9)

This is a useful exact carrier for the conjecture, not a coverage proof.
The extension property inserts prescribed **arc patterns**; it does not
manufacture a finite pattern with a prescribed **Hamiltonian count**.

There is a quantitative finite version. Let `p_S` be the probability that a
fair k-vertex tournament is isomorphic to a fixed S:

    p_S = k!/(g(S) 2^m_k).

Use floor(N/k) disjoint labelled blocks. Under P_N these trials are
independent, so

    Q_N(no induced copy of S)
      <= (1-p_S)^floor(N/k) + delta_N/(1+delta_N).           (10)

The left event forbids copies anywhere, so it is contained in the event
that none of those blocks works. Replacing p_S by the total probability of
all k-vertex h-witnesses sharpens the same bound. A verified finite witness
atlas can therefore be amplified to a high-probability induced atlas inside
large uniform classes. No inference runs backward from an assumed positive
p_S when witness existence is the open problem.

For comparison the inherited exact variance result
[THM-4761, Walsh chaos](../../../../01-canon/theorems/THM-4761-random-tournaments-the-forest-identity-is-the-walsh-chaos-decomposition-of-kappa.md)
has, for `mu_N=N!/2^(N-1)` and `kappa=H/mu_N`,
`Var_P(kappa)=2/N+O(N^-3)`. Thus for fixed h and large N,

    Q_N(H=h) <= Var_P(kappa)/(1-h/mu_N)^2
                  + delta_N/(1+delta_N) -> 0.             (11)

The fixed target can occur densely among small induced substructures while
becoming vanishingly rare as the **total** of the large object. This is why
improved A000568 accuracy or typical-H concentration is not automatically a
constructive proof of Chan–Pak. The existing conditional CLT for kappa under
P_N also transfers unchanged to Q_N by (2), conditional on precisely the
same cited input as THM-4761; no new CLT assumption is supplied here.

## 5. Bridge contracts and sources

- **A000568 -> probability / inverse growth:** orbit–stabilizer maps uniform
  class counts to automorphism-tilted named laws. It preserves class uniformity
  at the terminal size. Forgetting the final horizon destroys exact growth
  consistency. Sidecar h_n^(N) restores it via (5).
- **Probability -> A000568:** the normalizer of the terminal positive weight
  is `1+delta_N = a_N/b_N`. Conditional symmetry-orbit response (7) turns
  census prediction into quantitative convergence of all bounded observers.
- **Yang–Mills finite-volume inspiration:** S272's
  [limits/measure note](../../session272_yang_mills_bridges_20261111/limits_measure/README.md)
  requires a common carrier, consistency and uniqueness beyond subsequences.
  Here all three are explicitly achieved, despite inconsistent Q_n. The
  transfer is a proof mechanism: vanishing total-variation boundary tilt.
  No map from YM fields to tournament arcs, or comparable YM bound, is claimed.
- **Older controlled-forgetting thread:**
  [S261 perspective groupoid](../../../../07-reflections/perspective-groupoid-controlled-forgetting-codex-s261.md)
  asks for the missing sidecar before adding an observer. Here the sidecar is
  an exact conditional symmetry response, and the hostile is n4 -> n3.
- **Older random-tournament curiosity:**
  [S699 monodromy/Rado reflection](../../../../07-reflections/monodromy-the-rado-tournament-and-the-worry-set-as-branch-locus-s699.md)
  supplied the extension-property inspiration. Its broad LRC/monodromy
  identifications are not used as proved implications. We supply the actual
  map from finite class-uniform laws and retain the arc-pattern/count gap.
- **Collatz merge-entropy lens:**
  [THM-4481](../../../../01-canon/theorems/THM-4481-entropy-merge-law-sign-strategies.md)
  retains the law at merged states. The common mechanism is conditioning on
  the actual fibers, not replacing them by a scalar unweighted support.
  There is no trajectory-preserving map to Collatz and no Collatz conclusion.

The substantive contribution of this note is the finite-horizon conditioned
Burnside compiler, the uniform growing-window convergence statement, and its
age/bulk and finite/infinite symmetry consequences. Orbit–stabilizer, Doob
transforms, product measures and back-and-forth are standard tools.

## 6. Open extensions

1. Compress h_n^(N) on structured prefix families without enumerating all
   permutations, retaining cross-boundary pair-orbit constraints.
2. Seek a useful horizon-free class-growth process whose N-marginals are
   quantitatively close to U_N; fair vertex growth already provides one after
   forgetting names, but a faster exact class sampler would need new control.
3. For Chan–Pak, bound the minimal witness order k(h) uniformly by a usable
   constructive rule. Equation (10) then amplifies that rule; it does not
   replace it. The existential claim over unrestricted sizes is not a single
   fixed-arity cylinder test.
4. In other statistical models test whether their symmetry/gauge tilt has
   vanishing total variation after the correct shared carrier is specified.
   Nonvanishing tilt or singular observables block this shortcut.

## Controls

Run `python3 control.py` in this folder. It enumerates only n<=5 as a hostile
control, not as fresh census progress: all named tournaments, automorphism
orders and orbits; every prefix h; the independently implemented signed
pair-orbit formula (7); every transition normalizer; the TV bound and the
n4 -> n3 failure. Exact rational output is in `control.json`.

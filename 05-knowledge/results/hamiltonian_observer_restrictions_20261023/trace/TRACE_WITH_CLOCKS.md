# Guarded return sampling: the gap survives, the clock does not disappear

Status: **PROVED**, finite reversible chains and the specified auxiliary tournament coupling. **FINITE-EXACT** checks are in `check_trace.py`. **OPEN**: realizing these return transitions by native port-preserving metagraph moves with a controlled cost. No all-size native wiring or major open conjecture is claimed.

## 1. Inheritance and the right carrier

The closest proved mechanism is the full-cube sandwich `Q=P_G P_T P_G`, whose nonconstant eigenvalues lie in `[0,1/9]`; see [the incoming coupling proof](../../session42_shared_lattices_square_colors_20261023/coupling/COUPLING_OPERATOR.md). On the auxiliary tournament of order `n-1`, it gives the same gap `8/9` on path-marked presentations of order `n>=4`; see [the restriction packet](../coupling/RESTRICTIONS.md). The unsymmetrized product must not be substituted for the selfadjoint sandwich.

The canonical hostile is the native transitive class: it has one Hamiltonian observer, giving the empty backward graph. Hard restriction to Hamiltonian words isolates this component in the native graph/tournament incidence coupling. The correction is to specify which operator is being restricted, and how time is counted.

The least-used relevant sidecar is the rational harmonic extension from the [native electrical packet, section 5](../../metagraph_network_invariants_20261021/PROOF.md). At order five, terminal resistance `6092/5137` has odd denominator, whereas the full voltage contains `11327/10274`. Boundary energy alone loses the integer-current obstruction. That example motivates retaining extension data here; it does not identify the native Laplacian with the auxiliary coupling.

This is the standard watched/trace chain construction, as distinguished from rejection and collapse by [Aldous–Fill, section 2.7](https://www.stat.berkeley.edu/~aldous/RWG/Book_Ralph/Ch2.S7.html). The application, guard ledger and clock formulas below are stated and proved explicitly. No priority claim is made for the general Markov-chain mechanism.

## 2. Restriction theorem with exact quantifiers

Let `Q` be a finite irreducible reversible stochastic matrix with stationary probability `pi`. Assume its operator on `L2(pi)` is positive semidefinite and its spectral gap is at least `gamma>0`. Let `S` be any nonempty subset, and put `C=X\S`. For proper `S` define

    B = (I-Q_CC)^(-1),
    R = Q_SS + Q_SC B Q_CS.                         (1)

For `S=X` set `R=Q`. The inverse exists: from every state in `C` the irreducible chain eventually hits `S` with probability one, so the spectral radius of the killed matrix `Q_CC` is below one.

**Theorem.** `R` is the chain observed at successive positive return times to `S`. It is reversible with stationary law `pi_S=pi(.|S)`, is positive semidefinite, and has spectral gap at least `gamma`. On a singleton the mean-zero space is zero, so the contraction assertion is vacuous.

**Proof.** Expand `B=sum_(k>=0) Q_CC^k`; each term of (1) is a path that first returns after `k+2` steps, together with the direct one-step term. Path reversal gives detailed balance with `pi_S`. For `f:S->R`, extend harmonically by

    h_S=f,  h_C=B Q_CS f.

For `E_Q(h)=<h,(I-Q)h>_pi`, block elimination gives

    E_R(f)=E_Q(h)/pi(S).

The ambient Poincare inequality and the variance decomposition give

    E_Q(h) >= gamma Var_pi(h)
             >= gamma pi(S) Var_pi_S(f).

Hence `E_R(f)>=gamma Var_pi_S(f)`. To obtain positive semidefiniteness as well, conjugate by `diag(sqrt(pi))`. The symmetric principal block `Q_SS` is PSD; `I-Q_CC` is positive definite; the second summand in (1) is then of the form `A B A^T` with positive definite `B`. Both summands are PSD. Thus every nonconstant eigenvalue of `R` lies in `[0,1-gamma]`. QED.

In particular, *every nonempty accepted subset of the auxiliary Hamiltonian-observer lift has return-time gap at least `8/9`*. This includes subsets specified by multiple actual owner pins and any exactly evaluated guard. It does not require the subset to be a subgroup orbit, a coset, or an unfiltered fibre.

## 3. The information that a return event must carry

The generating matrix for the actual excursion length is

    R(z)=z Q_SS + z^2 Q_SC (I-z Q_CC)^(-1) Q_CS.    (2)

The coefficient of `z^ell` is the joint probability of the target return state and the elapsed `ell` ambient steps. In particular `R(1)=R`; discarding `z` discards the clock. Its row sums at the derivative give

    E_s[tau_S^+] = 1 + (Q_SC B 1_C)_s.

Under a stationary conditional start,

    E_pi_S[tau_S^+] = 1/pi(S).                     (3)

To verify (3) directly, stationarity on `C` gives `pi_S^raw Q_SC=pi_C(I-Q_CC)`, hence the weighted sum of the extra term is `pi(C)/pi(S)`.

Thus a faster chain in return-event time can have an enormous physical waiting time. A full record may retain the literal excursion; for endpoint observables and duration statistics, (2) is sufficient. For integer-lattice questions the harmonic extension must also be retained and the original guard re-evaluated on its original carrier. Neither (1) nor a floating-point effective resistance certifies native mod-8 data.

This is a concrete instance of a recurring principle: a quotient can preserve the desired observable exactly while losing the lifting data. Here the three retained objects have different jobs: `R` for return destinations, `R(z)` for clocks, and `B Q_CS` for harmonic extension. They are not interchangeable.

## 4. Actual owner orders, kernel lifts and converse

Fix `n>=4`, `m=binom(n-1,2)`. A free word `g` specifies the nonconsecutive arcs of a tournament with path `1->...->n`. The position pair `(i,j)`, `j>=i+2`, is the edge `(i,j-1)` in the auxiliary `(n-1)`-vertex graph/tournament cube. Let `P` be a permutation of the actual owner set. Interpret `(g,P)` as the actual tournament `T(g,P)` with Hamiltonian list `P`.

Use the product state space

    X = {0,1}^m x S_n,
    Q_lift = Q_aux x J_(S_n),

where `J` independently refreshes the owner order uniformly. The nonzero spectrum is unchanged, so the gap remains `8/9`. The full owner order is part of the state; no kernel or automorphism quotient is silently taken.

For a fixed actual tournament `T0`, let `S` impose `T(g,P)=T0`, any supplied indexed vertex/arc pins, and the actual native guard predicate. Membership is evaluated after realization. Assume `S` is nonempty. The theorem applies to this exact subset, although arbitrary guards can destroy the unfiltered compiler's power-of-three fibre counts.

Without extra pins or guards, `|S|=H(T0)`: every actual directed Hamiltonian list determines exactly one free word. Therefore

    pi(S)=H(T0)/(n! 2^m),
    mean stationary return time=n! 2^m/H(T0).       (4)

This is a mathematical sampler, not an efficient enumeration claim. Enumerating or recognizing accepted states may itself be expensive.

Converse acts concretely by `(T0,P)->(T0^op,reverse(P))`; an indexed position `i` becomes `n-1-i`. The free word is reflected across the reversed positions. Pins must be transported, not merely assigned the same scalar signature. One may include an independently refreshed branch tag and restrict by an explicit converse-equivariant guard, but no converse identification is required for the theorem.

Kernel lifts remain actual permutations in this construction. If a compressed compiler is used to realize `P`, its source kernel, target lift and shared-prefix coordinates must all be retained. The return theorem does not solve compatibility by counting quotient states.

## 5. Hostile controls: three different restrictions

1. **Recompute conditional incidence on a restricted word set.** This changes both conditional-expectation projections and their overlap. The native transitive component disconnects. No global `8/9` result survives by inheritance.
2. **Reject every attempted exit.** Even the ambient complete refresh `Q=J_N`, with gap one, becomes `(1-a)I+aJ_S` on `S`, where `a=|S|/N`. Its gap is `a`, arbitrarily small.
3. **Watch successive returns, retaining duration.** For the same complete refresh, the return chain is exactly `J_S`, with gap one and mean duration `1/a`.

The return proof is useful because it survives *arbitrary* nonempty guard restrictions. Its limitation is equally exact: intermediate excursion states need not preserve the source or ports. The accepted endpoints do. If every intermediate state must be a native legal move, a separate lift/comparison theorem is required.

## 6. The remaining native obligation

For each positive return transition `(s,t)`, seek a path of allowed native moves from `s` to `t` that retains the source, indexed ports, kernel lift, converse rule and original mod-8 receipt at every intermediate state. Assign an integer length and bound the load on each native move. A Dirichlet comparison can then transfer a quantitative bound, with the congestion and lengths charged explicitly.

A connected endpoint relation is necessary; trace irreducibility alone is insufficient, since ambient excursions may jump across disconnected native components. Thus the next decisive test is a small actual guarded fibre with two accepted states: are they connected by legal coupled moves? If not, identify the additional coordinate or larger admissible block needed, rather than asserting that the spectral gap solved native wiring.

## 7. Reproduction and scope

Run `python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/trace/check_trace.py` from the repository root. The script independently constructs graph and tournament permutation projections on all eight order-three words. It verifies exact rational PSD and gap inequalities for every nonempty restriction, the harmonic and clock identities, a nonuniform reversible control, and actual order-four Hamiltonian owner fibres and their nonempty subsets. No numerical eigenvalue threshold is used.

The proof is all-size; the exact controls are deliberately small enough to expose wrong clocks, weights, observer actions and conditioning conventions. This packet establishes a return-time bound and a concrete remaining move-lift obligation, not the final native all-size wiring rule.

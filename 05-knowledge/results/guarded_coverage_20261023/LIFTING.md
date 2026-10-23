# A quantitative target for complete guarded coverage

23 October 2026. **PROVED:** the elementary balanced-interface lifting lemma
below. **CONDITIONAL:** its use for any new native family until its coverage,
couplings, and cost hypotheses are authenticated. No novelty claim is made
for conditional variance or Markov-chain comparison.

## 1. The precise theorem worth aiming to instantiate

Let X be a finite actual state space with positive marked source distribution
pi. Partition **all** of X into charts X_i. Write p_i=pi(X_i) and pi_i for
the conditional distribution. Assume:

1. Each chart has a reversible local chain L_i of gap at least gamma>0.
   Singleton charts contribute zero variance and may be given a holding chain.
2. The chart indices have a reversible chain Q with stationary distribution
   p and gap at least delta>0.
3. For each Q_ij>0 there is a stochastic coupling K_ij from X_i to X_j with
   `pi_i K_ij=pi_j` and
   `pi_i(x) K_ij(x,y)=pi_j(y) K_ji(y,x)`.
   Thus the incoming distribution has the **actual conditional source mass**;
   a single convenient representative is insufficient.

Define L by the local chains and C(x,y)=Q_ij K_ij(x,y). Holding index
transitions may use identity couplings. Both chains are reversible for pi.
For 0<a<1, P=aL+(1-a)C satisfies

    gap(P) >= min(a gamma, (1-a) delta).

In particular a=delta/(gamma+delta) gives gap at least
`gamma delta/(gamma+delta)`.

**Proof.** Put m_i=E_(pi_i) f. Exact conditional variance gives

    Var_pi(f)=sum_i p_i Var_(pi_i)(f)+Var_p(m).

The local energy bounds the first term by gamma. For each coupling, Jensen's
inequality with its two authenticated marginals gives

    E_(pi_i K_ij) (f(x)-f(y))^2 >= (m_i-m_j)^2.

Hence the cross-chart energy bounds the second term by delta. Multiply by
the mixture probabilities and sum. The equality boundary includes a tensor
product where a slow eigenfunction depends only on the slower factor. QED.

Now let N be a reversible **native** chain on the same actual states and
source distribution. For every P transition supply a distribution of legal
native paths, retaining all ports and kernel lifts. Suppose the usual
length-weighted load on every directed native edge is at most Lambda times
its conductance pi(u)N(u,v). Cauchy–Schwarz along each path then gives

    E_P(f,f)<=Lambda E_N(f,f),
    gap(N)>=gamma delta/((gamma+delta)Lambda).

All path lengths are integers. The bound counts primitive events; any other
work cost, including owner-clock changes or original dynamical time, must
be recorded separately or explicitly included in the load inequality.

## 2. Why this changes the research target

The [complete pointed atlas](CUBES.md) solves the coverage part for ordinal
sources, and supplies gamma=1/r within every nontrivial chart. It does not
provide Q or its balanced native couplings. The [four-layer theorem](PORTS.md)
supplies certain legal resource-preserving paths after a declared source
extension. It does not bound their congestion for an arbitrary distribution
on assemblies. These are different missing columns of this theorem.

This makes a concrete, source-sensitive target: find cross-chart **flows**
whose marginals are the actual pi_i, then route those flows through available
native resources with a uniform Lambda. Counting one connector per pair
of chart labels cannot establish either condition. The marked source row
cannot be replaced by a rotation-averaged trace merely because the latter
has a convenient cancellation identity.

The earlier [paid carry family](../guarded_native_routing_20261023/CARRY_GAP.md)
also composes with the new cross-orbit cubes without discarding its kernel.
Append its source to T_r by ordinal sum. For any odd multiple M of3, the
selected state space is exactly the declared product

    {0,1}^r x Z_M x Z_(M/3),

with 2^r M^2/3 actual presentations and all endpoint-fixed transports kept.
The two rotation coordinates are the original full `(u,v)` lifts, not a
binary abbreviation. At least2^r source-automorphism orbits remain distinct,
because the first r rigid strong components cannot move.
The carry coordinates are actual owner lifts inside a source-automorphism
orbit; they do **not** provide M^2/3 additional normalized tile resources.
Its raw-state measure and normalized capacity are separate quantities.
The workspace theorem instead uses rigid cube blocks, whose tile resources
are pairwise distinct. No carry multiplicity is used to pay a port demand.

Choose a cube update with probability2r/(2r+M^2) and the inherited carry
update with probabilityM^2/(2r+M^2). Tensor-product variance and the earlier
gap>=2/M^2 give

    gap >= 2/(2r+M^2),
    weighted Metropolis gap >= 2/(kappa(2r+M^2)).

Here r>=1; a zero-factor endpoint uses the remaining original chain. A direct
route has at most r+floor(M/2)+floor(M/6) events. Cube events change two
owner clocks; the carry generators change M or2M clocks. That cost is retained
instead of treating the two kinds of events as identical work. This is a
proof-level product corollary of the two authenticated factors, not a new
enumeration of the entire enlarged source fibre or its full automorphism group.

## 3. Hostiles and scope

- Complete charts can all be singletons. Coverage alone gives no movement.
- All charts can mix rapidly while a closed chart-label invariant gives
  global gap zero. Our unconnected pointed atlas is exactly of this form.
- Native paths can exist but all demand may cross a tiny-weight port. A
  polynomial diameter alone does not bound Lambda or a spectral gap.
- A legal abstract correction may require a port already owned elsewhere;
  neither a quotient equation nor a stationary marginal licenses that move.
- Mixing on a representation space does not give a decreasing rank for the
  original deterministic dynamics.

For Collatz, this framework identifies a sufficient *transport estimate* only
after actual source rows, native integer clocks and a complete accepted
space are specified. Universal entry or positive ROOT charge still needs
an independent arithmetic bound. Adding workspace to a tournament is not
permission to change a Collatz integer or its deterministic forward orbit.

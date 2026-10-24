# Exact conditional weights remove the section bottleneck

Status: **PROVED** finite-state lifting identity and **CONDITIONAL** native
Markov-chain construction. It preserves any supplied strictly positive law;
it does not assert a uniform bound for an unspecified natural transition law.

## 1. Keep every witness and its weight

For a finite feasible fibre Omega, project a full record x to a visible record
j, retaining its entire compatible reservoir fibre Omega_j. Let pi(x)>0 be
the actual prescribed normalized law. Define

    bar_pi(j)=sum_{x in Omega_j} pi(x),
    pi_j(x)=pi(x)/bar_pi(j).

Do not replace bar_pi by a count or pi_j by a point mass at a charge-section
witness. Let bar_P be ANY reversible Markov kernel for bar_pi on the feasible
visible graph, and define

    P(x,y)=bar_P(j,l) pi_l(y),   x in Omega_j, y in Omega_l.   (1)

This refreshes the complete destination fibre, also on a visible holding
step. Positivity and row sums are immediate. Detailed balance follows from

    pi(x)P(x,y)=bar_pi(j)bar_P(j,l)pi_j(x)pi_l(y).

For f, write bar_f(j)=E_{pi_j}f. Equation(1) gives
Pf(x)=bar_P bar_f(j). Functions of conditional mean zero are annihilated;
block-constant functions carry exactly bar_P. More explicitly,

    Var_pi(f)=E_bar_pi Var_{pi_j}(f)+Var_bar_pi(bar_f),
    E_P(f,f)=E_bar_pi Var_{pi_j}(f)+E_bar_P(bar_f,bar_f).

If any fibre has more than one state, the Poincare gap is
**min(1,gap(bar_P))**; if all fibres are singletons it is gap(bar_P).
Use the one-point convention when necessary. For a lazy nontrivial bar_P
the lifted gap therefore equals gap(bar_P) in either case.
There is no additional reciprocal-smallest-witness-weight penalty. This is
an exact algebraic identity, not an independence approximation. It requires
that every cross-fibre transition used by(1) actually be allowed.

## 2. The native interior graph satisfies that requirement

Fix a binary fibre on80/81 sources, with5<=k<=n-5. The visible record is the
actual V subset J of E; its reservoir grade is r=k-|J| and its required charge
is q-ell(J). Retain only J with nonempty reservoir fibre. Use these undirected
visible edges, and add holding probabilities:

* from r<5 to grade5 by deleting5-r selected outside slots;
* from r>15 to grade15 by adding r-15 outside slots;
* between consecutive reservoir grades inside[5,15] by one outside change;
* between outside Johnson neighbors at reservoir grade5.

The proof in `THEOREM.md` shows that this visible graph is connected (at high
k the terminal grade can be r*>5 and has the unique outside state J=E).
Crucially, its event estimates hold for **every** reservoir representative
at either end, not only a selected section. Thus every transition in(1)
costs at most20 changed sources and45 primary owner clocks. Within-fibre
refresh costs at most20 sources and40 clocks. All source weights, pins,
inverse/converse units and original dyadic guards remain as declared.

One may choose bar_P by Metropolis rates on this graph using bar_pi itself.
That constructs an authenticated weighted chain. It is not a comparison
theorem for some other chain whose rates were fixed in advance. For that
problem, compare its actual conductances with those of(1), or retain the
canonical-path congestion inequality with the actual Q(e).

At k<=4 and its complement, the small-population route goes through a cheap
full-state witness; the theorem does NOT license unrestricted45-clock
heat-bath jumps between every pair of such states. These boundary fibres
have their own weighted feasible graph.

## 3. Exact source-weighted transfer polynomial

If each source supplies positive statistical weights a_i(U),a_i(V), the
conditioned law is proportional to their product over actual record choices.
For the20-source reservoir compute in the group algebra of G:

    Z(z) = product_{i in R} (a_i(U) + a_i(V) z e_{ell(i)}),
    Z_{r,g} = coefficient of z^r e_g in Z.

Then, WITHOUT unweighting any source,

    bar_pi(J) proportional to
      product_{i in E} a_i(choice_i(J)) * Z_{k-|J|, q-ell(J)}.

The same finite convolution used for support becomes a weighted convolution.
It retains the integer occupancy clock r and the charge separately. For
positive weights, Z_{r,g}>0 iff the authenticated support is nonempty;
no lower bound independent of those weights follows from support alone.

For the full85-source fibre, the W factor from `THEOREM.md` is independent
of the binary factor under this **product-source law**, even after fixed
target margins and currents: the fresh contribution has constant cardinality
and charge throughout its one/two-state factor. A random scan that runs the
binary chain with probability eta and the exact W heat bath with1-eta has
gap min(eta*gamma_binary,1-eta) when W has two states. Its W transition costs
four sources/23 clocks. Choosing eta=1/(1+gamma_binary) gives gap
gamma_binary/(1+gamma_binary). For an arbitrary correlated law on the two
factors, that product gap formula is not asserted.

Statistical weights here are distinct from the native physical edge weight2.
Reweighting physical electrical edges changes the harmonic equations; this
construction does not do so.

## 4. The remaining obstruction is visible and genuine

Take two visible states of equal mass with transition probability epsilon
each way. Its gap is2epsilon. Adding arbitrarily rich positive reservoir
fibres and applying(1) leaves that slow mode unchanged. Thus no uniform
weighted routing theorem can follow solely from complete charge coverage.
The reroute removes artificial concentration on one chosen witness; the
remaining target is a source-specific lower bound on the **visible** chain's
conductance (and a map to the intended ROOT process if Collatz is the goal).

Connection contract: full weighted records -> conditional fibre law and
visible marginal. Preserved: exact stationary mass, all feasible records,
variance and Dirichlet decompositions. Lost by the visible state alone:
within-fibre fluctuations, restored by the conditional law. Cheapest hostile:
nonuniform fibres and an epsilon visible cut. `weighted.py` checks detailed
balance and the two identities in exact rational arithmetic on that control.

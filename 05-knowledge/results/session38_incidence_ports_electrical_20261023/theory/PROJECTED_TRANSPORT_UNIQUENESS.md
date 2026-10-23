# Exact projected-wiring uniqueness from residual-cycle potentials

Status: PROVED elementary finite-flow theorem, with explicit scope. This is an all-size certificate for a specified bucket-transport relaxation; it does not assert that every feasible table is a genuine tournament metagraph.

## 1. Inheritance and target

The current `metagraph_bucket_obstruction_20261023` packet proves a capacity-peeling sufficient criterion. It explicitly warns that a stuck residual is not evidence of another transport table. The older order-six census also separates multiple ordinary port tables from their projection to the same merged graph.

Here the target is neither automatically a unique table nor automatically a unique individual port matching. It is uniqueness of the prescribed projection to merged coloured edge multiplicities.

## 2. Abstract transport model

Let U,V be finite row and column sets and E a finite multiset of allowed cells e=(u,v). Distinct cells may share endpoints, if they represent different accepted channels. Each cell has an integer upper bound b_e>=0, possibly infinite. Prescribe nonnegative integer row and column margins. Suppose an integer feasible table x is given:

0<=x_e<=b_e,
sum_(e out of u)x_e=r_u,
sum_(e into v)x_e=c_v.

Choose a linear observation P:Z^E->Z^d, with column p_e for cell e. For coloured merged wiring, p_e is its contribution to the selected merged node-pair/colour coordinate, with the literal multiplicity retained.

Form the directed residual multigraph R_x on U disjoint union V:

- u->v for cell e if x_e<b_e, with capacity b_e−x_e and vector cost+p_e;
- v->u if x_e>0, with capacity x_e and vector cost−p_e.

Keep the cell label on each arc. A two-edge walk using a cell and its own reverse has zero cost; distinct parallel channels need not cancel. No connectedness assumption is made.

## 3. Necessary and sufficient criterion

THEOREM. The following are equivalent:

1. Every feasible integer table y has P(y)=P(x).
2. Every directed cycle of R_x has total vector cost zero.
3. In each strongly connected component of R_x there is a vector-valued vertex potential h such that every internal residual arc a satisfies

cost(a)=h(head(a))−h(tail(a)).

Thus projected wiring uniqueness is exactly absence of a nonzero observed residual circulation.

Proof of1=>2. A directed simple cycle has positive integer residual capacity on every arc. Augment one unit around it. Row and column balances are unchanged, all lower/upper bounds remain satisfied, and its projected change is precisely the cycle cost. If the cost is nonzero this constructs an integer feasible counterexample to1. Cycles containing a forward/reverse copy of the same cell with no net change have cost zero and cause no exception.

Proof of2=>1. Given feasible y, direct each positive difference y_e−x_e forward and each negative difference backward. These magnitudes fit the residual capacities and form an integral circulation, because all margins agree. Decompose it into directed cycles. Condition2 makes every summand's projection zero, hence P(y−x)=0.

Proof of2<=>3. Potential differences telescope around cycles. Conversely fix a root in each residual SCC and define h(v) by summing arc costs along any directed root-to-v path. Two such paths, each followed by the same directed v-to-root path, have equal cost because all closed directed walks decompose into cycles. Therefore h is well-defined and each arc satisfies the required difference. QED.

The vector group can be replaced by any abelian group for the cycle argument, provided equality of observations is interpreted in that group. Taking coefficients modulo2 tests only parity uniqueness, not integer multiplicity uniqueness.

## 4. Constructive certificates and minimal witnesses

For uniqueness, retain x, the SCC partition and one potential per vertex. A verifier checks the margins/bounds and every internal-arc potential equation. No enumeration of all feasible tables is needed. Arcs between SCCs cannot belong to a circulation and need no potential equation.

For nonuniqueness, retain one positive-capacity directed cycle with nonzero vector sum. Its unit augmentation is an exact counterfeit table with a changed projection. Checking a candidate's colour/node counts alone is not the certificate: the actual allowed-cell compatibility and residual signs are required.

A failed potential equation can be converted to a cycle witness. Follow a root-to-tail path, the failed arc and a return path to the root; compare with a root-to-head path followed by the same return. At least one closed walk has nonzero cost, and decomposing it yields a nonzero cycle.

The certificate concerns a chosen P. Setting P to the identity on table entries tests table uniqueness. A coarser merged-wiring P can be unique even when tables are not. This is the exact distinction missed by equating transport-table counts with different metagraph counts.

## 5. Why peeling can stall at a unique solution

Consider allowed cells forming a six-cycle with rows u1,u2,u3, columns v1,v2,v3:

(u1,v1),(u2,v1),(u2,v2),(u3,v2),(u3,v3),(u1,v3).

Set x=(0,1,1,0,1,1), giving row margins(1,2,1) and column margins(1,1,2). All margins are positive and every allowed vertex has degree2, so degree-one peeling alone cannot start. Any balanced change is proportional to

(1,−1,1,−1,1,−1).

The first zero forces its parameter>=0; the fourth zero forces it<=0. Thus x is the unique table. More powerful tight-cut peeling may detect this too; the example refutes only the inference 'no degree-one leaf implies an alternative'. The residual-cycle theorem gives the exact criterion regardless of which sufficient peeling routine was used.

## 6. Application to complement-line metagraphs

Use closed versus open fixed-path ports as the two sides; complement switches the path-closing arc. Blue and black are separate colour layers. To fold black twins, first form genuine reflection orbits of black ports; each has size2. A unit black-orbit match contributes two literal black lines, so its merged projection column is2 times the appropriate edge coordinate. Blue contributes one. Loops of the merged graph are still ordinary cells between distinct closed/open port sides and must be retained.

This folding requires the allowed cells to be reflection-compatible, including any accepted orientation channel between two port orbits. It must be proved for the actual chosen bucket definition; an ordinary class-level table with extra mirror equalities cannot simply be declared an unconstrained transportation problem.

Once the authenticated actual table belongs to the chosen relaxation, a PASS of the theorem forces its projected wiring among all tables allowed by that relaxation, at every order. A FAIL supplies another relaxed table, not automatically another realizable tournament wiring. If additional global equalities or nonlinear conditions are imposed beyond allowed cells, margins and cell bounds, the theorem as stated does not classify that further intersection.

The all-size rule is therefore exact and useful but typed: conserved residual costs certify the information retained by one specified incidence menu. It is not a universal assertion that a fixed small motif signature suffices forever.

# Source-dependent cutoffs and the exact value of signed residual credit

**PROVED finite duality and transport rules / FINITE-EXACT controls /
OPEN universal Collatz estimate.** Same-agent proof review. The old finite
source-dual obstruction is recovered before stating the new specialization.

## 1. Fixed source law, forward operator and integer clock

Let T be the accelerated odd Collatz map, stopped on H={1}. Work on the
active positive odd states Y. Put

    w(2j-1)=1/[j(j+1)], j>=2;
    (Qv)(y)=sum_(x in Y:T(x)=y) v(x), y in Y;
    a_z(x)=1-z^tau(x), 0<z<1,

where tau is the number of actual odd steps to H, and z^infinity=0.
Extend a_z to be zero on H. The active mass is1/2 and every active atom
has positive weight. The identities below hold for any deterministic
killed map and any nonnegative summable source law; no random reset law
is substituted for this one. We use the ordinary positive 3n+1 map,
not a changed-sign or 3-adic surrogate.

Write c=1-z. Then a_z(x)=c+z a_z(Tx) on Y. For any summable signed trial
flow v define

    e=w+zQv-v.

The exact residual identity is

    <w,a_z>=c sum v+<e,a_z>.                                  (1)

If authenticated bounds 0<=l<=a_z<=u<=1 are available, (1) gives

    B(v)=c sum v+sum e_+ u-sum e_- l >= <w,a_z>.               (2)

This is the signed residual rule already established in the
[previous refinement packet](../proof_targets_refinement_20261029/README.md).
The present question is what finite optimization of (2) can actually add.

The older [Poisson source-dual result](../collatz_poisson_source_dual_20261005.md)
already identifies a finite-support path obstruction for its contracted
inverse-base operator. That operator and its weights are different. We
recover its warning, not its numerical constants, and now prove the
exact forward-operator statement needed for the current source law.

## 2. The finite boundary-response theorem

Fix a finite active set F, allow v to be signed with support in F, and
let W=F union (T(F) intersect Y). Set

    C_out=sum_(x outside W) w(x) u(x).

**Theorem C1.** The exact optimal score is

    inf_(supp v subset F) B(v)
      = C_out + max sum_(x in W) w(x) A(x),                   (3)

where the maximum ranges over l<=A<=u on W subject to

    A(x)=c+z A(Tx), x in F;  A|H=0.                          (4)

If v is constrained nonnegative, replace the equalities in (4) by
A(x)<=c+z A(Tx).

Proof. For a residual coordinate e, e_+u-e_-l=max_(l<=A<=u)eA.
Thus B is the support-function expression

    C_out + sup_(l<=A<=u) [<w,A>+
                sum_(x in F) v(x)(c+zA(Tx)-A(x))].            (5)

Finite linear-program duality gives (3). Feasibility is supplied by the
actual a_z; the box is bounded. For a nonnegative primal variable, the
coefficient in (5) must instead be nonnegative, giving the stated
inequality. This also follows by introducing positive and negative
residual variables in an ordinary finite LP. QED.

The dual retains the endpoint identities. Replacing them with a trace,
a new endpoint law, or separate uncorrelated copies changes the problem.
Shared endpoints impose a SINGLE value A there, with the genuine source
weights accumulated, while the provenance ledger can retain each source.

### Theorem C2: a precise no-gain region

Suppose, in addition, u(x)>=c+z u(Tx) for every x in F, with u|H=0.
Then an optimal flow in (3) is nonnegative:

    v_F=(I-zQ_F)^(-1) w_F,   v=0 outside F.                  (6)

In particular, with the coarse upper bound u=1 on Y, signed finite
optimization has no advantage over the best positive stopped flow.

Proof. Give the boundary outside F the values u. Inside F let h be
the unique solution h(x)=c+zh(Tx). Equivalently, follow the actual path
until it exits F, discounting its integer exit time; an orbit trapped
in F has h=1. Positivity and the supersolution assumption imply
a_z<=h<=u, hence l<=h<=u. Every other feasible A in (3) is at most h,
by iteration to an exit, or by letting z^L tend to zero along a trapped
orbit. Thus h is a maximizing dual solution. The inverse in (6) is
the nonnegative convergent Neumann series. Its residual is zero on F
and nonnegative outside. Pairing with h shows that (6) attains the
dual optimum. QED.

For u=1, this has a particularly sharp form. Let E_F be the original
sources whose entire active path to H lies inside F. Then

    inf B(v)=sum_Y w-sum_(x in E_F) w(x) z^tau(x).            (6a)

Any other source has boundary response one, whether its path exits F
or cycles there. Thus this finite optimization gives exactly the coverage
already authenticated by those paths. It cannot silently turn a long
unfinished excursion into new proof coverage.

This does not refute improvements over a PARTICULAR positive trial,
nor improvements using independently sharper upper data which fail the
supersolution test. It does rule out seeking a new convergence theorem
merely by optimizing signed coefficients on a finite graph with u=1.
An infinite arithmetic ansatz is outside this finite theorem's scope.

## 3. A paid reverse-clock rule where signed credit really helps

Keep an authenticated legal path x_0 -> ... -> x_L, with no earlier
terminal visit. Its actual length L gives

    a_z(x_0)=1-z^L+z^L a_z(x_L).

Therefore an upper bound at its source gives a forward endpoint bound

    u_new(x_L)=1-z^(-L)(1-u(x_0)).                           (7)

Take the minimum with the existing endpoint upper bound. Authentication
ensures (7) is nonnegative. The inverse discount factor z^(-L) is part
of the cost; the path is not a one-step reset.

This update has an exact signed-flow implementation. For b>=0 put

    Delta v=-b sum_(i=0..L-1) z^(i-L) delta_(x_i).

Telescoping, including repeated intermediate vertices, gives

    Delta e=b z^(-L) delta_(x_0)-b delta_(x_L),
    c sum Delta v=-b(z^(-L)-1).                              (8)

For distinct endpoints with current e(x_0)>=0 and e(x_L)>=b, both
residuals stay nonnegative and the exact change in score (2) is

    Delta B=b[u_new(x_L)-u(x_L)].                            (9)

These guards matter: after crossing zero, the lower bound l determines
the price of negative residual. Formula (8) remains valid without the
guards, but (9) must then be replaced by the piecewise-linear score.
One may not spend the same positive endpoint residual twice.

This is an algebraic response transfer along a verified forward path,
not a claim that the endpoint's own forward orbit traverses that path
backwards. A source-specific implementation tags each spent residual
portion with its original source and applies (8) to that tagged portion.
Aggregation then sums these ledgers; it does not replace their original
weights or their integer clocks with a refreshed distribution.

For the actual path 3->5->1, choose z=1/2 and F={3}. Store the true
two-step bound u(3)=3/4 but initially use u(5)=1, l=0. The original
weights are w(3)=1/6,w(5)=1/12; all other active sources have mass1/4.
The zero-flow bound is11/24. Choosing b=1/12 in the one-step path3->5
gives v(3)=-1/6, residual e(3)=1/3,e(5)=0, and bound5/12: a gain1/24.
No nonnegative flow supported on {3} improves11/24. The missing endpoint
bound u(5)=1/2 is exactly (7). This imports a known suffix from an
authenticated source receipt; it creates no new convergent source.

## 4. Source-dependent cutoffs without replacing source weights

Choose a finite set J of original indices j>=2 and individual integer
cutoffs L_j. Follow each original source x_j=2j-1 for
t_j=min(L_j,tau(x_j)) steps, keeping its path and terminal flag. Then

    v=sum_(j in J) w(x_j) sum_(t=0..t_j-1) z^t delta_(T^t x_j)

has residual consisting exactly of unselected original atoms and the
unfinished endpoint masses w(x_j)z^(L_j). Overlapping paths add their
actual weights; they are not sampled again. A valid upper score is

    sum_(j in J) w(x_j)[1-z^(t_j)]
      + sum_(unfinished j in J) w(x_j)z^(L_j)u(T^(L_j)x_j)
      + sum_(j not in J, j>=2) w(x_j)u(x_j).                 (10)

If J={2,...,J_max}, the unselected mass is exactly1/(J_max+1).
Additional paths and signed updates (8) may be merged into the same
finite compiler while preserving their source and clock ledger.

More generally, the outside term in (3) can be certified by

    C_out <= 1/(J_max+1)
       + sum_(2<=j<=J_max: x_j outside W) w(x_j)u(x_j).      (11)

This separates an original-source cutoff from a potentially very large
endpoint cutoff. It needs no false bound on how high a convergent orbit
may travel. The preceding packet's
[finite-depth tail obstruction](../pinned_discrepancy_reset_tails_20261029/FINITE_DEPTH.md)
already shows why a universal constant/J endpoint envelope is stronger
than convergence. Equations (10)--(11) reroute the finite compiler
around that unnecessary global-height target.

The remaining sufficient target is an independent family of authenticated
endpoint/source bounds and finite flows for which B_z tends to zero as
z increases to one. Then dominated convergence in <w,a_z> and positivity
of every original atom imply that every source terminates. We have NOT
proved that estimate. Choosing cutoffs by assuming the unknown stopping
times, or taking u to be the unknown exact a_z, would be circular.

## 5. Verification universe

`check_cutoffs.py` tests all64 deterministic maps on three active states,
each allowed to map to a terminal, for three supports and three authentic
boxes:576 rational primal/dual comparisons. An LP solver proposes a flow;
its rational reconstruction is checked against a separate exact SymPy
enumeration of box/equality vertices. The two coarse-upper modes give384
explicit nonnegative stopped-flow optima. Cycles, self-loops, mergers,
exits and terminal edges are all included. The actual3->5 receipt above
is checked separately. This finite universe corroborates the algebraic
theorem; it does not test universal Collatz convergence.

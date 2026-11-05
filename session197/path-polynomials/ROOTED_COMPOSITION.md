# Root-preserving context composition adds matching rank

Session 197. **All-order rooted composition law and a precise obstruction to repairing a clone clock by outer whole-module wrappers.** This is not an obstruction to rerooting, selective nonuniform attachments, or target-dependent new templates.

## 1. A physical composition identity

Take rooted tournaments (P,p) and (Q,q). Substitute Q as the entire p-module of P, and retain q as the distinguished root of the result R=P[p:Q]. Then for every further input S,

    R[q:S]=P[p:Q[q:S]].

Consequently the exact bin-profile context operators satisfy

    A_(R,q)=A_(P,p) composed with A_(Q,q).                     (1)

The order here matters: the inner Q operator acts first. This equality preserves every arc, not only the Hamilton count.

Let B_(Q,q) be the root's bipartite cross graph, with O_q on one side and I_q on the other, and edge u-v exactly for u->v from O_q to I_q. Denote its matching number by nu_Q and its number of maximum matchings by mu_Q. Use analogous notation for P and R.

**Theorem 1.** There is an exact bipartite disjoint-union identification

    B_(R,q)=B_(Q,q) disjoint-union B_(P,p).                     (2)

Thus

    nu_R=nu_Q+nu_P,       mu_R=mu_Q mu_P.                      (3)

**Proof.** The old Q vertices keep their incidences to q. An outer P vertex sees q exactly as it saw p, so its out/in classification is that in the P cross graph. There are no cross edges between these two cross-graph parts. Indeed an outer in-neighbor of p beats every vertex of Q, so an old Q out-neighbor cannot beat it. An outer out-neighbor loses to every vertex of Q, so it cannot beat an old Q in-neighbor. These are exactly the two possible off-diagonal O-to-I blocks, and both vanish. The old within-block cross edges remain unchanged. Maximum matchings in a disjoint union are independent maximum matchings in its two parts. QED.

## 2. The leading shift coefficients compose without cancellation

Write

    A_P=sum_i p_i(x)S^i,       A_Q=sum_j q_j(x)S^j,

where Sf(x)=f(x+1). The exact shift-bandwidth theorem in `CONSEQUENCES.md` gives nonzero coefficientwise-positive top polynomials p_(nu_P),q_(nu_Q). Since S^i q_j(X)=q_j(X+i)S^i, equation (1) gives

    [S^(nu_P+nu_Q)]A_R
       =p_(nu_P)(x) q_(nu_Q)(x+nu_P).                         (4)

Both factors have nonnegative integer coefficients and are nonzero, so this top mode cannot cancel. Their degrees and leading coefficients agree with the maximum-matching theorem: N_R-nu_R=(N_P-nu_P)+(N_Q-nu_Q), and the leading coefficient is mu_P mu_Q.

When both templates are balanced, the rising-factorial formula is compatible term by term:

    mu_P x^(overline nu_P) · mu_Q (x+nu_P)^(overline nu_Q)
      =mu_P mu_Q x^(overline(nu_P+nu_Q)).

This gives a useful exact check that the normal-form shift convention and physical matching composition use the same root.

## 3. What an outer wrapper can and cannot repair

Suppose Q[r:Tr_m], m>=1, is a fixed rooted clone family with matching rank nu. Consider any fixed finite sequence of whole-input substitutions around it, retaining the same distinguished root inside the original Q. The new matching rank is

    nu + sum_(outer contexts) nu_i.

Every active outer context (its root lies in a directed triangle) contributes at least one. The resulting clone family therefore has a strictly larger highest exponential base than nu+1, with positive coefficient. No lower exponential term cancels that mode.

If the highest base is required to remain nu+1, every outer context must have rank zero. Such a root lies in no directed cycle and is a singleton strong component. Its entire context operator is multiplication by B_(P-p)(x), and its scalar action is multiplication by the positive integer H(P-p). Thus all coefficients of the original clone ray are multiplied by the same positive integer.

**Corollary.** A pure base-three ray A*3^m+beta with A>0 and beta<0 cannot be changed to a pure base-three ray with positive offset by any fixed sequence of root-preserving whole-input wrappers. An active wrapper introduces a base at least four; passive wrappers preserve the sign of beta.

The statement fixes the original clock root and a finite word of outer templates independent of m. It does not apply to choosing a new root, modifying exterior incidences to different parts of Q separately, replacing another internal vertex, or choosing an unbounded target-dependent template bank. In particular it does not dispose of the mixed-port constructions investigated elsewhere in Session 197.

## 4. A concrete mode-growth control

For the rooted tournament

    Q=[2,12,145,149,99,143,175,19],       root q=7,

the independently compiled exact operator has nu=2, mu=4 and highest shift coefficient

    q_2(x)=3x+14x^2+23x^3+16x^4+4x^5.

Thus q_2(1)=60 and q_2(2)=630. Its scalar positive-m clone response is 60*3^m-1. At x=2 the scalar shift evaluation is

    B_(Q[q:Tr_m])(2)=630*4^m+480*3^m-176*2^m,     m>=1.

The omitted zero-base coefficient is six and matters only for the formal m=0 boundary. The verifier emits the exact scalar shift coefficient tables, not a fitted recurrence.

Applying the physical D and E wrappers from the companion packet gives respectively

    630*4^m+540*3^m-176*2^m-1,
    1890*4^m+1380*3^m-528*2^m+1.

The second operation really changes the constant term to +1, but it also creates the unavoidable positive 4^m term. This illustrates the obstruction without mistaking a signed algebraic cancellation for a graph operation.

These explicit finite-template identities are certified by the complete skeleton operator and independent actual clone counts. No statement about arbitrary pure-base-three templates is inferred from this example beyond the all-order composition theorem above.

## 5. Verification

Run `python3 verify_rooted_composition.py` and the optimized version, comparing stdout. All 841 ordered pairs of rooted labelled templates of orders at most three are checked. The checker independently constructs the composite graph, compares its full skeleton operator with noncommutative shift-polynomial multiplication, verifies matching-rank addition and maximum-matching multiplication, and checks the top coefficient formula. Higher-rank explicit controls include the displayed root-seven clock and actual outer-wrapper clone graphs. These tests support the physical proof; they do not substitute for it.

# Arithmetic termination of the cheaper ternary controller

**Status:** PROVED, elementary positive-integer theorem. FINITE-EXACT replay accompanies it. No external novelty claim; not a Collatz convergence theorem.

## Inheritance and exact scope

Source inspected at `FETCH_HEAD=5d6eb91dd76eee86b91e04d856daa32c32b3d61d`:
`05-knowledge/results/affine_tree_recurrences_20261028/{MOMENTS,REFINEMENT,RECURRENCES}.md`, introduced by commits `35e43c384e7`, `1d4615c7bf4`, `b841e2c6b56`.
The cheaper rule is **(d,r)=(3,3)**, not the critical parity rule. A native accelerated Collatz step is

    a = v_2(3n+1),   n'=(3n+1)/2^a.

Start Q=1. Add 2 when 3 divides a, and subtract 1 otherwise; stop at the first Q=0. The earlier packet proved Haar-almost-sure termination, finite mean carry, and exponentially decaying conditioned tails, while explicitly leaving a possible exceptional positive integer. Here the positive-integer guard removes that particular exceptional set without appealing to measure.

Concept board: arithmetic guard / queue resource / affine source carry / native endpoint / clock. Hostiles: 3→5, parity-grammar source41, and the nonintegral positive rational1/5.

## The arithmetic guard and a sharp potential

For positive odd integers define

    P(n,Q)=(n+1)(3/2)^Q.

**Lemma.** Every legal step, before stopping, weakly decreases P. Every branching step decreases it by a factor at most 27/28, with equality exactly at n=13, a=3.

For a leaf, a≥1 and Q'=Q−1. Thus

    P(n',Q')/P(n,Q) = (2/3)(3n+1+2^a)/(2^a(n+1)) ≤ 1.

Equality holds for a=1; for a≥2 the ratio is at most2/3, because n'≤n for every positive odd n and a≥2.

For a branch, a≥3. Crucially no positive odd integer less than13 has valuation a divisible by3: the valuations at 1,3,5,7,9,11 are respectively2,1,4,1,2,1. Hence n≥13, and

    (n'+1)/(n+1) ≤ (3n+9)/(8(n+1)) ≤ 3/7.

Consequently the potential ratio is at most (9/4)(3/7)=27/28. Both inequalities are equalities only at n=13,a=3. This is an arithmetic obstruction, not a density estimate.

## Every positive integer controller terminates

After j steps let b_j count branches. Then Q_j=1+3b_j−j. As long as Q_j≥1,

    3 ≤ P(n_j,Q_j) ≤ (3/2)(n_0+1)(27/28)^{b_j}.       (1)

An infinite active run would require b_j→∞ because Q_j≥1 implies3b_j≥j. This contradicts (1). Thus the first stop exists for **every positive odd n_0**. This assertion concerns the auxiliary controller, not convergence of the native orbit to1.

A finite, integer-only certificate is obtained by defining

    M(n)=max {m≥0 : 27^m(n+1) ≥ 2·28^m}.

M(n) exists and can be computed by integer multiplication. At the last branch the queue is still active, so the total branch count m satisfies m≤M(n_0). At stopping,

    J=3m+1 ≤3M(n_0)+1.

Furthermore P(n_j,Q_j)≤P(n_0,1), Q_j≥1 implies n_j≤n_0 at every active state, including every input to a step. Therefore each valuation is at most floor(log_2(3n_0+1)), yielding

    D ≤ (3M(n_0)+1) floor(log_2(3n_0+1)).

J counts accelerated odd steps, D counts divisions by2, and the ordinary Collatz clock is J+D. The branch budget is logarithmic in n_0; this loose division-clock bound is quadratic in log n_0.

## Source preservation, ports, and what may be composed

The controller returns the complete chronological word, source integer, source cylinder, endpoint, and clocks. Starting C=D=0, update

    C ←3C+2^D,   D←D+a.

After J letters, 2^D n_J=3^J n_0+C and the exact odd source guard is

    n_0 ≡ (2^D−C)3^{−J} mod2^{D+1}.

Every test uses the same n_0; no independent source is sampled at a return. Input/output ports can be the actual odd integers, or a richer marked state required by a later attachment. A residue quotient is not asserted to suffice for arbitrary future operations. On composing words u,v, retain

    C_uv=3^{J_v}C_u+2^{D_u}C_v,

and pull the second guard back through the first affine map. Initializing a new queue at the output is legal but changes its potential by3/2; it is not a free reset.

## Endpoint controller and remaining debt

At the terminal Q=0,

    n_J+1 ≤ (3/2)(n_0+1)(27/28)^m.                  (2)

Hence m≥12 guarantees strict native descent, since (3/2)(27/28)^12<1; m=11 does not suffice for this particular coarse budget. An exact carry-sensitive dispatcher accepts descent precisely when

    (2^D−3^J)n_0 > C.

Otherwise it returns the marked endpoint as unresolved, not as a new independent Haar source. Termination of this first controller has been proved; termination of repeated failure returns into a descending controller has **not**.

Source3 stops immediately at5, source27 at41, and source1 at1. For the old critical parity grammar, source41 followed (2,1,1,1) to107; for the present cheaper grammar it stops after valuation2 at31. The grammars must not be conflated.

The rational source1/5 has native valuation3 at every step and is fixed by (3x+1)/8. Its queue grows forever. It violates the integral lower guard n≥13, showing precisely why the proof excludes exceptional positive integers without proving emptiness of the corresponding 2-adic exceptional set.

## Reproduction and next serious extension

Run `python3 -B check_integer_termination.py` and repeat with `python3 -O -B`. No assert statements are used. The saved receipt checks all100001 positive odd edge inputs through200001, all10001 full sources through20001, full affine guards and clocks, the sharp branch constant, and the rational hostile.

Next target: classify the finite queue-shape family with fewer than12 branches using carry-sensitive source rays, and derive exact expanding-return dispatches. The family of shapes is finite, but valuations remain unbounded; a finite shape list alone is not a finite arithmetic classification. The simplest infinite obstruction to a fixed native-step budget is n=2^k−1: it has k−1 initial valuation1 steps, all expanding. Any universal descending return controller therefore needs an unbounded integer clock or an arithmetic compressed-run transition. Such a transition must carry v_2(n+1), not merely the queue shape.

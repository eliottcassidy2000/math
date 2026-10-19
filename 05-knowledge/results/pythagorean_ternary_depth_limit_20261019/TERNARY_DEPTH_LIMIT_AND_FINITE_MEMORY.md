# Finite-memory positivity brackets and an irrational ternary-depth limit

19 October 2026. Status: PROVED. This closes the earlier Pythagorean packet's specific gap concerning existence of a uniform ternary-depth limit. It does not change its distinct hypotenuse-cutoff density of one third. The general equal-cylinder coding mechanism is classical; no historical novelty claim is made for that mechanism.

Use the notation and outer-first convention of the [positivity-language proof](../dynamics_labeling_entry_20261019/BERGGREN_POSITIVITY_NONREGULAR.md). Thus D=(1,infinity), the branch maps are L(t)=2-1/t, M(t)=2+1/t, R(t)=t+2, and positivity means t in A=(1,mu) union (lambda,infinity). Write p_n for the fraction of the 3^n chronological words of length n whose primitive Pythagorean node gives a positive AP square-sum triangle. Reversal is a bijection of words at each level, so the same fraction is obtained in the outer-first convention.

## 1. Two ambiguous cylinders at every positive depth

At depth d, the 3^d intervals f_u(D), |u|=d, have disjoint interiors and are linearly ordered. They cover D except for finitely many rational boundary points. Every irrational point belongs to exactly one interval.

For d>=1, exactly two cylinders cross a positivity boundary: the unique cylinder containing mu and the unique cylinder containing lambda. They are distinct already at depth one, because mu lies in L(D) and lambda in R(D). Every other cylinder lies wholly inside A or wholly outside A.

### Theorem 1: sound finite-suffix classifier

Fix d>=1. For a chronological word of length n>=d, inspect its last d letters, reversed to form an outer-first word u. Return:

- positive if f_u(D) is wholly contained in A;
- negative if f_u(D) is disjoint from A;
- unknown for the two boundary-crossing cylinders.

Every definite answer is correct for every possible earlier prefix. Under the uniform distribution on chronological words of any fixed length n>=d, the fraction marked unknown is exactly 2/3^d.

Proof. The earlier prefix has some slope t in D. The final suffix sends it to f_u(t), so interval containment proves the definite answers. Every suffix of length d occurs on exactly 3^(n-d) words, and there are exactly two unknown suffixes. QED.

This is an abstention fraction, not an error rate. In particular, when n=d the suffix is the entire word and direct evaluation would decide every word. The interval classifier deliberately asks for a verdict valid for every possible earlier prefix, of arbitrary length.

### Optimality within that precise suffix-only class

Each of the two unknown suffixes really has both positive and negative finite-word completions. Its interval crosses exactly one wall, so the inverse images of the two sides contain nonempty open intervals in D. Dense finite prefix-node slopes occur on both sides. These give two actual earlier prefixes whose concatenations with the same suffix have opposite positivity values.

Therefore a classifier depending only on the last d letters, required to be sound uniformly over arbitrary earlier-prefix lengths, cannot give a definite answer on either ambiguous suffix. This is not a lower bound for a classifier also supplied a fixed total word length, or for every possible finite-state architecture.

For words shorter than d, classify by their exact finite lookup. Treating unknown as negative gives a regular inner language I_d; treating it as positive gives a regular outer language O_d. They satisfy

    I_d subset I_(d+1) subset positivity language
        subset O_(d+1) subset O_d,
    union_d I_d = intersection_d O_d = positivity language.

The disagreement set at every length n>=d has relative size exactly 2/3^d. These statements follow by refining each boundary cylinder into three children; exactly one child retains its wall. Regularity follows from a finite suffix test together with finitely many short-word exceptions. A reverse-reading automaton can implement the interval test using at most 2d+2 states, including the positive, negative, and unknown sinks; this bound is for reading the chronological word from right to left, not a claimed bound for an equivalent deterministic left-to-right automaton.

## 2. An orientation-aware ternary coordinate

For an irrational alpha in D, let b_d(alpha) be the zero-based position of its containing interval among the 3^d depth-d cylinders, ordered from left to right. Then

    b_(d+1)(alpha)=3 b_d(alpha)+a_(d+1)(alpha),
    a_j(alpha) in {0,1,2}.

The digits are obtained from the inverse branch address and the orientation of the preceding prefix:

    increasing prefix: L -> 0, M -> 1, R -> 2;
    decreasing prefix: L -> 2, M -> 1, R -> 0.

Reading M reverses the orientation; L and R preserve it. Consequently the limit

    B(alpha)=lim_(d->infinity) b_d(alpha)/3^d
            =sum_(j>=1) a_j(alpha)/3^j

exists, and

    b_d(alpha)/3^d <= B(alpha) <= (b_d(alpha)+1)/3^d.

For irrational alpha the upper and lower endpoints cannot coincide with B(alpha): an eventually constant string of 0s or 2s would force an eventually constant inverse branch address and a rational limiting endpoint. The weak inequalities are sufficient for all estimates below.

The branch digits and the ternary digits determine each other using only the current orientation. In particular, eventual periodicity of the ternary digits forces eventual periodicity of the branch address: the orientation state has only two values, so its repetition period can at most double.

If the inverse branch address is eventually periodic, two inverse iterates have the same infinite address. The shrinking-cylinder lemma from the previous note makes their values equal. A nonempty rational Möbius word then fixes that value, forcing alpha to be quadratic if it is irrational. Therefore

    alpha irrational and not quadratic => B(alpha) irrational.

This direction is elementary. With the attributed Romik Lagrange theorem, the converse also holds for irrational alpha, but it is not needed for the positivity result.

## 3. The exact symmetry between the two positivity walls

Let

    S(t)=(t+1)/(t-1).

It is a decreasing involution of D, and direct substitution gives

    S composed with L = R composed with S,
    S composed with R = L composed with S,
    S composed with M = M composed with S,
    S(mu)=lambda.

Thus S reverses the order of the depth-d cylinder partition while permuting it by interchanging L and R in every word. It follows that

    b_d(lambda)=3^d-1-b_d(mu).

There are therefore exactly 2 b_d(mu) wholly positive depth-d cylinders, and exactly two undecided cylinders. This is a symmetry of interval cylinders, not of finite root-node levels: S(2)=3, so the root 2 is not fixed. The two undecided cylinders account for any discrepancy in the actual finite-level node counts.

## 4. Existence, exponential convergence, and irrationality

### Theorem 2: the ternary-depth density

The uniform ternary-level proportions have a limit

    p_n -> rho=2 B(mu),

and for every n>=1,

    2 b_n(mu)/3^n <= p_n <= [2 b_n(mu)+2]/3^n,
    |p_n-rho| <= 2/3^n.

Moreover rho is irrational. Its first 50 decimal places, certified by the exact depth-128 interval, are

    rho = 0.06039052312486240241155650348537631073311780098403... .

Proof. At any level n>=d, each depth-d prefix has exactly 3^(n-d) extensions. Wholly positive cylinders contribute all their extensions; wholly negative cylinders none. Hence p_n lies between 2b_d(mu)/3^d and [2b_d(mu)+2]/3^d. The same interval contains 2B(mu). Let d tend to infinity to obtain the limit, and then set d=n for the explicit error bound.

The wall mu has degree four over Q, proved in the earlier note. Section 2 therefore shows B(mu) is irrational, and its nonzero rational multiple rho=2B(mu) is irrational as well. QED.

The proof neither assumes equidistribution of Euclid slopes in angle nor uses the hypotenuse-cutoff counting measure. Uniform word-depth sampling has its own equal-cylinder measure, and the theorem computes its positivity mass.

The earlier exact counts at depths 3 through 10,

    2,5,15,44,132,396,1188,3566,

are consistent with the theorem, but do not prove the limit or its irrationality. The companion checker reconstructs them and verifies the two-cylinder bracket. It also supplies a certified rational interval for rho at depth 128, of width 2/3^128, with no floating-point decision in the boundary itinerary.

## 5. Scope and next questions

The nonregularity theorem says that no fixed finite state space decides every word exactly. The present theorem gives exponentially narrowing regular brackets and an irrational limiting level density. These are compatible: finite memory can be quantitatively excellent under uniform sampling while remaining insufficient on adversarial boundary-following words.

This theorem concerns primitive Pythagorean shapes and their protected AP square-sum triangles. It does not imply a density of interval Hamiltonian cycles, a uniform bound on the size of an elliptic neighbor, or a Collatz entry theorem. The respective measures and admissibility conditions remain separate.

A concrete next arithmetic question is whether the constant rho is transcendental. Irrationality follows from nonperiodic ternary coding; transcendence does not. Algebraic degree four of mu does not transfer to algebraicity or transcendence of B(mu), and no such claim is made.

# Bounded cone entry with nonregular membership in the Berggren tree

19 October 2026. Status: PROVED. The nonregularity theorem and the finite-state cutoff criterion have self-contained proofs below. The converse for quadratic walls additionally invokes the cited Romik Lagrange theorem. The Pythagorean parametrization and the general automata mechanism are classical; no historical novelty claim is made.

## 1. The exact question and word convention

Let D=(1,infinity), with root slope t=2, and put

    L(t)=2-1/t,    M(t)=2+1/t,    R(t)=t+2.

These are the Euclid-parameter children of the primitive Pythagorean triple

    (a,b,c)=(n^2-m^2,2mn,n^2+m^2),   t=n/m>1.

A chronological word w=w_1...w_k acts by w_k composed with ... composed with w_1. Its node gives a positive arithmetic-progression square-sum triangle precisely when

    c^2>4ab,

or equivalently

    t in A=(1,mu) union (lambda,infinity),
    mu=2-sqrt(3)+sqrt(6)-sqrt(2),
    lambda=2+sqrt(3)+sqrt(6)+sqrt(2).

This is the existing exact Pythagorean triangle classification. The integer labels are 2(c^2-4ab),2c^2,2(c^2+4ab), up to a common square scale.

The previous result showed that appending LLLL or RRRR always enters A. Indeed

    L^4(t)=1+1/(4+1/(t-1))<5/4<mu,
    R^4(t)=t+8>9>lambda.

The new result is that the set of all admissible finite chronological words is nevertheless not a regular language. Thus no finite automaton, even one retaining a more sophisticated finite state than the last few letters, recognizes exact positivity for all words.

For the proof, reverse the word convention. For an outer-first word u=u_1...u_k define

    f_u=f_(u_1) composed with ... composed with f_(u_k).

Then f_(uv)=f_u composed with f_v. The outer-first positivity language is the reversal of the chronological language. Regular languages are closed under reversal, so proving nonregularity in this convention suffices. Word reversal is essential: a prefix in the following argument means an outermost composition, not an initial segment of a root-to-node trajectory.

## 2. Finite suffix nodes are dense

Every reduced rational n/m>1 with m+n odd occurs as f_v(2) for a finite word v. Here is a direct proof that fixes the suffix direction needed later.

For 1<n/m<2, remove its outer L by replacing (m,n) with (2m-n,m). For 2<n/m<3, remove its outer M by replacing (m,n) with (n-2m,m). For n/m>3, remove its outer R by replacing (m,n) with (m,n-2m). Each new pair is positive, increasing, coprime, and has opposite parity. Its sum is strictly smaller. The boundary n/m=3 is excluded by opposite parity; n/m=2 is exactly the root (1,2). Thus the procedure terminates at the root, and reversing the removal operations gives the required composition identity.

These rationals are dense in D: odd numerators over powers of two provide a dense subfamily. Consequently every nonempty open interval in D contains f_v(2) for an actual finite suffix v. No assertion about the distribution of infinite addresses is needed.

The branch images form the disjoint interval partition

    L(D)=(1,2),    M(D)=(2,3),    R(D)=(3,infinity),

apart from the two rational boundary points. Their orientations are respectively increasing, decreasing, increasing.

## 3. Exact finite-state criterion for a single irrational wall

For any irrational alpha>1, define

    C_alpha={u: f_u(2)<alpha}.

Define its inverse-branch map G by

    G(t)=1/(2-t)       for 1<t<2,
         1/(t-2)       for 2<t<3,
         t-2           for t>3.

All iterates alpha_j=G^j(alpha) are irrational and stay in D. Let u_j be the first j symbols of the unique inverse itinerary, so

    alpha=f_(u_j)(alpha_j).

### Theorem 1: the cutoff-state criterion

C_alpha is regular if and only if the set {G^j(alpha): j>=0} is finite.

Proof of necessity. The right residual after the finite prefix u_j is

    u_j^(-1) C_alpha={v: f_(u_j)(f_v(2))<alpha}.

If f_(u_j) increases, this is {v:f_v(2)<alpha_j}; if it decreases, it is {v:f_v(2)>alpha_j}. Two distinct cutoff values with the same orientation give distinct residuals: choose a finite suffix node strictly between the cutoffs by Section 2. Opposite orientations also give distinct residuals: choose a finite suffix node between 1 and the smaller cutoff. The lower-cut residual contains that suffix and the upper-cut residual does not. Thus a finite automaton, which has only finitely many residuals, forces only finitely many cutoff values.

Proof of sufficiency. Use states (beta,o), where beta ranges over the finite inverse orbit and o records whether the residual is a lower or upper cut, plus the two constant accept/reject sink states. At state (beta,o), read a symbol s. If beta is outside f_s(D), all possible finite suffix values f_s(f_v(2)) lie on one side of beta, so transition to the appropriate sink. If beta lies inside f_s(D), transition to (f_s^(-1)(beta),o'), reversing the inequality precisely when s=M. At the end of the word, accept if the stored inequality holds at the root 2. Because beta is irrational it never equals a rational branch boundary or 2. This is a finite deterministic automaton recognizing C_alpha. QED.

This proof makes the distinction between infinite addresses and finite words explicit. Infinite inverse addresses supply the test prefixes. Every separation of states is witnessed by a finite suffix whose value lies in an ordinary rational interval.

### Corollary 1.1: an elementary obstruction beyond quadratic degree

If C_alpha is regular, alpha is a real quadratic irrational.

Indeed some alpha_i=alpha_j with i<j. The nonempty branch word between these positions is a rational Möbius map h fixing alpha_i. It is not the identity, because it maps D into the proper interval determined by its first letter. A nonidentity rational Möbius fixed point satisfies a nonzero polynomial of degree at most two. Irrationality makes the degree exactly two. A rational invertible Möbius map then carries alpha_i back to alpha without changing its field degree.

This necessity proof also applies to transcendental alpha and does not invoke a Lagrange theorem.

### Corollary 1.2: complete irrational-wall classification, with an attributed input

For every irrational alpha>1,

    C_alpha is regular if and only if alpha is quadratic over Q.

For the converse use the conjugacy

    P(t)=((t^2-1)/(t^2+1), 2t/(t^2+1)),
    t=(1+x)/y.

Direct substitution identifies G with the Romik parent map, with branches R,M,L corresponding to Romik digits 1,2,3. Cha and Kim's Lagrange theorem says that a nonrational point on the unit quarter circle defined over a real quadratic field has an eventually periodic Romik orbit. Since P and its inverse are rational, a quadratic alpha has a finite G-orbit. Apply Theorem 1.

Primary source: Byungchul Cha and Dong Han Kim, “Number theoretical properties of Romik's dynamical system,” Section 3, https://arxiv.org/html/1902.03571 . The conjugacy and finite-word language deduction are supplied here. The more general dynamical background is Dan Romik, “The dynamics of Pythagorean triples,” https://arxiv.org/abs/math/0406512 .

## 4. The positivity language is nonregular

### Theorem 2

The finite chronological Berggren words producing positive AP square-sum triangles form a nonregular language.

Proof. Work first in the outer-first convention. The wall mu lies in L(D)=(1,2), whereas lambda>3. Therefore the residual after the one-letter prefix L sees only the first positivity boundary: for every suffix v,

    L(v) is positive iff L(f_v(2))<mu
                   iff f_v(2)<G(mu).

So this residual is exactly the single-wall language C_(G(mu)). The field Q(mu) is Q(sqrt(2),sqrt(3)). To see degree four directly, write

    mu=2-sqrt(3)+sqrt(6)-sqrt(2).

Its images under the four sign choices for sqrt(2),sqrt(3) are distinct. Since 1,sqrt(2),sqrt(3),sqrt(6) are linearly independent over Q, no nontrivial sign change fixes mu. Its minimal polynomial is

    X^4-8X^3+2X^2+8X+1.

The invertible rational Möbius map G on the initial L interval preserves degree, so G(mu) also has degree four. By Corollary 1.1, the residual cannot be regular. A residual of a regular language would be regular, contradiction. Reversal transfers the conclusion to chronological words. QED.

The theorem is stronger than failure of a fixed suffix-length rule. It excludes every exact finite-state recognizer on all finite words. It does not exclude exact integer arithmetic: positivity is decidable directly by c^2>4ab, or by maintaining the two Euclid parameters.

## 4A. Complete finite-boundary criterion

The argument extends to every piecewise-constant membership predicate on D with finitely many genuine jump boundaries. A genuine jump means that membership is different on the two neighboring open intervals. Values assigned at the finitely many boundary points can be arbitrary.

First, a rational wall has a regular cut language. Its reduced rational inverse-branch algorithm strictly decreases numerator plus denominator until it reaches 2 or 3. At wall 2, the L cylinder lies below and the M,R cylinders above; at wall 3, the L,M cylinders lie below and the R cylinder above. These are terminal decision states, with the empty suffix tested separately at root 2. The preceding finite-state construction therefore extends to rational walls.

Second, branch cylinders along any irrational wall alpha eventually isolate it from every other fixed finite collection of walls. This needs a proof beyond the density slogan. Every finite node q=f_v(2) is excluded from all open branch cylinders of depth greater than |v|: inverse decoding reaches the partition boundary 2 at that stage. Given any neighborhood of alpha, choose two finite nodes p<alpha<q inside that neighborhood using Section 2. At depths greater than both node depths, the cylinder containing alpha contains neither p nor q. Since it is an interval containing alpha, it must lie between p and q. Thus the nested cylinders shrink to alpha and isolate it from any finite set of distinct other walls.

### Theorem 3: regularity of finite-jump predicates

The finite-word language of such a predicate is regular if and only if every genuine jump boundary in D is rational or quadratic irrational.

For sufficiency, express membership on the open intervals as a finite Boolean combination of the regular lower-cut languages. Altering the finitely many boundary values affects at most finitely many words: an irrational boundary is never a rational finite node, and every rational finite node has a unique word by its deterministic parent decoder. Regularity survives finite changes.

For necessity, suppose a genuine jump alpha is neither rational nor quadratic. It is irrational. Isolate it in a cylinder f_u(D) as proved above. The residual after u is a lower or upper cut at f_u^(-1)(alpha), with no endpoint equality possible at a finite node. A regular language would give a regular residual and hence, by Theorem 1 and Corollary 1.1, a quadratic pulled-back wall. Invertible rational Möbius maps preserve the field Q(alpha), so alpha itself would be quadratic, a contradiction. QED.

This includes the positivity cones as a two-jump instance. It also explains why the degree of the actual Euclid-slope boundary matters: the quadratic leg-ratio threshold is not the coordinate being read by the ternary word automaton.

## 5. What this says about entry, and what it does not

There is no tension between bounded controllable entry and unbounded memory for exact recognition. Four chosen outer steps guarantee entry into the positivity cone from any slope. Deciding whether an arbitrary supplied word has already entered demands arbitrarily many distinct cutoff states.

This is a precise model for separating a local repair mechanism from its recognition or coverage problem. It does not identify the Collatz parity language with the Berggren language. In the Berggren tree the child letter can be chosen; in a Collatz trajectory the next parity is determined by the source. Appending LLLL is a total map on every positive slope. Appending a Collatz reset word instead restricts the source to its lawful cylinder, and need not be available at a given trajectory exit.

Nor does the theorem change the existing square-sum Hamiltonian threshold theorem. It concerns the classification of the protected AP seed shapes and their Pythagorean ancestry, not the coverage of every unused interval label.

## 6. Reproduction

The companion checker uses rational Möbius matrices and the exact quartic sign on (1,2) to generate 64 boundary prefixes and distinguish every pair by finite suffix nodes. Decimal approximations are used only to propose an ordering and candidates; all accepted witnesses are then checked with exact integer arithmetic. It also checks the direct Pythagorean positivity predicate against the quartic predicate on a finite word tree and independently decodes every suffix witness to a finite Berggren word.

Finite controls audit conventions, direction of composition, and witnesses. The nonregularity proof rests on infinitely many residuals, not extrapolation from the 64-state check.

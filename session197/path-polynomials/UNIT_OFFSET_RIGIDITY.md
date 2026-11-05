# Arbitrary-module extension and unit-offset rigidity

Session197. **PROVED physical five-module formula and complete unit-offset classification within its quotient.** This extends the hybrid lane's transitive identity; it does not classify all tournament clocks.

## 1. The quotient and inherited transitive identity

Use the seven-vertex quotient

    Q=[14,32,2,22,103,13,47],

in the order `(r,P,U,a,R,b,V)`. Put nonempty tournaments T,P,U,R,V into the five indicated modules, leaving a,b singleton. Between modules, every arc follows Q. Write h_A=H(A) and b_A=B_A(2), where B_A(d) counts path covers in d labelled bins, allowing empty bins. In particular B_A(1)=h_A.

The exact transitive-module identity is proved, with uniform run bounds and a complete48-point physical grid, in

`../session197-hybrid-constructions/positive-clock-family/PROOF.md`

SHA256 `be3c701188a4a354cee62f21b6ee76bf172f67351028d88a872f4c95c9b865c3`.

In its notation x=2^|P|, y=2^|R|, z=2^|U|, w=2^|V| and T=Tr_m, it gives

    H = 2zw(x-1)(y-1) 3^m
      + [w(y-1)(4-x)+z(x-1)(4-y)] 2^m
      + (x-3)(y-3).                                       (1)

We use this named theorem as an input, rather than pretending that its finite computation alone proves an unbounded formula. Its preceding cross-cover argument supplies the all-order dimension bound.

## 2. Extension to arbitrary physical modules

**Theorem1.** For all nonempty T,P,U,R,V, the actual substituted graph has

    H = A B_T(3) + B B_T(2) + C h_T,                        (2)

where

    A = 2 b_U b_V (b_P-h_P)(b_R-h_R),

    B = b_V (b_R-h_R)(4h_P-b_P) h_U
      + b_U (b_P-h_P)(4h_R-b_R) h_V,

    C = h_U h_V (b_P-3h_P)(b_R-3h_R).                       (3)

**Proof.** Fix all but one block. The full path-fragment kernel proves that Hamilton substitution is linear on the variable block's complete bin polynomial, as in `UNIVERSAL_CONTEXTS.md`. Transitive inputs supply every monomial x^n, n>=1, spanning all polynomials with zero constant term. Every nonempty tournament has such a bin polynomial.

Thus extend the transitive identity successively in all five modules. An exponential base2 becomes evaluation B_A(2), base3 becomes B_A(3), and a constant-in-positive-multiplicity factor1 becomes evaluation B_A(1)=h_A. In particular a variable absent from a term in (1) contributes its h_A factor; it must not be silently replaced by the scalar1 for a nontransitive block. Applying this rule to every term gives exactly (2)-(3). QED.

A is strictly positive because b_A>=2h_A>h_A. Therefore, when T=Tr_m varies, the largest exponential base is3. The coefficients B and C can have either sign; entrywise positivity of the universal central-factorial matrix does not make this exponential expansion coefficientwise positive.

## 3. A sharp elementary boundary statistic

**Lemma2.** If A has n>=1 vertices, then

    B_A(2) >= (n+1) H(A).                                  (4)

Moreover

    B_A(2)-3H(A)=+1  iff A is the transitive2-vertex graph;
    B_A(2)-3H(A)=-1  iff A is K1.                           (5)

**Proof.** Each Hamilton word has n-1 internal cuts. For each cut, assign its two nonempty path pieces to the two labelled bins in either order. This produces2(n-1)H(A) incidences. A fixed ordered pair of disjoint spanning paths has at most two preimages, since only the two concatenation orders are possible and each has a uniquely determined cut. Therefore at least(n-1)H(A) ordered two-bin path covers have both bins nonempty. The two choices of a single occupied bin contribute2H(A), proving (4).

For n>=4, subtracting3H gives at least(n-2)H>=2. For n=3, a value+1 would force H=1 by (4). An H1 tournament is transitive: a nontransitive tournament contains a directed triangle, whose three Hamilton paths inject into those of the full tournament by deterministic insertion. A transitive3-vertex graph instead has B(2)-3H=8-3=5. A value-1 is already excluded by (4). Finally K1 has2-3=-1, and every2-vertex tournament is transitive and has4-3=+1. This proves (5). QED.

The cut map is only used for the inequality; it is not claimed to be injective. The factor-of-two upper bound on its fibers is essential.

## 4. Complete pure-unit classification for arbitrary modules

Fix P,U,R,V and vary T=Tr_m over positive m. Suppose the count eventually has the form alpha*3^m+epsilon, with alpha>0 and epsilon in{+1,-1}. Linear independence of the three exponential bases in (2) forces B=0 and C=epsilon, so the identity then holds at every positive m.

Because all h values are positive integers, C=+/-1 forces h_U=h_V=1. Hence U,V are transitive, of some positive orders u,v. It also forces both factors `(b_P-3h_P)` and `(b_R-3h_R)` to be units. Lemma2 says that P,R must each be K1 or Tr2.

- For C=+1, either both are K1 or both are Tr2. If both are K1, then B=2(b_U+b_V)>0, so cancellation fails. If both are Tr2, then B=0 and

      alpha = 18*2^(u+v) = 72*2^(u+v-2).

- For C=-1, one is K1 and the other Tr2. Then B is either6b_U or6b_V, strictly positive. Cancellation again fails.

**Theorem3.** Within this quotient with arbitrary nonempty physical modules, the pure positive-unit ternary rays are exactly those with P=R=Tr2 and U,V transitive. Their amplitudes are precisely72 times a nonnegative integral power of2. There are no pure negative-unit ternary rays.

In particular this entire arbitrary-module family cannot realize180*3^m+1: the ratio180/72=5/2 is not an integral power of2. Thus replacing the leaves by nontransitive tournaments does not repair the missing amplitude180 route. This is a stronger statement than a transitive-grid failure, but remains a restriction on this specified quotient. Other physical quotients, arc edits, or port-changing moves are outside it.

## 5. Verification and provenance

`verify_unit_offset_rigidity.py` constructs the seven-block graphs directly, without importing the hybrid lane's verifier. It checks (2) on all32 input choices where each of the five modules is K1 or C3, giving genuinely nontransitive controls, and on twelve positive-unit graphs. It checks the cut inequality and the complete +/-1 boundary inventory on every labelled tournament through order5 (1,099 graphs). It also checks the algebraic coefficient classification over a grid of physical module profiles.

Run normally and with `-O`; the explicit exception gates and stdout receipts agree. The finite checks support the identities and hostile boundaries. The all-order extension and rigidity are proved above using the named transitive theorem and the elementary cut map.

The transitive five-parameter identity and physical four-port return belong to the hybrid lane. The new material here is the complete arbitrary-module coefficient formula, the boundary-statistic rigidity proof, and the resulting stronger all-order classification. No independent new claim is made for the already shared positive source.

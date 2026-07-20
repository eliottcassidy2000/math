        # Message: mac-mini-2026-07-20-S140: THM-1600 -- the LAPLACE-GMC(1) LAYER IS AN IDENTITY AT DEGREE 1 (truncated exponentials; L((v-1)^m) = the DERANGEMENT NUMBERS, so 1/e was the derangement asymptotic), CHARGE SPAN 2 ELIMINATED EXACTLY in three lines, and the ten Erdos problems are a CLEAN NEGATIVE

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 13:18

        ---

        OWNER: 'prove the one sided conjecture for bounded charge span by exact elimination and settle the laplace-GMC(1) layer ... show int_0^inf CT_w(P^m) e^{-s} ds != 0 for some m', plus ten Erdos problem numbers.

(A) THE LAPLACE LAYER AT DEGREE 1 IS AN IDENTITY, NOT AN ASYMPTOTIC.
THM-1520's saddle lemma was a Laplace-method sketch without error bounds (HYP-8350). At degree 1 it is not asymptotic at all. Since C(m,k)*k! = m!/(m-k)!, the binomial expansion reindexes into a PARTIAL EXPONENTIAL SUM:
    L((av+b)^m) = m! * a^m * e_m(b/a),    e_m(x) = sum_{j<=m} x^j/j!   (TRUNCATED EXPONENTIAL)
Verified for several (a,b), m = 1..8. Three rigorous consequences:
  * e_m(x) -> e^x, so the S135 saddle limit exp(a_{D-1}/(D a_D)) is EXACT here, = e^{b/a}.
  * |e_m(x) - e^x| <= (|x|^{m+1}/(m+1)!) e^{|x|}, so e_m(x) != 0 once that tail drops below e^{Re x} -- an EXPLICIT threshold m_0(x). No saddle point, no error analysis.
  * p = v - 1 gives L(p^m) = m! e_m(-1) = !m, THE DERANGEMENT NUMBERS 0, 1, 2, 9, 44, 265. Nonvanishing for m != 1 is the classical fact that DERANGEMENTS EXIST, and !m/m! -> 1/e IS the saddle limit. The constant 1/e I reported in S135 was the derangement asymptotic all along.
HYP-8350 IS CLOSED AT DEGREE 1. Degrees >= 2 remain -- see handoff.

(B) CHARGE SPAN 2 ELIMINATED EXACTLY.
Charge support {-1, 0, +1} is the smallest case NOT covered by THM-1540's two-charge theorem -- and is precisely the trichotomy of the S137 reflection. Write P = W s(V) + q(V) + Z r(V) with V = ZW. Charge-0 needs equally many +1 and -1 factors, so
    E[P^m] = sum_k [ m!/(k!^2 (m-2k)!) ] * L( v^k (rs)^k q^{m-2k} ).
With CONSTANT r,q,s this eliminates in THREE LINES:
    1. m=1 gives E[P] = q, so q = 0.
    2. With q=0 only m=2k survives: E[P^{2k}] = (2k)! (rs)^k / k!.
    3. That vanishing forces rs = 0, i.e. r = 0 or s = 0.
So P = Z r (charges {+1}) or P = W s (charges {-1}) -- ONE-SIDED. QED
Brute-force confirmed on 729 triples (r,q,s in [-4,4], m=1..10, zero exceptions), and extended to DEGREE-1 coefficients: 15624 triples, coefficients in [-2,2], m=1..8, ZERO exceptions -- verified but NOT proved there; a Groebner elimination was attempted and did not complete.

(C) THE REQUESTED INTEGRAL FORM.
THM-1540(A) gives E[P^m] = int_0^inf CT_w(H_sqrt(s)(w)^m) e^{-s} ds. Every two-sided case tested is nonzero at some m; the one-sided cases are identically zero -- they ARE in the nullcone, as the conjecture says. The {-1,+1} rows are nonzero exactly at EVEN m, matching THM-1540(B)'s prediction B'+C' = 2. Example: (r,q,s) = (3,0,-2) gives 0, -12, 0, 432, 0, -25920.

(D) THE TEN ERDOS PROBLEMS -- A CLEAN NEGATIVE, recorded so nobody repeats the search.
All ten (#1016, 506, 742, 19, 580, 547, 460, 556, 475, 848) were retrieved and read; no dead links. NONE of them involves polynomials, power sums, moments, constant terms, Laurent polynomials, truncated exponentials, zeros of partial sums, sum-product structure, or Mathieu-Zhao / Duistermaat-van der Kallen material. The DvdK circle needs a polynomial algebra, a linear functional, and a vanishing-on-all-powers hypothesis; NOT ONE of those three appears anywhere on the list. Six are extremal graph theory, one combinatorial geometry, three elementary number theory.
#475 IS A FALSE FRIEND: it is literally phrased with 'partial sums', but they are prefix sums of a permuted finite subset of F_p -- a Hall-Paige / sequenceability problem, nothing to do with truncating a power series. If the list was assembled by string-matching, that is very likely why it landed there.
BUT THE LIST IS NOT RANDOM, and the real axis is worth having: 8 of the 10 carry Bloom's DECIDABLE badge -- 'true for all sufficiently large n, small cases open'. That is an unusual concentration (the median Erdos problem is flatly OPEN), and it is EXACTLY THE SHAPE OF THE m >> 0 QUANTIFIER THAT DEFINES A MATHIEU-ZHAO SUBSPACE, and exactly the shape of (A)'s explicit threshold. THE RESONANCE IS OF PROOF SHAPE, NOT CONTENT -- almost certainly how the list was filtered -- and it should not be cited as a mathematical connection.

HANDOFF -- three:
(i) HYP-8350 is REDUCED, not closed: degree 1 exact, degrees >= 2 open. TWO UNTRIED ROUTES. (a) Find the degree-D analogue of the truncated-exponential collapse -- at D=1 the weights k! exactly cancelled C(m,k); does the same happen against multinomial coefficients at higher D? (b) THE ELEMENTARY ONE I would try first: integration by parts gives L(f') = L(f) - f(0), hence L(m p^{m-1} p') = L(p^m) - p(0)^m. So L(p^m) = 0 for all m forces L(p^{m-1} p') = -p(0)^m/m -- a recursion, completely elementary, and UNTRIED.
(ii) HYP-8430 is the live edge of the nullcone conjecture: span-2 with NON-CONSTANT coefficients (m=1 now gives L(q) = 0 rather than q = 0, so the whole tower is needed -- try inducting on deg q, or feeding (A)'s identity into the inner L(...) calls, which are exactly of Laplace type), and spans >= 3 (untouched; they need THM-1540(C)'s Newton-polygon top-edge reduction, whose missing input is THE SAME domination estimate as HYP-8350).
(iii) HYP-8435 records the Erdos negative. Do not re-run that search.

SCOPE I want carried: (A) closes degree 1 ONLY; (B) proves the CONSTANT-coefficient case only, with degree-1 coefficients being bounded-box evidence; (C) is a restatement of (B) plus THM-1540(A) on the tested cases, not an independent result; (D) is a negative and should be cited as one. The nullcone conjecture overall remains OPEN -- now proved for one-sided support, exactly two charges, and span-2-with-constant-coefficients.

Artifacts: THM-1600; 04-computation/laplace_layer_and_span_elimination_macmini_S140.py; 05-knowledge/results/span2_elimination_quick_macmini_S140.out; HYP-8350 reduced, HYP-8430 and HYP-8435 new.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*

# Higher-order bottleneck clocks and the lower digits they conceal

19 October 2026. A scoped extension of the near-cusp scaling discussion. No claim about prime occurrence, arbitrary Collatz entry, or historical priority is made.

## Results

Fix a=1/8. Start y_0=-a and iterate y_(j+1)=y_j+f(y_j), stopping at the first N with y_N>=a.

1. For any fixed even integer r>=2 and 0<epsilon<=a^r, take f(y)=epsilon+y^r. Then

   N = K_r epsilon^(-(r-1)/r) + O_r(1),
   K_r = (2 pi/r) csc(pi/r).

   An explicit valid absolute error bound is

   2 a^(1-r)/(r-1) + (r/2) log(5/4) + 18a.

2. For 0<=lambda<=1, 0<epsilon<=a^4, and f(y)=epsilon+y^4+lambda y^6, uniformly in lambda,

   N = (pi/sqrt(2)) epsilon^(-3/4)
       - (3 pi sqrt(2)/8) lambda epsilon^(-1/4) + O(1).

   One valid uniform absolute error bound is

   1 + 3 log(5/4) + 18a + 2/(3a^3) + 2/a + 2a.

These are full-passage, fixed-threshold statements. They do not describe entry into the bottleneck from an arbitrary orbit.

For every integer base b>=2, the first result at epsilon=b^(-rn) gives an integer count within a bounded distance of K_r b^((r-1)n). In the quartic family, epsilon=b^(-4n) gives

   N_lambda = (pi/sqrt(2)) b^(3n)
              - (3 pi sqrt(2)/8) lambda b^n + O(1).

Thus a common leading scaled real can coexist with an unbounded difference of order b^n in the integers being generated. Keeping the next term restores bounded absolute uncertainty. This is an explicit natural polynomial example of the precision obligation in a scaling argument.

## 1. An elementary corrected-clock lemma

Let f>0 be C² on J=[-a,5a/4]. Suppose every pre-exit point y lies in [-a,a), the step stays in J, and on each such step:

- f(y)<=1/32;
- sup_J |f'|<=5/16 and sup_J |f''|<=2;
- f'(y)²<=4f(y).

Set

   Psi(y)=integral from 0 to y of du/f(u) + (1/2)log f(y).

Then each step has

   |Psi(y+f(y))-Psi(y)-1| <= 8 f(y).

Proof. Write F=f(y), A=f'(y), q(t)=f(y+tF)/F for 0<=t<=1. Taylor's theorem gives q(t)=1+tA+R(t), |R(t)|<=F t², and the derivative bound gives |q(t)-1|<=5/16. The integral part of the clock increment is integral_0^1 dt/q(t).

Using 1/(1+s)=1-s+s²/(1+s),

   |integral_0^1 dt/q(t) - (1-A/2)|
   <= F/3 + [(2/3)A²+(2/5)F²]/(11/16)
   <= [1/3 +(8/3+1/80)/(11/16)] F.

Here (q(t)-1)² <= 2t² A²+2F²t^4. Also

   |log(1+s)-s| <= s²/[2(1-|s|)]

for |s|<1, by integration of -s/(1+s). Hence

   |(1/2)log q(1)-A/2|
   <= F/2 + [A²/2+F²/2]/(11/16)
   <= [1/2 +(2+1/64)/(11/16)]F.

The two constants sum to 5057/660<8. Summing the defects telescopes the actual displacement:

   sum_(j<N) f(y_j)=y_N+a<=9a/4.

Therefore

   |Psi(y_N)-Psi(-a)-N|<=18a.                 (1)

The stopping time is finite because positive continuous f has a positive minimum on the compact interval J (in the two applications, every increment is at least epsilon). The proof is uniform as epsilon decreases; it never treats the discrete iteration as exactly a differential equation.

## 2. Pure even powers

For f=epsilon+y^r, pre-exit F<=2a^r<=1/32 and y_N<=a+2a^r<=5a/4. On J,

   |f'|<=r(5/32)^(r-1)<=5/16,
   |f''|<=r(r-1)(5/32)^(r-2)<=2.

Both bounds hold at r=2 and decrease along even r>=2. At a pre-exit point,

   f'(y)² <= r² a^(r-2) y^r <=4 f(y).

Thus (1) applies. Since y_N>=a, the logarithmic increment is between zero and (r/2)log(5/4). Extending the passage integral from [-a,y_N] to the entire real line loses at most 2a^(1-r)/(r-1), because f(y)>=|y|^r outside the passage. Finally the Euler beta integral gives

   integral_R dy/(epsilon+y^r)
   = epsilon^(-(r-1)/r) integral_R dt/(1+t^r)
   = K_r epsilon^(-(r-1)/r).

This proves the first statement and its explicit bound. The special case r=2 recovers the pi coefficient. The logarithmic correction in the clock is necessary for a uniformly summable local defect, although the endpoint logarithms in a symmetric full passage remain bounded.

## 3. A quartic family with an unbounded next term

For f=epsilon+y^4+lambda y^6 with the stated ranges, the same clock hypotheses hold. For example, on J,

   |f'|<=4(5/32)^3+6(5/32)^5<5/16,
   |f''|<=12(5/32)^2+30(5/32)^4<2,

and at a pre-exit point

   f'(y)²<=(4+6a²)² a² y^4 <=4 f(y).

The endpoint overshoot satisfies y_N-a<=f(a), so its additional integral is at most one. The endpoint logarithm is at most 3 log(5/4), since f(y_N)/f(-a)<=(5/4)^6. It follows from (1) that

   N = integral_(-a)^a dy/(epsilon+y^4+lambda y^6) + E,
   |E|<=1+3log(5/4)+18a.                    (2)

Put A=epsilon+y^4. The exact identity

   1/(A+lambda y^6)
   =1/A -lambda y^6/A²
     +lambda² y^12/[A²(A+lambda y^6)]

has nonnegative last term bounded by lambda², since A>=y^4. Its integral is at most 2a. The two main full-line integrals are

   integral_R dy/(epsilon+y^4) = (pi/sqrt(2)) epsilon^(-3/4),
   integral_R y^6 dy/(epsilon+y^4)²
       =(3 pi sqrt(2)/8) epsilon^(-1/4).

The second constant is (1/2)B(7/4,1/4), using Gamma(7/4)=(3/4)Gamma(3/4) and Gamma(1/4)Gamma(3/4)=pi sqrt(2). Truncating these full-line integrals costs at most 2/(3a³) and 2/a respectively. Combining these bounds with (2) proves the second result, uniformly for lambda in [0,1].

## 4. What this adds to a large-integer argument

The real parameter remains in a small interval, but its accuracy and the terms retained in its clock matter. For the quartic family, omitting the next term changes the count by order b^n at scale b^(3n). In particular a relative asymptotic with error tending to zero does not determine the count's residue modulo two, much less its prime correlations. This is a demonstrated loss of information in the approximation, not a claim that a real encoding cannot preserve arithmetic.

The constructive approach is to retain enough terms to obtain a certified bounded interval of integer candidates, then keep exact threshold or residue information to select the integer. The quadratic escape-threshold work in the companion note provides one such exact selection mechanism. The general question is which clock coefficients, entry/exit phases, and integer guards a proposed transfer retains.

## Attribution and scope

The near-parabolic escape phenomenon and its O(1) refinement are established in the quadratic literature: Aaron Klebanoff, “pi in the Mandelbrot Set” (2001), and Brockmoeller–Scherz–Srkalovic, “Pi in the Mandelbrot set everywhere” (2025), arXiv:2505.07138, discuss the quadratic setting. Higher-order continuous-time bottleneck scaling is also classical; see Christian Kuehn, “Scaling of Saddle-Node Bifurcations: Degeneracies and Rapid Quantitative Changes” (2009), arXiv:0807.1546, Proposition 2.1 (the continuous-time exponent is -(r-1)/r). The present elementary estimates are given in full so the stated higher-order polynomial families and uniform error bounds can be used without assuming an unproved numerical pattern. No historical priority for the bottleneck method or these asymptotics is asserted.

- https://www.doc.ic.ac.uk/~jb/teaching/jmc/pi-in-mandelbrot.pdf
- https://arxiv.org/abs/2505.07138
- https://arxiv.org/abs/0807.1546
- Euler beta integral: https://dlmf.nist.gov/5.12

The accompanying checks are small controls of the algebra, bounds, and representative iteration counts. The proof above establishes the quantified statements for all allowed parameters.

# Independent audit of the anti-fixed R5 circulation series

Session 23b, 2026-10-21. This addendum independently proves the rational formula in the anti-fixed pentagon cone. It does not modify or replace the frozen primary packet, does not assume free primitive-cycle coordinates, and does not fit a numerator to finite coefficients.

## 1. Exact object being counted

Let R5 have vertices Z/5 and arcs i -> i+1 and i -> i+2. Write their nonnegative integer flow values as x_i and y_i. Balance at vertex i is

    x_i + y_i = x_(i-1) + y_(i-2).

Consequently x_i+y_i+y_(i-1) is independent of i. Calling this common integer t gives the integral bijection

    x_i = t-y_i-y_(i-1),
    y_i >= 0,  t >= y_i+y_(i-1).

The weight, total arc multiplicity, is 5t-sum_i y_i.

Reflection s(i)=-i is an anti-isomorphism. On flow coordinates it acts by reversing each arc after applying s. Thus y_i is sent to y_(-i-2). Its fixed vectors have

    y = (a,b,b,a,c),
    x = (t-a-c, t-a-b, t-2b, t-a-b, t-a-c).

The anti-fixed flow series is therefore exactly

    k(q) = sum_(a,b,c>=0, t>=max(a+b,2b,a+c))
                              q^(5t-2a-2b-c).                  (1)

For completeness, Aut(R5) is its cyclic rotation group: among the two out-neighbors of each vertex, the first successor is distinguished as the source of the arc between those neighbors. An automorphism is therefore determined by the image of a vertex and must be a rotation. The five anti-isomorphisms are rotation-conjugate because 2 is invertible modulo 5. The anti-coset Burnside average is consequently (1) itself. Thus k also counts ordinary flow-decoration orbits that are fixed by converse; it does not count all flow-decoration orbits over this self-converse underlying tournament.

## 2. Disjoint lattice decomposition

Put M=max(a+b,2b,a+c) and write t=M+s, where s>=0. The variable s contributes 1/(1-q^5) in each case below. All other parameters introduced below are nonnegative integers, except those explicitly declared positive. These cases partition the entire lattice domain, including all boundaries exactly once.

### I. a>=b and c<=b

Write b=c+v and a=b+u, with c,u,v>=0. Then M=a+b and

    5M-2a-2b-c = 5c+3u+6v.

The contribution is

    K_I = 1 / ((1-q^3)(1-q^5)^2(1-q^6)).          (2)

### II. a>=b and c>b

Write a=b+u and c=b+v, with b,u>=0 and v>=1. Then M=a+c and

    5M-2a-2b-c = 5b+3u+4v.

The contribution is

    K_II = q^4 / ((1-q^3)(1-q^4)(1-q^5)^2).       (3)

### III. a<b and c<=a

Write b=a+u and a=c+v, with c,v>=0 and u>=1. Then M=2b and

    5M-2a-2b-c = 5c+6v+8u.

The contribution is

    K_III = q^8 / ((1-q^5)^2(1-q^6)(1-q^8)).      (4)

### IV. a<b and a<c<=2b-a

Write b=a+u and c=a+v, with a>=0, u>=1 and 1<=v<=2u. Then M=2b and

    5M-2a-2b-c = 5a+8u-v.

There is a disjoint parity split:

- v=2w, u=w+h, with w>=1 and h>=0, giving 5a+6w+8h;
- v=2w-1, u=w+h, with w>=1 and h>=0, giving 5a+6w+8h+1.

Thus

    K_IV = (q^6+q^7) / ((1-q^5)^2(1-q^6)(1-q^8)). (5)

### V. a<b and c>2b-a

Write b=a+u and c=a+2u+v, with a>=0 and u,v>=1. Then M=a+c and

    5M-2a-2b-c = 5a+6u+4v.

The contribution is

    K_V = q^10 / ((1-q^4)(1-q^5)^2(1-q^6)).       (6)

Each parametrization is bijective on its stated region. Therefore the sum of (2)-(6) proves the generating function without any enumerative extrapolation.

## 3. Closed formula and exact algebra certificate

Define

    P(q)=1+q^4+q^5+q^6+q^7+q^8+q^12.

The result is

    k(q) = P(q) / ((1-q^3)(1-q^5)(1-q^6)(1-q^8)). (7)

One exact polynomial certificate for the simplification is

    (1-q^4)(1-q^8)
      + q^4(1-q^6)(1-q^8)
      + (q^6+q^7+q^8)(1-q^3)(1-q^4)
      + q^10(1-q^3)(1-q^8)
    = (1-q^4)(1-q^5) P(q).                        (8)

Indeed the common denominator of (2)-(6) is

    (1-q^3)(1-q^4)(1-q^5)^2(1-q^6)(1-q^8),

and (8) cancels exactly its extra factors (1-q^4)(1-q^5).

The denominator in (7) is fully reduced over Q[q]. Its cyclotomic factorization is

    (1-q)^4 Phi_2(q)^2 Phi_3(q)^2
          Phi_4(q) Phi_5(q) Phi_6(q) Phi_8(q).

The remainders of P modulo these irreducible cyclotomic polynomials are:

    index       1   2   3   4      5       6          8
    remainder   7   3   1   3    1+q^2     3    -q-q^2-q^3.

Every remainder is nonzero. In particular, cancellation does not reduce the genuine pole at any of these primitive root classes.

## 4. Recurrence, quasiperiod and growth

Let k_n=[q^n]k(q), let k_n=0 for n<0, and let p_n=1 for

    n in {0,4,5,6,7,8,12},

and p_n=0 otherwise. Expanding the reduced denominator gives

    D(q)=1-q^3-q^5-q^6+q^9+2q^11+q^13
                         -q^16-q^17-q^19+q^22.

Hence the exact recurrence is

    k_n = k_(n-3)+k_(n-5)+k_(n-6)
          -k_(n-9)-2k_(n-11)-k_(n-13)
          +k_(n-16)+k_(n-17)+k_(n-19)-k_(n-22)+p_n. (9)

It is homogeneous for n>=13, interpreting negative indices as zero; in the usual order-22 presentation it is homogeneous for every n>=22. The first 22 values are

    1,0,0,1,1,2,3,2,4,4,5,7,9,9,11,13,15,17,21,21,26,29.

Coprimality in Section 3 also proves that 22 is the minimal order of an eventual homogeneous constant-coefficient recurrence over Q.

The pole at 1 has order 4; the poles at -1 and at primitive cube roots have order 2; all remaining poles are simple. Partial fractions therefore show that k_n is an exact quasipolynomial in n>=0 of degree 3 and period dividing 120. Its minimal quasiperiod is 120: a quasipolynomial of period L can have poles only at L-th roots of unity, while (7) has genuine primitive third, fifth and eighth root poles, forcing 3,5,8 to divide L.

For the leading terms, write

    k(q)=H(q)/(1-q)^4,
    H(q)=P(q)/product_(d in {3,5,6,8})(1+q+...+q^(d-1)).

The exact values are

    H(1)=7/720,
    H'(1)=(42-7*(2+4+5+7)/2)/720=-7/240.

Thus the two highest singular terms at q=1 are

    7/(720(1-q)^4) + 7/(240(1-q)^3).

The remaining contribution is O(n), including the double poles away from 1. Consequently

    k_n = 7*n^3/4320 + 7*n^2/288 + O(n),           (10)

and in particular k_n ~ 7*n^3/4320. This is not merely a subsequence asymptotic.

## 5. Reciprocity and an intrinsic check on the shift

The numerator is palindromic of degree 12 and the denominator has degree 22 with four factors of the form 1-q^d. Therefore

    k(q^(-1)) = q^10 k(q).                        (11)

The shift 10 has a direct flow meaning. The relative interior of the anti-fixed cone consists of flows positive on all ten arcs. Subtracting one from every arc preserves balance, anti-invariance and integrality, and gives a bijection from its positive integral flows to all its nonnegative integral flows. This works because R5 has indegree and outdegree 2 at every vertex and the constant-one flow is anti-fixed. Thus the interior lattice-point series is q^10 k(q), in agreement with the reciprocal identity. This observation is a consistency check, not an assumption in the proof of (7).

## 6. Exact independent controls

The standalone standard-library checker is `check.py`. Every invariant uses an explicit failure check, not Python assert; normal and optimized runs have byte-identical output, saved as `CHECK.json` and `CHECK-OPTIMIZED.json`.

It verifies:

1. The exact polynomial identity (8), the expanded denominator, all cyclotomic remainders, and numerator/denominator coprimality by exact Euclidean division.
2. Direct enumeration of all 3,528,687 lattice tuples in (1) of total weight at most 300. Every coefficient agrees with (7), with (9), and with the sum of the five separately checked case series. The five cases contain respectively 854,691; 1,201,405; 298,266; 608,375; and 565,950 such tuples.
3. A structurally independent actual-arc calculation through weight 48. It constructs the six arc orbits of the reflection on the actual ten-arc R5, enumerates 217,035 feasible-weight assignments to its first five orbit variables, determines the final singleton's value from one vertex balance, and checks all five vertex balances, the literal arc total, and anti-invariance. All 3,324 resulting anti-fixed flows agree coefficientwise with (7).
4. The two leading asymptotic constants using exact rational arithmetic.

The cone cutoff is proved complete as follows. Adding the five inequalities y_i+y_(i-1)<=t gives 2 sum_i y_i<=5t; hence weight=5t-sum_i y_i>=5t/2. At total weight <=N one has t<=floor(2N/5), a,c<=t, and b<=t/2. These are the exact finite loop bounds used in the direct cone enumeration.

Selected exact coefficients are k_120=3170, k_240=23839, and k_300=45986. The finite controls check the proof; they are not its source.

## 7. Scope

This addendum is independent of the producer's geometric-sum derivation and its raw R5 triangulation. It proves only the anti-fixed series (7) and its stated consequences. The raw circulation series and ordinary flow-orbit series require their separate certificates. No formula here treats primitive directed cycle rays as algebraically independent semigroup generators.

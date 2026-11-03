# Right triangles, odd Collatz strands, and Fibonacci hypotenuses

Let T(x)=(3x+1)/2^{v_2(3x+1)} for positive odd x. Every result below is elementary and independent of Collatz convergence.

## 1. Odd strands and rigid gap-two triangles

For odd y, every odd T-preimage has the form x=(2^q y-1)/3, with q>=1. Such a preimage exists exactly when 2^q y=1 mod 3. Thus y=3 mod 6 has no odd preimages; y=5 mod 6 has exactly the odd q, and y=1 mod 6 has exactly the even q. These are inverse-exponent strands, not a claim of orbit convergence.

A positive right triangle cannot have three odd sides, as odd squares are 1 mod 4. In fact, an odd leg a and odd hypotenuse c cannot differ by 4: b^2=8(a+2) implies (b/2)^2=2(a+2)=2 mod 4. They differ by 2 if and only if (a,b,c)=(k^2-1,2k,k^2+1), k positive even. Indeed b^2=4(a+1) and b/2=k, and parity of a forces k even. These triples are primitive. For 3 not dividing k the adjacent residues (a,c) are (3,5) mod 6, while for 3|k they are (5,1) mod 6. Consequently the residues (1,3) in a 6n+1,6n+3,6n+5 block are never a gap-two odd-leg/hypotenuse pair, but (3,5) occur infinitely often. The full ordered progression (3h,3h+2,3h+4) satisfies the Pythagorean equation only at h=2, giving (6,8,10).

In this gap-two family, T(a)=c occurs only when k=2: 3(k^2-1)+1=2^q(k^2+1) forces q=1, then k=2. T(c)=a never occurs, since q>=2 and 3k^2+4=2^q(k^2-1): at q=2 it would require k^2=8; for q>=3 it would require k^2<=12/5, contrary to k>=2. Thus 3->5 is the unique accelerated Collatz edge in either direction between such a pair.

## 2. Complete forward Collatz classification

Consider any positive integer right triangle with odd leg a, even leg b and hypotenuse c and T(a)=c. Since c>a and a>1, 3a+1=2^q c forces q=1. Then a=4t-1, c=6t-1 for a positive integer t; write b=2y. The Pythagorean condition is exactly

(10t-1)^2 - 20y^2 = 1.

Every positive solution of X^2-20y^2=1 is X+y sqrt(20)=(9+2 sqrt(20))^r, r>=1. For an elementary descent, (9-2 sqrt(20))(X+y sqrt(20)) has nonnegative integer coefficients, with smaller y, since X/y<=9/2 for y>=2, and the only positive solution with y<2 is absent. Since X=9^r mod 10, t=(X+1)/10 is integral exactly for odd r. Thus the complete list is

X+y sqrt(20)=(9+2 sqrt(20))^(2j+1), t=(X+1)/10,
(a,b,c)=(4t-1,2y,6t-1), j>=0.

Conversely all these are positive right triangles with T(a)=c and q=1. They are primitive because gcd(a,c) divides 3a-2c=-1, and a=3, c=5 mod 6: X is 0 mod 3 for odd r, so t=1 mod 3. Initial examples are (3,4,5) and (1155,1292,1733).

## 3. Infinite reverse family

Let X+Y sqrt(7)=(381+144 sqrt(7))(127+48 sqrt(7))^j, j>=0, and let a=(X+4)/7, b=Y/3, c=(4a-1)/3. The norm is X^2-7Y^2=9. Modulo 21, X=381=3 mod 21, so a=1 mod 3; Y is divisible by 6, and X is odd, so a,c are odd and b even. A direct expansion yields a^2+b^2=c^2 and 3c+1=4a, hence T(c)=a. All are primitive since any common divisor of a,c divides 3c-4a=-1. Initial examples are (55,48,73), (13825,12192,18433). This is a subfamily, not a classification of all reverse edges: for example T(13)=5 has exponent 3.

## 4. One-digit Zeckendorf hypotenuses

With F_0=0,F_1=1, n>=2 and n !=1 mod 3, the coprime numbers F_n,F_(n+1) have opposite parity, so Euclid produces the primitive triple

(F_(n+1)^2-F_n^2, 2F_n F_(n+1), F_(2n+1)).

The last identity is F_n^2+F_(n+1)^2=F_(2n+1); the hypotenuse therefore has a single Fibonacci summand in its canonical Zeckendorf representation. The period-8 Fibonacci residues mod 3 show that this odd hypotenuse is 1 mod 6 for n=0 or 3 mod 4 and 5 mod 6 for n=1 or 2 mod 4. The inverse exponent q is therefore even in the first case and odd in the second. For n=3, the triangle is (5,12,13), and T(13)=5, giving an example simultaneously involving the Fibonacci hypotenuse and a reverse Collatz edge. No assertion is made that Zeckendorf digit operations commute with T.

## Prior context

The inverse-exponent strands and their return analysis appear in `arithmetic_braids_20260917_collatz.md` and `ternary_berggren_20260925.md`. THM-3339 already studies consecutive-Fibonacci Euclid parameters and Berggren rays; the Fibonacci hypotenuse identity is classical. The Pell intersection and residue comparison here are presented as a synthesis, without a literature-priority claim.

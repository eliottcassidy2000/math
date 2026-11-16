# From a decimal coincidence to a cubic orbit

Start with 16=10+6 and 1024=1000+24. The three different exact extensions below retain different features. None is a proof of Collatz or Chan–Pak.

## 1. Keep the Collatz arrow: a complete radix classification

On odd source n, the standard Collatz arrow is n→3n+1. Encode each vertex n by 2^n. The synchronized radix equation
`2^(3n+1)=B^n+3·2^n` forces B even and (B/2)^n=2·4^n−3.
This is above 4^n for n≥1, so B/2≥5. For n≥4, (5/4)^n>2, an impossibility. The three signals n=1,2,3 give 5,29,125 on the right, so n=1,3, with B=10, are the only solutions even when n is allowed to be even. The middle value 29 is not a square. With the sign reversed, `2^(3n−1)=B^n−3·2^n` forces (B/2)^n=4^n/2+3. For n≥2 this lies strictly between 3^n and 4^n (check n=2, then (4/3)^n>2 for n≥3). The only solution is n=1, B=10: 4=10−6, the 1→2 arrow of the positive 3n−1 map (equivalently signed Collatz at −1).

If instead we keep B=10 and let the normalized signed correction vary, the actual exponent-graph embedding gives K(n)=2·4^n−5^n on the odd branch. Its only positive values over all n≥1 are 3,7,3 (at 1,2,3); the 7 is an even-source ghost. It is strictly decreasing from n=2 onward, since K(n+1)−K(n)=6·4^n−4·5^n<0 for n≥2. More generally a collision Au^n−v^n=Au^(n+d)−v^(n+d) with 1<u<v satisfies (v/u)^(n+d)<A, by dividing the two geometric differences. Thus echoes between two exponentiated scales occur strictly before the sign crossover.

## 2. Keep the gap three: classify the entire enlarged family

Write α=∛2 and define F(x,y,z)=x³+2y³+4z³−6xyz, the field norm of x+yα+zα². Then
`1024−1000=8(128−125)`, and `F(−5,4,0)=3`.
Let u=α−1. The inverse u⁻¹=1+α+α² is integral, and u has norm one. Multiplication by u acts on coefficients as
`(x,y,z) → (2z−x,x−y,y−z)`. It gives
```
(1,1,0) → (−1,0,1) → (3,−1,−1) → (−5,4,0).
```
The opposite direction yields (3,2,2), (11,9,7), (43,34,27), … . That direction is eventually strictly positive with no zero-coordinate return: its matrix is ((1,2,2),(1,1,2),(1,1,1)).

**Completeness.** If β=x+yα+zα² has norm 3, its real embedding t is positive and its two conjugates have common magnitude √(3/t). Multiply β by an integral power of u to make 1≤t<u⁻¹<4. Inverting the three embeddings by the 3-point Fourier transform puts each of x,yα,zα² in
`[(t−2√(3/t))/3,(t+2√(3/t))/3]`.
The lower endpoint exceeds −1 since the numerator increases and at 1 is 1−2√3>−3. The upper endpoint is below 2 by convexity and checking the endpoints 1 and 4. Therefore x,y,z∈{0,1}, and the only such triple of norm 3 is (1,1,0). All integer solutions are the single unit orbit of 1+α (norm −3 is its negative). In particular the decimal identity is the third forward step of a complete recurrence, with a vanishing last coordinate and two remaining pure cubes. The modulo-M orbit, in contrast, is always periodic by determinant one. This gives a small exact model of the distinction between returning a physical register to zero and returning its residue to zero.

## 3. Keep the square–cube geometry: a movable chord

On y²=x³+24 the line y=3x+2 passes through the three integer points (−2,−4),(1,5),(10,32). Indeed x³+24−(3x+2)²=(x+2)(x−1)(x−10). Thus the elliptic group law gives (10,32)=(−2,4)−(1,5), and forces 32²=10³+24 from two tiny points. Do not assume elliptic addition preserves integral coordinates.

A nontrivial deformation keeps the first two abscissae and varies the chord slope a≠0:
`b=(a²+3)/(2a), c=(9a⁴+14a²+9)/(4a²)`.
Then x³+c−(ax+b)²=(x+2)(x−1)(x−a²−1). Positive integer a can have integer anchor ordinates only for a=1 or 3, since a|3. Both work: a=1 gives the small sibling 16=8+8, and a=3 gives the decimal identity. For every positive rational a, clearing denominators by weighted scaling gives integral examples; this deformation changes the ratio of the abscissae, not just the scale of the original triangle of points.

### Exact zero-coordinate returns, not just returns modulo a power

The cubic orbit also gives a second, radix-free classification of x³+2y³=3: its only integer solutions are (1,1) and (−5,4). Indeed negative k have three positive coordinates by the inverse matrix. Write u³=1+3v with v=α−α². For k≡1,2 mod 3 the last coordinate is respectively 1,−1 modulo 3. For k=3t, put a_j=[α²](1+α)v^j; the first five are 0,0,3,0,−18. The last coordinate is Σ 3^j a_j binom(t,j), vanishing for t=0,1. For t≥2 divide by t(t−1). The j=2 contribution has 3-adic valuation 3 (27/2); all j≥4 contributions have valuation at least 4, since the remaining falling factorial is integral, j−v₃(j!)≥4 for j≥5, and a₄=−18. This proves the exact valuation 3+v₃(t(t−1)) and excludes every further zero. For the original exponential equation the cubic slice classifies the case n≡0 mod3, k≡1 mod3; the full exponential theorem uses the previous modular argument. (Do not infer the full classification from the slice.)

### The swapped cubic: same recursion, different coordinate

For norm 5, normalize in the same interval [1,u⁻¹). The Fourier bounds give x∈{−1,0,1,2}, y,z∈{0,1}; in particular the lower endpoint is greater than −α and upper less than 2α. Of the 16 candidates only (1,0,1) has norm 5. So all norm-5 triples are generated by (1+α²)u^k. For k<0 all coordinates are positive. Modulo 3 the z residues at k=0,1,2 are 1,−1,2; y vanishes only at k=0,2 mod 3. In those two classes β=(1+α²)u^r has y=0 and [α]βv is respectively −1,−7. If k=3t+r with t>0, the leading term after dividing the binomial expansion by t has valuation 1, all higher terms at least 2, so y cannot vanish. Thus x³+4z³=5 has precisely (x,z)=(1,1),(−3,2), the second encoding 32=27+5. The first and second normalized cubic coincidences really are two coordinate slices of the same recurrence.

For norm 3, the complementary slice y=0 can occur only for k≡1 mod 3. At k=1 the triple is (−1,0,1). For k=3t+1 the first three positive-degree y coefficients in βv^j are −3,0,18. After division by t the leading term has valuation 2 and higher terms valuation at least 3 (j=3 by its coefficient); hence no new return. In particular for arbitrary positive b,m the full cubic-radix classifications are:
* 2^m=b³+3: (b,m)=(1,2),(5,7).
* 2^m=b³+5: (b,m)=(3,5).

To see this last implication, split m modulo 3. For m≡0 two distinct positive cubes differ by at least 7. For m≡1 use the z=0 slice and for m≡2 the y=0 slice of the corresponding norm.

## 4. Keep the same coefficient at two source sizes: two distinct growth laws

For the more general synchronized family A^(m n+d)=B^n+mA^n, with integer A,B≥2, m≥1, d integer, calibrate n=1 to be a hit. Then b=B/A=A^d R−m, R=A^(m−1), and f(n)=A^d R^n−b^n. The only nondegenerate repeated hit f(n)=m is A=2,B=10,m=3,d=1, at n=1 and 3. The degenerate case A=2,B=2,m=1,d=1 has every hit. Indeed m≤R. For d≤0, the case m=1 is impossible and otherwise b<R, so f(n+1)=Rf(n)+(R−b)b^n increases from m. For d≥2, b/R≥A^d−1>√(A^d), forcing f(n)<0 for all n≥2. This also works for d=1,A≥3. When A=2,d=1, m≥4 gives b/R≥3/2>√2; m=2 gives f(n)=2^n and m=3 is classified above.

Change instead the repeat positions to n=1,2 and drop the binary-power restriction on R. For K(n)=2R^n−b^n the repeated value forces b(b−1)=2R(R−1), equivalent to the negative Pell equation (2b−1)²−2(2R−1)²=−1. Multiplying 1+√2 repeatedly by 3+2√2 gives (R,b,K(1))=(3,4,2),(15,21,9),(85,120,50),… . Here b>R, and 2R^t−b^t has at most one interior turning point since its derivative divided by R^t is a constant minus an increasing exponential. There is no third hit. The decimal repeat at n=1,3 instead equates two *tetrahedral* numbers in a 2:1 ratio. This gives infinite Pell analogues in one deformation, while the Collatz synchronization is rigid in another.

# Exact labelled square-residue magic counts over 2-power and 3-power rings

Status: PROVED elementary formulas, independently verified finite controls. Counts concern ordered, labelled 3-by-3 arrays of nine pairwise distinct **square residues**, zero allowed. No quotient by the eight dihedral symmetries, no count of chosen roots, no integer-square witness, and no external priority claim is intended.

Write N_p(k) for this count over Z/p^k Z, and set N_p(0)=0. Then

- Sum_{k>=0} N_2(k) z^k = 1536 z^7 / [(1-z^2)(1-2z)(1-4z)(1-8z)].
- Sum_{k>=0} N_3(k) z^k = 216 z^3 (1+27z) / [(1-z^2)(1-3z)(1-9z)(1-27z)].
- N_2(k)/2^(3k) tends to 1/504.
- N_3(k)/3^(3k) tends to 27/728.

In particular N_2(7)=1536 and N_3(3)=216. Both generating functions have nonnegative coefficients and positive first coefficients, recovering the exact existence thresholds k>=7 and k>=3 respectively.

## 1. Parametrization over every commutative ring

Every magic array, with common line sum s and center c, has the form

M(c,a,b) = (a, 3c-a-b, b; c-a+b, c, c+a-b; 2c-b, a+b-c, 2c-a).

This does not divide by 2 or 3. Indeed, the two diagonals give bottom-right=s-c-a and bottom-left=s-c-b. The outer columns then give middle-left=c+b-a and middle-right=c+a-b. The middle row gives s=3c. The remaining entries follow from the first and last rows. Conversely direct substitution verifies all eight lines. Thus this is a bijection with triples (c,a,b), including in characteristics 2 and 3.

Put a=c+u, b=c+v. The nine offsets from c are 0, ±u, ±v, ±(u+v), ±(u-v). Consequently the number A_p(h) of labelled ordinary magic arrays with distinct entries over Z/p^h Z equals p^h times the number of offset pairs producing distinct offsets.

## 2. Ordinary counts for powers of 3

Let m=3^h with h>=1. Since 2 is invertible, a collision among the nine offsets occurs exactly on one of the eight lines

u=0; v=0; v=±u; v=±2u; v=±u/2.

This follows by taking all differences among the displayed offsets: after multiplying by units every difference is one of u, v, u±v, 2u±v, u±2v. Each line contains m pairs.

Group the lines by their reduction modulo 3: the vertical and horizontal lines are singleton groups, and the remaining two groups each contain three slopes. They are {1,-2,-1/2} and {-1,2,1/2}. Within either triple every difference of distinct slopes has 3-adic valuation exactly 1. Two lines therefore meet in exactly three pairs, and the three pairwise intersections coincide with their common triple intersection: u is a multiple of 3^(h-1), with v the common slope times u. This includes h=1, where all three lines in a group coincide. The union of a triple consequently has 3m-6 elements.

Lines from different groups meet only at (0,0), since their slope difference is a unit (and the vertical case is immediate). The forbidden union has

m + m + (3m-6) + (3m-6) - 3 = 8m-15

pairs. Hence

A_3(h)=m(m^2-8m+15)=m(m-3)(m-5), h>=1,

with A_3(0)=0 separately. The polynomial at m=1 must NOT be used.

## 3. Ordinary counts for powers of 2

Let m=2^h. Call a pair primitive if at least one coordinate is odd. For h>=3 and u odd, normalize by u and put t=v/u. Comparing the nine offsets shows that collision occurs exactly for

t in {0, 1, -1, 2, -2, m/2, m/2+1, m/2-1}.

These are eight distinct residues for m>=8. To check completeness without dividing by 2, the difference forms have coefficients among u, v, u±v, 2u, 2v, 2(u±v), 2u±v, u±2v. With u odd, the last two forms u±2v cannot vanish; solving the others gives precisely the displayed list. There are (m/2)(m-8) valid pairs with u odd.

For u even and v odd, use symmetry and normalize by v. Only the even forbidden ratios 0, 2, -2, m/2 occur, giving (m/2)(m/2-4) valid pairs. Thus the valid primitive count is

P(m)=3m(m-8)/4, m>=8.

It is zero for m=1,2,4 because nine distinct offsets are impossible. For a nonzero pair, divide both coordinates by their common exact power of 2. Equality of offsets before division is equivalent to equality after division in the reduced modulus. Valuation strata are disjoint and the division map is bijective on each stratum. Therefore, for h>=2, the total valid-pair count is

Sum_{j=3}^h (3/4)(4^j-8*2^j) = 4^h-12*2^h+32.

The empty sum at h=2 also equals zero; h=3 yields zero as well. It follows that

A_2(h)=m(m-4)(m-8), h>=2,

and A_2(0)=A_2(1)=0 separately. These exceptions are essential: the polynomial incorrectly gives positive values at m=1,2.

## 4. Unit/nonunit partition and descent

Opposite entries sum to twice the center. Square residues modulo 3 are 0 and 1. If the center is 0 modulo 3, each opposite pair is (0,0); if it is 1, each opposite pair is (1,1). Thus every square magic grid over a power of 3 has either all entries divisible by 9, or all entries congruent to 1 modulo 3.

Likewise square residues modulo 4 are 0 and 1. Every square magic grid over a power of 2 with k>=2 has either all entries divisible by 4, or all entries odd and hence congruent to 1 modulo 8 (when k>=3).

All residues congruent to 1 modulo 3 are squares modulo 3^k, and all residues congruent to 1 modulo 8 are squares modulo 2^k. This can be proved by elementary root lifting: for odd p=3, the derivative 2r is invertible; at powers of 2, changing an odd root r by 2^(j-1) changes its square by 2^j modulo 2^(j+1), for j>=3. The modulus-8 base case starts the induction.

The maps x -> 1+3x and x -> 1+8x give bijections from ordinary distinct magic arrays over moduli 3^(k-1) and 2^(k-3), respectively, onto the all-unit square grids. Each line has exactly three entries, so equality of line sums is equivalent under these maps. Consequently

U_3(k)=A_3(k-1), k>=1;
U_2(k)=A_2(k-3), k>=3,

with smaller U_2 values zero. The nonunit grids are bijective, by division by p^2, with square grids over modulus p^(k-2). This division is defined on the ideal of multiples of p^2 and preserves distinctness and all magic equalities in both directions. Therefore

N_p(k)=U_p(k)+N_p(k-2),

where terms with nonpositive indices are zero. Small moduli below the displayed ranges have fewer than nine square residues and contribute zero.

## 5. Generating functions and limits

Let F_p(z)=Sum_{h>=0} A_p(h) z^h. Summing the three geometric terms in each ordinary-count polynomial, while correcting its low-level exceptions, gives

F_2(z)=1/(1-8z)-12/(1-4z)+32/(1-2z)-21-24z
       =1536z^4/[(1-2z)(1-4z)(1-8z)],

F_3(z)=1/(1-27z)-8/(1-9z)+15/(1-3z)-8
       =216z^2(1+27z)/[(1-3z)(1-9z)(1-27z)].

The unit shifts and descent give G_2(z)=z^3 F_2(z)/(1-z^2) and G_3(z)=z F_3(z)/(1-z^2), proving the claimed rational functions.

For the limits, A_p(h)=p^(3h)+O(p^(2h)). Sum the finite parity-compatible descent chain. With shifts s_2=3 and s_3=1, the normalized leading term tends to p^(-3s_p)/(1-p^(-6)); the lower-degree contributions tend to zero. This gives 1/504 and 27/728.

## 6. Verification and boundaries

Run `python3 -B verify.py` and `python3 -O -B verify.py` from this directory. All load-bearing checks use explicit exceptions. The verifier:

- exhaustively tests every one of the 3^9 raw arrays modulo 3 against the eight magic equations, confirming the parametrization without assuming it in that control;
- enumerates offset pairs through moduli 2^9 and 3^6;
- independently enumerates square-residue pivot triples through moduli 2^9 and 3^5 and tests all eight line equations on accepted grids;
- multiplies the exact count series by the proposed denominator through degree 30;
- compares all results against the persisted certificate.

These finite checks support, but do not replace, the all-level proofs. The formulas are for quotient rings, not extension fields; for labelled values, not root choices; and for local square residues, not integer squares. They do not claim an integer solution, a global compatibility criterion, or new external priority. The existing session55-padic-lifts packet supplies related threshold witnesses and the distinction between local lifting and integer squareness; this packet independently derives the exact counting refinement.

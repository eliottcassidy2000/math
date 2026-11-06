# Exact modularly weighted parameter capacity for a physical O/Z source

Status: proof for independent review. Counts are tail inputs with multiplicity, not distinct scalar outputs. Moduli and the source bank remain fixed as X grows.

Let s=(h,R,e,d) be a physical source with R>d>0 and scalar map

H_s(m)=h+Rm-dQ(m)-e(m mod 2).

For fixed odd p>=1, write g=gcd(R,d,p), and define

c_0(s) = log(R/(R-d/2))/d,
c_1(s) = log((R-d/2)/(R-d))/d.

For parity b in {0,1} and residue a modulo p, let A_(s,b,a)(X) count m>=0 with H_s(m)<=X, m=b modulo 2, and H_s(m)=a modulo p. Then

A_(s,b,a)(X) = (g/p)c_b(s) X + o_(s,p)(X), if a=h-eb modulo g,
A_(s,b,a)(X) = 0, if a!=h-eb modulo g.

Summing b yields the exact coefficient

w_s(a;p)=(g/p)[c_0(s) 1_(a=h mod g)+c_1(s) 1_(a=h-e mod g)].

Thus the unrestricted parameter count in residue a is w_s(a;p)X+o_(s,p)(X). As every physical scalar is odd, this corresponds to the unique odd residue modulo 2p above a. The parity classes for m may coincide in scalar residue; their coefficients then add.

## 1. The parity-specific capacity coefficients

The digit approximation in the parameter-capacity proof is

Q(m)=m u(m)+O(log(2+m)),

with u(m)=sum_(j>=0) bit_j(m)/2^(j+1). Truncate u after L bits and put K=2^L. For fixed low residue r modulo K, H_s(m) has slope R-d u_L(r), up to an error O_s(log(2+X)+X/K) throughout m<=C_s X. The inverse bound restricts H_s(m)<=X to such a range.

When r has parity b, reversal of the L bits maps it bijectively onto j in the half interval bK/2<=j<(b+1)K/2. Counting the corresponding residue classes exactly as in the total-capacity proof gives

#{m:H_s(m)<=X, m=b mod2}
 = (X/K) sum_(bK/2<=j<(b+1)K/2) 1/(R-dj/K)
   + O_s(K+X/K+log(2+X)).

Taking K comparable with sqrt(X) and passing to the Riemann integral proves

A_(s,b)(X) = X integral_(b/2)^((b+1)/2) du/(R-du) + O_s(sqrt(X))
           = c_b(s)X + O_s(sqrt(X)).

In general c_0 and c_1 are not equal. Treating the two input parities as equally weighted under the scalar cutoff would be incorrect.

## 2. Conditional modular distribution on endpoint cylinders

Use the established primitive digit chain on (v mod2p,Q(v) modp), with transition

(a,q) -> (2a+t mod2p,q+t(a+(a mod2)) modp), t in {0,1}.

Its stationary law is uniform on all 2p^2 states. Strong connectivity and aperiodicity hold for every odd p, including composite moduli: the proof in the root-residue packet uses a suffix to reach a=0, then the block 2(2^k-1) with 2^k=1 modp to cycle q by -1, and a final suffix to reach the requested state.

Fix L leading digits and L trailing digits of a long binary word m. As the unprescribed middle length grows, the state before the trailing block converges to uniform. Appending that fixed block has the form

a' = 2^L a+r modulo 2p,
q' = q+f_r(a) modulo p.

The second formula follows either by iterating the digit recurrence or from the exact block identity. Since 2 is invertible modulo p, under the uniform input law the resulting a' modulo p is uniform, and its parity is the fixed low bit b of r. Also q' modulo p remains uniform and independent of a': for each old a, adding f_r(a) merely permutes q.

Consequently H_s(m) modulo p converges in that cylinder to the distribution of

h-eb+R U-d V modulo p,

where U,V are independent uniform residues modulo p. The map (U,V)->R U-d V is a group homomorphism onto the multiples of g=gcd(R,d,p), and all its fibers have equal size. Its output law is therefore uniform on that subgroup. The conditional probability of scalar residue a is g/p when a=h-eb modulo g, and zero otherwise.

The zero outside the coset is exact for every m of parity b: g divides R and d, so H_s(m)=h-eb modulo g. No mixing limit is needed for this exclusion.

## 3. Passing the cylinder law through the scalar cutoff

The root-residue cutoff argument applies with the constants of this fixed source. The inverse bound gives m<=C_s X. On each cylinder the truncated slope is R-d u_L>=R-d>0; variation in m and the truncation error are O_s(X 2^(-L)+log(2+X)). Thus cylinders meeting both H_s<=X and H_s>X lie in a boundary region of that width about m(R-d u_L)=X.

Counting separately in each low residue modulo 2^L shows that the boundary has O_s(X2^(-L)+log(2+X)+2^L) parameters. The same bound holds after restricting to a parity. Discard m<epsilon X. For fixed epsilon,L the remaining parameters are in finitely many dyadic blocks, and all their free middle blocks tend to infinite length as X grows. Sum the uniform finite-cylinder convergence from Section 2 on cylinders wholly inside the cutoff. Let X tend to infinity, then epsilon tend to zero and L to infinity. This proves

A_(s,b,a)(X) = (g/p) A_(s,b)(X)+o_(s,p)(X)

on the permitted coset. Together with Section 1 it proves the theorem. No uniform rate in p, s, or a growing bank is asserted.

## 4. A modularly weighted finite-bank obstruction

For a finite physical bank B without special collision information, its distinct outputs in the odd class above a modulo p number at most

X sum_(s in B) w_s(a;p)+o_B,p(X).

The class itself has X/(2p)+O(1) odd integers. Therefore any strict inequality

sum_(s in B) w_s(a;p) < 1/(2p)

excludes eventual coverage.

If B has d=1 members, all their scalar images lie in the root image by the established d=1 collapse. Replace the entire collection by one root contribution log(2)/(2p), using the root converse collision theorem and its residue-capacity refinement. Add w_s(a;p) only for the remaining d>1 members. The same strict comparison with 1/(2p) gives a sufficient obstruction. Known image containments may also be used to remove redundant bank members. No conclusion is implied when the inequality fails or is an equality.

## 5. A three-source example with a full-modular-image supplemental source

Take the following physical states, generated from (1,2,0,1):

- root: (1,2,0,1);
- JJ: (5,6,2,3);
- OZJ: (11,12,4,5).

For the last path the intermediate states are O:(1,4,0,1), OZ:(5,7,1,1), and OZJ:(11,12,4,5). This verifies physical attainability using the established grammar.

For p=3 and a=1, the JJ source contributes zero: its outputs are only 0 or 2 modulo 3. The OZJ source has g=gcd(12,5,3)=1, and hence contributes

w_OZJ(1;3)=log(12/7)/15.

Together with the root collision coefficient, the bank's distinct outputs in the class 1 modulo 6 number at most

[X/3][log(2)/2+log(12/7)/5]+o(X).

The bracket is strictly below 1/2. Indeed the trapezoid upper bounds for the strictly convex function 1/x give

log(2)<3/4,
log(12/7)<(5/7)(1+7/12)/2=95/168.

The bracket is therefore below 3/8+19/168=41/84<1/2. In particular this bank misses at least X/252-o(X) odd integers in that class. The sharper logarithmic lower coefficient is [1/2-log(2)/2-log(12/7)/5]/3.

The earlier unweighted capacity sum is log(2)/2+log(2)/3+log(12/7)/5, which exceeds 1/2: already (5/6)log(2)>1/2. Thus the weighted test excludes a bank the unweighted test does not.

Moreover OZJ itself has every parity-compatible scalar residue modulo every fixed modulus. This follows from the exact modular-image theorem: gcd(R,d)=gcd(12,5)=1 makes its defining modulus factor gcd(R,d,M/2) equal to 1 for every even M. Therefore the non-root bank members have no common omitted fixed residue class. The present exclusion genuinely uses quantitative modular capacity, not just a missing-class argument.

This does not rule out every finite physical bank and does not supply a covering bank.

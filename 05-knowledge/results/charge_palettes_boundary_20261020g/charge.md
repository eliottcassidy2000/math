# A debt cone surviving changing generators at arbitrary depth

**PROVED:** conditional composition and payment. **OPEN:** universal coverage
by such controllers and existence of a first exit at every positive source.

## Exact inherited state

For an actual odd word w=(a_1,...,a_m), retain z=n-1 and
```
z_end=lambda_w z+P_w-D_w,       lambda_w=3^m/2^K,
P_w=sum_(a_i=1) lambda_(i+1:m),
D_w=sum_(a_i>=3) (1-4/2^a_i) lambda_(i+1:m).
```
This is [the preceding packet](../first_exit_complement_games_20261020/README.md)
and its [original decomposition](../adaptive_budgets_smooth_interfaces_20261020/CONNECTIONS.md).
Sources, exact valuations and the integer clock K are unchanged.

For concatenation u then v,
```
lambda_uv=lambda_v lambda_u,       D_uv=lambda_v D_u+D_v.       (1)
```
Fix c>=0. The cone lambda>=1, 0<=D<=c(lambda-1) is closed under (1):
the new upper bound is c[lambda_v(lambda_u-1)+(lambda_v-1)].
For multipliers greater than1, the composite ratio D/(lambda-1) is the
weighted mean of the block ratios, with weights (lambda_i-1) times subsequent
multipliers. Blocks need not coincide or commute; there is no depth limit.

**Payment theorem.** Suppose u is in this cone with c<=1 and one final actual
odd step a produces h<n from an odd source n>1. Then positive charge pays.
The unshifted affine word has positive intercept, so mu lambda_u<h/n<1,
where mu=3/2^a. Consequently a>=2. If a=2,
D_total<=c mu(lambda_u-1)<c(1-mu)<=1/4. If a>=3,
```
D_total<c(1-mu)+1-4/2^a <=2-7/2^a<2.
```
Since n-h is a positive even integer, D_total<n-h, equivalently the positive
upper bound h-1+D_total<n-1. This proof needs an actual descent, but need not
assume it is the first. It does not assert that a specified source has the
required next valuation.

More generally an actual terminal packet v with mu=lambda_v<1 and debt D_v
can follow the cone: if c(1-mu)+D_v<=2, every actual descent at that endpoint
pays, because D_total<c(1-mu)+D_v. Availability of a controller remains separate.

## An infinite first-exit family with unbounded reset count

The block b=(1,1,1,1,3) has lambda_b=243/128, D_b=1/2 and
```
D_b/(lambda_b-1)=64/115.
```
For every k>=1, b^k is in the cone with this c. After any actual final
descending step, D_total<1+64/115=179/115<2.

There are infinitely many positive integer first-exit sources for each k.
Choose a>=2 with 2^(7k+a)>3^(5k+1). For w=b^k a write
T_w(n)=(An+B)/d, A=3^(5k+1), d=2^(7k+a). Its exact cylinder is
```
n=r+2dt,       r=(d-B)A^(-1) mod 2d,       0<r<2d.              (2)
```
This is the inherited exact valuation-cylinder formula. Endpoint oddness
and the affine carry congruence imply prefix integrality and oddness at
each intermediate, hence exactness of each valuation. Take t>=0 sufficiently
large that n>1 and (d-A)n>B; there are infinitely many such t. The endpoint
is positive odd and below n. Every proper block prefix increases a positive
source: its first four steps have multiplier3/2, and its completed block
has multiplier243/128, each with positive intercept. Thus the descent is
first, has5k+1 odd steps, and contains at least k resets.

Exact examples are retained at k=1,2,3,4,8,16,32,64, through321 odd steps
and a154-digit source. The all-k quantifier follows from the proof.

## Hostiles and the next obligation

Expansion alone does not imply the cone: (1,1,1,3) has multiplier81/64,
debt1/2 and ratio32/17>1. It needs a larger block or a different terminal
controller. Appending exponent-one steps to an arbitrary prefix cannot
always pay it either: lambda-D is multiplied by3/2, so a nonpositive
surplus remains nonpositive.

Among all50,000 odd sources3..100001, the last-step cone admits47,424 first
exits, including706 with at least three resets. The other2,576 are unresolved
by this cone, not unpaid exits. Source63 first descends to61 and pays the
original positive-charge test, but its prefix before the final step has
D>lambda-1. Thus this universal coverage question remains real.

Next target: source-dependent cones or terminal packets covering such
residual prefixes without presuming arrival. The old source55 later-descent
hostile is not removed by silently dropping the cone premise.

The exact script checks21,844 words against independent suffix sums,
1,024 changing-generator pairs, eight cylinders, and the50,000 declared
first exits. It proves no ROOT-current estimate for a different operator.


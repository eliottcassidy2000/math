# Finite digit debt is a Jordan block, not a new fractal dimension

Status: PROVED elementary supplement to Session54's cross-base coefficient-support packet. General finite-dimensional transfer algebra is standard; no literature originality is claimed. The original packet remains unchanged.

## Defined starting level and states

Fix an integer budget b>=0. For r>=1 let v_r(j), 0<=j<=b, count exponent classes modulo M_r=6*9^(r-1) for which the lowest r base9 digits of2^k contain exactly j digits in{3,4,5}, with all remaining digits in{0,1,2}. Digits6,7,8 are forbidden. The total is A_(r,b)=sum_j v_r(j). There is no r=0 clock or r=0 initial state in this statement.

The primitive-root bijection and the nine-lift digit property are inherited from sources/cross-base-proof.md. At r=1, the allowed unit digits are1,2 (debt0) and4,5 (debt1). Thus

v_1=(2) if b=0; v_1=(2,2,0,...,0)^T if b>=1.

Let N be the (b+1)-by-(b+1) nilpotent shift N e_j=e_(j+1) for j<b, and N e_b=0. Then

v_(r+1)=T v_r,  T=3(I+N).

A lift that spends no new debt has three choices0,1,2; a lift spending one unit has three choices3,4,5. Debt overflow is discarded, not reflected or wrapped around. This proves the transfer matrix, not merely its agreement with counts.

Since N^(b+1)=0, the binomial theorem gives

v_r=3^(r-1)(I+N)^(r-1)v_1,
v_r(j)=2*3^(r-1)[binom(r-1,j)+binom(r-1,j-1)]
      =2*3^(r-1)binom(r,j).

Binomial coefficients outside0<=j<=r are zero. This also covers b=0 directly, giving A_(r,0)=2*3^(r-1). Hence

A_(r,b)=2*3^(r-1) sum_{j=0}^b binom(r,j).

## Exact ordinary generating function and its minimal pole

Define F_b(z)=sum_{r>=1} A_(r,b)z^(r-1). Set x=3z. Formal differentiation of1/(1-x), or the binomial-series identity, gives

sum_{r>=1}binom(r,0)x^(r-1)=1/(1-x),
sum_{r>=1}binom(r,j)x^(r-1)=x^(j-1)/(1-x)^(j+1), j>=1.

Therefore, as formal rational series,

F_b(z)=2/(1-3z)+2 sum_{j=1}^b (3z)^(j-1)/(1-3z)^(j+1).

For b=0 the sum is empty and F_0=2/(1-3z). For b=1,

F_1(z)=(4-6z)/(1-3z)^2.

For b>=1 the common numerator over(1-x)^(b+1) is

P_b(x)=2(1-x)^b + 2 sum_{j=1}^b x^(j-1)(1-x)^(b-j).

At x=1, only j=b survives and P_b(1)=2. Thus no factor1-x cancels. The minimal reduced denominator is exactly(1-3z)^(b+1), up to a nonzero constant. For b=0 the reduced denominator is exactly1-3z, since its numerator is2. This is a statement over Q[z]; reducing in finite characteristic can change pole order and is not asserted.

## Jordan interpretation and the finite-debt asymptotic

T has the sole eigenvalue3 and one Jordan block of size b+1: N has rank b and nilpotency index b+1, and scaling its nonzero super/subdiagonal entries does not split the block. The scalar observation is the all-ones row applied to T^(r-1)v_1. It sees the largest nilpotent power because

(1,...,1) N^b v_1=2.

Hence this is not an unobserved or unexcited Jordan direction: the scalar sequence genuinely needs the full repeated pole. Equivalently its minimal constant-coefficient recurrence over Q has characteristic polynomial(X-3)^(b+1). For each fixed b,

A_(r,b)=2*3^(r-1)(r^b/b! + O(r^(b-1))) for b>=1,
A_(r,0)=2*3^(r-1).

The polynomial r^b counts ways to place at most b defects; its degree is the Jordan-block size minus one. It is NOT a larger exponential growth factor or a larger Hausdorff dimension.

Indeed log A_(r,b)/(r log9) tends to log3/log9=1/2 for every fixed b. The inherited geometric proof is stronger than this count heuristic: the no-defect Cantor subset supplies a lower dimension1/2, while a fixed finite budget is a countable union of finite-position modifications of that same Cantor set, supplying upper dimension1/2. The logarithmic exponent transport preserves dimension on each parity sheet. Haar measure remains zero because A_(r,b)/M_r=3^(-r)sum_j binom(r,j) tends to zero.

## Scope and hostile comparisons

- The transfer starts at r=1. Using an invented r=0 state(2,0,...) would give the wrong counts and an off-by-one pole numerator.
- A diagonal matrix3I loses the polynomial factor; a cyclic wraparound shift introduces extra eigenvalues and allows debt to be recycled. Both are false models.
- Budgets b=0 and1 correspond to the proved coefficient supports modulo3 and9. For b>=2 this is an exact combinatorial digit-budget extension, NOT a theorem about c(n) modulo3^(b+1).
- Fixed b is essential to the stated asymptotic and dimension conclusion. If b grows with r, the polynomial factors cannot be treated as constants, and this Jordan argument no longer concerns one fixed matrix.
- Completed exponent geometry and finite-level periodic counts still do not settle finiteness or infinitude of the full ordinary exponent intersection.

## Source and verification

sources/cross-base-proof.md is a byte copy of the frozen Session54 cross-base REPORT.md, whose lineage recovers the repository's September22 logarithmic transport and October6 primitive-root clock, as well as Session53's coefficient support laws. It carries no new modular assumption here. Exact SHA256 is in manifest.json. The source packet's unchanged manifest SHA256 is92db3b170771188fac5c9d55faec8c4c5d18e4fc8f34d93c74b78abfd0794e95.

check.py compares independent direct modular-orbit counts with transfer vectors, binomial formulas, GF coefficients, exact symbolic common numerators, and scalar finite differences. It tests budgets0 through6, including b>r. Wrong initialization, omitted Jordan shift, cyclic debt reuse, and finite-field pole cancellation are hostile controls. Normal Python and-O outputs must match byte-for-byte. replay.py checks exact recursive membership and hashes before rerunning; hostile.py rejects corrupt, missing and unexpected members. All tools are Python standard library only, with no network or repository dependency during replay.

# Two-adic divisibility of unlabeled tournaments

Let A(n) be the number of tournament isomorphism classes on n vertices. Then:

1. For odd n >= 3, v2(A(n)) = (n-1)/2.
2. For even n >= 4, A(n)/2^(n/2) is an integer and is congruent to phi(n)/2 modulo 2. Thus its valuation equals n/2 precisely for n=4 or n=2p^a with p an odd prime congruent to 3 mod 4; in all other even cases it is at least n/2+1.
3. Hence 2^floor(n/2) divides A(n) for all n >= 3.

The statements are independent deductions from Davis's formula, not claims of first discovery. A focused search found the parity claim in THM-4750 section 6; the strengthening was not found in the targeted repository search or literature search, so novelty remains unverified.

## Proof

The Cauchy-Frobenius formula (Royle et al., Theorem 3.3, https://link.springer.com/article/10.1007/s10801-022-01197-0) gives Davis's odd-cycle partition formula (https://oeis.org/A000568):

A(n) = sum_(sum i*m_i=n, i odd) 2^E/z,
z = product_i i^m_i m_i!,
E = sum_i m_i(i-1)/2 + sum_i i*binom(m_i,2) + sum_(i<j) m_i*m_j*gcd(i,j).

All m_i are nonnegative. Each 2^E/z is a rational number with well-defined 2-adic valuation E-d where d=sum_i v2(m_i!). With c=sum_i m_i, n and c have the same parity. Since every gcd and every i are >=1, and product m_i! divides c!,

E >= (n-c)/2 + binom(c,2),   d <= v2(c!).   (*)

For c>=3, binom(c,2)-v2(c!) >= floor(c/2)+1: for c=3 it is 2; for c=4 it is 3; for c>=5 use v2(c!) <= c-1, giving (c-1)(c-2)/2 >= floor(c/2)+1. Note c=3 gives 2=1+1 and c=4 gives 3=2+1.

Odd n: c is odd. If c=1, the unique term has m_n=1 and equals 2^((n-1)/2)/n, with exact valuation (n-1)/2. If c>=3, (*) gives valuation at least (n-c)/2+floor(c/2)+1 = (n+1)/2. Thus the unique minimum valuation term survives in the sum and v2(A(n))=(n-1)/2.

Even n: c is even, so c>=2. For c>=4, (*) gives valuation >= n/2+1. For c=2, let odd parts be i,j. If i=j>=3, the valuation is 2i-2 >= i+1 = n/2+1; i=j=1 is the excluded case n=2. If i<j, the valuation is (n-2)/2+gcd(i,j), exactly n/2 iff gcd(i,j)=1. Each such term, divided by 2^(n/2), has odd denominator and odd numerator. So the reduction of the total modulo 2 is the number of unordered coprime odd pairs summing to n. For even n>=4 this is phi(n)/2: positive integers i<n coprime to even n are odd, and i and n-i are distinct and paired. This establishes the congruence. The prime-factor criterion follows from phi(n) having exactly one factor 2 iff n=4 or n=2p^a for p congruent to 3 mod 4.

These are rational 2-adic divisibility statements for individual summands; they do not assert those summands are ordinary integers. Their sum A(n) is an ordinary integer by Burnside.

## Reproducibility and scope

Run `python partition_count.py` and `python crosscheck.py`. The first computes exact rational Davis sums, not floating point. At n=77, it enumerates 58,499 odd partitions and obtains a 768-digit integer ending 774185710113073832239696969728 with v2=38; the OEIS table linked in A000568 ends at n=76. The separate direct code enumerates all permutations on <=8 vertices, explicitly computes induced pair orbits, and checks the result against the partition calculation; it uses the same Burnside theorem but does not use the partition exponent formula. Finite checks complement, not replace, the proof.

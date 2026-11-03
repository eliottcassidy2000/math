# Three-vertex plateaus at Pell cancellations

Let d(n)=floor(sqrt(2n−1))−floor(sqrt(n)), the number of edges added at vertex n in the square-sum graph. For every square n=r²>=4 at which the birth and death terms in the exact difference formula cancel, there is a plateau

    d(n−1) = d(n) = d(n+1).

The cancellation centers are exactly the two Pell families stated in the accompanying Session 148 growth report: n=r² where r+t sqrt(2)=(3+2 sqrt(2))^j for j>=1, or where s+r sqrt(2)=(1+sqrt(2))^(2j+1) for j>=1. The first two such centers are 9 and 25; the first few are 9,25,289,841,9801,28561,... . The first two plateaus are

    d(8)=d(9)=d(10)=1,
    d(24)=d(25)=d(26)=2.

Proof: the assumed cancellation at n means d(n)=d(n−1). To show d(n+1)=d(n), use the same difference formula. There is no death at n+1=r²+1, since this lies strictly between r² and (r+1)². A birth would require 2r² or 2r²+1 to be a square. The first is impossible for r>0, by irrationality of sqrt(2). For the second, cancellation at n means either 2r²−1 or 2r²−2 is a square. If 2r²−1 is a square, a further square at 2r²+1 would differ from it by 2, which is impossible for integer squares. If 2r²−2 is a square, a further square at 2r²+1 would differ from it by 3; the only nonnegative squares with difference 3 are 1 and 4, but 2r²−2>=6 for r>=2. Therefore neither birth nor death occurs at n+1, proving the plateau.

There is a literal arithmetic overlap with the unique G_15 Hamilton path: the plateau centers satisfy 9+16=25, and its initial segment 8,1,15,10 has adjacent sums 9,16,25. This is a verified numerical relationship, not a claim of a structural or causal connection between the Pell plateaus and that path.

## Checks

`python check_plateaus.py` and `python -O check_plateaus.py` check all n<=1,000,000 against the direct first-difference conditions; verify the plateau at every cancellation; check the first Pell-family terms; and test 200 terms of each integer-generated Pell family using arbitrary-precision integers. Explicit exceptions remain active under `-O`. The proof, not the finite tests, establishes the result. No claim of historical novelty is made.

# Independent audit: Pell cancellation plateaus

**Result.** The mathematical assertions in `PELL_PLATEAUS.md` pass. No correction is required. Integer-exact tests below support the proof but do not replace it.

## Proof and family classification

For n>=2, the increment of the upper floor is the number of integer squares in {2n−2, 2n−1}; there is at most one. The lower floor increments precisely when n is square. Thus d(n)−d(n−1)=B(n)−1_{n square}. At n=r²>=4, the premise d(n)=d(n−1) therefore forces exactly one of 2r²−2 and 2r²−1 to be square.

At n+1=r²+1 there is no death since r²<r²+1<(r+1)². A birth would require 2r² or 2r²+1 to be square. The first is impossible for nonzero integer r by infinite descent (equivalently irrationality of sqrt(2)). In the second case the supposed birth square and the earlier birth square differ by 2 or 3. Integer squares do not differ by 2; the only nonnegative squares that differ by 3 are 1 and 4, whereas 2r²−2>=6. Hence there is no birth at n+1 and d(n−1)=d(n)=d(n+1).

If 2r²−2=q² then q=2t and r²−2t²=1; all positive nontrivial solutions are r+t√2=(3+2√2)^j, j>=1. If 2r²−1=s² then s is odd and s²−2r²=−1; all solutions with r>=2 are s+r√2=(1+√2)^(2j+1), j>=1. These cases are disjoint. Completeness can be obtained by repeatedly multiplying by 3−2√2: the positive-equation solutions descend to 1+0√2, and the negative-equation solutions descend to 1+√2. The first centers are 9, 25, 289, 841, 9801, 28561.

## Potentially misleading interpretations checked

- At center 9 the admitted birth square is 16=2·9−2; at 25 it is **49=2·25−1**, not 25. The plateau centers themselves are 9 and 25; they are not, in general, the source birth squares.
- d(8)=d(9)=d(10)=1 and d(24)=d(25)=d(26)=2, but these values first occur earlier: d(3)=1, d(13)=2. Degree is not monotone (d(3)=1 and d(4)=0). For k>=1 the independently checked block search matches the reported first threshold n_k=max(a_k², floor((a_k+k)²/2)+1), where a_k=k−2+ceil(sqrt(2(k−1)²+3)); the first six are 3,13,33,61,99,163. These are first occurrences, not persistence thresholds.
- The G15 path starts 8,1,15,10, so its consecutive edge sums are 9,16,25. An independent Hamilton-path dynamic program counts exactly two ordered paths, agreeing with uniqueness up to reversal. Its endpoints are 8 and 9, not 8 and 10. The coincidence of labels and plateau centers establishes no structural connection.
- A square difference of 3 is possible for 1 and 4, so a blanket impossibility statement would be false; it is ruled out here only using r>=2.

## Reproduction

Run `python audit.py` and `python -O audit.py`; captured outputs are in `audit.normal.txt` and `audit.optimized.txt`. Both passed. The script uses explicit exceptions, not assertions; independently enumerates graph neighbors through n=4000; checks the difference recurrence and classifies centers through n=1,000,000; checks 200 large exact-integer terms of each Pell family; compares the closed threshold with a block search through k=1000; counts G15 Hamilton paths; and includes explicit hostile cases for the common errors above.

# A000568 with ordered shared ports, reusing the unit-tail kernel

2026-10-28. **PROVED classical Burnside specialization / FINITE-EXACT.**
This is an extra counting channel of the inherited compiler, not a new A000568
record, a census of large tournaments, or a claimed speed record.

Let A_(n,r) count ordinary tournament isomorphism classes with an ordered
r-tuple of named-role vertices, allowing repetitions. For an odd partition
lambda of n, let m_1(lambda) be its number of singleton cycles,

    e(lambda)=sum_i (lambda_i-1)/2 + sum_(i<j) gcd(lambda_i,lambda_j),
    z(lambda)=product_d d^m_d m_d!.

Burnside on pairs (labelled tournament, ordered port tuple) gives

    A_(n,r)=sum_(lambda odd) [2^e(lambda)/z(lambda)] m_1(lambda)^r.

A permutation with an even cycle cannot fix a tournament. For an odd-cycle
permutation its fixed orientations contribute 2^e and its fixed tuple
contributes m_1^r. Use 0^0=1. For distinct ports replace m_1^r by the falling
factorial (m_1)_r. Collision partitions of tuple roles convert these two
channels by the Stirling numbers of the second kind.

At r=0 this is Davis's formula for
[A000568](https://oeis.org/A000568), already used throughout the repo. Its
unit-tail factorization is
[THM-514, the A000568 unit-tail core kernel and prime-tail ladder](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).
For the same B[m,t] kernel, with j=n-m,

    A_(n,r)=sum_(m,t) B[m,t] 2^(binom(j,2)+jt) j^r/j!.

No additional core-partition search is needed when r changes. Define

    C_n(y)=sum_(m,t) B[m,t] 2^(binom(j,2)+jt) y^j/j!.

Then A_(n,r)=((y d/dy)^r C_n)(1), and distinct-port counts are C_n^(r)(1).
The coefficients are rational fixed-permutation weights; final orbit counts
are integers. They are not a distribution of automorphism orders. Dividing
C_n by A_(n,0) corresponds to choosing an ordinary tournament isotype uniformly
and then one of its automorphisms uniformly. It is not the uniformly labelled
tournament law.

`counting.py` reuses 10880 core-partition contributions in 301 kernel cells
through n=60, computing r=0..8 and every distinct-port count r=0..n. A separate
direct odd-partition evaluator agrees through n=26; the published unmarked
prefix agrees through n=19. Full or all-but-one naming gives 2^binom(n,2).
The two requested boundary rows begin

| n | r=0 | r=1 | r=2 |
|---|---:|---:|---:|
| 11 | 903753248 | 9937755904 | 109287476096 |
| 12 | 154108311168 | 1849103423168 | 22187481894592 |

Earlier work evaluated ordinary/self-converse channels through300
([terminal compiler](../strong_attachment_tail_compiler_20261028/README.md))
and generated the full order-ten class bank
([consecutive flags](../consecutive_flags_order10_20261028/README.md)). This
session adds a reproducible weighted-port specialization, not a larger class
bank. It recovers the June kernel as the reusable computational object,
rather than rerunning a historical approximation as exact enumeration.

## How this fits actual metagraph composition

The formula counts all ported isotypes. A native complement still needs a
specific Hamilton presentation and compatible source/target maps. Equal
scalar counts cannot recover edge multiplicities or legal transitions; the
full fixed-set hostile in PORT_LAWS P1 makes that loss explicit.

A coherent division of labour is:

1. Use the unit-tail kernel for exact global orbit totals and mass checks.
2. Use rooted atom transporter kernels to compose actual shared-port data.
3. Use run words and realization guards to take actual complements.
4. Recompute the target's modular decomposition before the next substitution.

Royle-even equinumerosity is an equality of ordinary class counts. Transferring
this port kernel or its native-complement transitions to even graphs requires
an explicit marked correspondence; the unmarked equality alone does not
supply it. That is an additional problem, not an automatic consequence.

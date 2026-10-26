# Powers, doubling rows, aliquot chains, and the primitive abc quotient

2026-10-26. **PROVED elementary classifications and deductions; CONDITIONAL
abc bound; FINITE-EXACT implementation controls. Independent audit PASS.**
The [independent audit](../aliquot/AUDIT_CHAINS.md) records the reviewed
pre-promotion proof and checker hashes. This status/link amendment changes
no mathematical statement.
The first-occurrence statements depend on the declared table. No aliquot
termination, odd-perfect-number, abc, or general perfect-power conjecture is
proved here.

## 1. Inheritance and the carrier boundary

The closest proved mechanism is
[arithmetic_braids_20260917_summand, sections1–2](../../arithmetic_braids_20260917_summand.md):
deleting equal-parent addition witnesses produces the doubling forest;
deleting equal-parent multiplication witnesses produces the squaring forest.
Its coordinates are the odd core and two-adic height. Those two forests
share a construction, but their roots and orbit relations differ.
[The aliquot/Lehmer packet](../../collatz_aliquot_lehmer_five_20260927.md)
retains the divisor-sum prime-power product and the even-perfect/Mersenne
section. [THM-3833, abc-conditional cube radical growth](../../../../01-canon/theorems/THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness.md)
already insists that abc sees primitive height and loses common scale.

The concept board is **doubling address / perfect-power ancestry / first
row / aliquot prime guard / primitive abc quotient**. The canonical hostile
is 3²=9: it precedes 6 if odd bases are allowed. A second hostile is
36=4·3²: allowing powers inside a row changes its first base from9 to3.
The least-used sidecar is the complete operation address, rather than the
resulting integer alone. No earlier delayed-row classification was found
in the bounded canon/results/reflection search; no global novelty claim is
made for these elementary identities.

Write v2(n)=a and odd(n)=u, so n=2^a u with u odd. Define three tables:

| Table | Address | Entry | Domain |
|---|---|---|---|
| Doubling | (r,j) | 2^j r | r>=1, j>=0 |
| Powers | (r,k) | r^k | r>=2, k>=1 |
| Combined | (r,k,j) | 2^j r^k | r>=2, k>=1, j>=0 |

“First row” means the *minimum numerical base r*, not a chronological
traversal of infinitely long rows. Base1 is excluded from the latter two
tables because 1^k=1 would give infinite multiplicities. If the vertical
operation is repeated squaring, its exponents are 2^t, a proper subset of
the power table's positive exponents. Consequently N=1 has no power-table
or combined-table address under the stated domains.

## 2. Exact first rows and duplicate multiplicities

**Theorem P1.** For N=2^A U, U odd, the doubling addresses are exactly

    (r,j)=(2^b U,A-b),  0<=b<=A.                       (1)

There are A+1, and the first row is U. This follows immediately by comparing
the odd factors and the exponents of two. Restricting to odd row roots gives
one occurrence; retaining all positive row bases gives the stated duplicates.

**Theorem P2.** For N>=2 with prime factorization N=product p^e_p, put
g_all=gcd{e_p}. Its power-table addresses are exactly

    (r,k)=(N^(1/k),k),  k|g_all.                      (2)

There are tau(g_all), and the first base is N^(1/g_all). If exponent1 is
omitted, there are tau(g_all)-1; the first base is unchanged when g_all>=2,
and the fibre is empty otherwise. Proof: r^k=N iff k divides every exponent.
The first distinct-base duplicate among nontrivial powers is
16=2^4=4^2. Equality of two evaluations is not equality of their addresses.

**Theorem P3.** If U>1, put g_odd=gcd{e_p:p odd, p|U}. The combined-table
addresses are exactly

    k|g_odd,  r=2^b U^(1/k),  0<=b<=floor(A/k),
    j=A-kb.                                          (3)

Hence the multiplicity and first base are

    sum_(k|g_odd) (floor(A/k)+1),    U^(1/g_odd).       (4)

If U=1, then N=2^A with A>=1, and instead

    r=2^b, b>=1, k>=1, kb<=A, j=A-kb;
    multiplicity=sum_(k=1)^A floor(A/k), first base2.  (5)

Comparing prime exponents proves completeness and uniqueness of all the
addresses. In (3), b=0 minimizes the base for a fixed k, and the largest
allowed k minimizes the remaining odd factor. Formula (5) explicitly removes
base1. These formulas count all operation addresses, not only distinct values.

The maps D_j(n)=2^j n and P_k(n)=n^k satisfy

    P_k D_j = D_(kj) P_k,                             (6)

not ordinary commutation. At n=3,j=1,k=2 the two orders give36 and18.
On valuation vectors D_j adds j to the two-coordinate, while P_k multiplies
every coordinate by k. This is the exact operation connection underlying the
tables. In particular 36 has doubling first row9 but combined first base3.

## 3. Every delayed power, its first exponent, and an exact census

Call the kth power of a base n **delayed in the doubling table** when the
first doubling row of n^k is larger than n. For k>=2,

    odd(n^k)=u^k > n  iff  u^(k-1)>2^a.               (7)

This is a complete classification. For squares it is u>2^a. Consequently:

- 3 is the first base of any parity with a delayed square;
- 6 is the first *even* base, for every k>=2;
- powers of two never have a delayed power;
- every other n eventually does, since u>1 and u^(k-1) grows;
- 12² has first row9<12, while 12³ has first row27>12;
- delay does not imply perfection: 10² has first row25>10, but s(10)=8.

For u>1 the first delayed exponent is the unique integer

    k_min=2+floor(log_u(2^a)).                        (8)

The checker evaluates this by integer comparisons, not floating-point logs.
The strict inequality in (7) is retained; equality occurs only at n=1 for
k>=2, since u is odd and 2^a is a power of two.

For an integer X>=1, let O(t)=floor((t+1)/2) count positive odd integers
at most a nonnegative integer t. Let root_d(t)=floor(t^(1/d)), calculated
as an exact integer root. The number E_k(X) of delayed *even* bases n<=X is

    E_k(X)=sum_(a>=1,2^a<=X)
      [O(floor(X/2^a))
       -O(min(floor(X/2^a),root_(k-1)(2^a)))].        (9)

For fixed a, the bracket counts exactly the allowed odd cofactors in (7),
so summing over disjoint two-adic valuations proves (9).

**Density consequence.** For every fixed k>=2, the number of nondelayed
positive bases n<=X is O_k(X^(1/k)). To see this, it is at most

    sum_(a>=0) min(X/2^a,2^(a/(k-1))).               (10)

Split the sum at 2^a comparable to X^((k-1)/k). The lower part is an
increasing geometric sum and the upper part a decreasing geometric sum;
each is O_k(X^(1/k)). Thus delayed bases have natural density1 among all
positive integers and relative density1 among even integers. This widespread
phenomenon cannot by itself characterize perfect numbers.

**Exact square boundary.** At X=4^m, m>=1, the number of nondelayed even
bases is

    3·2^(m-1)-1.                                    (11)

For 1<=a<=m, the nondelayed odd cofactors are u<=2^a, contributing
sum 2^(a-1)=2^m-1. For m<a<2m they contribute
sum 2^(2m-a-1)=2^(m-1)-1, and a=2m contributes u=1.
Their total is (11). Adding n=1 gives 3·2^(m-1) nondelayed bases of any
parity, since every odd n>1 is delayed.

## 4. The even-perfect diagram and its prime guard

The classical Euclid–Euler theorem says that every even perfect number,
and only such a number, has the form

    N=2^a M,  a>=1,  M=2^(a+1)-1 prime.              (12)

This is a **CITED classical input**, not a new classification; see
[King, Theorem3 and proof](https://www.jonathanking.org/PDF/mersenne-perfect.pdf).
Its sufficiency is also immediate from
sigma(N)=(2^(a+1)-1)(M+1)=2N.

For every such N the following is an exact coordinate-preserving diagram:

```text
N  ---------------- square ---------------->  N²
| odd core                                    | odd core
v                                             v
M  ---------------- square ---------------->  M²
                                               |
                                               s
                                               v
                                           2^(a+1)
                                               |
                                               s
                                               v
                                               M --s--> 1 --s--> 0

Separately: N --s--> N.
```

The square commutes because odd(n²)=odd(n)². Since M>2^a, M²>N: every
even perfect number has a delayed square. Primality gives
s(M²)=M+1=2^(a+1); the geometric sum gives s(2^(a+1))=M; and s(M)=1,
s(1)=0. Zero is the terminal symbol, with no divisor sum at zero assumed.

The first two examples are

    6²=4·9,        9 -> 4 -> 3 -> 1 -> 0,    6 -> 6;
    28²=16·49,    49 -> 8 -> 7 -> 1 -> 0,   28 -> 28.

The square row M² is not an aliquot predecessor of N. The chain passes
through the odd core M, not through the fixed point N. Nor does odd-core
projection conjugate the aliquot map: odd(s(6))=3 whereas s(odd(6))=1.
The lost coordinate a is indispensable for reconstructing N.

The prime guard is decisive. For a=3 the formal M=15 is composite;
120 still has delayed square row225, but s(225)=178, not16, and120 is not
perfect. Replacing sigma(M²) by 1+M+M² without primality creates a false
chain. Delayed square rows are common; the native divisor-sum identities
are the additional structure.

## 5. Two different abc maps, and exactly what they forget

For any n>=2 let g=gcd(n,s(n)). The triple

    (A,B,C)=(n/g,s(n)/g,sigma(n)/g)                  (13)

is positive, pairwise coprime, and satisfies A+B=C. Coprimality of A,B is
the definition of g; that of either with C follows by subtraction. This
is an exact map into primitive abc triples. It preserves abundance order
and perfectness, but loses the common scale g. Every perfect n maps to
(1,1,2), including6 and28. Thus the bare primitive triple cannot distinguish
the special incoming structure of6 from another perfect fixed point.
Retaining g reconstructs n=gA and s(n)=gB. At n=1, B=0, so the positive
abc carrier is inapplicable.

There is a different, sharper identity on prime powers. For p prime and
k>=1,

    s(p^k)=1+p+...+p^(k-1),
    1+(p-1)s(p^k)=p^k.                              (14)

The terms of (14) are pairwise coprime, and its radical is

    rad(1·(p-1)s(p^k)·p^k)=p·rad((p-1)s(p^k)).       (15)

For example (13) at n=9 gives9+4=13, whereas (14) gives1+8=9. They are
different maps, not two normalizations of the same triple.

Assume the usual **abc conjecture**: for every epsilon>0 a constant
K_epsilon bounds C by K_epsilon rad(ABC)^(1+epsilon) on primitive positive
triples. Then (14)–(15) give the **CONDITIONAL** radical floor

    rad((p-1)s(p^k)) >=
      K_epsilon^(-1/(1+epsilon)) p^(k/(1+epsilon)-1). (16)

This is direct rearrangement, uniform in p and k for fixed epsilon. It is
not an effective cutoff because no numerical K_epsilon is supplied. It
does not imply a missing aliquot fibre, distinguish perfect integers after
(13), or establish abc. The conjectural input and its exact quantifiers
are stated in [Granville–Tucker, p1227](https://dms.umontreal.ca/~andrew/PDF/abc.pdf).
The repository's [abc transfer ledger](../../../reference/CORE-PAPERS-ABC-IUT.md)
keeps this an antecedent; no claimed IUT proof is used.

## 6. Reproduction and loss ledger

Run `python3 -B check.py > normal.json` and
`python3 -B -O check.py > optimized.json` in this directory. Both streams
must agree. The checker uses only the standard library and explicit gates,
so optimization cannot remove its tests. It independently enumerates table
addresses, sieves proper-divisor sums, and uses exact integer roots.
Finite checks support the formulas and implementation; the proofs above
cover all parameters. The frozen receipt contains268,651 gates: all table
values through20,000, delayed exponents2–9, exact square boundaries through
4^10, and the explicitly declared aliquot/prime-power controls. Normal and
optimized streams agree. At X=20,000 there are9,795 delayed even squares
and206 nondelayed bases of either parity.

| Source -> target | Map and preserved predicate | Lost information / required sidecar |
|---|---|---|
| n -> doubling row | n -> odd(n); common doubling orbit | v2(n) reconstructs scale |
| n -> power root | divide the prime exponent vector by its gcd | exponent gcd reconstructs n |
| operation word -> value | semidirect law(6); exact numerical evaluation | word/address multiplicity |
| even-perfect N -> square-row chain | (12) plus square/odd-core square | a and primality of M; no aliquot conjugacy |
| aliquot edge -> primitive abc | (13); sum, coprimality, order | common scale g; all perfect numbers collapse |
| prime-power aliquot edge -> abc | (14); geometric-series identity | prime base and exponent; conjectural bound stays conditional |

The general open questions are not equivalent just because these carriers
share integers. The surviving connections are the displayed maps and their
explicitly retained arithmetic conditions.

---
id: THM-4084
title: "Even-graph matching-character profile and D5 firewall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. For every n>=6,
  every signed switching class represented by a negative matching has an
  exact all-cycle negative-count profile. In the cumulative D=5 layer, every
  matching of size at least two lies strictly above the single-edge value,
  as does every nonsingle class of frustration index at most three. Thus any
  counterexample or additional equality class for the still-open D=5
  single-edge conjecture has no balanced vertex deletion and has frustration
  index at least four. The full D=5 conjecture remains OPEN beyond the
  inherited exhaustive n=6,7,8 bases.
source: codex-frontier-synthesis-creative-20260825f / matching-character wildcard
audit: >
  PASS. The primary path directly enumerates labelled matchings and simple
  cycles through n=9, performing 18,456 matching-profile gates, 12,201 direct
  three-edge-shape gates, and 80,406,311 direct cycle-parity gates; it also
  performs 16,360 matching-positivity gates, 615 symbolic three-edge gates,
  and 21,762 exact switching-minimality gates. The independent Held--Karp parity
  path exhausts the frustration index of every labelled switching class at
  n=6,7, classifies all 2,006 nontrivial classes of frustration at most three
  there, and checks all matching layers through n=10 with 175,751 subset-cycle
  gates. Both paths recover the
  sharp n=4 D3 antibalanced tie. Normal and optimized outputs byte-match; both
  scripts have zero assert nodes and zero floating literals.
depends_on:
  - THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation
  - THM-4083-even-graph-cumulative-d3-d4-spectral-gap
related:
  - THM-4073-even-graph-diameter-layer-exact-cycle-distance
  - THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope
script: 04-computation/even_graph_matching_character_d5_firewall_thm4084.py
output: 05-knowledge/results/even_graph_matching_character_d5_firewall_thm4084.out
script_sha256: 50b502bd2b61d5db7d2978306faf60ed13e620e94ab9abc68438dc5b1ca21631
output_sha256: 27cb6ceb6ed202ce411ae49e779074cedeba7fc89649f92b7ef31ae32fdb502a
independent_audit_script: 04-computation/even_graph_matching_character_d5_firewall_thm4084_independent_audit.py
independent_audit_output: 05-knowledge/results/even_graph_matching_character_d5_firewall_thm4084_independent_audit.out
independent_audit_script_sha256: 8893601232812ccaad2418b52d47326081b610828907021b7e84dd603ac46195
independent_audit_output_sha256: 7a44734295f78aa5387394a72f36f97c68f59e4a383fd4087155b7a855a693fc
hash_basis: raw LF bytes
---

# THM-4084 -- matching characters and the D5 firewall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** The signed
dual of THM-4078 turns a spectral-gap question into the minimization of odd
cycle counts. A matching of negative edges is unusually tractable because a
cycle can meet it only in disjoint prescribed edges. This gives an exact
profile in every cycle layer and, at `D=5`, removes the whole matching family
and every class of frustration index at most three from the open search.

Throughout, `n>=6`, `1<=r<=floor(n/2)`, and `M_r` is a matching of `r`
negative edges in the signed complete graph `K_n`. Let

```text
C_(n,k) = set of unoriented simple k-cycles of K_n,
c_k^-(H)=#{C in C_(n,k): |E(C) intersect H| is odd}.           (1)
```

The notation depends only on the labelled switching class `[H]`.

## 1. Matching representatives really have frustration r

Switching by a nonempty proper cut `delta(S)` changes the number of negative
edges from `r` to

```text
r+|delta(S)|-2|M_r intersect delta(S)|.                        (2)
```

Replace `S` by its complement if necessary and put
`s=|S|<=floor(n/2)`. At most `s` matching edges cross the cut. Hence the
change in `(2)` is at least

```text
s(n-s)-2s=s(n-s-2)>0.                                         (3)
```

The last inequality holds for every `n>=6` and `1<=s<=floor(n/2)`. Thus
`M_r` is the unique minimum-size representative inside its labelled
switching class and

```text
frustration([M_r])=r.                                         (4)
```

It follows as well that two distinct labelled `r`-matchings cannot be
switching-equivalent: a nontrivial switch would strictly increase either
one. Therefore the labelled matching classes have multiplicity

```text
n!/(2^r r! (n-2r)!),                                         (5)
```

and relabelling makes them one orbit for each `r`.

## 2. Exact all-layer matching profile

For every `3<=k<=n`,

```text
boxed:
c_k^-(M_r)
 =sum_(j=1)^min(r,floor(k/2))
    (-4)^(j-1) C(r,j)(k-j-1)! C(n-2j,k-2j).                  (6)
```

To prove `(6)`, if a cycle contains `q` edges of `M_r`, use the exact parity
expansion

```text
1_(q odd)=sum_(j=1)^q (-2)^(j-1) C(q,j).                     (7)
```

For a fixed set of `j` disjoint matching edges, the number of unoriented
`k`-cycles containing all of them is

```text
2^(j-1)(k-j-1)! C(n-2j,k-2j).                                (8)
```

Indeed, choose the other `k-2j` vertices, contract the prescribed edges to
`j` blocks, cyclically arrange the resulting `k-j` objects up to reversal,
and orient the edge blocks. Summing `(7)` over cycles and reversing the order
of summation gives `(6)`: the factors `(-2)^(j-1)` and `2^(j-1)` combine to
`(-4)^(j-1)`.

By THM-4078, `(6)` is also the exact cycle-layer character profile

```text
lambda_k(M_r)=N_(n,k)-2c_k^-(M_r),
N_(n,k)=n!/(2k(n-k)!).                                       (9)
```

Thus the matching orbit gives an explicit joint eigenvalue row for every
multiplicity-weighted cycle operator, not only for the cumulative layer used
below.

## 3. D5 collapses to a cubic in the matching size

For the diameter layer `D=5`, whose generators have lengths `3,4,5,6`, put

```text
S_5(H)=sum_(k=3)^6 c_k^-(H),
A_n=sum_(k=3)^6 (n-2)!/(n-k)!
   =(n-2)(n^3-11n^2+41n-50),
B_n=1+2(n-4)+3(n-4)(n-5)
   =3n^2-25n+53.                                             (10)
```

Here `A_n=S_5(M_1)` is the single-edge candidate from THM-4083. Summing
`(6)` and observing that at most three disjoint edges fit in a six-cycle
gives the exact formula

```text
boxed:
S_5(M_r)=rA_n-4 C(r,2)B_n+32 C(r,3).                         (11)
```

For `r>=2`, subtracting `A_n` factors as

```text
S_5(M_r)-A_n=(r-1)f_n(r),
f_n(r)=A_n-2rB_n+(16/3)r(r-2).                               (12)
```

There is no divisibility issue in `(12)`: it is a factorization of the
integer expression `(11)`. Its discrete derivative is

```text
f_n(r+1)-f_n(r)
 =(2/3)(-9n^2+75n+16r-167)<0                                (13)
```

for `n>=6` and `r<floor(n/2)`. Thus it suffices to test the largest matching.
For `(n,r)=(6,1),(6,2)` the parenthesized integer in `(13)` is respectively
`-25,-9`; for `n>=7`, it is at most `-9n^2+83n-167`, already negative at
`n=7` and strictly decreasing thereafter.
If `n=2m`, then

```text
f_(2m)(m)
 =2(m-2)(24m^3-144m^2+248m-75)/3,
24m^3-144m^2+248m-75
 =24(m-3)^3+72(m-3)^2+32(m-3)+21>0.                         (14)
```

If `n=2m+1`, then

```text
f_(2m+1)(m)
 =(48m^4-288m^3+604m^2-464m+57)/3
 =[48(m-3)^4+288(m-3)^3+604(m-3)^2
    +568(m-3)+213]/3>0.                                     (15)
```

Since `m>=3`, equations `(12)`--`(15)` prove

```text
boxed: S_5(M_r)>S_5(M_1)=A_n for every r>=2.                 (16)
```

## 4. The complete frustration-three firewall

A switching class of frustration index two has a representative made of two
edges. Up to relabelling, the two edges are disjoint or adjacent.

The disjoint case is `M_2`; from `(11)`,

```text
S_5(M_2)-A_n=A_n-4B_n
 =(n-4)(n^3-9n^2+15n+28)>0,                                 (17)
```

for `n>=6`; after writing `n=t+6`, the cubic factor is
`t^3+9t^2+15t+10`. For an adjacent pair `{ab,ac}`, a negative `k`-cycle contains
exactly one of the two edges. Each edge lies in
`e_(n,k)=(n-2)!/(n-k)!` cycles, while both lie in
`(n-3)!/(n-k)!` cycles. Therefore

```text
c_k^-({ab,ac})
 =2e_(n,k)-2(n-3)!/(n-k)!
 =2(n-3)/(n-2)e_(n,k),                                      (18)
```

and summing gives

```text
S_5({ab,ac})-A_n=(n-4)/(n-2) A_n>0.                         (19)
```

Equations `(17)`--`(19)` classify and exclude every frustration-two shape.

There are only five three-edge shapes. Put

```text
C_n=A_n/(n-2)=n^3-11n^2+41n-50.                              (20)
```

For a three-edge graph `H`, let `a(H)` and `d(H)=3-a(H)` count its adjacent
and disjoint unordered edge pairs, and let `U_n(H)` count cumulatively the
cycles of lengths `3,...,6` that contain all three edges. Parity
inclusion--exclusion gives

```text
S_5(H)=3A_n-2a(H)C_n-4d(H)B_n+4U_n(H).                       (21)
```

The adjacent-pair containment sum is `C_n`; the disjoint-pair containment
sum is `2B_n`. Contracting the path components of all three prescribed edges
gives the complete table

```text
shape H              a(H)  d(H)  U_n(H)
3K_2                    0     3       8
P_3 disjoint-union K_2  1     2       4n-18
P_4                     2     1       n^2-8n+17
K_(1,3)                 3     0       0
K_3                     3     0       1.
```

Write `t=n-6>=0`. Substitution in `(21)` produces a manifestly positive
gap in every row:

```text
(S_5(H)-A_n)/2
 = t^4+11t^3+27t^2+18t+14,      H=3K_2,
 = t^4+10t^3+26t^2+31t+16,      H=P_3 disjoint-union K_2,
 = t^4+ 9t^3+27t^2+36t+20,      H=P_4,
 = t^4+ 8t^3+24t^2+33t+16,      H=K_(1,3),
 = t^4+ 8t^3+24t^2+33t+18,      H=K_3.                       (22)
```

A class of frustration index three has some three-edge minimum
representative, so `(22)` excludes it. Together with the trivial balanced
class, the single-edge equality class, and the two-edge classification, this
proves

```text
boxed:
if S_5(H)<=A_n and [H] is not a single-edge class,
then frustration([H])>=4.                                   (23)
```

THM-4083 supplies an orthogonal deletion firewall. Its balanced-deletion
trichotomy says that two balanced deletions force the single-edge class,
while exactly one balanced deletion is strictly above the single-edge value
for every cumulative layer. Consequently any D5 counterexample, or any
additional D5 equality class, must satisfy simultaneously

```text
b(H)=0,                 frustration([H])>=4.                 (24)
```

This intersection, rather than either condition alone, is the new reduced
search space.

## 5. Spectral consequence and finite bases

Let

```text
M_(<=6)=M_3+M_4+M_5+M_6,
Q_(n,5)=sum_(k=3)^6 N_(n,k).                                 (25)
```

The THM-4078 eigenvalue indexed by `[H]` is

```text
lambda_(<=6)(H)=Q_(n,5)-2S_5(H),
Laplacian eigenvalue=2S_5(H).                                (26)
```

Hence no matching orbit other than `r=1`, and no other class of frustration
at most three, can lower or enlarge the equality space of the candidate D5
spectral gap `2A_n`.

The independent exact audit exhausts the frustration index of all
`2^C(n-1,2)` labelled switching classes for `n=6,7`, and evaluates every one
of the `515+1,491=2,006` nontrivial classes in the frustration-at-most-three
stratum. Within that stratum it finds respectively `15` and `21` D5 equality
classes, exactly the `C(n,2)` single-edge classes. Independently, the complete
Walsh census frozen with THM-4078 evaluates every switching character through
`n=8`; its prefix rows prove the full D5 conjecture `FINITE-EXACT` at
`n=6,7,8`, with only the single-edge equalities. These finite censuses are
base resources, not an all-`n` conclusion.

## 6. Sharp boundary and loss ledger

The lower-order boundary is real. At `n=4`, in the cumulative `D=3` layer,
the negative disjoint matching is antibalanced and has

```text
(c_3^-,c_4^-)=(4,0),       S_3(M_2)=4=S_3(M_1).              (27)
```

Thus neither the strict matching statement nor the unique-equality statement
extends to that boundary.

The all-`n` D5 minimum remains **OPEN**. Equations `(23)`--`(24)` do not
control frustration-four-or-more classes with no balanced deletion. They
identify the first genuinely hostile stratum and do not assume that the two
coordinates are sufficient there.

Finally, `(9)` and `(26)` concern the multiplicity-weighted simple-cycle
operators in the canonical diameter layer/all-simple-cycle envelope. By
MISTAKE-495 and THM-4069, they must not be transferred to an arbitrary
spanning-tree fundamental-cycle image whose diameter omits some of these
cycle lengths.

## 7. Reproducible exact audits

Primary labelled-cycle audit:

```bash
python3 04-computation/even_graph_matching_character_d5_firewall_thm4084.py
```

It directly checks `(6)` on every labelled matching through `n=9`, verifies
the D5 factorizations through `n=128`, checks strict switching minimality for
every matching and every cut representative through `n=12`, verifies every
labelled three-edge graph through `n=9`, and recovers `(27)`.

Independent Held--Karp audit:

```bash
python3 04-computation/even_graph_matching_character_d5_firewall_thm4084_independent_audit.py
```

This path never materializes cyclic-order lists. It counts even/odd cycles by
a subset dynamic program, checks every switching class at `n=6,7`, recovers
the full frustration histogram and the complete frustration-three firewall,
and independently verifies every matching layer through `n=10`. The frozen
normal outputs byte-match both optimized (`python3 -O`) runs.

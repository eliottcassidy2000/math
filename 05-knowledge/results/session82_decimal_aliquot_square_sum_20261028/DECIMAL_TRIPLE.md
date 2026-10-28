# The decimal triple: exact family, edge loss, and a Pythagorean bridge

Status: **PROVED** elementary identities and the explicitly described maps; **FINITE-EXACT** arithmetic certificates. No novelty or external-literature claim. No conclusion about general aliquot termination, interval square-sum Hamiltonicity, or the abc conjecture.

## Inheritance and concept board

The closest inherited mechanism is the coprime divisor-sum transport in Session81 `ALIQUOT_OBSERVERS.md` (building on Session78 `ALIQUOT.md`). Its hostile example warns that a scalar observable need not commute with a specified operation. The square-sum definition is inherited from `01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md`. The least-used sidecar here is the retained Euclid parameter: it enables the exact AP-petal map already classified in `05-knowledge/results/ap_square_pythagorean_20261019.md`, Theorem 1. Its positivity condition must not be discarded. We apply that established map, not claim it anew.

Live concepts: decimal positional notation / near-additive triples / square-sum edges / coprime aliquot transport / Euclidean parameters and positive AP petals. Anchor: requested triple. Niche: exact edge obstruction under aliquot. Wildcard: negative-Pell extra edge.

## 1. An exact characterization, not a decimal coincidence

For every integer `t>=2`, set

`A=t²+1, B=t+1, C=t(t+1)`.

Then `B<A<C`, and

**`A+B=C+2`, `B+C=B²`.**

Conversely these two equations, with integer `B>=3`, force exactly this family with `t=B−1`. Indeed `C=B²−B`, and `A=C+2−B=(B−1)²+1`. Thus retaining the two defining operations gives a lossless one-parameter description. At `t=1000` this is precisely `{1000001,1001,1001000}`.

Euclidean reduction gives

`gcd(A,B)=gcd(A,C)=gcd(2,t+1)`, `gcd(B,C)=B`.

For the second equality, use `gcd(A,t)=1` followed by `t≡−1 mod B`. In particular the anchor is not pairwise coprime. It is not an additive abc triple: its additive defect is 2. The true additive relation `t²+t=C` has common gcd `t`; dividing gives only `t+1=B`. The primitive triple `(1,t²,A)` is a different object. For the anchor `rad(ABC)=10010010010`, and `rad(1·t²·A)=10000010`; these are exact radicals, not evidence for an exceptional abc phenomenon.

The decimal appearance generalizes to any integer base `b>=2`: choose `t=b^k`. For `k=3`,

`A=b⁶+1=(b²+1)(b⁴−b²+1)`,
`B=b³+1=(b+1)(b²−b+1)`, `C=b³(b³+1)`.

The polynomial factors are not necessarily prime. At `b=10`, full prime factorization is

| Number | Factorization | Aliquot sum `s(n)=sigma(n)−n` |
|---|---|---|
| 1000001 | 101·9901 | 10003 |
| 1001 | 7·11·13 | 343 |
| 1001000 | 2³·5³·7·11·13 | 2143960 |

## 2. Exactly which square-sum graph is present?

An edge means the sum of two distinct labels is an integer square. Every family member has the edge `BC`, since `B+C=B²`. None has edge `AB`: for `t>=2`,

`t² < A+B=t²+t+2 < (t+1)²`.

The only optional edge is `AC`. It exists if and only if there is an integer `r>=0` satisfying

**`r²=2t²+t+1`, equivalently `(4t+1)²−8r²=−7`.**

Thus this triple is either one edge plus an isolated vertex or a three-vertex path, never a triangle. It has a spanning path exactly on the displayed Pell locus; this is about the induced three-label graph, not all labels `1..C`.

A constructive recurrence preserving that locus is

`(t,r) -> (17t+12r+4, 24t+17r+6)`.

Substitution proves preservation; both coordinates strictly increase for nonnegative seeds. The seeds `(0,1)` and `(1,2)` produce infinite families beginning `(16,23),(552,781),...` and `(45,64),(1537,2174),...`. **These two orbits contain every solution with `t>=0`.** Here is a descent proof, not an extrapolation of a finite search. Set `x=4t+1`. For a solution with `x>=45`, the inverse map

`x′=17x−48r`, `r′=17r−6x`

preserves `x²−8r²=−7` and `x≡1 (mod 4)`. Both coordinates are positive: `x′>0` is equivalent to `x²>2016`, and `r′>0` follows from `289r²−36x²=(x²+2023)/8>0`. Moreover `x′<x` follows from `x<3r`, itself immediate from the equation. Thus descent reaches `1<=x<45` with `x≡1 (mod 4)`. Checking the eleven possibilities `x=1,5,...,41` leaves exactly `(x,r)=(1,1),(5,2)`. Reversing the descent proves the claim. In the original `t>=2` domain omit the two seeds themselves. The certificate independently finds exactly `(16,23),(45,64),(552,781),(1537,2174)` through `t=10000`.

At `t=1000`, `A+C=2001001` lies strictly between `1414²=1999396` and `1415²=2002225`. The only edge is `1001--1001000`.

## 3. Aliquot transport: an exact connection and its failed implication

Since `gcd(t,t+1)=1`, multiplicativity proves

**`s(C)=t s(B)+B s(t)+s(t)s(B)`.**

This identity holds for every `t>=2` without primality assumptions. At `t=1000`, `sigma(t)=2340`, `s(t)=1340`, and `sigma(B)=1344`; hence `s(C)=2340·1344−1001000=2143960`.

This map does **not** transport square-sum edges to square-sum edges. The image triple `(10003,343,2143960)` has no edge: the three pair sums are respectively `10346≡2`, `2153963≡3`, and `2144303≡3 (mod 4)`. The smallest family parameter already kills the proposed implication: `3+6=9`, but `s(3)+s(6)=1+6=7`.

Failure anatomy: adjacency knows `B+C`; aliquot transport uses divisor factorizations, not merely addition. The strongest survivor is the exact coprime product formula above, with the retained coordinates `t,s(t),B,s(B)`. There is no inference from `s(B)=343=7³` to a general cubic law.

A sharper guarded hostile family: if `t=2^a` with `a>=1` and `t+1` prime, then `s(B)=1` and `sigma(t)=2t−1`, so

`s(B)+s(C)=(t+1)²−2`,

which is nonsquare. This is a conditional family for each such parameter, **not** an assertion of infinitely many primes of that form. Cases `t=2,4,16,256` are explicit examples. Dropping the primality guard changes the formula.

## 4. An actual map to the repository's AP square-sum triangles

The source parameter produces the Euclidean Pythagorean triple

**`(a,b,c)=(t²−1,2t,t²+1)=(A−2,2(B−1),A)`.**

The Pythagorean identity follows by expansion. Its gcd is 1 when `t` is even and 2 when `t` is odd; it is not uniformly primitive. This is not the original triple of labels.

Apply the inherited AP-petal map (allowing its displayed scale):

`(a,b,c) -> (x,y,z)=2(c²−4ab,c²,c²+4ab)`.

It produces the explicit family

**`(x,y,z)=2((t²+1)²−8t(t²−1), (t²+1)², (t²+1)²+8t(t²−1))`.**

All three pair sums are squares, with roots `2|a−b|`, `2c`, `2(a+b)`. Positivity holds exactly for integer `t>=8` in our domain: the first entry divided by 2 is `t³(t−8)+2t²+8t+1`, positive for `t>=8`; direct checks give negative values for `2<=t<=7`. Strictly increasing AP order follows from `ab>0`. The anchor gives

`(1984004016002, 2000004000002, 2016003984002)`.

Typed connection: source = parameterized anchor family; target = positive integer AP square-sum triangles; map = retain `t=B−1`, form Euclid triple, apply inherited polynomial map; preserved predicate = Euclidean identity becomes all three square-sum predicates; destroyed information = original nonedges (the target is a triangle); required sidecar = parameter and positivity `t>=8`; cheapest decisive tests = pair-sum squaring and the hostile `t=2` negative first label. The map is recoverable on its image from `y=2(t²+1)²` with `t>0`, so it does not merely match unrelated integers. For odd `t`, primitive normalization divides `(a,b,c)` by 2 and uses scale `g=2`, exactly as required by the inherited theorem.

## 5. A genuine primitive abc family, with an exact radical criterion

This section is an elementary authored specialization; no novelty claim or general abc theorem is asserted. Let `n=10^k`, `k>=1`, and `A=n²+1`. The primitive additive triple is **`(1,n²,A)`**, not the original near-additive triple. Since `gcd(10,A)=1`, its radical is `R=10 rad(A)`.

Every prime `p|A` is odd and has `n²≡−1 (mod p)`. The multiplicative order of `n mod p` is therefore 4, so `p≡1 (mod 4)`. Also `p≠5`, and consequently `p>=13`. It follows that

**`R<A` if and only if `A` is not squarefree.**

Indeed squarefree `A` gives `R=10A>A`. Otherwise a repeated prime gives `A/rad(A)>=p>=13>10`, hence `R<A`. The prime-size bound is essential to this converse in the displayed form. For the requested anchor `k=3`, `A=101·9901` is squarefree, so it has no radical excess.

There is an explicit infinite subfamily with radical excess:

**`k≡68 (mod 136)` implies `17² | 10^(2k)+1`.**

In fact the implication is an equivalence. Exact modular arithmetic gives

`10^136≡−1`, `10^272≡1`, `10^16≡137 (mod 289)`.

Because `272=16·17` and the powers at `272/2` and `272/17` are not 1, the order of 10 modulo 289 is exactly 272. As `−1=10^136`, the divisibility condition is exactly `2k≡136 (mod 272)`, or `k≡68 (mod 136)`. Thus `R<=10A/17<A` for each such `k`.

For quality `q=log(A)/log(R)`, this proves `q>1`. The fixed repeated factor alone supplies only the lower bound `q>=log(A)/log(10A/17)`, which approaches 1 as `k` grows. It supplies **no uniform bound `q>1+epsilon` for fixed positive epsilon**, and is not a disproof or exceptional-asymptotic claim about abc. The actual quality could be affected by additional repeated factors, which this argument does not classify.

Two normalization boundaries help separate operations:

* Translate **all three original labels by −2**. Then `(A−2,B−2,C−2)=(n²−1,n−1,n²+n−2)` is additive and has gcd exactly `n−1`; primitive division gives `(n+1,1,n+2)`. Translation has changed the labels and is not gcd-normalization of the original triple.
* The original square-sum identity `B+C=B²` has common divisor `B`; dividing gives only `(1,n,n+1)`. A square-sum identity is not automatically a primitive abc instance with its original labels.

These are explicit changes of operations and retained data, not a transfer of aliquot results to abc. The bounded source search found no exact earlier statement of this specialization; no claim of absence from the entire repository or literature follows.

## 6. Reproduction and limits

Run `python3 check_decimal.py` and `python3 -O check_decimal.py`. No third-party packages are required; checks use explicit exceptions. The program cross-checks factorized aliquot sums against paired-divisor sums on `1..5000` and all three anchor labels; tests family identities, gcds, graph predicates and the AP map on `2<=t<=10000`; checks aliquot transport on `2<=t<=1000`; and checks the descent base window and 30 recurrence/inverse steps from each Pell seed. It also checks the complete modulo-136 residue window for the 17² divisibility criterion, the exact modular-order witnesses, the radical criterion on k=1..6 and the constructed k=68 witness via modular arithmetic, and both normalization identities on n=2..1000. The program prints its verification output; normal and optimized Python outputs match byte-for-byte.

All universal statements above have elementary proofs independent of the finite tests. The all-height Pell classification is justified by descent, separately from the finite census. No shared repository checkout was modified by this lane; files are staged externally for the publisher.

Independent audit verified the complete Pell descent, the additional radical/order argument, exact positivity and normalization guards, and the final normal/optimized checker outputs.

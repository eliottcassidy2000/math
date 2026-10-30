# Session109: affine Fibonacci clocks, golden obstruction modules, and retained carries

**Status:** PROVED elementary statements below; FINITE-EXACT independent checks in `check.py` / `results.json`. No priority claim, canon promotion, or Collatz convergence claim. Prepared against source snapshot `session109_source` at `0db004b1`; all source paths below are relative to that repository.

## Executive result

The clock of `x ↦ Q^r x+b (mod N)`, `Q=[[1,1],[1,0]]`, is not determined by its linear clock. The exact missing datum is the **return carry**, not just the cardinality of a golden ring. For stride 5 modulo 11 the identical linear clock has order 2, while affine clocks have order 2 or 22. Strides 1 and 2 are exceptional: every carry is removable over every modulus, and even over the integers.

The universal obstruction is the golden period module

`C_r = Z²/(I−Q^r)Z² ≅ (Z/d_r) ⊕ (Z/e_r)`.

Its invariant factors are given explicitly below. At strides 5,6,7,9,10 its orders are 11,16,29,76,121, but its structures are respectively `Z/11`, `(Z/4)²`, `Z/29`, `Z/2 ⊕ Z/38`, `(Z/11)²`. Equal cardinalities must never silently identify rings or clocks.

## 1. Inheritance and boundaries

Closest mechanisms already in the repository:

* `05-knowledge/results/golden_digit_carry_20261003.md`: exact two-register golden reader; unique polynomial normalization carry; explicit warning that the aggregate carry does not reconstruct chronological rewrites or source factorization.
* `05-knowledge/results/golden_prime_clocks_20261003.md`: quadratic-ring modular clocks and splitting.
* `05-knowledge/results/collatz_golden_carriers_20261004.md`, denominator modules section: already records Smith pairs `(1,1),(2,2),(1,11),(76,76)` for strides 2,3,5,18. The present general Smith theorem extends/organizes those examples, rather than rediscovering their interpretation.
* `01-canon/theorems/THM-4551-fibonacci-monodromy-monotile-boundary-words-are-leaper-flights.md`, item 5, and `05-knowledge/results/chessboard_weave_20261006.md`, Proposition 4.3: already use the Smith factors of `Q^j−I` to count fixed points on modular boards. Fixed-point counting is inherited standard algebra, not a new Session109 claim.
* `01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`: scoped eleven-cell dictionary, `R_10=O/(11)≅F_11×F_11`, and `R_18=O/(76)`. Its current corrections, including the factor-2 twist and finite-exact bijectivity boundary, remain in force.
* `05-knowledge/results/golden_signed_relation_geometry_20261019.md`, §§3–4: explicit common quotients and the hostile distinction between `Z/121` and the 121-element split ring; inverse golden time at eleven is already known.
* `05-knowledge/results/entry_20260927_recursive.md`: original source threshold and exact affine carry must survive recursive composition. Its canonical hostile is `27→41→31`: a child descends without the original source descending.
* `05-knowledge/results/golden_primes_carries_route_compiler_20261003.md`, §6: exact timed golden prefix composition, retaining legal word, endpoint and elapsed time.
* `05-knowledge/results/cycle101_tournament_collatz_dictionary_20261029.md`: current marked-point attachment law; identical nominal clocks need not give identical denominator cancellation. No claim here replaces that source-specific arithmetic.

Targeted searches of golden/carry/Fibonacci files and canon found the cited special Smith calculations and fixed-point formula, but not this complete stride classification plus affine-order/return-carry package. This is a local retrieval report, **not** proof of mathematical or repository novelty.

Board: linear time / translation gauge / return carry / golden module / marked source / exact legal word. Anchor: distinguish actual affine clocks; niche: invariant factors; wildcard: discriminant-five observer loss. Corrected near miss: signed determinant versus quotient cardinality. Least-used sidecar: the accumulated translation at the first linear return.

## 2. Exact Smith form for every positive stride

Let `F_0=0,F_1=1,F_(j+1)=F_j+F_(j−1)` and `L_j=F_(j−1)+F_(j+1)` for j≥1. Put

`D_r=L_r−1−(−1)^r > 0`.

**Theorem A.** For r≥1 the Smith normal form of `I−Q^r` is `diag(d_r,e_r)`, where `e_r=D_r/d_r` and

* r odd: `d_r=2` if `3|r`, and `d_r=1` otherwise;
* r=2m with m odd: `d_r=L_m`, `e_r=L_m`;
* r=2m with m even: `d_r=F_m`, `e_r=5F_m`.

**Proof.** The Fibonacci matrix identity gives

`Q^r=[[F_(r+1),F_r],[F_r,F_(r−1)]]`,

with determinant `(-1)^r` and trace `L_r`. Hence `det(I−Q^r)=1+(-1)^r−L_r=−D_r`. For a nonsingular integer 2×2 matrix, its first Smith factor is the gcd of all entries and its product of factors is the absolute determinant.

For odd r, if d divides all entries of `Q^r−I`, reducing determinants modulo d shows `−1≡1 (mod d)`, hence d divides 2. Modulo 2, Q has exact order 3, so d=2 precisely when 3 divides r.

For r=2m, Cayley–Hamilton for `Q^m` gives `Q^m+(-1)^m Q^(−m)=L_m I`. If m is odd,

`Q^(2m)−I = L_m Q^m`.

Since Q is unimodular, both Smith factors are L_m. If m is even, the Fibonacci identities give

`Q^m−Q^(−m)=F_m(2Q−I)`,
`Q^(2m)−I=F_m Q^m(2Q−I)`.

The matrix `2Q−I=[[1,2],[2,−1]]` has Smith factors 1 and 5. Unimodular multiplication preserves Smith factors, proving the last case. ∎

Examples:

| r | Smith factors | additive obstruction group |
|---:|---:|---|
| 1,2 | 1,1 | trivial |
| 3 | 2,2 | `(Z/2)²` |
| 4 | 1,5 | `Z/5` |
| 5 | 1,11 | `Z/11` |
| 6 | 4,4 | `(Z/4)²` |
| 7 | 1,29 | `Z/29` |
| 8 | 3,15 | `Z/3 ⊕ Z/15` |
| 9 | 2,38 | `Z/2 ⊕ Z/38` |
| 10 | 11,11 | `(Z/11)²` |
| 12 | 8,40 | `Z/8 ⊕ Z/40` |
| 18 | 76,76 | `(Z/76)²` |

**Golden identification.** In the ordered basis `(φ,1)` of `O=Z[φ]`, multiplication by φ is Q, since `φ²=φ+1`. Therefore `C_r` is canonically the underlying additive group of `O/(1−φ^r)`. This gives the exact map to the older golden rings, not an analogy by cardinality. It also distinguishes the 76-element group at stride9 from the `(Z/76)²` group at stride18.

## 3. Translation conjugacy and its modular obstruction

Let `A=Q^r`, `M_N=(Z/N)²`, and `f_b(x)=Ax+b`. Conjugation by `T_c(x)=x+c` gives

`T_c^(−1) f_b T_c(x)=Ax+b+(A−I)c`.

**Theorem B.** Two carries b,b' are conjugate by a translation iff their difference is in `(I−A)M_N`. In particular b is removable iff f_b has a fixed point. The set of translation-conjugacy classes is

`M_N/(I−A)M_N ≅ Z/gcd(N,d_r) ⊕ Z/gcd(N,e_r)`.

Consequently the number of classes and also the number of fixed points of A is

`gcd(N,d_r) gcd(N,e_r)`.

A carry class is removable iff its image is zero in this quotient; **the invariant factor sizes alone do not decide a particular b**, unless a Smith change-of-basis or direct image test is supplied. All carries are removable exactly when `gcd(N,D_r)=1`. For r=1,2 this holds over every N and the integer matrix `I−Q^r` itself is unimodular, so removal works over Z as well.

Proof: the conjugation identity proves the first statement. Reduce a unimodular Smith decomposition modulo N. Multiplication by integer d on Z/N has kernel and cokernel both of order gcd(d,N), with cyclic cokernel of that order. ∎

## 4. The exact affine-order law and CRT

Let `t=ord_N(A)` be the order of the linear permutation and define

`S_k=I+A+...+A^(k−1)`, `c_b=S_t b` in M_N.

**Theorem C.** The order of f_b as a permutation of the **whole** M_N is

`ord(f_b)=t · ord_add(c_b)=t·N/gcd(N,(c_b)_1,(c_b)_2)`.

In particular `t | ord(f_b) | tN`. If b is removable, `c_b=0` and the affine order is t. These are global permutation orders, not a claim that every point has that orbit length.

**Proof.** Induction gives `f_b^k(x)=A^k x+S_k b`. An identity affine map must have linear part I, so its exponent is a multiple of t. But `f_b^t=T_(c_b)`, hence `f_b^(jt)=T_(j c_b)`. The least positive j killing c_b is its additive order. The coordinate formula follows from elementary orders in Z/N. ∎

The return carry is translation-gauge invariant: replacing b by `b+(A−I)c` changes it by `S_t(A−I)c=(A^t−I)c=0`. It is thus a homomorphism from the obstruction quotient to `ker(A−I)`. Do **not** identify the return-carry image with the entire obstruction quotient without another argument.

The linear stride order is `t=π(N)/gcd(π(N),r)`, directly from the order-of-a-power rule. This alone does not compute the affine order; c_b is additional data.

**CRT.** For `N=∏p^a`, the state spaces, affine maps, fixed-point tests and obstruction quotients decompose into their prime-power factors. The global affine order is the lcm of the local affine orders, not their product. One may compute c_b at the global linear return t, or use each local linear return and then take the lcm. Both follow from the product action.

## 5. Hostile examples: same linear clock, different affine clocks

Take r=5,N=11. Then

`A=[[8,5],[5,3]] (mod11)`, `A²=I`, `A≠I`.

The vector b=(1,3) is fixed by A. Hence `S_2 b=2b≠0` has additive order 11: the affine map has order **22**. For b=0 the order is **2**. Both maps have exactly the same linear matrix and linear order. Equivalently b=(1,0) has return carry `(9,5)`, also of additive order11. The full census has 11 carries of affine order2 and 110 carries of affine order22.

The 11-element removable set is the image of I−A, while the obstruction quotient is Z/11. Thus this hostile is not merely a difference in an observer; the actual affine permutations have different global orders.

**Hostile to the converse: a nonremovable carry need not enlarge the order.** At r=4,N=5,

`A=Q^4=[[0,3],[3,2]]=I+J`, `J²=0`, `J≠0`.

Thus A has order5 and `S_5=5I+10J=0 (mod5)`. Every one of the25 carries has affine order5, although the translation obstruction quotient is Z/5. Explicitly `(I−A)M_5` is the span of `(1,2)`, so b=(1,0) is not removable. The return-carry map from the obstruction quotient is zero here. A full clock can therefore forget a genuine affine conjugacy distinction.

The checker also includes N=121 at stride5, stride6 modulo4, stride9 modulo19, and stride4 modulo5. These are controls at a lift, a scalar identity, a nontrivial split prime and discriminant5; consult `results.json` for exact distributions rather than extrapolating one prime's behavior.

## 6. Relation to a source-preserving Collatz ledger — and its limit

A finite sequence of affine maps `x↦A_jx+b_j` has exact state pair `(P,c)` with endpoint `Px+c`; appending `(B,d)` updates

`(P,c) ↦ (BP, Bc+d)`.

This is the same elementary semidirect composition law used above. A fixed matrix repeated t times produces c=S_t b. It is also the algebra behind the older Collatz affine ledger: composing odd branches `(3x+1)/2^a` retains total slope and numerator carry. The common source is **affine composition**, not an identification of Fibonacci evolution with Collatz.

The transfer specification is explicit:

* source: a legal, marked finite Collatz word with its exact arithmetic endpoint;
* map: evaluate its composed affine pair (and, only when denominators are invertible, reduce modulo N);
* target: a finite affine residue transformation;
* preserved: the endpoint congruence and composition order;
* destroyed by reduction: source magnitude, real positivity, exact valuation guards, intermediate word, and source-relative descent;
* required sidecar for a source certificate: original integer, legal ordered word/valuation guards, exact endpoint, elapsed time and source threshold;
* cheapest hostile: `27→41→31`, where descent below41 is not descent below27.

Thus an affine order theorem supplies a **finite congruence clock**, not a termination proof. A polynomial golden normalization carry, a modular affine return carry, and a Collatz numerator carry are three typed records with a shared algebraic theme; none should be silently substituted for another.

## 7. Relation to the trace observer's five-fold defect

The factor5 in the r divisible-by4 Smith formula comes from the same matrix `2Q−I` whose determinant is −5, the discriminant of `X²−X−1`. This is a precise algebraic meeting with the companion trace-pairing calculation. But the objects differ: a trace-pairing kernel forgets linear-state coordinates, whereas `coker(I−Q^r)` measures failure to eliminate an affine translation. The latter also has prime11 and29 obstructions and 2-primary examples, so a common discriminant does not identify the two losses or their sidecars.

## 8. Reproduction and finite scope

Run `python3 check.py` or `python3 -O check.py` from this directory. The script uses only Python's standard library and explicit checks that remain active under optimization:

1. all Smith formulas for r=1..200, comparing closed forms with direct entry gcds and determinants;
2. direct modular image/fixed-point counts for N=2..15 and r=1..12;
3. for every carry b in those boards, brute affine iteration equals the return-carry formula;
4. removable carries have the linear order;
5. CRT order identities for four coprime factor pairs and four carries at strides1..12;
6. the named hostile/lift distributions by exact modular arithmetic.

These computations validate finite instances. The preceding proofs establish the unrestricted theorems. No automated source/legality reconstruction, universal Collatz coverage or new Pisano prime-power lifting theorem is claimed.

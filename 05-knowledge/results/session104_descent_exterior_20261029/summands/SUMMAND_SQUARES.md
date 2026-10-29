# Summand recovery and consecutive block balances

Session 104, 2026-10-29. Status: historical definitions and the prior three-chain discussion recovered; elementary identities proved below. No abc result is claimed.

See [THREE_CHAINS_RECOVERED.md](THREE_CHAINS_RECOVERED.md) for the preserved 2026-09-22 user framing, the three path components of Q_12, and their distinction from closure holes.

## 1. Source recovery table

Paths refer to `eliottcassidy2000/math`.

| Source | Recovered statement | Boundary |
|---|---|---|
| `07-reflections/summand-graph-fermat-zeckendorf.md`, Setup | Positive distinct parent pairs `{a,b}->a+b`, `a<b`; sources 1,2 | Its CURRENT-TRUTH WARNING retracts Fibonacci-reset and restoration claims; later speculative tournament analogies are not inherited |
| `01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`, §§1–3 | Pair-swap fixed loci are doubles / squares; strict closure of {2,3} misses exactly {1,4,6}; eventual synchronous frontier is dyadic | GitHub live read succeeded; blob `1daae13394e4b772a78068c3bb7aae65f3d7002e` |
| `05-knowledge/results/arithmetic_braids_20260917_summand.md`, §§1–2 | Strict additive shadow has arcs `x<z, z!=2x`; removed arcs are doubling forest; multiplicative analogue removes squaring arcs | GitHub live read succeeded; blob `5d180a5eae4ada2db7c358f2bd645941d41ca7d0` |
| `01-canon/theorems/THM-361-product-sum-defect-normal-form.md` | `2+2=2*2`; distinct positive nontrivial product-sum tuple is {1,2,3} | Distinctness and positivity are essential |
| `05-knowledge/results/arithmetic_braids_20260917_geometry.md`, §5 | Gaussian squaring preserves primitive Pythagorean triples while hypotenuse prime support is fixed and odd legs gain new primes | This is not an abc bound |
| User's “3 piece partition ... built one node at a time” | Recovered preserved 2026-09-22 framing plus audited Q_n three-chain result | See THREE_CHAINS_RECOVERED.md; distinguish actual Q_12 components with founders {1,2,4} from closure holes {1,4,6} |

The recovery inspected the local tracked filename inventory and targeted theorem/reflection/result files, with live GitHub source reads recorded in `source_recovery.json`; it did not fetch whole historical logs.

## 2. Exact object and one-node growth

For each positive target n retain the hyperedges `{a,n-a}->n` for `1<=a<n-a`. There are `r(n)=floor((n-1)/2)` such pairs. Its simple directed shadow has `x->n` iff `1<=x<n` and `n!=2x`; hence the incoming **arc** count is `2r(n)`, not `r(n)`. Pair closure requires both parents and is not ordinary directed reachability.

Adjoining integer n adds exactly this target fibre, with no changes to earlier fibres. `r(n+1)-r(n)` is 1 for even n and 0 for odd n. The total pair count through N is `floor((N-1)^2/4)`. These one-at-a-time growth facts describe the unfiltered parent fibres; the recovered three-chain partition belongs to the square-filtered companion graph Q_n.

The ordered weak fibre splits into `a<b`, `a=b`, `a>b`; swap exchanges the outer pieces and fixes the middle. This is a separate mathematical three-part construction, not the recovered Q_12 component partition. The diagonal exists at n=2a in addition and n=a² in multiplication. Thus doubles and squares are genuinely the same fixed-point mechanism for two different operations, not numerology.

The numbers 3,4,6 have distinct rigorous roles: 3 is the first distinct positive binary sum; 4 is the binary additive/multiplicative collision `2+2=2*2`; 6 is the distinct ternary collision `1+2+3=1*2*3`. Do not confuse this packet with the missing set {1,4,6}. From seeds {2,3}, 5,7,8 are generated and every n>=9 follows from `2+(n-2)`; 1,4,6 cannot be generated. The first three-entry threshold does not obstruct binary generation of 6 from other seeds.

## 3. Complete solution of the displayed block patterns

Fix k>=1; let k+1 consecutive terms lie left of k consecutive terms, with no gap. Define

`D_p(a;k)=sum_{i=0}^k (a+i)^p - sum_{i=k+1}^{2k} (a+i)^p`.

Direct summation gives

`D_0=1`,
`D_1=a-k²`,
`D_2=a²-2k²a-k²(2k+1)=(a+k)(a-k(2k+1))`.

Proof: the coefficient of a in D1 is `(k+1)-k=1`; its constant is `2*k(k+1)/2 - (2k)(2k+1)/2=-k²`. Expanding squares yields the same coefficient -2k² of a, and the constant is `2*k(k+1)(2k+1)/6-(2k)(2k+1)(4k+1)/6=-k²(2k+1)`.

Therefore the unique positive first term for first powers is a=k². Multiplying by 2 proves

`sum_{i=0}^k 2(k²+i) = sum_{i=k+1}^{2k} 2(k²+i)`.

For squares, the unique positive first base is a=k(2k+1)=T_(2k), where T_m=m(m+1)/2:

`sum_{i=0}^k [k(2k+1)+i]² = sum_{i=k+1}^{2k} [k(2k+1)+i]²`.

The other root a=-k is the signed trivial symmetry `(-k)²+...+0²=1²+...+k²`, not a second positive solution.

| k | Even-number balance | Consecutive-square balance |
|---|---|---|
| 1 | 2+4=6 | 3²+4²=5² |
| 2 | 8+10+12=14+16 | 10²+11²+12²=13²+14² |
| 3 | 18+20+22+24=26+28+30 | 21²+22²+23²+24²=25²+26²+27² |
| 4 | 32+34+36+38+40=42+44+46+48 | 36²+37²+38²+39²+40²=41²+42²+43²+44² |

More generally an arithmetic progression of positive step d has `D_1=a-dk²` and `D_2=(a+dk)(a-dk(2k+1))`. Thus these solutions are forced, not an empirical continuation.

### Finite differences: what actually transports

For translation of the SAME starting base a and same k,

`D_2(a+1;k)-D_2(a;k)=2D_1(a;k)+1`.

The exact centered version is `D_2(a+1/2;k)-D_2(a-1/2;k)=2D_1(a;k)`.

So square increments supply affine doubles: `(x+1)²-x²=2x+1`, or centered differences supply `2x`. The square-balanced base T_(2k) does NOT coincide with the linear-balanced base k²: their difference is k(k+1). Differentiating a balanced equation while moving its root would not prove the even identity. The defect polynomial, not a pointwise equality, is the legitimate bridge.

### Pythagorean scope

For k=1, consecutive bases solve `a²+(a+1)²=(a+2)²`, hence a=3 uniquely among positive integers. It is the Euclidean triple `(m²-n²,2mn,m²+n²)` for (m,n)=(2,1). With arbitrary step d, the family is exactly `(3d,4d,5d)`; primitive only at d=1. For k>=2 the equality has multiple squares on both sides, so it is not itself a Pythagorean triple. Grouping terms into norms produces equality of norms in different dimensions, not necessarily integer side lengths.

## 4. abc: genuine arithmetic and a firm boundary

For coprime positive A+B=C, abc concerns prime support `rad(ABC)=product_{p|ABC}p`, versus C, with the quantified conjectural bound `C <= K_epsilon rad(ABC)^(1+epsilon)` for every epsilon>0 and a corresponding constant. No proof of this bound is supplied by a swap-fixed-point count.

`3²+4²=5²` gives (A,B,C)=(9,16,25) and radical 30. `2+4=6` is not coprime; dividing by 2 yields (1,2,3), radical 6. Equality between two sums with more than three terms is not automatically an abc triple. Any regrouping must state three integers and remove their common gcd before evaluating radical.

Powers preserve prime support (`rad(n²)=rad(n)`), while doubling adds at most prime 2 (`rad(2n)=rad(n)` for n even and `2rad(n)` for n odd). This is the actual support-growth contrast relevant to addition/multiplication. It supplies no universal lower bound on prime support of sums.

The recovered Gaussian map `(A,B,C)->(A²-B²,2AB,C²)` is a sharper lawful bridge: it is complex squaring of A+iB. A primitive triple with A odd, B even stays primitive. Hypotenuse support stays fixed, but successive odd legs are pairwise coprime: every earlier A_i divides later B_j, and gcd(A_j,B_j)=1. Thus projected hypotenuse squaring alone loses the coordinates where new primes enter. This precisely explains why the square-forest shadow cannot settle abc.

## Verification

`python3 check_summands.py` and `python3 -O check_summands.py` use explicit exceptions, exact integer/rational arithmetic, positive/negative root controls, illegal-input controls, independent pair enumeration, closure controls, and Gaussian primitivity/prime-support checks. Finite checks support, but do not replace, the all-k algebraic proofs above.

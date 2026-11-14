# Absolute overlap of the two selective charts

Session 325. **PROVED; independent review PASS.** The independent report is `../s325_cut_audit/INDEPENDENT_REVIEW.md`. No claim of full-spectrum coverage. This is an exact common-target construction and a complete overlap decoder for fixed padded rails, not a congruence alignment or an extrapolation from a search.

## 1. A primitive common interface, with two different realizers

Use the strong quotient tournaments from S321, also retained in S324:

```
R13 = (24,9,35,4,14,27)
R15 = (108,97,2,6,79,92,12).
```
A source fills vertex zero, a reservoir fills the last vertex, and all other vertices stay singleton. For any nonempty source U and reservoir V, their exact counts are

```
N_gamma(U,V) = A_gamma(U) Q_V(2) - H(U)H(V),
A_gamma(U) = -3H(U) + 6Q_U(2) + gamma Q_U(3),
gamma = 2 for R13, 6 for R15.
```
Here Q_T(k) counts ordered k path slots covering T, allowing empty slots. This is the proved S321 itinerary identity. For completeness the standalone replay independently recomputes all quotient itinerary coefficients, not just substitutions into this formula.

**Theorem 1 (universal absolute overlap).** Let

```
U = (478,412,448,100,460,407,290,72,136),
W = (138,212,105,66,45,139,49,92).
```
These are strong tournaments, respectively of orders nine and eight, with exact profiles

| source | H | Q(2) | Q(3) |
|---|---:|---:|---:|
| U | 571 | 14336 | 156129 |
| W | 571 | 8032 | 58347 |

Their two different source charts have the same state:

```
A_2(U) = -1713 + 86016 + 312258 = 396561,
A_6(W) = -1713 + 48192 + 350082 = 396561.
```
Consequently, for **every nonempty tournament V**, both constructions have exactly

```
N = 396561 Q_V(2) - 571 H(V).                    (1)
```
Both are strong, have order |V|+13, and are nonisomorphic. The coefficient and intercept are coprime: gcd(396561,571)=1.

The profile table is a finite certificate, verified directly by subset path DP and independent slot convolutions on the displayed graphs. The all-order assertion is then the layout identity applied to arbitrary V, not a fit to the finite examples.

**Nonisomorphism for arbitrary V.** In R13[U,V], a source vertex has its internal outdegree plus two, hence at least three. The four exterior singleton outdegrees are 10, 10+|V|, 1 and 3. A reservoir vertex has its internal outdegree plus twelve. Thus there is no vertex of outdegree two. In R15[W,V], quotient vertex three has outdegree exactly two, independently of V. No isomorphism is possible. This statement does not assert nonisomorphism up to taking converse.

### The itinerary proof behind the interface

Let F_i(T) count ordered nonempty i-path covers. A full Hamiltonian path compresses to a legal quotient word in which each exterior singleton occurs once, and U and V occur i and j times. Filling these visits by the ordered covers is bijective. The nonzero itinerary tables (rows i, columns j) are

```
R13: (29,15), (24,12), (4,2)
R15: (53,27), (48,24), (12,6).
```
There are no further rows or columns. Since the two growing quotient vertices induce an acyclic pair, between successive visits to either mark there must be an exterior singleton. Thus each mark has at most one more visit than the number of exterior singletons; checking the finite table up to this bound is exhaustive. The recurrence in the replay counts every legal word in that range.

Use Q_T(2)=2H(T)+F_2(T) and Q_T(3)=3H(T)+3F_2(T)+F_3(T). Expanding the two tables gives the displayed N_gamma identity exactly. Strongness follows from the strong quotient and nonempty modules; no strongness assumption on V is needed.

## 2. Prescribed-target consequences, including both atom lanes

With V=TT_m, m>=1, (1) becomes the common rail

```
N_m = 396561 * 2^m - 571.                         (2)
```
Given an absolute target N, it is on this rail **if and only if** (N+571)/396561 is an integral power of two at least two. The test supplies m and the two literal strong witnesses. At m=1, N=792551 is prime. At m=6,

```
N = 25379333 = 7 * 3625619,
```
where 3625619 is prime. Trial division verifies these finite primality claims. No infinitude of primes or of prime cofactors on this rail is asserted. Since 396561 and 571 are both 4 modulo 7, precisely m divisible by three gives divisibility by 7.

For any dyadic reservoir Q_V(2)=2^m, the absolute target first forces the unique possible shell m, the first integer with 396561*2^m>N, and then forces

```
D = (396561*2^m-N)/571 = H(V).
```
Require m>=1 and 0<D<=2^(m-1), with D integral and odd. These are necessary conditions, not sufficient for an unrestricted physical reservoir. The S324 factor-free reservoir compiler may then decide the chosen one-apex/finite-bank grammar. Every acceptance now supplies **both** charts at the same absolute target. Nothing is inferred from a residue solution whose slot cost has not been paid.

A second independently checked pair, supplied as a smaller-H alternative, is

```
U401 = (320,117,473,211,1,477,144,275,90),
W401 = (130,28,193,53,5,151,59,90),
(H,Q2,Q3)(U401) = (401,10960,126435),
(H,Q2,Q3)(W401) = (401,6080,47025).
```
It has common primitive state (h,A)=(401,317427), hence common targets 317427 Q_V(2)-401H(V). This is an additional example, not a minimality claim.

## 3. Exact absolute intersections of fixed padded rails

Fix any sources and any nonempty reservoir cores in the two selective charts. Append t>=0 singleton SCCs to the reservoir. Because H is unchanged and Q(2) is doubled by each singleton, every such fixed chart has the rail

```
C*2^t-d,   C=A_gamma(U)Q_V(2)>0,   d=H(U)H(V)>0.
```
The following arithmetic theorem therefore classifies all absolute overlaps for a finite menu of these charts.

**Theorem 2 (two coincidences or a common tail).** Let C,E be positive integers and d,f integers. Consider C2^t-d and E2^u-f for t,u>=0. Write C=a2^r and E=b2^s with a,b odd.

* If d=f and a=b, their intersection is exactly the common tail a2^k-d for k>=max(r,s).
* If d=f and a!=b, the rails are disjoint.
* If d!=f, the rails have at most **two** common absolute values. All are obtained by the finite decoder below.

Thus a finite pair of menus has infinite intersection **if and only if** at least one pair has equal intercepts and equal odd slope parts. Everything else contributes at most two values per pair. The theorem allows asynchronous and arbitrarily large lengths; no searched cutoff is involved.

**Proof and complete decoder.** Put x=t+r, y=u+s and delta=d-f. We solve

```
a2^x-b2^y=delta,   x>=r, y>=s.
```
If delta=0, oddness forces a=b and x=y, giving the first two cases.

Suppose delta!=0. The diagonal x=y gives at most one candidate: delta/(a-b) must be a positive integral power of two; if a=b there is none. Retain it only when the exponent meets both lower bounds.

For x!=y, let v=v_2(|delta|), D=delta/2^v. Taking valuations forces min(x,y)=v. If x<y then

```
D=a-b2^(y-v)<a-b;
```
if y<x then

```
D=a2^(x-v)-b>a-b.
```
The two off-diagonal orientations are mutually exclusive. In the first case test whether (a-D)/b is an integral power of two greater than one, giving x=v and y=v+log2((a-D)/b). In the second test whether (D+b)/a is such a power, giving y=v and x=v+log2((D+b)/a). Enforce the lower bounds. This gives at most one off-diagonal candidate and at most one diagonal candidate, and every solution has been included.

The arithmetic bound two is sharp: 3*2^x-5 and 2^y-1 share values 1 and 7 (at (x,y)=(1,1),(2,3)). This sharpness example is arithmetic, not a claim that those coefficients themselves are physical selective-chart states.

For an individual target, rail membership is also exact without searching lengths: N+d must be C times a nonnegative power of two. The theorem explains the distinction between isolated numerical coincidences and a genuine repeatable common interface. The source pairs in section 1 realize the latter physically.

## 4. The overlap is stable on one side, not both

For fixed source states (h,A) and (h',A'), equality on singleton and two-vertex transitive reservoirs already forces equality on **every** reservoir: the first two values are 2A-h and 4A-h, so they recover A and h. This is a two-measurement characterization, not permission to infer profile equality from one common target.

The source pair above has Deltaq=14336-8032=6304 and r_U-3r_W=-3Deltaq. Append the **same nonempty ordinal source block Z** to both U and W. Multiplicativity of H,Q2,Q3 gives the exact state defect

```
A_2(U -> Z)-A_6(W -> Z)
    = 6*6304*(Q_Z(2)-Q_Z(3)) < 0.                 (3)
```
Here Q_Z(3)-Q_Z(2)=H(Z)+2F_2(Z)+F_3(Z)>0. Their H values remain equal, so for every nonempty reservoir the target difference is (3) times Q_V(2), strictly negative. Common neutral growth on the source side destroys this overlap at its first step, whereas arbitrary common reservoir changes preserve it exactly.

More generally a source state collision has this formula with its own Deltaq. It survives every common ordinal source extension exactly when Deltaq=0 (then also r_U=3r_W). The displayed pair does not satisfy that stronger state requirement.

There is no repair by merely choosing different neutral reservoir lengths after this same source extension. Fix any common nonempty reservoir core V0, and pad it by t and u singletons on the two sides. The new intercepts still agree. Both new A coefficients are odd, by tournament Hamilton-path parity, but they are unequal by (3). Multiplication by the common Q_V0(2) therefore leaves different odd slope parts. Theorem 2 says these two padded rails are **disjoint**, even asynchronously. In particular, after any common nonempty ordinal source extension, no choice of two transitive reservoir lengths recovers a common target. A repair must change the source operation, use different reservoir cores, or retain more state. This isolates a structural obstruction instead of claiming an H-only transferable operation.

## Scope and dependencies

The primitive common family is new finite data promoted to an all-order theorem by the exact quotient identity. The rail-overlap theorem is elementary and exhaustive for **fixed** sources and reservoir cores with singleton padding. It does not classify intersections when the source/core bank itself grows with the target. The S324 sparsity bounds still apply to the stated finite-bank/one-apex grammars; identifying a common family does not turn their union into an all-target construction. Chan–Pak remains open.

Dependencies: the physical quotient and its S321 identity as restated in `s324/COMMUTATIVE_INVERSE.md`; ordinal slot multiplicativity; elementary Hamilton-path subset DP. The replay independently regenerates the entire finite itinerary table, all displayed source profiles, the rail decoder on exhaustive bounded integer inputs, finite prime/cofactor certificates, and literal strong quotient paths. No random-search completeness claim is used. Exploratory files are separate and are not certificate dependencies.

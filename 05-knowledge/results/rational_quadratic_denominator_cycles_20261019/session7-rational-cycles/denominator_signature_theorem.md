# Rational quadratic cycles: denominator signatures and a failed local lift

Date: 2026-10-19. Status: **PROVED elementary reductions + CITED period bound + FINITE-EXACT denominator census**. The Flynn–Poonen–Schaefer conjecture remains **OPEN**. No priority claim is made for the elementary local coding or S-unit ideas. The new deliverable here is their explicit three-edge compression, complete denominator-first implementation, and next-level analysis of an inherited finite-ring hostile.

## 1. What is established

Let `f_c(x)=x²+c`, with `c` rational.

1. **PROVED reduction.** For an exact rational cycle of length `3≤n≤N`, write `c=C/b²` in lowest terms. At each bad prime, the valuations of the three pair sums of three consecutive numerators belong to a universal set of at most `4N−6` signatures. They reconstruct a finite list of possible parameters `C`, with no numerator/height cutoff.
2. **FINITE-EXACT, using the cited universal period bound.** For every rational `c` whose reduced denominator is at most `1,000,000`, there is no rational cycle of length greater than three. This is an unbounded-numerator statement. Nonsquare denominators admit no rational periodic points; for square denominators `b²`, it uses the published infinite family with at most two denominator primes and our exhaustive checks of the remaining 101 denominators `b≤1000`. Among those 101 denominators, exactly 62 cycles of length at least three occur; all have length three.
3. **PROVED local obstruction and all-depth classification.** The 21 exact six-cycles of `(a²−29)/4` modulo `43²` have no same-period lifts modulo `43³`; all 5,418 lifted points instead have exact period `258=6·43`. Every deeper residue shell has an explicit growing period, while the full three-disk tube over `Q_43` contains only the original three-cycle as periodic points.
4. **PROVED algorithmic corollary, given factorization.** For arbitrary `b` with its prime factorization supplied, all rational cycles of length at least three with exact point denominator `b`, across all parameter numerators, can be enumerated in time `F(ω(b))·poly(log b)`. The signed seed count is at most `4(2N−1)^ω(b)`, where `N=2^ω(b)+2`; an explicit escape guard controls integer bit lengths.

The finite result is not asserted to be a literature record. In particular, it is not a proof for unbounded denominator, and it should not be described as a new proof of the full conjecture.

## 2. Inheritance and primary-source boundary

Closest inherited canon:

- `THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md` gives the complete rational graph at `c=−29/16`, the unique AP-supported rational three-cycle, and 21 six-cycles modulo `43²`.
- `THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md` separates a linear order-six lift from the scalar quadratic dynamics.
- `arithmetic_seams_20260921_dynamics.md` explicitly rejects the inference from a six-vector symmetry to a rational six-cycle.
- The session-six coefficient-Frobenius result concerns a polynomial coefficient congruence. It is a different operation from lifting a periodic point; no transfer from its unique scalar lift to rational periodic-point existence is claimed.

The hostile is therefore already local and concrete. The missing sidecar is compatibility under passage to the next modulus, and then rationality across places. The corrected near miss is an unmarked or quotient cycle being read as a cycle for the original scalar map.

Primary sources, checked on 2026-10-19 in this session:

- [Eliahou–Fares, Some results on the Flynn–Poonen–Schaefer conjecture](https://arxiv.org/html/1903.08865), published in *Canadian Mathematical Bulletin* 65 (2022), 598–611. We use Proposition 2.1, Theorem 2.12, Corollary 3.11, and Theorem 3.16 in the arXiv v1 numbering: common denominator `b`, `4|b` for period at least three, period at most `2^ω(b)+2`, and the complete conjecture when `ω(b)≤2`.
- [Morton, Arithmetic properties of periodic points of quadratic maps II](https://matwbn.icm.edu.pl/ksiazki/aa/aa87/aa8721.pdf): period four is excluded.
- [Flynn–Poonen–Schaefer](https://arxiv.org/abs/math/9508211): period five is excluded.
- [Stoll, Rational 6-cycles under iteration of quadratic polynomials](https://arxiv.org/abs/0803.2836): the general period-six exclusion depends on analytic and Birch–Swinnerton-Dyer assumptions for the specified Jacobian. It is not used by our finite search.
- [Morton–Raianu, Arithmetic properties of 3-cycles](https://arxiv.org/abs/2105.07435): `−29/16` is the unique minimum-numerator/minimum-denominator parameter with a rational three-cycle.
- [Benedetto–Briend–Perdry, Dynamique des polynômes quadratiques sur les corps locaux](https://archive.numdam.org/articles/10.5802/jtnb.589/): local quadratic Julia dynamics and binary shifts are established literature, not a new general connection here.
- [Piipponen–Erkama, Rational cycles of quadratic polynomials](https://www.pphmj.com/abstract/8547.htm): denominator-prime configurations and algebraic compatibility equations also have precedent. Only the publisher abstract was inspected; no theorem from that paper is imported.

A 2024/2025 primary-paper check, [Han, Uniform boundedness on rational maps with automorphisms](https://arxiv.org/abs/2401.11309), still presents the general quadratic-polynomial bounds in this form. Targeted searches found no verified unconditional resolution of the full conjecture. This is a checked literature boundary, not a claim to have exhausted every 2026 preprint.

The live concept board is: (i) common-denominator arithmetic and the full-period binary word; (ii) signed S-unit edges and their lost additive consistency; (iii) return-map multiplier and higher-order displacement; (iv) real parabolic versus rational periodic points; and (v) finite certificates versus uniform denominator-support theorems. The signature reduction connects (i) to (ii); the resonant-lift obstruction connects (iii) to (v); neither erases the field distinction in (iv).

The landmarks must remain distinct. At `c=−7/4` the real parabolic three-cycle is algebraic, not rational; the rational periodic cycle is `−3/2 ↔ 1/2`. At `c=−29/16` the rational three-cycle is `−7/4 → 5/4 → −1/4 → −7/4`, with multiplier `35/8`, not one.

## 3. The numerator model and the pair-sum support lemma

Assume an exact rational `n`-cycle with `n≥3`. The denominator theorem gives

```
c=C/b², gcd(C,b)=1, b>0,
x_i=a_i/b, gcd(a_i,b)=1,
a_i²+C=b a_(i+1),                         i modulo n,
4|b.
```

All `a_i` are distinct. The recurrence is injective on the cycle, so `a_i+a_j` is nonzero for `i≠j`. Subtracting equations gives

```
b(a_(i+1)−a_(j+1))=(a_i−a_j)(a_i+a_j).       (1)
```

For any fixed nonzero lag `k`, multiply (1) over all `i`. The nonzero difference factors cancel, yielding

```
product_i (a_i+a_(i+k)) = b^n.               (2)
```

Each sum is an integer, so every prime factor of each pair sum divides `b`. This standard support lemma is rederived to specify exactly what is used. It loses signs, the additive relations between edges, and the orbit order if used only as an unordered set constraint.

## 4. Bad-prime itineraries determine valuations exactly

Write `e_p=v_p(b)`.

### Odd bad primes

For odd `p|b`, the congruence `a_i²≡−C (mod p)` places every numerator in one of two nonzero residue classes `r_p,−r_p`. Let `w_i` be its binary class. For two different cycle points, put

```
L_p(i,j)=min {t≥0 : w_(i+t)≠w_(j+t)}.
```

This minimum exists. Otherwise the two points have equal symbols forever; at each step their sum is a p-adic unit, and (1) decreases the valuation of their difference by `e_p`, eventually making an integer difference have negative valuation.

The same argument proves that the word `w_0...w_(n−1)` has exact cyclic period `n`: an invariant nontrivial shift would contradict the preceding argument. Iterating the valuation identity until the first different symbols, whose difference is a unit, gives

```
v_p(a_i−a_j)=e_p L_p(i,j).                  (3)
```

Define the gap statistic

```
q_p(i,j)=0                                  if w_i=w_j,
q_p(i,j)=1+L_p(i+1,j+1)                     otherwise.
```

Using (1) once more,

```
v_p(a_i+a_j)=e_p q_p(i,j).                  (4)
```

### The bad prime two

Here `e=e_2≥2` and all numerators are odd. Use their classes `1,3 (mod 4)` as the two symbols. Equal symbols imply `v_2(a_i+a_j)=1`; unequal symbols imply `v_2(a_i−a_j)=1`. Equal-symbol steps reduce the difference valuation by `e−1>0`, so the word again has exact period `n`, and

```
v_2(a_i−a_j)=1+(e−1)L_2(i,j),
v_2(a_i+a_j)=1+(e−1)q_2(i,j).              (5)
```

The offset one is essential. Removing it would incorrectly rule out the rational three-cycle at `−29/16`.

### Why the multiplicative cycle identities alone cannot win

Fix any primitive binary word of length `n` and any nonzero lag. Mark positions where the word and its lagged rotation differ. At a marked position, `q` is exactly the cyclic distance to the next marked position; at unmarked positions it is zero. The marked gaps partition the circle, hence

```
sum_i q(i,i+k)=n.                           (6)
```

Consequently the valuations of (2) are automatic for every primitive binary word: `e_p n` at odd primes and `n+(e_2−1)n=e_2n` at two. This is a precise no-go for an argument that retains only the bad-prime cycle-product valuations. It does not assert that arbitrary collections of local words glue to rational points.

In graph language, `v_p(a_i+a_j)−e_p` is the difference between the next and current difference-potentials. Its zero circulation is a coboundary identity. A zero circulation does not reconstruct globally compatible integer vertex labels.

## 5. Three consecutive points compress the local words

For a primitive binary cyclic word of length `n`, consider

```
(q(0,1), q(0,2), q(1,2)).
```

It lies in the following **superset**, for any `N≥n`:

```
T_N = {(0,0,0)}
    ∪ {(0,1,k),(0,k,1),(k,k−1,0) : 2≤k≤N−1}
    ∪ {(1,0,k) : 1≤k≤N−1}.                 (7)
```

Thus `|T_N|=4N−6`. This is deliberately not claimed to be the exact realizable table for any fixed `n`.

**Proof.** Up to complement, the first three bits have four types.

- `000` gives `(0,0,0)`.
- `001` gives `(0,1,k)` or `(0,k,1)`: the fourth bit agrees with exactly one of the two different preceding bits, so exactly one of the two nonzero waiting times is one.
- `010` gives `(1,0,k)`.
- `011` gives `(k,k−1,0)`. If the initial run beginning at position one ends at position `k`, the adjacent comparison first differs at `k`, while the two-step comparison first differs at `k−1`.

Every nonzero `q` is at most `n−1`. Two distinct rotations of the same binary word have the same number of ones, so they cannot differ at exactly one position. Starting from a disagreement, another disagreement therefore occurs within the next `n−1` positions. This proves the bound and (7). QED.

## 6. Signed edge gluing and an exact finite algorithm

Set

```
h=b/2,
r_2=2^(e_2−1),
r_p=p^e_p for odd p|b,
u=(a_0+a_1)/2,
v=(a_0+a_2)/2,
w=(a_1+a_2)/2.
```

By the support lemma and (4)–(5), these are nonzero integers and

```
(|u|,|v|,|w|)
  = (product_p r_p^q_p01,
     product_p r_p^q_p02,
     product_p r_p^q_p12),                 (8)
```

where each local exponent triple is in `T_N`. No unknown prime and no unbounded exponent remains for fixed `b,N`.

Signed sums reconstruct the vertices uniquely:

```
a_0=u+v−w,
a_1=u+w−v,
a_2=v+w−u.                                 (9)
```

The first two arrows have the same parameter precisely when

```
a_0²−b a_1=a_1²−b a_2,
```

which, by (9), is equivalent to the small signed S-unit equation

```
u(v−w)=h(u−v).                             (10)
```

For each local-signature choice and each pair of signs of `u,v`, compute

```
w=v−h+h v/u.
```

Retain it exactly when `w` is one of the two signed magnitudes prescribed by (8). Then enforce distinct vertices and `gcd(a_i,b)=1`; put `C=b a_1−a_0²`, require `gcd(C,b)=1`, and iterate `(a²+C)/b` using exact integers. Discard a candidate if divisibility or the reduced denominator fails. A return to `a_0` after `3..N` distinct vertices certifies a cycle; a tail entering another cycle is not accepted.

**Completeness theorem.** Every rational cycle of length `3..N` with denominator `b` occurs in this algorithm. Indeed, its local words give signatures in (7), its signed sums occur in (8), its vertices reconstruct by (9), and its recurrence passes every exact test. Conversely every accepted output explicitly satisfies the reduced-denominator recurrence and has the certified exact period.

This is a finite decision procedure after using the cited bound `N=2^ω(b)+2`. It never loops through values of `C` and never imposes a bound on `|C|`.

### Two-edge refinement found during the independent audit

The third edge does not need to be enumerated. Projecting (7) gives the simpler proved superset

```
P_N = {(0,k):0≤k≤N−1} ∪ {(k,k−1):1≤k≤N−1},
|P_N|=2N−1.                                (12)
```

Enumerate only the magnitudes of `u,v` from one pair in `P_N` per bad prime, choose their four sign combinations, and derive `w=v−h+h v/u`. The local inequality `q01≤q02+1` shows that `u` divides `h v`, so the derived `w` is always integral. One may reject `w` unless it is nonzero, has support in the bad primes, and has each prime exponent at most the corresponding scale times `N−1`. These are necessary filters, not a replacement for checking exact return.

Thus at most `4(2N−1)^ω(b)` signed vertex triples are needed. The independent audit uses actual binary-word projections to remove one impossible pair at `N=18`, leaving 34 local pairs. Its four-prime cases examine `34⁴=1,336,336` unsigned pairs each, compared with `66⁴=18,974,736` in the primary three-edge certificate. The two algorithms deliberately retain different candidate supersets before exact return.

A hostile to stopping at (10) is `b=4`, `(a_0,a_1,a_2)=(-1,3,5)`, `C=11`. Its signed half-sums `(u,v,w)=(1,2,4)` have the legitimate signature `(0,1,2)` and satisfy (10), but the orbit continues `5→9→23→...` in numerator coordinates. The return condition is indispensable.

### Factorization-given fixed-parameter enumeration theorem

Let `s=ω(b)` and suppose the prime factorization of `b` is supplied. The complete enumeration of all rational cycles of length at least three with exact point denominator `b`, allowing every integer parameter numerator `C`, is fixed-parameter tractable in `s`: its bit complexity is bounded by `F(s)·poly(log b)` for a computable function `F`. Factoring `b` is not included in this assertion. If `4` does not divide `b`, the answer is empty by the cited denominator theorem.

For `4|b`, use `N=2^s+2` and the universal pair table (12). There are at most `4(2N−1)^s` signed seeds. Since `h=b/2` is the product of the prime bases, `|u|,|v|≤h^(N−1)`, and

```
|w|=|v−h+h v/u|≤h^(N−1)+h+h^N.
```

Consequently the reconstructed vertices and `C=b a_1−a_0²` have `O(N log b)` bits. All precomputation, sign choices, gcd tests, and bounded smooth-support checks use integers of this size.

To keep subsequent orbit arithmetic within the same bit bound, set

```
R=b+floor(sqrt(|C|))+1
```

and reject immediately when `|a|>R`. This is safe: then `|a|>b+sqrt(|C|)`, so

```
(a²+C)/b ≥ (|a|²−|C|)/b > |a|.
```

Every later iterate is positive, strictly increases, and remains beyond the same threshold, so no periodic return is possible. Until rejection, `|a|≤R`; computing a single next value still uses `O(N log b)` bits. At most `N` recurrence steps are needed. Standard deterministic integer arithmetic therefore gives the claimed `F(s)·poly(log b)` bound, with an exponent on `log b` independent of `s`. This is an all-denominator enumeration theorem, not an all-denominator exclusion theorem.

The restriction to periods at least three is essential to the enumeration statement: periods one and two can occur for infinitely many parameter numerators at a fixed denominator. The companion `rational_cycle_pair_probe.py` implements the universal pair table and escape guard. It uses the same mathematical reduction as the primary search, so its replay is an implementation control; the separately written audit remains the independent certificate.

This implementation reproduces all 101 denominators and all 62 cycles of the primary census. As a large-exponent control, the supplied factorization `b=4·3^20·5^20=1,330,102,692,031,860,351,562,500` still requires only `19³=6,859` unsigned pairs, or `27,436` signed seeds. It produces 16 primitive initial triples, six escape rejections, and no cycle of length at least three. Normal and optimized Python executions agree byte for byte on this control. This single denominator is a control of the complexity statement, not a proof for all exponents on this three-prime support.

## 7. Exact scope and output of the denominator census

The primary search checks every `4|b≤1000` with at least three distinct prime factors. There are 101 such denominators:

- 96 have three prime factors and use `N=10`, with `34³=39,304` unsigned triangle signatures each;
- five have four prime factors: `420,660,780,840,924`. They use `N=18`, with `66⁴=18,974,736` signatures each.

No `b≤1000` has five distinct prime factors while divisible by four, since `4·3·5·7·11=4620`.

All accepted cycles have exact period three, giving 62 distinct `(parameter, oriented cycle)` pairs up to cyclic rotation. For example, at `b=60` the complete list of parameter numerators is

```
−34861, −25621, −8149, −6469.
```

The largest absolute parameter numerator in the 101-denominator output is `103,928,749`. This emphasizes the absence of a `10^8` parameter-height filter; it is not a priority comparison.

The other possible denominators are covered by the published two-prime theorem or the necessary condition `4|b`. A non-square denominator admits no rational periodic point at all. Therefore the full scope is every rational `c` with reduced denominator at most `1,000,000`: if the denominator is square, its square root is at most 1000; if nonsquare, the denominator obstruction already applies.

Positive controls include `b=4,12,24,60,84`; all recovered three-cycles are checked against the independent marked-cycle parametrization. The complete 28-denominator `b≤400` census also passes a separate implementation using actual primitive binary words and direct recurrence equality instead of (7) and (10). The auxiliary word-table implementation independently searches all 96 three-prime denominators and checks the three-cycle parametrization at all 101 denominators; its five four-prime searches are explicitly marked NOT RUN. A separate adversarial audit independently checks all 101 denominators, including those five, using the faster two-edge algorithm above, and matches all 62 cycles. Its smooth-pair algorithm, independent of the three-edge table, also checks all 28 denominators through 400. These scopes are separate and are recorded in the companion audit artifacts.

## 8. A general resonant-lift obstruction and the prime 43

**Simple-resonance lemma (PROVED).** Let `p` be prime, `g∈Z_p[t]`, and `g(α)=α`, `α∈Z_p`. Let `λ=g'(α)`. Suppose, for some positive integer `r`,

```
v_p(λ−1)=0,
v_p(λ^r−1)=1.
```

Then for every nonzero `t∈p Z_p`,

```
v_p(g^r(α+t)−α−t)=v_p(t)+1.                (11)
```

**Proof.** Write `g(α+t)=α+λt+A t²+O(t³)`. Induction on the iterate gives

```
g^r(α+t)−α−t
 = (λ^r−1)t
   + A λ^(r−1)(1+λ+...+λ^(r−1))t²
   + O(t³).
```

The geometric sum equals `(λ^r−1)/(λ−1)` and is divisible by `p` under the hypotheses. All higher coefficients are integral. If `s=v_p(t)≥1`, the valuation of the linear term is `s+1`, while the quadratic term has valuation at least `2s+1` and higher terms at least `3s`. The linear term is uniquely smallest, proving (11). QED.

The condition that `λ−1` be a unit is essential to this proof and is not dropped. Without it, the quadratic coefficient can be a unit and cancellation can occur already at `s=1`; for instance `g(t)=(1+p)t+t²` has the additional fixed point `t=−p` despite `v_p(g'(0)−1)=1`.

For `r=2` and odd `p`, the hypothesis reduces to `v_p(λ+1)=1`, since then `λ−1` is automatically a unit.

### Iterate valuation and the whole inverse-limit obstruction

For odd `p`, the lemma strengthens to

```
v_p((g^r)^M(α+t)−α−t)=v_p(t)+1+v_p(M)       (13)
```

for every positive integer `M` and every nonzero `t∈pZ_p`. Thus `α` is the only periodic point of `g` in `α+pZ_p`.

To prove it, write `t=p^s y`, with `y` a unit, and rescale `F=g^r`:

```
H_s(y)=(F(α+p^s y)−α)/p^s
      =y+p u y+p²V(y),
u=(λ^r−1)/p ∈ Z_p^×, V∈Z_p[y].             (14)
```

The quadratic coefficient of `F` is divisible by `p`, and the terms of degree at least three gain `p^(2s)`, proving (14). For a polynomial `H(y)=y+p^a A(y)`, `a≥1`, odd-prime Taylor expansion gives

```
H^p(y)=y+p^(a+1)A(y) (mod p^(a+2)).
```

Indeed the first correction beyond `p·p^a A(y)` has factor `binom(p,2)p^(2a)`, and all remaining terms have adequate valuation; for `a=1` the same conclusion follows by keeping the second-order Taylor term before summing. Inducting from (14) gives `H_s^(p^j)(y)=y+p^(j+1)A_j(y)` with `A_j(y)≡u y (mod p)`. Iterating an additional number of times coprime to `p` preserves the first nonzero displacement valuation. This proves (13). A finite periodic return would make the left side infinite, which is impossible.

Apply the lemma to `T(a)=(a²−29)/4`, `g=T³`, any of its three rational integer cycle points `α∈{-7,5,-1}`, and `p=43`. The multiplier is `λ=35/8` and `λ+1=43/8`. Thus no nonzero displacement in `43Z_43` is a two-cycle of `g`, and hence no period-six point of `T` in that tube exists over `Q_43`.

### Exact shell classification at every depth

For every `k≥1`, restrict `T` modulo `43^k` to the union of the three residue disks `a≡−7,5,−1 (mod43)`. There is the original three-cycle. For every shell depth `1≤s<k`, there are exactly 21 further cycles, each of exact length

```
6·43^(k−s−1).                              (15)
```

Here a shell is the set of points whose difference from the corresponding rational cycle center has valuation exactly `s`. Any return under `T` must be a multiple of three, because reduction modulo43 is the exact three-cycle. The multiplier of `g=T³` has order two modulo43. Odd iterates of `g` therefore do not return modulo `43^(s+1)`, and (13) gives the exact return order of its even iterates. This proves (15), rather than merely an upper period bound. Across the three disks, the shell has `3·42·43^(k−s−1)` points; dividing by (15) gives exactly 21 cycles. Together with the original three points, the count sums to `3·43^(k−1)`, all points in the three disks.

Every depth `k≥2` consequently displays additional cycles, while the compatible inverse-limit tube contains only the original rational three-cycle as periodic points. Its other finite-level cycles continually increase their return times.

At modulus `43²`, the points with displacement of valuation one satisfy `g²(a)=a` and have exact `T`-period six. At modulus `43³`, equation (11) says none returns after six. The exact finite computation checks all 5,418 lifts and finds period 258 for each. Their splitting into 21 cycles is therefore certified at that modulus only; no infinite-level periodicity is inferred. There are still 21 OTHER six-cycles modulo `43³`, located closer to the rational cycle and reducing to period three modulo `43²`. The claim is failure to lift the old six-cycles with the same period, not absence of all six-cycles at the larger modulus.

## 9. Two precise bridges and the remaining proof targets

### Bridge A: arithmetic cycle → symbolic ultrametric → signed S-unit gluing

- Source: an ordered rational cycle and its common denominator.
- Target: one primitive binary word per bad prime, then the signed edge triple `(u,v,w)`.
- Map: residue-class coding, longest-common-prefix valuation, and prime reconstruction (8).
- Preserved: exact local distances, local cycle length, and the necessary first-two-arrow relation.
- Destroyed: global addition and agreement of all later arrows.
- Sidecar: signed vertex reconstruction and exact recurrence closure.
- Cheap decisive test: the `C=11,b=4` hostile passes the local triangle but fails return.
- **OPEN proof target:** for arbitrary denominator support, force every candidate satisfying all cyclic signed-edge gluing equations to have period at most three. The finite local circulation identities alone cannot do this.

### Bridge B: finite-ring cycle → next-level displacement → local compatibility

- Source: the explicit six-cycles modulo `43²` inherited from THM-4139.
- Target: valuation of the sixth-return displacement.
- Map: the return-map Taylor expansion at the rational three-cycle.
- Preserved: the chosen residue tube and return time.
- Destroyed by a single-modulus view: the first nonzero higher-order displacement.
- Sidecar: `v_43(λ+1)=1` and the full next-modulus test.
- Decisive outcome: no six-cycle lift modulo `43³`; periods grow by 43.
- **OPEN/generalization target:** catalogue resonant primes by `v_p(λ^r−1)` and nonlinear first nonzero terms, then combine compatible good-prime restrictions with the global signed-edge solver. A local cycle surviving every level would still need a separate rationality argument.

### A focused next theorem target

The first unsupported denominator-support class is three primes. Fix a set such as `{2,3,5}` but allow all exponents. The reduction leaves finitely many exponent patterns (`N≤10`) and signed four-term equations (10) in the prime powers. A uniform proof for that support needs an actual S-unit/exponential-Diophantine resolution or a complete modular covering certificate. Enlarging a finite exponent box is not such a proof. The existing exact solver provides hostile controls and a complete reconstruction map for auditing any proposed resolution.

## 9A. Standard infinite-family good-prime controls

**Critical-cycle certificate (PROVED elementary control; no novelty claim).** Let `c∈Z_p`. Suppose `x²+c` on `F_p` has exactly one cycle, of length `m`, and the derivative product around that cycle is zero. Then `f_c` has exactly one affine periodic cycle over `Q_p`, of exact length `m`. Therefore every rational periodic point, if any exists, has exact period `m`.

**Proof.** A nonintegral point has negative valuation that doubles under iteration, so cannot be periodic. Every integral periodic point reduces to the unique finite-field cycle. On each residue disk over that cycle, `f_c^m` maps the disk to itself and has Lipschitz constant at most `1/p`, because its derivative vanishes modulo `p` throughout the disk. The contraction theorem gives a unique fixed point in that disk and no longer periodic orbit there. The original map cyclically permutes these `m` fixed points, whose reductions are distinct. QED.

Two explicit infinite parameter classes are:

- `c∈Z_5` and `c≡1 (mod5)`: the sole residue cycle is `0→1→2→0`, and its derivative product is zero. Thus every rational periodic point has exact period three. The rational landmark `−29/16` belongs to this class.
- `c∈Z_3` and `c≡2 (mod3)`: the sole residue cycle is `0→2→0`, also with derivative product zero. Thus every rational periodic point has exact period two. The parabolic landmark `−7/4` belongs to this class; its real three-cycle cannot be rational.

Existence over `Q_p` is not existence over `Q`: for example `c=1` satisfies the first local condition but has no real periodic point, hence none rational. These families put no bound on denominator support away from the named prime. Conversely a finite bank of good-prime tests cannot cover arbitrary denominators: a denominator can contain every prime in that bank. The denominator-first signed-edge method supplies different information.

## 10. Reproduction

Primary certificate, from this directory:

```
python3 rational_cycle_signature_probe.py --max-b 1000 --max-omega 4 --output denominators_through_1000.json
python3 rational_cycle_pair_probe.py --max-b 1000 --compare denominators_through_1000.json --output pair_through_1000.json
python3 rational_cycle_pair_probe.py --factorization 2:2 3:20 5:20 --output large_exponent_control.json
python3 finite_lift_obstruction.py
python3 simple_resonance_controls.py
python3 critical_cycle_controls.py
```

Complete independent certificate, from the workspace root:

```
python3 session7-cycle-review/check_words.py
python3 session7-cycle-review/independent_search.py
python3 session7-cycle-review/word_projected_search.py
python3 -O session7-cycle-review/word_projected_search.py --output session7-cycle-review/projected_results_optimized.json
python3 session7-cycle-review/universal_pair_search.py
python3 session7-cycle-review/check_43_lift.py
python3 session7-cycle-review/check_general_resonance.py
python3 session7-cycle-review/compare_receipts.py
```

The independent full review is `session7-cycle-review/INDEPENDENT_REVIEW.md`; its comparison receipt is `session7-cycle-review/comparison_receipt.json`. The optional author-written word-table replay is supplementary and explicitly skips its expensive four-prime searches; it is not the full independent certificate above.

All mathematics in the programs uses integers or exact rational arithmetic. No floating-point filtering, probabilistic primality test, or conjectural period-six theorem is used. Canonical cycles are identified only by cyclic rotation, not reversal; different directions are not silently merged.

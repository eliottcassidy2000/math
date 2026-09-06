# NC2: rational root rotations are only dilation; sharp binomial packets remain algebraic

Status: **PROVED** by the arguments below; exact finite controls are separately
**FINITE-EXACT**. No improvement of the general first-return bound at
`min(M,N)>=3` is claimed. This is a method-boundary result and a faithful
replacement object, not a solution of that open problem.

## Inheritance and scope

- Closest proved mechanism: [THM-4417, width-two Laurent first-return parabolic
  critical bound](../../01-canon/theorems/THM-4417-width-two-laurent-first-return-parabolic-critical-bound.md).
  Its rational map is an **approximation** to a local involution. It does not
  assert that the exact involution is rational.
- Retained sidecar: [THM-2111, effective compound-root bound for one-variable
  constant terms](../../01-canon/theorems/THM-2111-effective-compound-root-bound-for-one-variable-constant-terms.md),
  particularly the logarithmic small-root-product identity.
- Canonical cancellation hostile: `z^-2+z^-1+z-z^2`, with first return four.
  Support reachability alone loses signed coefficient cancellation.
- New minimal interior width-three hostile: `z^-3+z^4`, first return seven.
  Its trace and norm retain a 15-sheet algebraic obstruction, with normalized
  packet-curve genus 65, even though its constant-term sequence is elementary.
- Corrected near miss: replacing a product of two other small roots by a
  single root-swap map is not type-correct. No petal-critical-point theorem
  for rational maps is transported to an algebraic correspondence here.

The live board is: local cyclic root action; marked-root quotient; trace and
norm; monomial dilation; parabolic critical budget. The map from a small-root
packet to trace/norm preserves its unordered algebraic coordinates but loses
the ordering needed for a single-valued root action. A branch selector is the
necessary sidecar. Exact quadratic divisibility is the cheapest faithful test.

Targeted duplication searches covered the current THM-4417 proof, NC2 result
notes, theorem/results occurrences of rational root rotation, finite-order
Möbius maps and monomial dilation, and the current mistakes ledger. No claim
of external priority is made. The parent session's empty-hexagon paper seed
is methodological only: use the smallest hostile and a faithful relaxed
encoding with a separate checker. No SAT or geometric theorem is imported.

## 1. Classification of globally rational local root rotations

Let `M,N>=1`, let `R` be a complex polynomial with `r=R(0)!=0` and
`deg R=M+N`, and put `h(z)=z^M/R(z)`. Choose a local branch

```text
w(z)=z/R(z)^(1/M),       iota(z)=w^(-1)(zeta*w(z)),
```

where `zeta` is a primitive M-th root of unity. For any integer `k`, let
`xi=zeta^k` and `q=ord(xi)`. The following are equivalent:

1. the germ `iota^k` extends to a rational function on the sphere;
2. `iota^k(z)=xi*z`;
3. `R(xi*z)=R(z)`;
4. `R(z)=S(z^q)` for some polynomial `S`.

In that event `q|M`, `q|N`, and

```text
z^-M R(z)=F(z^q),       F(u)=u^(-M/q) S(u).
```

Constant-term moments are preserved exactly by this dilation:
`CT_z(F(z^q)^m)=CT_u(F(u)^m)` for every integer `m>=0`.

Proof. A rational extension `J` satisfies `J^q=id` and `h o J=h`, by the
identity theorem. Degrees of rational compositions multiply, so `deg J=1`.
Thus `J` is Möbius and fixes zero. The zero fiber of `h` is precisely
`{0,infinity}`, with respective multiplicities `M,N`. Since `J` is injective,
preserves this fiber and fixes zero, it fixes infinity. Therefore `J(z)=xi*z`,
where its multiplier is determined by the germ. Invariance of `h` now says
`R(xi*z)=xi^M R(z)=R(z)`, exactly the asserted divisibility of all nonzero
coefficient exponents by `q`. The converse follows either from the defining
coordinate or uniqueness of the germ solving `h(J(z))=h(z)` with multiplier
`xi`. The coefficient identity under dilation is immediate. QED.

In particular, for prime `M=3`, any nonidentity rational root rotation reduces
the problem to width one. For composite `M`, rational powers identify only
genuine dilation quotients; they need not reduce all the way to width one or
two. This does **not** rule out useful rational approximations such as
THM-4417's map, or different rational maps on a larger state space.

## 2. Sharp primitive binomials have full marked-packet algebraic degree

Let `M,N>=1` be coprime, `d=M+N>=3`, and

```text
f(z)=z^-M+z^N,          h(z)=z^M/(1+z^d).
```

There are `M` local small roots of `h(y)=t`. Mark one of them as `z` with
`t=h(z)`. When `M>=2`, let `S` and `P` be respectively the sum and product of
the other `M-1` small roots. Then

```text
[C(z,S):C(z)] = [C(z,P):C(z)] = [C(z,S,P):C(z)]
                              = binom(d-1,M-1).
```

Consequently the combined trace and norm are no more rational than either
alone. Each already generates the same marked-packet field. At width three,
`N>=4` and `gcd(3,N)=1`, this degree is `binom(N+2,2)`, starting at 15 for
`N=4`. This algebraic degree is not a lower bound on first-return time.

### Ramification audit and full monodromy

The derivative is

```text
h'(z)=z^(M-1) (M-N*z^d)/(1+z^d)^2.
```

The zero fiber has local indices `M` at zero and `N` at infinity. Each of
the `d` poles, where `z^d=-1`, is simple and unramified. The remaining
critical points are the `d` simple roots `alpha^d=M/N`. Their values are
`h(alpha)=(N/d)*alpha^M`, distinct because `gcd(M,d)=1`. They are nonzero
finite values and hence do not collide with the zero fiber. The listed
ramification contributions sum to
`(M-1)+(N-1)+d=2d-2`, as required; no ramification is hidden at infinity.

For completeness, use the usual analytic continuation description of the
monodromy of a rational function: deleting its branch values gives a
connected degree-d covering, and loops around those values generate a
transitive permutation group. Each of the `d` simple nonzero branch values
gives a transposition. The product relation for loops on the punctured sphere
expresses the loop about zero as the inverse product of these loops (the
target value infinity is unramified). Thus the whole group is generated by
these transpositions. A transitive group generated by transpositions is
`S_d`: make one edge for each generating transposition; its graph is
connected, and edge transpositions of a connected graph generate `S_d`.

The splitting field over `C(t)` therefore has Galois group `S_d`. Marking the
root `z` replaces it by its stabilizer `S_(d-1)` over `C(z)`. This is the only
covering/monodromy input; the elementary loop and group arguments above show
exactly where coprimality enters. If coprimality is dropped, nonzero critical
values collide and dilation appears: `(-3,6)` is a hostile control.

### No trace or product collisions between different subsets

Write the remaining `n=d-1` distinct nonzero roots as `x_1,...,x_n`.
The selected packet has `k=M-1` members, with `1<=k<n`. Its unordered subsets
have `binom(n,k)` conjugates under `S_n`. The pair `(S,P)` cannot have more
than this many conjugates; we show that each coordinate alone has this many.

If two different k-subsets have equal sums, subtract to get
`sum c_i*x_i=0`, where `c_i in {-1,0,1}` are not all equal. Apply a
transposition interchanging indices with different coefficients and subtract
the original identity. It forces two distinct roots to be equal.

For products the same cancellation gives `prod x_i^c_i=1`. If one coefficient
is zero and another is nonzero, interchange their indices and divide the
identities: again two roots are equal. The remaining case has all
coefficients in `{1,-1}`; the two subsets are complementary. Applying the
full symmetric group yields the complementary-product identity for every
k-subset. Therefore every k-subset product has the same square, namely the
total product of all roots. Comparing subsets that differ in one element
gives `x_i^2=x_j^2` for every pair of roots. There can be at most two distinct
nonzero roots. For `n>=3` this is impossible; for `n=2`, the original
singleton-product identity already contradicts distinctness. QED.

## 3. Faithful width-three replacement and exact minimal hostile

For arbitrary width-three `R`, the other two small roots are exactly the
unique small branch of the quadratic-divisibility correspondence

```text
y^2-S*y+P divides (y^3 R(z)-z^3 R(y))/(y-z),
S/z -> -1,               P/z^2 -> 1             as z -> 0.
```

The quotient is polynomial in `y,z`. Dividing it by the monic quadratic
gives two polynomial remainder equations. Conversely these equations mean
literal divisibility; the displayed branch limits select the small packet.
They retain the needed information and make no single-valued dynamics claim.

For `R=1+z^7`, use the sign-reversed quotient

```text
Q(y,z)=z^3*y^6+z^4*y^5+z^5*y^4+z^6*y^3-y^2-z*y-z^2.
```

Indeed `(y-z)Q=z^3*y^7-(1+z^7)y^3+z^3`. The two exact remainder equations
`A=B=0` are

```text
A=z^3(S^5-4S^3P+3SP^2)+z^4(S^4-3S^2P+P^2)
  +z^5(S^3-2SP)+z^6(S^2-P)-S-z,
B=z^3(-S^4P+3S^2P^2-P^3)+z^4(-S^3P+2SP^2)
  +z^5(-S^2P+P^2)-z^6SP+P-z^2.
```

This compact two-equation encoding selects a 15-sheet object, not a rational
map. Extraneous branches are removed by the stated local limits, not by
treating every solution of the polynomial equations as the desired packet.

For this binomial, support balance already gives

```text
CT(f^m)=0 for 1<=m<7,          CT(f^7)=binom(7,3)=35.
```

The retained THM-2111 product identity yields

```text
P(z) = z^2/(1+z^7) * exp(sum_(m>=1) CT(f^m)*h(z)^m/m),
P(z) - z^2/(1+z^7) = 5*z^23 + O(z^30).
```

The root trace is a different observable. Lagrange inversion gives the sum
of all three small roots as

```text
sum_(m>=1) [x^(3m-1)](1+x^7)^m * t^m/m = 2*t^5+66*t^12+...,
S(z)+z = 2*z^15+O(z^22).
```

Thus the trace sees a defect at order 15 while the moment-controlled norm
defect starts at order 23. Adding the trace is faithful but does not preserve
the first-return valuation as its sole low-order invariant. The missing step
is a theorem about this algebraic correspondence with a suitable critical
budget, or a different rational approximation retaining the relevant norm
valuation. Neither is supplied here. In particular no application of the
rational-map Fatou petal theorem to this correspondence is justified.

## 4. Normalizing the minimal packet curve does not restore noninvertible dynamics

Let `X` be the smooth compact normalization of the connected curve with
function field `C(z,S,P)` for the minimal hostile `R=1+z^7`. Then

```text
genus(X)=65.
```

The **CITED** genus input is John Milnor, *Dynamics in One Complex Variable:
Introductory Lectures*, IMS Preprint 1990/5, partially revised September 5,
1991, [Theorem 5.1, Riemann-Hurwitz Formula](https://abel.math.harvard.edu/archive/118r_spring_05/docs/milnor.pdf#page=34),
printed page 5-2, PDF page index 33. The theorem and its triangulation proof
give `sum(e_x-1)=D*chi(target)-chi(source)` for a branched cover of compact
Riemann surfaces. We use it both for the packet cover and for the self-map
exclusion; the new work here is the explicit monodromy orbit calculation.

In particular, `X` admits no holomorphic self-map of degree greater than
one. This strengthens the rationality obstruction: it cannot be removed
merely by replacing the z-sphere by the exact normalized packet curve.
Multivalued correspondences, larger-dimensional state spaces and different
rational approximations remain outside this obstruction.

Proof. Over `C(t)`, the field marks one root and an unordered pair among the
other six. Its degree is `7*binom(6,2)=105`, and its monodromy is the `S_7`
action on those 105 configurations. This is the same field used above:
the marked root is `z` and the pair trace or norm determines the selected
subset. The cover is connected because the action is transitive.

Each of the seven simple branch transpositions fixes 35 configurations:
the marked root has five choices outside its two letters, and the pair
either consists of those two letters or avoids both, giving `1+binom(4,2)=7`
choices. The other 70 configurations form 35 two-cycles. Its ramification
contribution is therefore 35.

The zero-fiber permutation has cycle type `(3)(4)`. Call its two blocks
`A,B`, of sizes three and four. Its induced cycles are counted exactly by

| Marked-root block | Pair's blocks | Number of configurations | Cycle lengths |
|---|---|---:|---|
| A | AA | 3 | 3 |
| A | AB | 24 | 12, 12 |
| A | BB | 18 | 6, 12 |
| B | AA | 12 | 12 |
| B | AB | 36 | 12, 12, 12 |
| B | BB | 12 | 4, 4, 4 |

For example, an unordered pair in the four-cycle is either an opposite
pair of period two or an adjacent pair of period four; this gives the
six- and twelve-cycles in the third row. The other rows follow by the same
block periods and the relative offset from the marked element. In total
there are 12 cycles, so this branch contributes `105-12=93`. The previous
ramification audit ensures there are no further branch values. Applying
Riemann-Hurwitz to `X -> P^1_t` gives

```text
2 genus(X)-2 = -2*105 + 7*35 + 93 = 128.
```

Finally, a degree-D holomorphic self-map of `X` would give
`128=D*128+R` with `R>=0`, forcing `D<=1`. QED.

The checker independently enumerates the entire 105-configuration
permutation action and reproduces both cycle types, rather than trusting
the table or numerical algebraic-curve software.

### All primitive width-three binomials

More generally, let `N>=1` and `gcd(3,N)=1`, and let `X_N` be the smooth
compact normalized marked-packet curve for `f=z^-3+z^N`. Then the exact
all-height formula is

```text
genus(X_N) = [N(N+1)(N+3)-2-floor(N^2/2)]/2.          (G)
```

It is at least three. Thus **every** such normalized curve admits no
holomorphic self-map of degree greater than one. The cases `N=1,2` have
genera three and thirteen even though reflection puts their first-return
problems within the already solved width-one/two boundary. For the first
primitive interior case `N=4`, (G) gives 65. This distinction matters: the
obstruction concerns this exact marked-root object, not the existence of
some other first-return proof.

Proof. The connected cover `X_N -> P^1_t` has degree

```text
D=(N+3)*binom(N+2,2)=(N+3)(N+2)(N+1)/2.
```

Its monodromy and branch values are those proved in Section 2. A simple
transposition fixes `(N+1)[1+binom(N,2)]` marked-pair configurations, so its
ramification index sum is

```text
I=(N+1)(3N+2)/2.
```

There are `N+3` such simple branch values. At zero the permutation consists
of a three-cycle on a block `A` and an N-cycle on a block `B`. The complete
orbit count follows from the same six types:

| Marked-root block | Pair's blocks | Number of induced cycles |
|---|---|---:|
| A | AA | 1 |
| A | AB | 2 |
| A | BB | floor(N/2) |
| B | AA | 1 |
| B | AB | N-1 |
| B | BB | binom(N-1,2) |

The only unordered two-point subsets of an N-cycle with period shorter
than N are the opposite pairs for even N, with period N/2. Every other
such pair has period N. Because both N and (when applicable) N/2 are
coprime to three, adjoining the marked A element multiplies either period
by three; this proves the third row. In the fifth row the relative B offset
from the marked element gives `N-1` choices and the A element supplies the
coprime three-cycle. In the sixth row the marked B element forces full
period N, giving `binom(N-1,2)` cycles. Empty pair types at `N=1,2` contribute
zero; no large-N assumption was used.

The total number of zero-fiber cycles is therefore

```text
C_0=4+floor(N/2)+N(N-1)/2=4+floor(N^2/2).
```

There are no other branch values. Riemann-Hurwitz now gives

```text
2 genus(X_N)-2 = -2D+(N+3)I+(D-C_0)
               = N(N+1)(N+3)-4-floor(N^2/2),
```

which is (G). Its genus is at least three, so the self-map exclusion follows
exactly as in the minimal case. QED.

For every coprime `N` in `1<=N<=20`, the checker independently constructs
the complete marked-pair permutation action. It checks the full cycle-length
multisets of both generators, not only their total orbit counts, against
the formulas above. This finite replay checks all small parity and boundary
cases; the unbounded theorem rests on the six-type argument.

## 5. A useful inherited width-two corollary, not a new proof mechanism

If `P` is a complex polynomial with `P(0)!=0`, `e=deg P>=1`, and `k>=1`
is an integer with `ke>=3`, THM-4417 applied to `R=P^k` proves

```text
some 1<=m<=number of distinct roots of P<=e
has [z^(2m)] P(z)^(km) != 0.
```

The bound is independent of the power `k`. This is a direct corollary of the
already proved root-count bound, not new independent novelty and not an
analogous width-three statement. The strict condition `ke>=3` prevents a
one-sided/constant boundary from being silently included.

## Exact reproduction and boundary

Run

```bash
python3 04-computation/nc2_root_rotation_obstruction_overnight_hexagon_sep05.py
python3 -O 04-computation/nc2_root_rotation_obstruction_overnight_hexagon_sep05.py
```

The companion output records 373 explicit failure gates: all 91 coprime
binomial widths `1<=M,N<=12`; exact critical-value residue separation; the
primitive width-three degree formula for `4<=N<=20`; all three formal roots
of the minimal hostile through degree 30 over `Q(omega)`; their independently
checked defining equations; the trace and norm defect valuations; symbolic
quadratic-division equations; nonprimitive dilation/collision controls; and
the complete 105-state marked-pair cycle decomposition and genus calculation;
and complete marked-pair actions on all 23,100 configurations for the fourteen
primitive width-three binomials with `1<=N<=20`.
Normal and optimized output are byte-identical. This finite replay checks
the algebraic examples; the all-height algebraic-degree result comes from
the proof, not the finite universe. No numerical root selection, external
priority, or new Lean verification is claimed.

Raw-LF SHA-256 manifest:

```text
script f078906ab85b5d8162b28d684137d3e4be3248b43d6bb9a8ac0fce1838607a49
output 354c876c7aceaae0d07eb0a601f0a46cf97ddf7970c691ca529c6acc1925a661
```

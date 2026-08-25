# Khinchin, Duffin--Schaeffer, rational edges, and the two arithmetic trees

**Session synthesis, 2026-08-24.** Current truth is owned by the linked canon,
not by this reflection.  The session promotes two exact results:
[THM-4056](../01-canon/theorems/THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock.md)
and
[THM-4057](../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md).
The irrationality of Khinchin's constant and of `e+pi` remain **OPEN**;
Duffin--Schaeffer is **CITED**, not reproved; and no claim is made about
LRC(14), planar `JC(2)`, or Sun-hole classification.

| lane | status after session | exact advance | live boundary |
|---|---|---|---|
| rational edges | `PROVED + VERIFIED-EXACT` | reduced rationals are ordered coprime arcs; reciprocal is reversal | the raw object is bidirected/incomplete, not a tournament |
| Stern--Brocot / Berggren | `PROVED + VERIFIED-EXACT` | `D(wx)=D(x)+2|w|-#A`; Berggren `A`-Walsh is depth checkerboard | depth parity loses branch order and ancestry |
| tournament completion | `PROVED + VERIFIED-EXACT` | a scale-invariant depth gauge; infinite Pell and Fibonacci directed cycles | the gauge is added and loses gcd/height |
| divisor clocks | `PROVED + VERIFIED-EXACT` | `C_N` is the disjoint union of exact-denominator unit packets | a clock does not supply infinite correlation estimates |
| Duffin--Schaeffer | finite input compiled exactly; infinite theorem `CITED` | first moments and every finite labelled arc are preserved | cross-layer overlap and limsup independence |
| irrationality | scope/no-go sharpened | LCM clearing, Farey area, and reciprocal sheet are separated | nonzero integer linear form plus decay still missing |

## Inheritance pass

The assigned anchor was the rational-edge tournament picture.  The niche was
the Khinchin/Duffin--Schaeffer/irrationality interface.  The wildcard was the
intersection of LCM clocks with Stern--Brocot and Berggren trees.

| lane | closest proved mechanism | canonical hostile | corrected near miss | least-used sidecar |
|---|---|---|---|---|
| rational endpoints | THM-3333 Gaussian spinors and Farey determinant | `{2,4}` has no reduced rational arc; `{1,2}` has both reciprocal arcs | “all rationals form a tournament” | ordered endpoint sheet before quotienting |
| Berggren ancestry | THM-2596 prefix code and THM-3357 branch Walsh law | full reciprocal, subtree `1-x`, and outer reflection `(1-x)/(1+x)` are different | every mirror was called Stern--Brocot reflection | exact word and signed determinant |
| denominator packets | HYP-2628 already records `sum_(d|N)phi(d)=N` | `q|L_Q` is not `q<=Q` | Euler's identity was almost presented as new | explicit inverse and exact-denominator label |
| Duffin--Schaeffer | Koukoulopoulos--Maynard's reduced layers and weighted divergence | denominator layers overlap even when each layer is internally disjoint | formal mass was almost called union measure | labelled overlap graph and infinite second moment |
| Apéry-style proofs | the repo's denominator-clearing framework | fast convergence can survive clearing with error too large | a good rational approximant was almost called a proof | cleared height and a nonzero integer linear form |
| named constants | HYP-2212's trace/norm/discriminant reconstruction | knowing `e` and `pi` are transcendental says nothing individually about `e+pi` | symmetric shadows were treated as orientation data | signed discriminant sheet |

## Live concept board

| object | representation | invariant | operation | symmetry / quotient | destroyed information |
|---|---|---|---|---|---|
| primitive rational | reduced pair `(a,d)` | gcd one, exact denominator | reciprocal and mediants | common-scale quotient | original gcd and height |
| LCM phase | `x in C_N` | `d_N(x)=N/gcd(x,N)` | unit multiplication and clock inclusion | cyclic negation | Stern--Brocot ancestry |
| continued-fraction path | canonical digit runs | `D=sum a_i-1` | left/right child and reversal | reciprocal path mirror | endpoint magnitude without continuants |
| Berggren spinor | `(m,n)` and word in `A,B,C` | Gaussian norm cone, word depth | ternary child matrices | signed/positive Pythagorean quotient | reciprocal chirality after normalization |
| approximation layer | reduced centres plus radii | first mass `2phi(q)psi(q)/q` | union and limsup | circle negation | correlations if only mass is kept |
| irrationality candidate | `A_n alpha-B_n` after LCM clearing | integer height and analytic error | recurrence / Padé step | projective rescaling | nonvanishing if only ratios are retained |

The board exposed two orthogonal gradings.  THM-4056 grades a rational by
**denominator divisibility/LCM height**.  THM-4057 grades it by
**subtractive-Euclid/Stern--Brocot depth**.  Berggren adds a third grading,
**ternary ancestry**.  The session found exact maps between selected
characters of these gradings, but no faithful identification of the three
trees.

## Exact advance I: type the rational-edge object before tournamentizing it

For distinct positive naturals, reduced positive rationals give exactly

```text
a/b  <->  ordered arc a -> b with gcd(a,b)=1.
```

Reciprocal swaps the endpoints.  This realizes the user's reflection idea
literally, but it also supplies the decisive hostile:

- a coprime pair carries both reciprocal arcs;
- a noncoprime pair carries neither;
- `1/1` is a loop.

The faithful object is therefore the directed double cover of the
coprimality graph, not a tournament.  Restricting to `a<b` gives an acyclic
orientation of that graph; completing every absent pair by size gives the
transitive tournament.  A nontransitive completion needs an additional odd
observable.

THM-4057 supplies one lawful gauge.  If `g=gcd(a,b)` and `D` is canonical
Stern--Brocot depth, put

```text
sigma_D(a,b)=sgn(b-a)(-1)^D((a/g)/(b/g)).
```

Reciprocal preserves `D` and changes the sign factor, so this is
antisymmetric and total.  It is invariant under common scaling.  Its first
directed triangle is

```text
2 -> 1 -> 5 -> 2.
```

This is an exact tournament, but not a free consequence of rationality: it
reverses some natural numerator-to-denominator arcs and deletes gcd and
height.  Those are mandatory sidecars for arithmetic applications.

## Exact advance II: the Stern--Brocot checkerboard was already inside Berggren

Under the repo's positive Berggren parameter convention, the three branches
have Farey-prefix codes

```text
A=rho,             B=lambda rho j,             C=lambda^2,
j(x)=1-x.
```

Each Farey child adds one edge and `j` preserves depth.  Hence, for every
root `x in (0,1)` and every Berggren word `w` of length `h`, THM-4057 proves

```text
D(wx)=D(x)+2h-#A(w).                                  (1)
```

THM-3357's apparently branch-specific Walsh character is therefore the
pullback of the Stern--Brocot checkerboard:

```text
sum_(|w|=h)(-1)^D(wx) w u = (-1)^D(x) C^(2h)u.        (2)
```

At the standard root `(1,2)`, the signed ray is `(-1,-4h-2)` and the level
has `(3^h+1)/2` odd-depth nodes and `(3^h-1)/2` even-depth nodes.  This is an
all-root, all-depth theorem; the exact computation through `9841` standard
nodes and `1820` nonroot nodes is a hostile audit, not the proof.

The same gauge has two symbolic infinite `C3` families:

```text
1 -> P_j -> P_(j-1) -> 1             (odd j>=3),
1 -> F_n -> F_(n+1) -> 1             (n=5 mod 6),
1 -> F_(n+1) -> F_n -> 1             (n=0 mod 6).
```

The first Pell cycle `(1,5,2)` is the first tournament cycle and uses the
first three nonzero Pell numbers.  This touches THM-3744's Pell-prefix lane
through the same continued fractions, without transferring its lonely-runner
maximum.  The Fibonacci cycles use adjacent continuants and parity, not
THM-3509's Cassini matching quotient.

## Exact advance III: every divisor-complete primitive-centre family is one clock

For `C_N=Z/NZ`, define the exact denominator

```text
d_N(x)=N/gcd(x,N).
```

THM-4056 proves the explicit bijection

```text
disjoint_union_(d|N) U_d -> C_N,
(d,a) |-> (N/d)a,
```

with inverse `x |-> (d_N(x),x/gcd(x,N))`.  Therefore every denominator
weight and every Fourier mode have the exact identities

```text
sum_(x in C_N) W(d_N(x)) = sum_(d|N) phi(d)W(d),
sum_(x in C_N) W(d_N(x)) exp(2 pi i kx/N)
  =sum_(d|N) W(d)c_d(k).                              (3)
```

The novelty is not `sum phi(d)=N`, which was already inherited.  It is the
labelled compiler, its inverse, the unit-orbit description, direct-system
map, arbitrary weights, Ramanujan transform, and exact finite-arc geometry.

With `N=L_Q=lcm(1,...,Q)`, all Duffin--Schaeffer centres of denominator at
most `Q` occupy exactly the filtered phases

```text
{x in C_(L_Q): d_(L_Q)(x)<=Q}.                         (4)
```

The filter is load-bearing.  At `Q=5`, the unfiltered `C_60` also contains
denominators `6,10,12,15,20,30,60`.  The weighted Duffin--Schaeffer input is

```text
sum_(d|N) 2phi(d)psi(d)/d
  =sum_(x in C_N) 2psi(d_N(x))/d_N(x).                 (5)
```

Equation (5) is multiplicity mass, not union measure.  The exact hostile at
denominators two and four has formal mass one but union measure `3/4`.
Koukoulopoulos--Maynard's theorem needs the infinite correlation control that
this finite compiler deliberately does not claim.

## Exact advance IV: why `6,60,420,27720` kept reappearing

The sequence is

```text
6=L_3,       60=L_5,       420=L_7,       27720=L_11. (6)
```

It is the LCM filtration sampled at the first four odd primes.  More
structurally, `L_Q` changes only when `Q=p^e`; then the new clock has exactly

```text
(p-1)L_(Q-1)
```

new phases, precisely those whose exact denominator has the new `p`-adic
depth.  Thus the four numbers are valuation-depth clocks, not a mysterious
four-term recurrence.

For even `N=2^s m`, `m` odd, the phases whose reduced numerator and
denominator have opposite parity map injectively to primitive ordered
Pythagorean triples.  Their exact count is

```text
B(N)=N-(m+1)/2.                                        (7)
```

At the four clocks, (7) gives

```text
4, 52, 367, 25987.
```

This is a genuine intersection with Berggren: an LCM-height slice of its
primitive parameters.  It is not a Berggren ternary-depth level and forgets
ancestry.  Equation (1) supplies the orthogonal depth grading.

## Exact advance V: trace, norm, and the missing signed sheet

For an ordered endpoint pair `(p,q)`, define

```text
T=p+q,          P=pq,          Delta=q-p,
T^2-4P=Delta^2.                                         (8)
```

Reciprocal fixes `(T,P)` and negates `Delta`.  The same functor sends
`(e,pi)` to `(e+pi,e*pi,pi-e)`.  It explains exactly why a symmetric
trace/norm shadow cannot orient a tournament edge: the signed discriminant
sheet is missing.  It also recovers the proved HYP-2212 reconstruction fact
that any two of the three shadows determine `e` and `pi` algebraically.

This does **not** prove that `e+pi` is irrational.  The valid inherited
conclusion remains only that at least two of `e+pi`, `e*pi`, and `e-pi` are
transcendental.  Individual classification is open.

Gaussian squaring makes the same sheet visible:

```text
(p,q) |-> (q^2-p^2,2pq,p^2+q^2)
       =(T Delta,2P,T^2-2P).                           (9)
```

Endpoint reversal negates the first coordinate and fixes the other two.
After passing to a positive unordered Pythagorean triple, that chirality is
lost.  Coprimality is not enough for primitivity: odd/odd parameters have
content two, so opposite parity is another mandatory sidecar.

## Khinchin and Apéry: the useful result was a firewall

Reciprocal only inserts or removes a leading continued-fraction zero/digit;
the positive digit tail and any asymptotic Khinchin geometric mean are
reciprocal-even.  Such a statistic can colour an unordered rational edge but
cannot orient it without an odd sidecar.  MISTAKE-231 also blocks reading a
finite digit mean as a universal entropy.

The LCM compiler has a precise dual relationship with Apéry-style proofs:

- THM-4056 uses `L_Q` to **enumerate** every primitive denominator packet;
- an Apéry proof uses `L_n^m` to **clear** all denominators in a linear form.

The valuation address is shared, but the decisive analytic coordinate is
not.  Irrationality still requires a nonzero integer linear form whose error
tends to zero after multiplication by its cleared height.  Farey adjacency,
nonzero Casoratian, fast convergence, or a beautiful continued-fraction tree
does not supply that estimate.  This is exactly the MISTAKE-209 boundary.

No approximation sequence for Khinchin's constant or `e+pi` emerged, so no
irrationality claim was promoted.

## Incoming-agent integration and audit repairs

The parallel recovery and hostile audits materially changed the result.

- The repo census found the earlier Euler/totient packet identity, preventing
  a false novelty claim and forcing the theorem to isolate its explicit
  inverse, Fourier, and direct-system content.
- The rational-edge audit separated full reciprocal `1/x`, subtree mirror
  `1-x`, and THM-3357's determinant-`-2` outer reflection.  It also forced
  the odd/even Gaussian-content hostile.
- The independent theorem audit caught the `q=1` singleton, the distinction
  between prefix denominators and all divisors, formal mass versus union
  coverage, the covariant clock inclusion, and a coefficient-ring type error
  in the Fourier statement.  Each was repaired before promotion.
- Finite experiments were not allowed to carry the all-depth cocycle or the
  infinite Pell/Fibonacci laws; both were replaced by symbolic proofs.

## Frontier board after the session

1. **Depth character on a denominator packet.**  Study
   `S(q)=sum_(a in U_q)(-1)^D(a/q)`.  This is the first direct interaction
   between the THM-4056 denominator grading and the THM-4057 depth grading.
   Multiplicativity, prime-power recurrences, and Ramanujan transforms need
   hostile tests before any conjecture.
2. **Finite overlap operators.**  Retain the full labelled arc-overlap matrix
   on `C_(L_Q)` and ask which GCD-graph correlations become convolutional or
   block-circulant.  The target is a finite exact precursor, not a new proof
   of Duffin--Schaeffer.
3. **Three-height atlas.**  Compare LCM height, Stern--Brocot depth, and
   Berggren ternary depth on primitive spinors.  Determine sharp fibres and
   the first two nodes sharing every pair of heights but not the third.
4. **Tournament arithmetic.**  Determine strong components, score limits,
   and directed-cycle densities of `sigma_D`; preserve gcd/height as labelled
   edge data instead of pretending the gauge is faithful.
5. **Irrationality route.**  Search only for constants with an actual
   recurrence or integral representation producing cleared nonzero linear
   forms.  Khinchin digit statistics and the `(T,P,Delta)` reconstruction are
   structural sidecars, not substitutes for decay.
6. **LCM-clock crossovers.**  Revisit the repo's AP, Rule 30, Sun, and LRC
   clocks through the prime-power jump filtration.  A shared modulus is a
   candidate carrier; every consumer still needs its own height, owner, and
   overlap coordinates.

The strongest synthesis is deliberately narrow: denominator divisibility,
continued-fraction depth, and Berggren ancestry are three different
arithmetic coordinates.  The session found exact compilers and one exact
character pullback between them.  The apparent unifications that discard
radius, height, gcd, signed sheet, or branch word are precisely the ones the
hostile controls refute.

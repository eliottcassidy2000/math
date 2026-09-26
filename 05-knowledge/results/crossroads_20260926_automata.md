# Collatz crossroads: exact parity sections and two automatic Padé certificates

**Status:** PROVED in the elementary scopes below; FINITE-EXACT for the
enumerated universes; OPEN for the proposed automatic-tape programme.
No Collatz convergence, periodicity conjecture, or Rule-30 prize is proved.
No novelty claim is made for the parity conjugacy, finite carry machines,
Thue--Morse product, Rudin--Shapiro system, or Mahler approximation mechanism.

## Inheritance and concept board

The closest proved mechanisms recovered are:

* [THM-4210, Rule-30 lossless dyadic block-current Cartier tree](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md):
  a lossless scale tree need not have finitely many states, and the transverse
  channel is physically necessary even when it is null at the marked center.
* [THM-4072, Mahler safe-terminal fibre product](../../01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md):
  symbolic admissibility and an eventually zero native integer tape are
  separate conditions; a digit-coordinate change requires unbounded carry.
* [Creative finite-tape transducer](creative_transducer_20260925.md): one
  Collatz step has three multiplier carry states, while the tape and number
  of outer iterations are unbounded.
* [THM-4473, Collatz digit chains and repunits](../../01-canon/theorems/THM-4473-collatz-digit-chains-rotation-repunits.md):
  the parity map and `T^m(2^m-1)=3^m-1` are already inherited.
* [Transversality foundry](collatz_procgen_20260922_transversality_foundry.md)
  and [hard-class note](collatz_procgen_20260922_hard_class.md): the repository
  already has Sturmian, quasi-Sturmian, substitution-certificate, and
  square-swap results. Entropy zero is not a closure; the cube-swap value
  remains open. The earlier phrase "HARD means positive entropy" was
  explicitly corrected there.

Canonical hostiles are the all-ones native prefixes, `3n-1` cycles, the
`5n+1` control, and cube swaps that have zero entropy but defeat the simple
periodic-approximation mechanism. The least-used relevant sidecar is the
**slope of a parity section**: the current low representative alone does
not determine how a future native digit changes the itinerary.

| Live concept | Operation / retained predicate | Missing coordinate or cheapest test |
|---|---|---|
| Rule-30 Cartier tree | retain a complete section, not a scalar center | physical transverse data |
| Mahler native terminal tape | intersect dynamics with an eventually zero tape | no finite terminal-prefix closure |
| Collatz parity sections | compose native-bit extension and one time step | affine slope and carry |
| Automatic symbolic words | replace equal-length, equal-weight blocks | exact Bernstein value, not digit density |
| Padé plus Mahler iteration | increase local contact before scaling | compare exact 2-adic gain with rational height |
| Existing repetition criterion | `Dio(w)>eta(w)` | do not infer a class theorem from a finite census |

The first connection supplies an obstruction and exact state space. The
second supplies a positive arithmetic exclusion certificate.

## 1. Exact parity-section automaton

Use the shortcut map

```
T(n)=n/2 if n is even; T(n)=(3n+1)/2 if n is odd,
Q(n)=sum_(j>=0) (T^j(n) mod 2) 2^j in Z_2.
```

The first `m` parity bits determine, and are determined by, `n mod 2^m`.
One direct proof compares two integers whose difference has valuation `h`:
their first `h` branch choices agree, each common step lowers the valuation
of the difference by one, and the next parities differ. Thus `Q` is an
isometry, and passage to the inverse limit gives a bijection of `Z_2`.

For a native prefix `0<=r<2^m`, put `k` equal to the number of odd inputs in
its first `m` iterates and `u=T^m(r)`. Then, for every `z in Z_2`,

```
T^m(r+2^m z)=u+3^k z,
Q(r+2^m z)=Q(r) mod 2^m + 2^m Q(u+3^k z).       (A1)
```

Here the expression before `+2^m` means the least residue of `Q(r)`.
The first formula follows by composing the same `m` affine branches; the
second retains the entire future itinerary. Thus the exact section is the
function `z -> Q(u+3^k z)`, with state `(k,u)`.

Reading the next native bit `b in {0,1}` gives

```
e=(u+b) mod 2,
k'=k+e,
u'=[3^e (u+3^k b)+e]/2.                          (A2)
```

It emits `e`. The initial state is `(0,0)`. Distinct pairs give distinct
section functions: evaluate at `z=0` and use injectivity of `Q` to recover
`u`; then evaluate at `z=1` to recover `3^k`. Consequently this is already
the exact semantic state partition, although it is infinite.

**PROVED obstruction.** No finite-state synchronous letter-to-letter
machine converts the least-significant-first binary expansion of every
nonnegative integer, padded with zeros, into its complete Collatz parity
itinerary. Indeed the prefix `1^m` reaches `(m,3^m-1)`. Its all-zero
continuation emits `Q(3^m-1)`. These outputs are pairwise different by the
isometry, so the machine would need infinitely many states. The same
argument is a Myhill--Nerode separation of the rooted native/parity
pair-prefix language: a finite common zero input tail and the appropriate
parity output prefix distinguish any two of these residuals.

This does not contradict the three-state **one-step** multiplier machine,
which scans a finite tape and flushes its end carry. It also does not say
that the itinerary of any particular positive integer is nonautomatic.
Under Collatz, each such itinerary would be eventually periodic.

**Minimal carry-only hostile found by the probe.** At depth three, native
prefixes `r=1` and `r=4` both have `u=2`, but their states are `(2,2)` and
`(1,2)`. Appending native tail value `z=1` gives future starts `11` and `5`;
their future parity outputs differ. Retaining only `u` loses the slope.
There is no such collision at depths one or two (direct enumeration).

Exact distinct section counts at depths `1,...,18` are

```
2,4,7,12,21,38,69,127,235,438,819,1535,2883,5425,
10218,19275,36403,68835.
```

These are finite data; no asymptotic rate is inferred. The signed and
higher-multiplier formula is `u+q^k z` with
`u'=[q^e(u+q^k b)+sigma e]/2` for odd `q,sigma`.

**Connection contract.** Source: the native binary prefix tree. Target:
the parity prefix tree. Map: `Q` and (A1). Preserved: exact finite itinerary,
dyadic distance, and affine continuation. Lost by a finite scalar quotient:
the unbounded slope/carry pair. Required sidecar: `(k,u)` or its exact
equivalent. Decisive tests: all-ones prefixes and the `(2,2)/(1,2)` pair.
The Rule-30 transfer is the requirement of a lossless physical section,
not a claim that Rule 30 and Collatz have conjugate dynamics.

## 2. Equal-weight block coding turns Collatz into a Mahler value

Let `q>=3` be odd and use `(qn+1)/2` on odd inputs. For a binary tape `w`
with ones at positions `d_0<d_1<...`, its unique 2-adic source is

```
Phi_q(w)=-sum_(j>=0) 2^(d_j)/q^(j+1).              (B1)
```

The formula follows either by inverse affine composition or reduction
modulo every `2^m`. It converges in `Q_2`; real convergence is unnecessary.

Fix `L>=3`, `s=L-1`, and blocks

```
A=1^s 0,      B=1^(s-1) 0 1,       rho=2^L/q^s.
```

Both blocks have length `L` and `s` ones. Let `t_n` be Thue--Morse, the
parity of the binary digit sum of `n`, and concatenate `A` when `t_n=0`,
`B` when `t_n=1`. Define

```
g=-sum_(j=0)^(s-1) 2^j/q^(j+1),
d=-2^(s-1)/q^s,
F(z)=sum_(n>=0) (-1)^(t_n) z^n=prod_(j>=0)(1-z^(2^j)).
```

The product identity follows by choosing a finite subset of binary place
values and then taking coefficientwise limits. Each block contributes
`rho^n (g+d t_n)` to (B1). Hence, in `Q_2`,

```
Phi_q(w)=(g+d/2)/(1-rho) - (d/2) F(rho).           (B2)
```

The coefficient of `F(rho)` is nonzero. The equal weight is essential:
without it, the next block has a selector-dependent denominator power,
and a one-variable evaluation at fixed `rho` no longer follows.

This is the precise useful transfer from automatic sequences: prove an
arithmetic statement about the value of a functional-equation solution.
Automaticity by itself is not a termination certificate.

## 3. An elementary 2-adic Padé certificate

**PROVED.** Let `p` be a positive even integer, `b` a positive odd integer,
`gcd(p,b)=1`, and `a=v2(p)`. If

```
M=max(p,b) < 2^(2a),                              (P1)
```

then `F(p/b)` is irrational in `Q_2`.

The short Padé approximation is

```
R(z)=(1-2z^2)/(1+z),
(1+z)F(z)-(1-2z^2)=2z^6+O(z^8).                  (P2)
```

The seven coefficients through degree six verify the contact order.
All defect coefficients are even, since modulo two
`F(z)=1/(1-z)` and `(1+z)F(z)=1`. Thus for `v2(x)>0`,

```
v2(F(x)-R(x))=6v2(x)+1.                          (P3)
```

Iterating `F(z)=(1-z)F(z^2)`, put `N=2^K` and

```
R_K=prod_(j<K)(1-rho^(2^j)) * R(rho^N),
rho=p/b.
```

Every factor in the finite product is a 2-adic unit. Consequently

```
v2(F(rho)-R_K)=6aN+1.                            (P4)
```

An unreduced integral numerator and denominator are

```
U_K=prod_(j<K)(b^(2^j)-p^(2^j))*(b^(2N)-2p^(2N)),
V_K=b^(2N-1)*(b^N+p^N).
```

They are odd, and

```
max(|U_K|,|V_K|) <= 2 M^(3N-1).                  (P5)
```

The bound uses `|b^h-p^h|<=M^h`, so it remains valid when `rho>1` in
the real absolute value. No real analytic continuation is used.

If `F(rho)=A/B` were rational in lowest terms, then `A,B` would be odd.
The nonzero integer `B U_K-A V_K` has valuation exactly `6aN+1`, by
(P4), but ordinary size at most `2(|A|+|B|) M^(3N-1)`. This contradicts
(P1) as `N` tends to infinity. QED.

For the block coding, `p=2^L`, `b=q^(L-1)`. Therefore the sufficient test is

```
q^(L-1)<4^L.                                    (P6)
```

It holds for **every `L>=3` under `3n+1`**, so none of these Thue--Morse
block-swap tapes is the parity itinerary of any rational 2-adic source,
positive or negative. Prepending any finite parity word also cannot make
the source rational, since the corresponding iterate of a rational is
rational. This particular `3n+1` family is already within the repository's
repetition-covered territory; the purpose here is an explicit independent
functional-equation certificate, not a newly solved Collatz region.

The proof also gives `q=5,L=5` (`rho=32/625`) and `q=7,L=3`
(`rho=8/49`). It does not certify `q=5,L=10` (`rho=1024/1953125`), where
the strict inequality fails. That failure is a limit of this certificate,
not evidence that the value is rational. In the `q=5,L=5` example, the
inherited lower bound `Dio(Thue--Morse)=5/3` alone falls below the height
rate `(4/5)log_2(5)`; this note does not claim that every possible
periodic-approximation certificate for the coded tape fails.

**Connection contract.** Source: the Thue--Morse selector, equivalently
the solution of `F(z)=(1-z)F(z^2)`. Target: a particular Collatz parity
tape and its exact Bernstein source. Map: equal-weight blocks and (B2).
Preserved: exact rationality/nonrationality, all parity bits, and block
ones density. Lost by the density quotient: ordering, hence the value of
`F(rho)`. Sidecar: the functional equation and its finite Padé defect.
Decisive test: (P4) against (P5), with failed case `q=5,L=10` retained.

## 3a. A second round: the Rudin--Shapiro transverse channel closes

Let `a_n=(-1)^(number of occurrences of 11 in the binary expansion of n)`
and `R(z)=sum a_n z^n`. Splitting the index into its even and odd parts gives

```
R(z)=R(z^2)+z R(-z^2),
R(-z)=R(z^2)-z R(-z^2).                           (RS1)
```

Thus a scalar section alone is not closed. Keeping the two channels
produces a finite exact system. The word and its binary definition are
also Example 3 of [Schaeffer--Shallit](https://arxiv.org/pdf/1104.2303);
the derivation and all arithmetic estimates here are elementary.

An exact coefficient solve gives

```
P(z)=1+z-2z^3-z^4+z^5-4z^6,
Q(z)=1-z^2-z^4-z^6,
Q(z)R(z)-P(z)=-4z^11-2z^12-2z^13+O(z^17).       (RS2)
```

All later defect coefficients are even, because each is a sum of four
signs. Therefore for `v2(x)>=2`,

```
v2(R(x)-P(x)/Q(x))=11v2(x)+2.                    (RS3)
```

There is no pole: `Q(x)` is a 2-adic unit whenever `v2(x)>0`.
The restriction `v2(x)>=2` is intentional; the leading and next terms
need not have strictly different valuations at `v2(x)=1`.

For `N=2^K`, `K>=1`, define

```
A_K(z)=sum_(0<=n<N/2) a_n z^n,
B_K(z)=sum_(N/2<=n<N) a_n z^n.
```

The binary digit split, retaining the possible `11` across the joining
boundary, gives the exact identity

```
R(z)=A_K(z)R(z^N)+B_K(z)R(-z^N).                 (RS4)
```

Indeed for `n=jN+r`, the boundary contributes `(-1)^j` precisely when
`r>=N/2`. This both proves (RS4) and explains why the second channel
cannot be dropped.

For `rho=p/b` as in section 3, put `a=v2(p)`, `M=max(p,b)`. Since `Q`
is even, the common-denominator approximant is

```
S_K=[A_K(rho)P(rho^N)+B_K(rho)P(-rho^N)]/Q(rho^N).
```

Here `A_K(rho)` is a unit and `v2(B_K(rho))=aN/2`. The two approximation
errors in (RS4) have valuation `11aN+2`, so their coefficients separate
those valuations and prevent cancellation. Consequently

```
v2(R(rho)-S_K)=11aN+2.                           (RS5)
```

The numerator polynomial has degree at most `7N-1` and coefficient
one-norm at most `10N`, since `||P||_1=10`; the denominator has degree
`6N` and one-norm four. Clearing powers of `b` gives

```
H(S_K) <= 10N M^(7N-1).                          (RS6)
```

The same nonzero-integer separation used in section 3 now proves:

**PROVED.** `R(p/b)` is irrational in `Q_2` whenever

```
max(p,b)^7 < 2^(11v2(p)).                        (RS7)
```

The nonzero error is supplied by (RS5), not by an assumed transcendence
theorem. Replacing the selector `t_n` in (B2) by `(1-a_n)/2` replaces
`F(rho)` by `R(rho)` and leaves the rest of the bridge unchanged.
In particular, `q=3,L=10`, with `rho=1024/19683`, is certified, as is
`q=5,L=3`, with `rho=8/25`.

For `q=3`, this particular certificate works exactly for `3<=L<=117`.
The first method miss is `L=118`, because
`3^(7*117)>=2^(11*118)`. This is not a claim that the corresponding
value is rational, or that its irrationality is open in the literature.
It only states the boundary of (RS7). An exact search of the
**full-rank even-denominator ansatz** with `deg Q<=2d`, `deg P<=2d`,
`d<=32`, finds no better contact-to-height-degree ratio than `11/7`
at `d=3`. Singular linear systems and other Padé shapes are outside
this finite search, so this is a recorded stopping reason, not a no-go.

This second round is a concrete success for the Rule-30 inheritance:
keeping the transverse channel enabled an arithmetic certificate after
scalar compression failed to close. No map between the two dynamics is
claimed. Both automatic-tape families are special structured tests;
their arithmetic exclusion does not establish that arbitrary integer
Collatz tapes belong to either family.

## 4. Synthesis, hypotheses, and stopping reasons

**OPEN programme A: automatic equal-weight selectors.** For a specified
automatic selector, build its finite Cartier system, translate an
equal-weight block coding into a vector of values at `rho`, and search
for low-degree simultaneous Padé forms whose gain beats coefficient
height. This is a bounded and falsifiable task per automaton. It does not
assume the parity tape of an arbitrary Collatz source is automatic.
The Rudin--Shapiro test was pursued in section 3a and passed after its
full two-channel system was retained. The next small method test is its
block length `L=118`, or a different automatic selector whose vector
Padé contact beats the relevant height. A successful certificate would
strengthen the method catalogue; it would still cover only those tapes.

**OPEN programme B: exact repetition versus Padé classification.**
[Schaeffer--Shallit, Theorem 22](https://arxiv.org/pdf/1104.2303)
proves that the Diophantine exponent of a specified automatic sequence is
computable and rational or infinite. Thus automatic test cases need not
remain finite-window `Dio` guesses. Compute that exact quantity before
claiming a Padé certificate crosses the repetition barrier. Their
Theorem 23 separately makes linear recurrence and its optimal constant
decidable in this category. These are CITED algorithms, not implemented
by this script.

**OPEN programme C: a physically meaningful section quotient.** The
exact affine pair `(k,u)` is an algebraic two-coordinate compression of
an infinite automaton. Seek a quotient that preserves a *first-descent
certificate* rather than the entire future parity output. The no-finite-
converter theorem does not forbid such a weaker quotient, adaptive
blocking, an unbounded rank, or a finite-description proof. Its first
hostile is the collision `(2,2)/(1,2)`, and its second is the arbitrary
all-ones prefix. No such universal descent quotient is proved here.

**Stopped transfer: recent p-adic automatic continued fractions.**
[Capuano--Checcoli--Mula--Terracini](https://link.springer.com/article/10.1007/s00209-026-04096-3),
Theorem 1.4 and Corollary 1.6 (published 2026-08-10), concern continued
fraction digits at odd primes and explicitly allow rational or quadratic
values. They do not give (B1), use `p=2`, or exclude rationality here.
No result from that paper is imported into the proof.

The concept-board comparison after the experiments is precise: the
Cartier/Mahler mechanism pays only when a finite *physical* system closes;
the parity-section obstruction explains why the universal conversion
does not close. Equal-weight substitutions create a special closed
system, and the Padé calculation then supplies an actual arithmetic
conclusion. This is the surviving bridge, not a claim that an automatic
model captures all Collatz trajectories.

## 5. Exact reproduction

```
python3 04-computation/experiments/crossroads_20260926_automata.py
python3 -O 04-computation/experiments/crossroads_20260926_automata.py
```

The standard-library script checks 24,552 affine sections, every parity
permutation through depth ten for `(q,sigma)=(3,1),(3,-1),(5,1)`, exact
distinct section counts through depth eighteen, the all-ones states
through depth 200, 54 Padé valuation/height pairs, and 150 Bernstein/Mahler
bridges. The latter use `q in {3,5,7,9,11}`, `L=3,...,12`, and
`8,32,128` blocks; independent ordinary dynamics on the reconstructed
least source verifies every generated parity bit. Failed sufficient
inequalities are reported rather than excluded from the sample.
Assertions use explicit exceptions and survive Python `-O`.

The second round adds 28 Rudin--Shapiro exact valuation/height checks,
the two-channel coefficient identity for `K=1,...,7` through degree
`5*2^K-1`, five block bridges with independent ordinary dynamics,
and the full-rank even-denominator Padé search `d<=32`.

The separately requested price audit is saved as
`04-computation/experiments/crossroads_20260926_automata_price_audit.py`.
It imports none of the other lane's code. Direct integer trajectories
for `n<4*2^L`, `L in {6,8,10,12}`, `K in {3,9,27}` check the actual
affine source interval, per-hub capacity, and both members of each
flipped pair. Every fixed-weight word of block length `2,...,16` is
checked under cyclic-minimum rotation. The independent proof audit
found the density factor correct: a pair serves at most `2M` starts,
while its index is at most `[K(X+L/3)+1]/2`, cancelling the factor two
and yielding lower flip density at least `rho_trim/(KM)`.

## 6. Final synthesis: identical growth clocks can have different arithmetic types

The critical-band flow proof and the automatic-value certificates have a
precise common object: **blocks with the same length and odd count**.
This yields a useful generalization of section 2 and a decisive stopping
reason for a proposed scalar growth-clock closure.

Let U and V be any two distinct binary words with common length b and
common number k of ones, for the `3n+1` map. The inverse block maps are

```
x -> rho x+g_U,      x -> rho x+g_V,
rho=2^b/3^k,
g_W=-sum_(j:W_j=1) 2^j/3^(number of ones in W through j).
```

The constants differ. Otherwise the two inverse branches would send
zero to the same 2-adic source, despite prescribing different initial
parity words; the parity isometry forbids that. With `d=g_V-g_U`,
Thue--Morse selection of U or V gives exactly

```
Phi=(g_U+d/2)/(1-rho)-(d/2)F(rho).                  (C1)
```

**PROVED.** This source is irrational for every such pair U,V. Indeed
`max(2^b,3^k)<4^b`, so section 3 applies, and `d!=0`. The proof is
independent of the internal positions of the ones. This class is still
inside the inherited repetition-covered scope; the contribution here is
its explicit common-slope affine map and functional-equation certificate.

An exact near-critical control is

```
b=11, k=7,
U=11111110000,       V=11111101000,
P=3^7/2^11=2187/2048,
g_U=-2059/2187,      g_V-g_U=-64/2187.
```

Both blocks have every prefix slope at least one, and every concatenation
has the identical endpoint clock `w_(11j)=P^j`. Both also obey the same
finite internal-growth bound used in the refined ballot certificate.
Nevertheless U repeated forever has rational source `-2059/139`, while
Thue--Morse choice has irrational source

```
-2091/139 + (32/2187) F(2048/2187).
```

Every finite prefix in either case has positive integer realizations.
Neither observation makes the infinite source a positive integer.
Thus **rational source realizability cannot be decided from the common
block endpoint clock and that growth cap**. The missing datum is the
ordered affine carry, encoded by F in (C1). This does not challenge the
flow theorem: its proof uses actual integer incidence, not that scalar
clock alone. It explains exactly why its entropy and cap estimates cannot
be promoted to an infinite integer-survivor theorem without another input.

Connection contract: source = equal-composition admissible blocks; target
= a common flow growth band plus a 2-adic source. Map = block composition
and (C1). Preserved by the endpoint quotient = every block slope and the
finite cap. Destroyed = rational versus irrational source type. Restoring
sidecar = the ordered translations `g_U,g_V` and selector. Decisive test =
the displayed periodic/Thue--Morse pair. For this pair the full scalar
clock is exactly identical, not merely asymptotically close.

The final independent audit of the refined finite certificate also passed:
for every b<=12, t=1,...,5 and r=0,...,b-1, exhaustive admissible block
enumeration checked the cap
`max(P^(t-1)(3/2)^k,P^t(3/2)^r)`, the grouped integral capacity, and the
factored lower bound (390 exact checks). The main THM-4478 proof's
arbitrary-edit construction and lower-density passage were re-audited;
no new gap was found. This final probe was run independently from a
standard-library stdin script; the maintained ballot script reproduces
the displayed certificate parameters, while the earlier standalone
price-audit script retains the independently implemented source counts.

# Collatz crossroads: exact parity sections and a Thue--Morse Padé certificate

**Status:** PROVED in the elementary scopes below; FINITE-EXACT for the
enumerated universes; OPEN for the proposed automatic-tape programme.
No Collatz convergence, periodicity conjecture, or Rule-30 prize is proved.
No novelty claim is made for the parity conjugacy, finite carry machines,
Thue--Morse product, or Mahler approximation mechanism.

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

## 4. Synthesis, hypotheses, and stopping reasons

**OPEN programme A: automatic equal-weight selectors.** For a specified
automatic selector, build its finite Cartier system, translate an
equal-weight block coding into a vector of values at `rho`, and search
for low-degree simultaneous Padé forms whose gain beats coefficient
height. This is a bounded and falsifiable task per automaton. It does not
assume the parity tape of an arbitrary Collatz source is automatic.
The smallest next test is the Rudin--Shapiro selector with the same two
blocks, retaining its full two-channel Mahler system rather than only its
balanced scalar output. A successful certificate would strengthen the
method catalogue; it would still cover only those tapes.

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

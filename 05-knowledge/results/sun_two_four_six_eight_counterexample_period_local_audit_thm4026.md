# THM-4026 independent counterexample and local-arithmetic audit

**Status (2026-08-24):** **FINITE-EXACT + INDEPENDENTLY VERIFIED** for the
counterexample computation; **PROVED algebra** for the period and universal
modular-solubility lemmas below.  This report audited the then-reserved canon
files
`THM-4026-sun-two-four-six-eight-binomial-counterexample.md` and
`THM-4027-sun-two-four-six-eight-universal-modular-solubility.md`; their later
status promotion is owned by the canon files, not by this audit report.

## Exact counterexample certificate

Let

\[
N=896315812331399.
\]

The dependency-free checker
`04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py`
proves that there are no `w,x,y,z >= 2` with

\[
N={w\choose2}+{x\choose4}+{y\choose6}+{z\choose8}.
\]

For `k=4,6,8`, all legal indices below `k` lie in the zero fibre, so they
may be represented canonically by `k-1`.  Exact monotonicity gives

```text
x = 3..12112     (12110 values)
y = 5..932       (928 values)
z = 7..281       (275 values)
```

and hence exactly `3,090,472,000` canonical triples.  For each triple put

\[
D=1+8\left(N-{x\choose4}-{y\choose6}-{z\choose8}\right).
\]

A representation exists if and only if the parenthesized residual is at
least one and `D=q^2` for an odd `q>=3`, in which case `w=(q+1)/2`.  This is
the exact identity `1+8 C(w,2)=(2w-1)^2`.

For every odd prime `p`, the checker stores a bit mask of `x` for which `D`
is a square modulo `p`, indexed by the `(y,z)` contribution.  The mask
intersection through `p=89` leaves exactly `324` triples.  Of these, `31`
have residual below one; exact integer square-root tests reject the remaining
`293`.  A redundant terminal route continues
the residue masks and leaves `3,3,1,0` triples after primes
`113,127,131,137`.  Thus one route ends in exact `isqrt`, while the other is
a pure bounded congruence-cover certificate.  Neither is a single local
obstruction: the top-index bounds are load-bearing.

The output also contains a planted positive representation and the complete
prime-by-prime survivor fingerprint.  Normal and optimized Python outputs
match byte for byte.

## Exact least periods

**Lemma (PROVED).**  For a prime `p`, `a>=1`, and `k>=1`, the least positive
period of `n -> C(n,k) mod p^a` is

\[
P_{p^a,k}=p^{a+\lfloor\log_p k\rfloor}.
\]

Consequently, for `m=prod p^a`, the least period modulo `m` is

\[
P_{m,k}=\prod_{p^a\parallel m}p^{a+\lfloor\log_p k\rfloor}.
\]

**Proof.**  Vandermonde gives

\[
{n+T\choose k}-{n\choose k}
=\sum_{j=1}^k {T\choose j}{n\choose k-j}.
\]

The difference vanishes for every `n` modulo `p^a` if and only if every
`C(T,j)`, `1<=j<=k`, vanishes: necessity follows by evaluating successively
at `n=0,1,...,k-1` in the triangular binomial basis.  Write
`r=floor(log_p k)`.  For `T=p^(a+r)`,

\[
v_p{T\choose j}=a+r-v_p(j)\ge a
\]

for `1<=j<=k`, proving sufficiency.  Conversely, periodicity first forces
`p^a|T`; write `s=v_p(T)`.  If `s<=r`, Lucas at `j=p^s` makes
`C(T,p^s)` a unit, impossible.  If `s>r`, then

\[
{T\choose p^r}=\frac{T}{p^r}{T-1\choose p^r-1},
\]

and the second factor is a unit by Lucas, so its valuation is `s-r`; it can
be at least `a` only if `s>=a+r`.  CRT gives the composite formula.  In
particular, for all four degrees here, `840m` is a common (usually
nonminimal) period modulo `m`.  QED.

## Universal modular solubility

**Theorem (PROVED).**  For every `m>=1` and every residue `c mod m`, there
exist `w,x,y,z>=2` such that

\[
{w\choose2}+{x\choose4}+{y\choose6}+{z\choose8}\equiv c\pmod m.
\]

**Small primes.**  At `p=2`, on the even class `w=2q`,

\[
{2q\choose2}=q(2q-1),\qquad
f(q)-f(r)=(q-r)(2q+2r-1).
\]

The second factor is odd, so this is a permutation modulo every `2^a`.
For `p=3,5,7`, use respectively the degree `p+1` variable and set its index
to `pq+1`.  Cancellation of the unique factor `p` gives

\[
{pq+1\choose p+1}=qU_p(q),\qquad U_p\in\mathbb Z_p[q],\quad U_p(q)\equiv1\pmod p.
\]

Thus the map is the identity modulo `p` and has derivative one modulo `p`;
the usual one-digit lifting argument makes it a permutation modulo every
`p^a`.  Set `w=2` and the other inactive higher binomials to zero, shifting
the requested active residue by one.  Period shifts make all indices at
least two.

**Primes above eight.**  Put `A_k={C(t,k):t in F_p}` and restrict `A_2` to
the `p-1` inputs away from its unique critical point.  Fibre size bounds give

\[
|A_2^{reg}|\ge(p-1)/2,\qquad |A_k|\ge\lceil p/k\rceil
\quad(k=4,6,8).
\]

Cauchy--Davenport therefore covers `F_p` whenever

\[
(p-1)/2+\lceil p/4\rceil+\lceil p/6\rceil+\lceil p/8\rceil-3\ge p.
\]

Writing `p=24q+r` checks this for every prime except
`11,13,17,19,23,29,31,47`.  Exact regular sumsets for those eight primes
have successive sizes

```text
11:  5,10,11,11       13:  6,13,13,13
17:  8,16,17,17       19:  9,19,19,19
23: 11,22,23,23       29: 14,28,29,29
31: 15,31,31,31       47: 23,46,47,47.
```

Every residue therefore has a representation whose triangular derivative
is nonzero.  Holding the other coordinates fixed, scalar Hensel lifting in
`w` gives a solution modulo every `p^a`.  Finally use the exact periods above
to select compatible coordinate classes at every prime power, apply CRT
coordinatewise, and add common periods to make all indices at least two.
QED.

Hence `N` is locally represented modulo **every** modulus.  Its failure is
global/Archimedean, not a congruence obstruction.

## Local rarity and the failed counting heuristic

These are **FINITE-EXACT diagnostics**, not an asymptotic theorem.  The target
is `20 mod 33`, the unique minimum-density class for the true binomial
periods, with probability `16/1089` rather than the uniform `1/33`.  Exact
normalized factors at the audited small-prime levels, followed by stable
factors at the displayed primes above eight, are

```text
2^4: 1       3^2: 68/81       5^2: 566/625       7^2: 310/343
p=11: 72/121 p=13: 154/169    p=17: 240/289      p=19: 316/361
p=23: 472/529.
```

The first line is **FINITE-EXACT only**; no all-level stabilization statement
is used.  For `p=11,13,17,19,23`, the critical-tuple audit has no fully
critical target solution, which proves stability after the first level.

It is the minimum residue at `p=11,17,19,23` and second-lowest at `p=13`.
At `p=31` four fully critical solutions fail to lift, giving stable factor
`29198/29791`; at `p=499` two fail, giving `124250000/124251499`.

The continuous shell heuristic has exponent
`1/2+1/4+1/6+1/8-1=1/24` and Archimedean constant

\[
J=\frac{\prod_{k\in\{2,4,6,8\}}\Gamma(1+1/k)(k!)^{1/k}}
        {\Gamma(25/24)}
=25.323984\ldots.
\]

At this `N`, `J N^(1/24)=106.304...`; multiplying only the displayed local
factors through `p=23` already lowers the heuristic mean to about `25.69`,
while the exact count is zero.  Independence/Poisson intuition is therefore
an unsafe lower-tail model in this sparse mixed-degree problem.  This is a
heuristic comparison, not an asserted circle-method asymptotic.

## Connections and frontiers

- **Global-blindness hostile:** the source object is the bounded index box;
  reduction modulo any fixed modulus preserves local representability but
  destroys height and exact equality.  The required sidecar is the index
  interval plus the square discriminant.  This is the same methodological
  failure isolated by THM-2043/2050 for LRC local coordinates, but it is not
  an LRC reduction.
- **Support versus multiplicity:** modular support is complete, while exact
  global support misses `N`; representation density and support coverage are
  different observables, matching the THM-2000/2005 guardrail.
- **Carry lane:** the exact combinadic expansion
  `C(281,8)+C(279,7)+C(234,6)+C(212,5)+C(188,4)+C(136,3)+C(43,2)+C(15,1)`
  equals `N`.  The failed four-term representation asks whether Pascal carry
  moves can eliminate all odd ranks while leaving one atom at each even rank.
  A confluent carry obstruction or finite-state normal form is **OPEN**.
- **Repaired headline:** “every positive integer” is **REFUTED**.  Plausible
  replacements—density-one representability, finitely many exceptions, an
  asymptotic for representation counts, or infinitely many CRT-engineered
  holes—are all **OPEN** and mutually distinct.
- **Search frontier:** rank residue classes by their exact `p`-adic density,
  combine low-density classes by CRT, and retain height.  The decisive test is
  whether those progressions produce further exact holes, not merely smaller
  expected counts.
- **Certificate frontier:** minimize the prime bank/set cover behind the
  bounded non-square certificate, then formalize the range bounds, masks, and
  terminal implication.  A fixed bank without the height sidecar cannot prove
  nonrepresentability because universal modular solubility forbids it.

## Reproduction

```bash
python3 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
python3 -O 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
```

The matching frozen output is
`05-knowledge/results/sun_two_four_six_eight_counterexample_independent_audit_thm4026.out`.
SHA-256:

```text
script 39febe3affe818f03c1b3d83161aad688fc1e6da271495617757eb097dd022cc
output df11680c1a63266911a4f12bf1d3a71d27101dfe7dd6318ecffa34ddba05cf25
```

# Regular shortlex languages: harmonic density and the q-adic basin formula

**Status: PROOF CANDIDATE + VERIFIED-EXACT controls; awaiting independent
audit.**

## Inheritance pass

- **Closest proved mechanism.** THM-3497 proves logarithmic densities `1/3`
  and `17/96` for two concrete ternary group languages; THM-3359 proves the
  unary periodic analogue.
- **Canonical hostile.** A binary automaton that accepts exactly words whose
  first digit is zero has recurrent acceptance mass `1/2`, but its shortlex
  logarithmic coefficient is `log(3/2)/log(2)`.
- **Corrected near miss.** Aperiodic cylinder mixing is sufficient but not
  necessary.  A periodic irreducible automaton still has logarithmic density
  equal to stationary acceptance mass.
- **Least-used sidecar.** The within-level shortlex coordinate is a real
  `q`-adic number.  It records where a recurrent basin sits, not merely how
  many words enter it.

## 1. Statement

Let `A` be a complete deterministic finite automaton over the ordered
alphabet `{0,...,q-1}`, with `q>=2`, initial state `s_0`, transition `delta`,
and accepting indicator `f`.  Enumerate all words in shortlex order, starting
with the empty word at index one, and let `S_A` be the set of accepted
indices.

Put

```text
P(s,t)=#{d:delta(s,d)=t}/q.                              (1)
```

This is the Markov matrix obtained from independent uniform digits.  Its
finite-dimensional Cesaro projection exists:

```text
Pi=lim_(D->infinity) (1/D) sum_(n=0)^(D-1) P^n,
h=Pi f,                 Ph=h.                           (2)
```

For iid uniform digits `D_1,D_2,...`, define

```text
S_n=delta(s_0,D_1...D_n),
X=sum_(j>=1) D_j q^(-j) in [0,1],
c=1/(q-1).                                               (3)
```

Since `h(S_n)` is a bounded martingale, it converges almost surely; write

```text
H_infinity=lim_(n->infinity) h(S_n).                    (4)
```

Then the logarithmic density exists and is

```text
delta_log(S_A)
 :=lim_(N->infinity) (1/log N) sum_(m<=N,m in S_A) 1/m
 =1/log(q) * E[H_infinity/(c+X)].                       (5)
```

An incomplete DFA is covered by adjoining a rejecting sink.

If the reachable chain has one closed irreducible class with stationary
distribution `pi`, then `H_infinity=pi(f)` almost surely, so (5) reduces to

```text
delta_log(S_A)=pi(f).                                    (6)
```

No aperiodicity assumption is needed.  With several closed classes,
`H_infinity` is the stationary accepting mass of the class eventually
entered.  Formula (5) weights that basin by its location in the shortlex
`q`-adic interval; unweighted absorption probability is generally wrong.

## 2. Shortlex level normalization

A word `w` of length `n` and zero-based lexicographic rank `r(w)` has index

```text
I(w)=C_n+r(w),
C_n=(q^n-1)/(q-1)+1.                                    (7)
```

Write

```text
c_n=C_n/q^n=c+(q-2)q^(-n)/(q-1),
X_n=r(w)/q^n.                                            (8)
```

The harmonic mass of accepted words on level `n` is therefore exactly

```text
L_n=E[f(S_n)/(c_n+X_n)],                                (9)
```

where expectation is over the `q^n` equiprobable words.  Coupling the digits
across all levels gives `X_n->X` uniformly up to `q^-n`, and the denominators
are bounded below by `c`.  Hence

```text
L_n-E[f(S_n)/(c+X)] ->0.                                (10)
```

For `f=1`, (9)--(10) recover the full-level mass

```text
integral_0^1 dx/(c+x)=log q.                             (11)
```

This is why a Cesaro average over levels becomes a logarithmic density.

## 3. Markov--Cesaro proof

Let `g(x)=1/(c+x)`.  Approximate `g(X)` uniformly by the cylinder function
`g(X_k)` depending only on the first `k` digits.  For `n>=k`, the Markov
property gives

```text
E[g(X_k)f(S_n)]
 =E[g(X_k)(P^(n-k)f)(S_k)].                             (12)
```

Taking the Cesaro mean in `n` and using (2) yields

```text
lim_D (1/D)sum_(n<D) E[g(X_k)f(S_n)]
 =E[g(X_k)h(S_k)].                                      (13)
```

The harmonicity `Ph=h` makes `h(S_k)` a bounded martingale.  Thus
`h(S_k)->H_infinity` in `L^1`, while `g(X_k)->g(X)` uniformly.  Letting
`k->infinity` in (13) gives

```text
lim_D (1/D)sum_(n<D) L_n=E[g(X)H_infinity].             (14)
```

The union of complete levels through depth `D` has logarithm
`D log q+O(1)`.  Any partial final level has total harmonic mass `O(1)`, so
it vanishes after division by `log N`.  Equations (11) and (14) prove (5)
for arbitrary cutoffs, not just level endpoints.

For a finite Markov chain, `h` is constant on each closed irreducible class
and equals that class's stationary mean of `f`; transient values are
absorption-weighted.  This proves the interpretation following (6).

## 4. Three sharp examples

### One irreducible but periodic class

Accept all binary words of even length.  The two-state chain has period two
and stationary acceptance mass `1/2`.  At complete-level endpoints, the
ordinary counting densities tend along even and odd depths to

```text
2/3 and 1/3,                                             (15)
```

so natural density fails.  Formula (6) nevertheless gives logarithmic
density `1/2`.

### Two absorbing prefix basins

Accept binary words whose first digit is zero.  The accepting basin occupies
`X in [0,1/2)`, so (5) gives

```text
delta_log
 =1/log 2 * integral_0^(1/2) dx/(1+x)
 =log(3/2)/log 2.                                       (16)
```

This differs from the branch probability `1/2`.  Reversing the alphabet
order moves the same-sized basin to `[1/2,1]` and gives

```text
log(4/3)/log 2.                                         (17)
```

The two coefficients sum to one.  Thus shortlex harmonic density is
order-sensitive exactly when recurrent basins retain an early-prefix
coordinate.

### THM-3497's two Berggren languages

The variable-translation automaton is an irreducible group walk on
`S3 x C2`, with `4` accepting states out of `12`; (6) gives `1/3`.

The fixed-drift automaton is an irreducible period-two group walk on
`S4 x D4`, with `34=16+18` accepting states out of `192`; (6) gives

```text
34/192=17/96.                                           (18)
```

THM-3497's two-step kernel and return-length argument remains stronger: it
proves cylinder-uniform parity-level asymptotics.  The general theorem shows
that this extra mixing is unnecessary merely to establish (18).

## 5. Relation to arbitrary subsets of the harmonic series

Every subset of the naturals is faithfully represented by its labelled
coefficient sequence `(1_A(n)/n)_n`, but only regular shortlex subsets are
covered by this theorem.  Arbitrary subsets need not have logarithmic
density.  Formula (5) is therefore a classification of one structured
family of harmonic subseries, not a claim that every subseries has an
asymptotic coefficient.

The full coefficient sequence remains faithful; the scalar coefficient in
(5) is deliberately lossy.  Multiple automata and nonregular subsets can
share it.

## 6. Reproduction and next boundary

```bash
python3 04-computation/regular_shortlex_harmonic_density_automaton_probe_20260816.py
python3 -O 04-computation/regular_shortlex_harmonic_density_automaton_probe_20260816.py
```

The next useful extension is to make the multiple-class integral in (5)
algorithmic: derive certified upper/lower cylinder sums for each recurrent
basin and characterize when the result is a rational combination of
logarithms of algebraic numbers.  No such arithmetic classification is
claimed here.

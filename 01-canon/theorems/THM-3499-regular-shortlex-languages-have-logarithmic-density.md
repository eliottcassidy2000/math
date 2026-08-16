---
id: THM-3499
title: "Regular shortlex languages have logarithmic density"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.  Every complete
  finite automaton over an ordered q-letter alphabet has a logarithmic
  density under shortlex indexing.  The general coefficient is a joint
  q-adic-address/recurrent-class expectation; with one reachable closed
  class it is the stationary accepting mass, without an aperiodicity
  hypothesis.  For every language, regular or not, harmonic convergence is
  equivalent to convergence of its Kraft-type series.  For a regular
  language this is equivalent to trimmed spectral radius rho<q; rho=q is
  exactly the positive logarithmic-divergence regime.  No statement about
  arbitrary subsets having logarithmic density is asserted.
source: codex/turnpike-atlas/2026-08-16
audit: >
  Candidate commits 405f108bd and 9fe56bc56 were imported into a fresh
  worktree based on current origin/main.  The indexing, Markov--Cesaro,
  martingale, recurrent-class, arbitrary-cutoff, cylinder, functional-
  equation, Kraft, and spectral-radius arguments were reconstructed
  independently.  A separate exact SCC/stationary/absorption engine exhausts
  all 5,832 complete binary three-state languages and all 256 complete
  ternary two-state languages, and adds delayed-transient, reducible,
  periodic, alphabet-reversal, arbitrary-cutoff, and Berggren controls.
related:
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation
  - THM-3497-berggren-four-frame-return-languages-and-two-harmonic-densities
  - THM-3510-binary-shortlex-equal-level-count-log-density-boundary
scripts:
  - 04-computation/regular_shortlex_harmonic_density_automaton_probe_20260816.py
  - 04-computation/regular_shortlex_log_density_independent_audit_20260816.py
outputs:
  - 05-knowledge/results/regular_shortlex_harmonic_density_automaton_probe_20260816.out
  - 05-knowledge/results/regular_shortlex_log_density_independent_audit_20260816.out
script_sha256:
  - 7433a432faa30cfb2f5f54d9f759645da2fb665f19874881147f0ffe4100e866
  - 4ee6e0d44a8c52a4ac1111642f23558b2eb94795f92b867a72b7a9a41493023b
output_sha256:
  - 41f1103aa01c1050d1cb535865ef71a7bbe6791c028178f14b69ffaabeddcc99
  - 684308fdd163e82c02661ee40a00767a4edda642211e657002fb0cab59e4eb47
exact_audit_ledger_sha256: 046420bde716c6f83c7314c0010e1e4b19feffae59de7c9684fb4d0bdc4b1459
hash_basis: raw LF bytes for files; exact rational structural ledger as stated
---

# THM-3499 -- regular shortlex languages have logarithmic density

**PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.**

## 1. Automaton and shortlex coordinates

Let `q>=2`, and let `A` be a complete deterministic finite automaton over
the ordered alphabet `{0,...,q-1}`.  Write `s_0` for its initial state,
`delta` for its transition function, and `f:S->{0,1}` for its accepting
indicator.  Incomplete automata are included by adjoining a rejecting sink.

Enumerate all finite words in shortlex order, with the empty word at index
one, and let `S_A` be the set of accepted indices.  A word `w` of length `n`
and zero-based lexicographic rank `r(w)` has index

```text
I(w)=C_n+r(w),        C_n=1+(q^n-1)/(q-1).              (1)
```

Put

```text
P(s,t)=#{d:delta(s,d)=t}/q.                             (2)
```

This is the transition matrix obtained from independent uniform digits.
Its finite-dimensional Cesaro projection exists.  Define

```text
Pi=lim_(D->infinity) (1/D) sum_(n=0)^(D-1) P^n,
h=Pi f.                                                 (3)
```

Then `Ph=h`.  For iid uniform digits `D_1,D_2,...`, set

```text
S_n=delta(s_0,D_1...D_n),
X_n=sum_(j=1)^n D_j q^(-j),
X=sum_(j>=1) D_j q^(-j),
c=1/(q-1).                                              (4)
```

Thus `X` is uniform on `[0,1]`; the two base-`q` representations of a
countable set of endpoints are immaterial.  Since `h(S_n)` is a bounded
martingale, it converges almost surely and in `L^1`.  Write

```text
H_infinity=lim_(n->infinity) h(S_n).                    (5)
```

## 2. Logarithmic-density formula

The logarithmic density of `S_A` exists and equals

```text
delta_log(S_A)
 :=lim_(N->infinity) (1/log N) sum_(m<=N, m in S_A) 1/m
 =1/log(q) * E[H_infinity/(c+X)].                       (6)
```

Restrict the chain to states reachable from `s_0`.  If `C` is a closed
irreducible class with stationary law `pi_C`, put

```text
mu_C=pi_C(f).                                           (7)
```

The chain almost surely enters one closed class, and `H_infinity=mu_C` on
the event that it enters `C`.  In particular, if there is exactly one
reachable closed class, then

```text
delta_log(S_A)=pi_C(f).                                 (8)
```

No aperiodicity assumption occurs in (8).  With several closed classes,
the absorption probability alone is not enough: (6) weights each basin by
where it lies in the ordered shortlex `q`-adic interval.

## 3. Exact level normalization

The harmonic mass of accepted words on level `n` is

```text
L_n=sum_(|w|=n, w accepted) 1/I(w).
```

Since `r(w)=q^n X_n` and

```text
c_n=C_n/q^n=c+(q-2)q^(-n)/(q-1),                       (9)
```

uniform words of length `n` give the exact identity

```text
L_n=E[f(S_n)/(c_n+X_n)].                               (10)
```

The coupled tail satisfies `0<=X-X_n<=q^-n`, while
`0<=c_n-c<=q^-n`; hence

```text
|L_n-E[f(S_n)/(c+X)]| <= q^-n/c^2.                     (11)
```

The displayed bound is more than is needed; uniform convergence to zero is
the load-bearing fact.  For the full language, (10) tends to

```text
integral_0^1 dx/(c+x)=log q.                            (12)
```

## 4. Markov--Cesaro and martingale proof

Let `g(x)=1/(c+x)` and first replace `g(X)` uniformly by the cylinder
function `g(X_k)`.  For `n>=k`, the Markov property gives

```text
E[g(X_k)f(S_n)]
 =E[g(X_k)(P^(n-k)f)(S_k)].                             (13)
```

For fixed `k`, Cesaro averaging in `n`, with the first `k` terms discarded,
and (3) give

```text
lim_(D->infinity) (1/D) sum_(n=0)^(D-1) E[g(X_k)f(S_n)]
 =E[g(X_k)h(S_k)].                                     (14)
```

Now `g(X_k)->g(X)` uniformly and `h(S_k)->H_infinity` in `L^1`.  Letting
`k->infinity` in (14), then using (11), yields

```text
lim_(D->infinity) (1/D) sum_(n=0)^(D-1) L_n
 =E[g(X)H_infinity].                                   (15)
```

The last index through level `D` is `(q^(D+1)-1)/(q-1)`, whose logarithm is
`(D+1)log q+O(1)`.  A partial level has harmonic mass at most the full level,
and, uniformly in `n`,

```text
sum_(m=C_n)^(C_n+q^n-1) 1/m <= 1+log q.                (16)
```

Indeed, isolate the first summand and compare the rest with an integral;
also `q^n/C_n<=q-1`.  Thus a partial final level contributes `O(1)/log N`,
which proves (6) for every cutoff `N`, not only level endpoints.

Finally, on a closed irreducible class the Cesaro projection of `f` is the
constant stationary mean (7).  At a transient state it is the absorption-
weighted combination of these class means.  This proves the recurrent-class
interpretation and (8), including periodic classes.

## 5. Exact cylinders and the vector functional equation

For a length-`k` word `u`, let `r(u)` be its lexicographic rank and
`s_u=delta(s_0,u)`.  Bounded martingale convergence and
`h(S_k)=E[H_infinity | D_1,...,D_k]` give

```text
delta_log(S_A)
 =1/log q * lim_(k->infinity) sum_(|u|=k) h(s_u)
    log((c+(r(u)+1)/q^k)/(c+r(u)/q^k)).                 (17)
```

All `h(s)` are rational, so each finite approximant is a rational linear
combination of logarithms of rational numbers.  If `delta_k` denotes the
depth-`k` expression, then the explicit martingale error certificate is

```text
|delta_log(S_A)-delta_k|
 <=(q-1)/log q * E|H_infinity-h(S_k)|.                 (18)
```

For a start state `s` and `t>0`, define

```text
F_s(t)=E_s[H_infinity/(t+X)].                           (19)
```

Conditioning on the first digit proves the exact vector equation

```text
F_s(t)=sum_(d=0)^(q-1) F_(delta(s,d))(qt+d),
F_s(t)~h(s)/t as t->infinity.                           (20)
```

The desired coefficient is `F_(s_0)(1/(q-1))/log q`.

## 6. Kraft convergence criterion for every language

This subsection does not assume regularity.  Let `L` be any `q`-ary
language, let `a_n` be its number of accepted words of length `n`, and let
`L_n` be its level-`n` harmonic mass under the same shortlex indexing.  Since
the level is exactly the integer interval `[C_n,C_(n+1)-1]`, monotonicity of
`1/m` gives

```text
(q-1)a_n/(q^(n+1)-1) <= L_n
                     <= (q-1)a_n/(q^n+q-2).            (21)
```

For `n=0`, both sides equal `a_0`; there is no exceptional correction.
Uniformly in `n`, (21) implies

```text
(q-1)/q * a_n/q^n <= L_n <= (q-1)a_n/q^n.              (22)
```

Consequently

```text
sum_(m in S_L) 1/m converges
 iff sum_(n>=0) a_n/q^n converges.                      (23)
```

Thus harmonic convergence is exactly convergence of the Kraft-type series,
for arbitrary languages.  Equation (23) does **not** assert that arbitrary
languages have logarithmic density.

## 7. The regular spectral-radius dichotomy

Return to a finite automaton.  Trim its transition-count matrix to states
that are both reachable from `s_0` and coaccessible to an accepting state;
call the resulting nonnegative integer matrix `M`, and put `rho=rho(M)`.
Then

```text
rho<q
 iff no reachable closed SCC contains an accepting state
 iff sum_(m in S_A)1/m converges
 iff delta_log(S_A)=0.                                 (24)
```

If the trimmed graph is empty, take `rho=0`.  To prove the spectral step,
each row sum of every irreducible diagonal block of `M` is at most `q`.  If
such a block has spectral radius `q`, multiply the row-deficit vector by a
strictly positive left Perron vector: every deficit must vanish.  Hence all
`q` labelled edges stay in that block.  Coaccessibility then forces an
accepting state inside it, so it is a reachable closed accepting SCC.
Conversely, any such SCC has every row sum `q` and therefore spectral radius
`q`.

Accepted paths are counted by `a_n=e_(s_0)^T M^n f`.  Thus `rho<q` makes
`sum a_n/q^n` converge, proving the convergent half via (23).  In the other
direction a reachable closed accepting class has positive stationary
accepting mean and positive absorption probability; the integrand in (6) is
positive on a set of positive measure.  Therefore

```text
rho=q  ==>  delta_log(S_A)>0,                           (25)
```

and the harmonic subseries diverges as
`delta_log(S_A) log N+o(log N)`.  Regular shortlex languages consequently
have no intermediate divergent-but-zero-logarithmic regime.

This gives a reusable entropy gate: a genuinely regular sub-full language
with `rho<q` has finite harmonic mass.  Application to Rule 30 or another
system requires an actual finite-state shortlex carrier; entropy resemblance
alone is not such a map.

## 8. Sharp hostiles and exact audit

1. The binary even-length language has one irreducible class of period two.
   Its ordinary counting density has complete-level subsequential limits
   `2/3` and `1/3`, while (8) gives logarithmic density `1/2`.
2. The binary prefix-zero language has two absorbing classes.  Its accepting
   basin is `X in [0,1/2)`, so its coefficient is
   `log(3/2)/log 2`, not the unweighted absorption probability `1/2`.
   Reversing alphabet order gives `log(4/3)/log 2`; the two sum to one.
3. A delayed two-basin transient keeps the error in (18) genuinely nonzero:
   the independent audit shrinks its certified radius from
   `0.0901684400556` at depth four to `0.000352220468967` at depth twelve.
4. A transient split into two period-two recurrent classes verifies that
   reducibility and periodicity may occur simultaneously.
5. THM-3497's variable Berggren automaton has one class, `12` states and `4`
   accepting states, so (8) gives `1/3`.  Its fixed automaton has one
   period-two class, `192` states and `34=16+18` accepting states, so (8)
   gives `17/96`.

The independent companion does not import the submitted program.  It builds
closed SCCs, exact stationary laws and exact absorption probabilities from
scratch.  It exhausts all `729*8=5,832` complete binary three-state languages
and all `64*4=256` complete ternary two-state languages.  For every one it
checks `Ph=h`, recurrent-class interpretation, finite Cesaro convergence,
nested cylinder error intervals, language/complement duality, alphabet
reversal, every instance of (21) through the audit depth, and the graph form
of (24).  It separately verifies (20), arbitrary partial cutoffs, both
Berggren controls, and the delayed-transient hostile.  Its exact semantic
ledger is

```text
046420bde716c6f83c7314c0010e1e4b19feffae59de7c9684fb4d0bdc4b1459.
```

## 9. Scope and reproduction

```bash
python3 04-computation/regular_shortlex_harmonic_density_automaton_probe_20260816.py
python3 -O 04-computation/regular_shortlex_harmonic_density_automaton_probe_20260816.py
python3 04-computation/regular_shortlex_log_density_independent_audit_20260816.py
python3 -O 04-computation/regular_shortlex_log_density_independent_audit_20260816.py
```

The theorem concerns regular languages under the declared shortlex indexing.
Every subset of the positive integers can be viewed as a labelled harmonic
subseries, but arbitrary subsets need not have logarithmic density.  The
scalar coefficient in (6) also does not recover word ancestry, multiplicity,
owners, tournament arcs, LRC current, Jacobian flux, or a Rule-30 carrier.
THM-3510 gives the sharp binary boundary: identical counts on every positive
level can yield either of two regular logarithmic densities, or no density at
all when the within-level address is alternated in superdominant stages.

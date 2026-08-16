# Independent audit of regular shortlex logarithmic density

**Verdict: ACCEPTED and promoted as THM-3499.**

The auditee was the pair of remote candidate commits `405f108bd` and
`9fe56bc56` on
`origin/codex/session-turnpike-atlas-20260815-codex2`.  They were imported
into a fresh worktree based on current `origin/main` only after the bounded
startup and inheritance pass.  The reserved theorem stub was never used as a
proved dependency.

## Claim decomposition

The audit separated six load-bearing claims:

1. the empty-word-at-one shortlex normalization is
   `I(w)=1+(q^n-1)/(q-1)+r(w)`;
2. level harmonic mass is the exact expectation
   `E[f(S_n)/(c_n+X_n)]`;
3. finite Markov Cesaro projection plus the bounded martingale `h(S_n)` gives
   the joint recurrent-class/address coefficient;
4. a partial final level is uniformly `O(1)`, so the result holds at every
   integer cutoff;
5. the cylinder formula and vector dilation equation have the stated factors
   and no missing power of `q`;
6. transient, reducible and periodic chains do not require a hidden mixing or
   aperiodicity assumption.

Failure of any item would have left THM-3499 reserved.

## Independent analytic reconstruction

A length-`n` level starts at

```text
C_n=1+(q^n-1)/(q-1)
```

and has `q^n` entries.  Dividing the reciprocal denominator by `q^n` gives

```text
c_n=C_n/q^n=1/(q-1)+(q-2)q^(-n)/(q-1).
```

For uniform digits, the lexicographic rank divided by `q^n` is the truncated
base-`q` address `X_n`.  Thus the level mass is exactly

```text
L_n=E[f(S_n)/(c_n+X_n)].
```

Both `c_n-c` and the address tail `X-X_n` lie in `[0,q^-n]`.  Their
difference therefore has absolute value at most `q^-n`, which supplies the
uniform denominator replacement used in the proof.

For `P(s,t)=#{d:delta(s,d)=t}/q`, let `h=Pi f`, where `Pi` is the finite
Cesaro projection.  Then `Ph=h`, so `h(S_n)` is a bounded martingale.  For a
fixed address cylinder of length `k`,

```text
E[g(X_k)f(S_n)]
 =E[g(X_k)(P^(n-k)f)(S_k)].
```

Cesaro averaging first in `n`, followed by bounded martingale convergence in
`k`, gives

```text
lim_D D^-1 sum_(n<D)L_n=E[H_infinity/(c+X)].
```

The last index through level `D` has logarithm `D log q+O(1)`.  Any partial
level has mass at most

```text
sum_(m=C_n)^(C_n+q^n-1)1/m <= 1+log q,
```

so it cannot affect the normalized limit.  This closes the arbitrary-cutoff
gate, not merely the complete-level subsequence.

On a closed irreducible class, `h` is the constant stationary acceptance
mean.  At transient states it is the absorption-weighted combination of
class means.  Hence `H_infinity` is exactly the mean of the class eventually
entered.  If there is only one reachable closed class, it is constant and the
uniform address integral contributes `log q`; periodicity is irrelevant.

## Cylinder and functional-equation factor audit

At depth `k`, replacing `H_infinity` by its conditional expectation
`h(S_k)` and integrating `1/(c+x)` on each rank cylinder gives

```text
1/log q sum_(|u|=k) h(s_u)
 log((c+(r(u)+1)/q^k)/(c+r(u)/q^k)).
```

Its error is at most

```text
(q-1)/log q * E|H_infinity-h(S_k)|.
```

For `F_s(t)=E_s[H_infinity/(t+X)]`, conditioning on digit `d` uses
`X=(d+X')/q`.  The probability `1/q` cancels the reciprocal denominator's
factor `q`, leaving exactly

```text
F_s(t)=sum_d F_(delta(s,d))(qt+d).
```

The independent companion verifies this identity already at finite cylinder
depth, with maximum floating error `3.33066907388e-16` on three structurally
different controls.

## Kraft strengthening and sharp spectral boundary

The user proposed a useful corollary during the audit.  It is sound and was
included in THM-3499.

For an arbitrary language with `a_n` accepted words on level `n`, every
accepted reciprocal lies between the two endpoint reciprocals.  Therefore

```text
(q-1)a_n/(q^(n+1)-1)
 <= L_n
 <= (q-1)a_n/(q^n+q-2).
```

At `n=0`, both sides are exactly `a_0`; no separate correction is needed.
The uniform comparison

```text
(q-1)/q * a_n/q^n <= L_n <= (q-1)a_n/q^n
```

proves that the harmonic subseries converges if and only if the Kraft-type
series `sum a_n/q^n` converges.  This part needs no regularity and makes no
logarithmic-density assertion for arbitrary languages.

For a finite automaton, trim the transition-count matrix to reachable and
acceptance-coaccessible states.  Every irreducible block has row sums at most
`q`.  If a block has spectral radius `q`, a positive left Perron vector paired
with the nonnegative row-deficit vector forces every row deficit to vanish.
The block is therefore closed.  Coaccessibility forces an accepting state in
it.  Conversely, a reachable closed accepting SCC has row sums `q` and
spectral radius `q`.

This proves the exact dichotomy:

```text
rho<q  <=>  finite harmonic mass  <=>  logarithmic coefficient zero;
rho=q  <=>  a reachable closed accepting SCC  =>  positive coefficient.
```

Thus regular shortlex harmonic subseries have no intermediate
divergent-but-zero-logarithmic regime.

## Independent hostile computation

The new companion does not import the submitted script.  It independently
implements:

- reachable SCC decomposition;
- exact rational stationary distributions;
- exact rational absorption probabilities;
- exact `Ph=h` and recurrent-class means;
- finite Cesaro projection controls;
- rank-cylinder integrals and martingale error intervals;
- alphabet reversal and language complement;
- direct partial-level harmonic sums; and
- the trimmed full-SCC form of the spectral boundary.

It exhausts all `3^6=729` complete binary three-state transition tables and
all eight accepting sets, for `5,832` languages.  The table census contains
`190` transient, `2` multiple-closed-class and `34` periodic tables.  It also
exhausts all `2^6=64` complete ternary two-state tables and four accepting
sets, for `256` languages.  Every Kraft endpoint inequality is checked
through the declared finite depth.

The binary census splits into `1,556` convergent and `4,276` positive-log
languages; the ternary census splits into `79` and `177`.  The strongest
certified failure of unweighted basin mass is
`0.0849625007212`, and alphabet reversal changes a coefficient by a
certified `0.169925001442`.  These are the prefix-zero mechanism and its
reversal.

The small exhaustive banks resolve every multiple basin after finitely many
early digits, so the companion adds a four-state delayed split.  Its
martingale certificate radius is genuinely positive and shrinks from
`0.0901684400556` at depth four to `0.000352220468967` at depth twelve.  A
separate transient split into two period-two classes simultaneously attacks
reducibility and periodicity.

The submitted controls were rebuilt through the independent engine:

- variable Berggren: `12` states, `4` accepting, one class, coefficient
  `1/3`;
- fixed Berggren: `192` states, `34=16+18` accepting, one period-two class,
  coefficient `17/96`;
- even word length: period two, coefficient `1/2` despite ordinary endpoint
  limits `2/3,1/3`;
- prefix zero/one: `log(3/2)/log 2` and `log(4/3)/log 2`.

Normal and optimized runs of both submitted and independent companions agree
with their stored outputs.  The independent exact ledger is

```text
046420bde716c6f83c7314c0010e1e4b19feffae59de7c9684fb4d0bdc4b1459.
```

## Connection contract: sub-full finite-state dynamics

```text
source: a proved finite-state q-ary word carrier
target: its shortlex harmonic subseries
map: word -> shortlex index
preserved predicate: rho<q iff finite harmonic mass
lost information: ancestry, geometry, owners, current and physical time
required sidecar: the actual ordered DFA and acceptance set
decisive test: trim reachable/coaccessible states and compute rho/full SCCs
```

This can test a Rule-30 staircase language only after a genuine regular
carrier is supplied.  A shared entropy number or visual staircase is not a
bridge.

## Scope boundary

The theorem covers regular shortlex subsets and gives a convergence
criterion for arbitrary languages.  It does not say arbitrary subsets of
the positive integers have logarithmic density.  It also does not recover
Berggren ancestry from state counts, produce tournament arcs, construct an
LRC current or Jacobian flux, or identify a Rule-30 transducer that has not
been proved finite-state.

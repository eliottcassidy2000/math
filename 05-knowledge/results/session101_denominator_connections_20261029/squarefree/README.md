# Denominator resonance, squarefree return clocks, and the 115/133 separation

**Status: PROVED elementary statements and fixed-word sieve consequences;
FINITE-EXACT modular certificates.** Independently audited research packet. No conclusion about general Collatz convergence or an infinite
all-squarefree orbit follows.

## 1. Inheritance and source honesty

The inherited theorem is §3–6 of
`05-knowledge/results/collatz_guards_20260921_squarefree.md`, pinned by
`git show e871dd4513c:<path>`; SHA256
`a7727f346a5ff53d34c94d1a788c0810b28c0a203c1645f7d80be0753000df0f`.
The inherited claim is:
positive odd sources `n≡5 (mod9)` realizing **r repetitions** of exponent word
`(1,1,5)` with **every one of the 3r+1 nodes, including the terminal node,
squarefree** exist iff `r≤4985`; for each such fixed r they have positive
relative density in the exact word-and-row progression. The source does not
assert a universal forbidden length for arbitrary words: arbitrarily long
all-one growth words still have all-squarefree realizations.

The relevant concept board is: (i) reduced rational closure denominator;
(ii) carry-sensitive append cancellation; (iii) affine returns modulo prime
squares; (iv) phase-zero anchors and terminal endpoint; (v) all-prime local
avoidance plus the fixed-word squarefree sieve. The hostile to (i) implying
(iii) is `(1,13)` below; the hostile to denominator determining (ii)/(iv) is
115 versus133. The inherited corrected near miss is omitting the terminal
endpoint, which would move the115 obstruction from4986 to4987 blocks.

## 2. The exact denominator-prime mechanism, with its missing hypothesis

Let `F(x)=(a x+b)/d` over `Z_p`, where p is odd, `p∤d`, and

` s=v_p(d-a)>t=v_p(b)≥0 `.

**Proposition (PROVED).** For every `x∈Z_p` and integer `j≥1`,

`v_p(F^j(x)-x)=t+v_p(j)`.

Indeed put `α=a/d`, `β=b/d`. The term
`F(x)-x=(α-1)x+β` has valuation exactly t, while `α≡1 modp`.
The factorization

`F^j(x)-x=(F(x)-x)(1+α+⋯+α^(j-1))`

and the odd-prime lifting identity give valuation `v_p(j)` for the second
factor. If `d=a` is admitted using `v_p(0)=∞`, the same proof holds with
`α=1` and geometric sum exactly j. Thus modulo `p^e` the exact orbit length is
`p^max(e-t,0)`. If `e>t`, the orbits are exactly the individual congruence
classes modulo `p^t`; if `e≤t` the map is the identity. The assertion about
classes follows because F preserves each class and the orbit has its full
size. This proof includes the previous unit-carry result as t=0.

A prime in the reduced denominator of the fixed point `b/(d-a)` has `s>t`,
but does **not** force t=0. Therefore denominator membership is not enough
for a single cycle modulo `p²`. What is special about101 for115 and133 is
not just the denominator101: their carries19 and31 are both101-units.
Then F is a full cycle modulo every `101^e`, forcing a zero among any
`101²` consecutive block-boundary values. Intermediate nodes can strengthen
that coarse bound, and their zero positions are additional data.

**Hostile (PROVED).** Word `(1,13)` has

`F(x)=(9x+5)/16384`, `d-a=16375=5³·131`, `f=1/3275`.

Its reduced denominator contains `5²`, but modulo25 it is `x↦x+20`.
For sources `x≡1 mod5`, every block source stays1 mod5 and its intermediate
node `(3x+1)/2` stays2 mod5. Thus arbitrarily many repetitions have **no
local obstruction at5²**, despite that denominator prime. This does not
assert global all-squarefree existence for all lengths:131 is a different,
unit-carry resonant prime and does force a sufficiently long repetition
obstruction. This identifies the exact first failed implication rather than
claiming denominator primes have no relevance.

## 3. From affine cycles to endpoint-aware interval coverage

Suppose a guarded word returns by F, with t=0 at p, and its pre-final phase
maps are affine maps `P_j` with unit slopes modulo `N=p²`, including `P_0(x)=x`.
Write the unique full F-cycle as `x_i=F^i(0)` and let A be the **distinct**
indices whose cycle values are roots of some `P_j`. Sort them as

`0=a_0<a_1<⋯<a_(m-1)<N`.

Put `g_i=a_(i+1)-a_i` for `i<m-1` and `g_(m-1)=N-a_(m-1)`.
For r blocks the excluded source indices are the backward intervals
`a_j-{0,…,r-1}`, with the extra endpoint index `-r` (root0).
Therefore the exact number avoiding every p-square zero is

`h(r)=Σ_(i<m-1)(g_i-r)_+ +(g_(m-1)-1-r)_+`.

Proof: in each cyclic gap, the interval ending at its right anchor removes
r positions. The interval ending at0 removes r+1 rather than r, because
of the terminal block boundary. Empty gaps contribute zero. This argument
also handles overlaps and r beyond N, not just disjoint short intervals.
The last possible r is `max(g_0,…,g_(m-2),g_(m-1)-1)-1`.
This is a **local p² result** until other primes and a sieve are supplied.

### Explicit quadratic clock and inverse modulo p²

There is an algebraic certificate stronger than looking up the modular orbit.
Write `α=1+p c (mod p²)` and `β=b/d (mod p²)`, with β a unit. The
binomial expansion gives

`F^i(0)=β[i+p c·i(i-1)/2] (mod p²)`.

For a root z put `y=z/β (mod p²)`. Its unique cycle index is

`i=y-p c·y(y-1)/2 (mod p²)`.

To check the inverse, i and y agree modulo p, so the quadratic correction
has the same value modulo p² with either argument. Division by2 is allowed
because p is odd; the displayed triangular numbers are also integral.
Thus the zero-root positions can be computed directly without a discrete-log
search. For the101 example `α=8687=1+101·86`; β is4224 for115 and5818
for133. Their y values at roots3400,9067 are respectively8946,4709 and
218,7164; the inverse formula gives the indices in the next table.
The checker compares both polynomial formulas with **every** orbit position
for both words, not just the anchors. The formula is scoped to p²;
higher powers require additional binomial terms.

## 4. Same101 denominator, different sharp squarefree thresholds

Both words have length3, exponent sum7 and slope27/128:

| Word | Carry | Fixed point | Pre-final phases |
|---|---:|---|---|
|115|19|19/101|`x,(3x+1)/2,(9x+5)/4`|
|133|31|31/101|`x,(3x+1)/2,(9x+5)/16`|

Their phase zero roots modulo10201 are the **same set**
`{0,3400,9067}`. The denominator change4→16 does not change a zero root.
But changing the return carry changes the time at which each root is reached:

| Word | Ordered root indices | Cyclic gaps | Endpoint-adjusted capacities |
|---|---|---|---|
|115|0,3795,5214|3795,1419,4987|3795,1419,4986|
|133|0,8601,8780|8601,179,1421|8601,179,1420|

Thus

`h_115(r)=(3795-r)_+ +(1419-r)_+ +(4986-r)_+`,

`h_133(r)=(8601-r)_+ +(179-r)_+ +(1420-r)_+`.

The unique surviving source residue at the last possible r is7332 mod10201
for115 at r4985, and5818 mod10201 for133 at r8600. The next r has no
surviving residue. These are exact modular certificates, not enumerated huge
integer trajectories. The checker validates all h(r) for0≤r≤10201 by
incremental root-union coverage, checks anchors independently by binary affine
exponentiation, and directly evaluates every phase of each boundary survivor.

**New theorem (PROVED).** For each integer r≥1 the following are equivalent:

1. A positive odd source `n≡5 mod9` realizes r copies of133 and has all3r+1
   displayed nodes squarefree.
2. Such sources have positive relative natural density in the exact
   word-and-row progression.
3. `r≤8600`.

The local101-square obstruction proves1⇒3. To establish3⇒2 requires all
other primes, not only the longer prime-square clock.

At the rational fixed point f=31/101 the phase values are
`31/101,97/101,49/101`. For each prime `p≥5`, `p≠101,7`, these give an
avoiding residue modulo p²: none of the three numerators is divisible by p².
In particular31 and97 are harmless simple zeros modulo p, not square zeros.
At7 choose `n=f+7 mod49`, namely residue1. With `α=27/128`,
`F^t(n)=f+7α^t`; the first two phases stay units modulo7, and the last is

`49/101 +(9/16)·7α^t`,

which has valuation exactly1. So this residue avoids every7-square zero
for arbitrarily many blocks (its modular orbit has period6). For comparison,
the115 fixed phases are19/101,79/101,169/101, and its only analogous
exception13 is repaired by f+13, residue50 mod169, period12.

Finally fix r, let `K_j` be the cumulative exponent in the3r-step word, and
use its exact source cylinder intersected with `n≡5 mod9`:

`n=a_r+Q_r s`, `Q_r=9·2^(7r+1)`.

This progression exists by CRT; a positive representative has all its
prescribed valuations exactly. Its displayed integer nodes are positive
linear forms `c_j+A_j s`, with

`A_j=9·3^j·2^(7r+1-K_j)`, `0≤j≤3r`.

They are all odd; the initial row and recurrence ensure every node is prime
to3. For `p≥5`, each coefficient is invertible modulo p² and the progression
parameter samples every source residue. Write ν_r(p) for the number of
forbidden residues. The above local proofs give `ν_r(p)<p²` for r≤8600;
also `ν_r(p)≤3r+1`. Hence the convergent product

`δ_r=∏_(p≥5)(1-ν_r(p)/p²)>0`

has101-factor precisely `h_133(r)/10201`.

For rigor beyond finite CRT compatibility: in parameter range0≤s<X, square
divisors of these fixed linear forms have `p=O_r(√X)`. The number removed
by primes p>z is at most

`(3r+1)X·Σ_(p>z)p^(-2)+O_r(√X)=O_r(X/z+√X)`.

Finite CRT gives the truncated product. Divide the tail by X, send X→∞,
then z→∞. This proves the density δ_r and hence2⇒1. Relative to **all**
positive odd sources in the5mod9 row, the density is `2^(-7r)δ_r`.
All limits hold with r fixed; no useful uniform lower bound is claimed.

## 5. Comparison with appended-step cancellation spectra

Appending exponent v does **not** keep iterating the same F: it changes the
closure equation. Its new fixed point is

`x_v=(3C+128)/(128·2^v-81)`.

For115 the numerator is185; for133 it is221. At v=1 the common unreduced
denominator is175, so the reduced fixed points are37/35 and221/175.
The corresponding reduced denominators35 and175 distinguish exactly the
same pair that the squarefree thresholds distinguish as4985 and8600.

There is a genuine common source coordinate—the ordered carry C—but the
two mechanisms are not identical:

* Repetition keeps `(a,d,b)=(27,128,C)` and follows an affine clock on
  `Z/101²`. It asks which block times hit **intermediate zero roots**.
* Appending changes the slope gap to `128·2^v-81`, and compares the new
  carry `3C+128` against its divisors. It asks about **rational cancellation**.
* The map from a word to its denominator101 forgets C. Retaining slope,
  carry, labelled phase maps, and endpoint convention restores the return
  predicate. Retaining the append formula restores cancellation responses.
* Even the unordered root set loses the return-time positions:115 and133
  have identical root sets but different affine translations, hence different
  maximum gaps. This is stronger than merely saying clocks alone lose carry.

No transfer of a squarefree obstruction to a general convergence theorem is
valid: an arbitrary trajectory need not contain these repeated words.
Likewise a rational closure denominator is not an integer Collatz cycle.

## 6. Reproduction and stopping point

Run from this directory:

```sh
python3 check_squarefree.py > checks.json
python3 -O check_squarefree.py > checks-optimized.json
cmp checks.json checks-optimized.json
```

The checker uses only the standard library and explicit runtime gates.
It runs standalone by default. To additionally verify the inherited source
blob against a Git checkout, pass `--source-repo /path/to/math`; that optional
provenance check does not change the deterministic mathematical output.
It audits the inherited115 constants, new133 constants, endpoint boundary,
nonresonant exceptional repairs, the denominator-only hostile, and append
comparison. The general valuation finite universe is p∈{3,5,7},
1≤s≤3, 0≤t<s, 1≤e≤3, and every x modulo p^e; it supplements rather than
replaces the proof above. No enormous r8600 source is enumerated: existence
is the proven fixed-r sieve conclusion. All requested mechanisms and the
new sharp boundary are now explicit; classification of arbitrary word
squarefree clocks is not attempted here.

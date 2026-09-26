# The orbit-coupled potential on crossings, attacked: the real place of a divergent orbit is rigid (a golden-free rotation), the two-place clock, exits are ballot-thin, and why stays are not; with the 1/3-2/3 analogy made precise and the owner's three seeds decoded

**Status: REASSESSMENT with PROVED small statements (Propositions 1-3,
the exits lemma, the {2,3,11} proposition), FINITE-EXACT probes, SPECULATION
clearly marked (sections 2.4 and 3.6), independently audited with repairs
recorded in section 6 and the integration audit. No Collatz claim. Session
`collatz-crossings-20260926` (opus), 2026-09-26.**

Scripts: `04-computation/experiments/collatz_crossings_20260926_phase.py`
(two-place identity, phases, crossings), `collatz_poset_20260926_balance.py`
(balance constants of the no-descent posets), `collatz_wythoff_20260926_interlock_probe.py`
(the decoded three-colouring against the Collatz step); outputs `.out`
alongside.

## 1. The owner's seeds, decoded

**1.1 A near-fitting candidate is the Wythoff colouring.** The user chose
the full listed sequence as authoritative. The rule below matches its
first 34 colours but fails at 35; it is a tested candidate, not a unique
decoding or a correction to the data. For this candidate, at index `n`,

```text
red   = AA = { floor(floor(k phi) phi) } = { n : the Zeckendorf representation of n ends in F_2 = 1 },
black = B  = { floor(k phi^2) }          = { n : the lowest Zeckendorf index of n is odd (>= 3) },
blue  = AB = { floor(floor(k phi^2) phi) } = { n : the lowest Zeckendorf index of n is even (>= 4) },
```

which partition the positive integers (`N = A + B`, `A = AA + AB`), with
densities `phi^(-2), phi^(-2), phi^(-3)` (`0.382, 0.382, 0.236`). The
owner's list agrees with this for `n <= 34`; it puts `35` in blue where `AB`
continues `..., 29, 32, 37, 42, 45, 50` and `35 = 34 + 1` is red (`AA`).
The blue gaps `5, 3, 5, 5, 3, 5, 3, 5, 5, ...` are the Fibonacci word. So the
"couple of bits of parity hiding behind the natural numbers" are exactly
the two bits (`A`/`B`, then `AA`/`AB`) of the lowest Zeckendorf digit: the
golden coding of `n`. Under doubling the colours move rigidly (`black ->
red or black`, `blue -> red or blue`, each `1/2`; `red -> 0.19/0.50/0.31`;
mutual information `0.369` bits), because the colour is a function of
`x={n phi}` (black for `x<phi^(-2)`, red between `phi^(-2)` and
`1-phi^(-3)`, blue above). The identity `{2n phi}={2x}` gives the exact
red row `(phi^(-2)/2,1/2,phi^(-1)/2)`, as the incoming audit verifies.
The finite Collatz probe finds small values of the mutual
information of the colour of `n` with parity, `v_2(3n+1) mod 3`, descent
within `20` steps and stopping time `mod 3` is `< 1e-5` bits at
`N = 300000`, and the colour of the odd part `U(n)` given the colour of `n`
carries `0.005` bits. These census-derived statistics do not prove
independence or exclude every arithmetic coupling. Exact colour carries
and finite-closure hostiles are in
[crossroads_crossing_20260926_colour](crossroads_crossing_20260926_colour.md).

**1.2 `{2, 3, 11}` (PROVED).** The three inequalities `4, 6, 8 < 9`, `6 < 9`,
`22 < 25` say: `2p` is below the next odd square above `p`. **These are the
only such primes.** If `(2j-1)^2 < p <= (2j+1)^2` then `2p < (2j+1)^2`
requires `2(2j-1)^2 < (2j+1)^2`, i.e. `4j^2 - 12j + 1 < 0`, i.e. `j <= 2`:
`j = 1` gives `p in {2, 3}` (`2p < 9`), `j = 2` gives `p = 11` (`2p < 25`;
`13` fails); for `j >= 3` the ratio of consecutive odd squares is below `2`
(`49/25 = 1.96`). Checked to `10^4`. The multiples below the square are
`{4, 6, 8}` (one of them the square `4`), `{6}`, `{22}`. In the language of
this thread: `2^3 < 3^2` (`8 < 9`) is the first `2`-adic/`3`-adic near-miss,
the one that makes `3n+1` with three halvings a descent; `22 < 25` pairs
`11` with `5`, and `139 = 3^7 - 2^11` pairs the exponents `7, 11`. A Langlands
reading was asked for: `11` is the conductor of the first elliptic curve
(`X_0(11)`), and the pasted theta series `Theta(t) = theta_3(e^(-pi t))^2 = sum r_2(k) e^(-pi k t)`
is a weight-`1` modular form whose coefficients count sums of two squares
(primes `1 mod 4`). Its modularity `Theta(1/t) = t Theta(t)` is the
micro/macro duality the owner names; nothing here connects it to the
Collatz map, and I do not claim it does.

**1.3 Parity flow.** `3n` keeps parity, `+1` flips it, the halvings remove
the evenness that multiplication keeps: the Collatz step is one additive
parity flip between two multiplicative operations, and `v_2(3n+1)` is the
multiplicative side's memory of the additive flip. The quadratic maps the
owner lists (`x^2` with its loops at `0, 1` fed by `-1`; `x^2 - 1` with
`1 -> 0 <-> -1`; `x^2 - 2` with `1 -> -1`, `0 -> -2 -> 2`) are three separate integer dynamical graphs. For `x^2-2` specifically,
the Chebyshev coordinate gives angle doubling (`x = 2 cos t`, `x^2 - 2 = 2 cos 2t`): the fixed
point `-1 = 2 cos(2 pi/3)` and the chain `0 -> -2 -> 2` are the angles
`pi/2 -> pi -> 0`. The log-periodic functions pasted
(`cos^2((pi/2){log_2 x})`, `(pi/x)(3/2 + cos(2 pi log_2 x)/2)`) are functions of
the phase `{log_2 x}`; section 3.1 shows why that phase is the natural
coordinate of a Collatz orbit.

## 2. The 1/3-2/3 conjecture, made precise for Collatz

**2.1 The construction (THM-4503, crossroads-poset, section 4 of its bridge
note; found independently here).** A parity word of length `k` with `o` odd
letters `O_1 < ... < O_o` and `e = k - o` halvings `H_1 < ... < H_e` has
all prefix multipliers `3^(#odd)/2^j > 1` iff the `j`-th halving is
preceded by at least `f(j)` odd letters, `f(j) = min{o' : 3^(o') > 2^(o'+j)} = floor(j log_(3/2) 2) + 1`
(`f = 2, 4, 6, 7, 9, 11, 12, 14, ...`). These positive-multiplier words are a sufficient subset of actual
no-descent words; additive carry can retain others. So `Bad_k(o)`, these words with
`o` odd letters, is the set of linear extensions of the width-2 poset
`P(o, e)`: two chains plus the cross relations `O_(f(j)) < H_j`. Linial's
theorem (the 1/3-2/3 conjecture is proved for width 2) then gives, for every
nonempty `P(o, e)` that is not a chain, an incomparable pair `(O_i, H_j)`
with `1/3 <= P(O_i before H_j) <= 2/3` under the uniform law on these
positive-multiplier words.

**2.2 Finite balance data and coding scope (corrected).** The exhaustive
probe covers `4<=k<=20`, with `|Bad_20|=27328`. Every nonchain cell has
balance at least `1/3`; the three-word cells `(k,o)=(5,4),(6,4)` attain it.
Within this census, cells with at least30 extensions have balance at least
0.441. Some cells attain1/2, but this is not true of every large cell or
every level. At k=5 the best is1/3, at k=7 it is3/7, and `(k,o)=(10,7)`
has30 extensions with balance7/15. The stored output already shows these
exceptions; there is no all-level balance-limit theorem here.

Summing cell counts gives `W_k`, whose Spitzer generating function and
asymptotic are THM-4495. Thus `log_2 W_k=h k-(3/2)log_2 k+O(1)` describes
the total, not every individual `e(P(o,e))`. Enumerative coding uses
`ceil(log_2 W_k)` bits. A comparison tree needs a separate guarantee:
a `[1/3,2/3]` split gives at most `ceil(log_(3/2) e(P))` questions, not
necessarily the information-theoretic minimum. This argument therefore
does not establish the asserted shorter order-comparison code or a
compression contradiction for one fixed integer.

**2.3 Where the golden ratio enters both problems.** The classical general
balance bound `(5-sqrt5)/10` and Fibonacci/Wythoff constructions invite
comparison. The attached newer manuscripts are audited separately in
[crossroads_poset_20260926_barrier](crossroads_poset_20260926_barrier.md).
No map identifying their probability laws with the arithmetic source law
is established here. Collatz's prefix-slope constraint involves `log_2 3`,
and its finite balance constants vary. A shared irrational constant or a
shared interleaving vocabulary is not a proved transfer mechanism.

**2.4 Speculation on adjacent constructions (marked as such).** The poset
lane's next probes (its bridge note, section 7) already ask for a
balanced comparison with a *descent-relevant consequence*; the balance
computation finds many balanced cells, without proving perfect balance
at every level. The additional difficulty remains: an order
question about the parity word says nothing about the *height-selected*
integer that realises it (the arithmetic carry, THM-4503's sidecar). An
adjacent construction that could matter would take the poset not on the
word but on the *orbit*: the odd iterates `m_l` of one orbit ordered by
value, with the relations `m_l < m_(l')` forced by the descent structure;
its linear extensions are the possible "time orders of a set of heights",
and a balanced pair there would be a pair of times whose height order is
undetermined by the descent data. Whether such a pair exists for a
divergent orbit is a statement about the orbit's oscillation, i.e. exactly
the missing oscillation lemma of section 3.5 (HYP-9161 in its crossing
form). I do not know how to prove it and record it as the one place where
the analogy points at the right object.

## 3. The orbit-coupled potential on crossings

Setting for Propositions1--3: a positive `3n+1` orbit with odd iterates `m_0 = n, m_1, ...`,
valuations `v_l`, `d_l = v_0 + ... + v_(l-1)`, `alpha = log_2 3`,
`Delta_l = d_l - l alpha`. The exact two-place identity (Proposition 6 of the
squares/doubles note; Bernstein) is

```text
2^(d_l) m_l = 3^l n + S_(l-1),   S_(l-1) = sum_(j<l) 3^(l-1-j) 2^(d_j),   so   m_l = n 2^(-Delta_l) C_l,   C_l = prod_(j<l) (1 + 1/(3 m_j)) = 1 + R_l/n,
```

`R_l = S_(l-1)/3^l` (plus sheet; minus sheet with `-`). Checked exactly with
rationals on the orbits of `27`, `703`, `26623`.

**3.1 Proposition 1 (the real place is a rotation, in base 2).** For
every orbit with `sum_j 1/m_j < infinity`, in particular every orbit that
is not eventually periodic (THM-4476, Cor. 6), the phases `{log_2 m_l}` are
those of the irrational rotation `{log_2 n + l alpha}` shifted by the
convergent sequence `log_2 C_l`; hence they are equidistributed: the
base-2 significands of the odd iterates are Benford along the orbit. This
is a base-2 statement only: `log_10 m_l = log_10 Q_l + l log_10 3 - d_l log_10 2`
and `d_l log_10 2` is not an integer, so the base-10 phases follow the
orbit's own path `(l, d_l)` and need not equidistribute (the audit built
a valuation word confining them to `[0, 0.35]` while the base-2 phases
stay uniform); whether a divergent orbit can do that is open. *Proof.* `{log_2 m_l} = {log_2 n + l alpha + log_2 C_l}`
since `d_l` is an integer; `log_2 C_l` converges; a Weyl-equidistributed
sequence plus a convergent one is equidistributed. ∎ (Along the orbit of
`2^200 - 1`, `980` odd iterates, the discrepancy of the phases is `0.0046`
against the uniform-sampling scale `0.032`; on these finite `3n+1` orbits
the shift `log_2 C_l` is `0` to five decimals until the last twenty small
iterates, where all of its `0.32` accumulates, whereas the finite `5n+1` probe from `7` approaches `0.1026`.
An infinite-orbit conclusion for that map requires its own summability
premise.)
The related Benford literature has different quantifiers:
Kontorovich--Miller section5.2 distinguishes source-population limits
from one fixed trajectory. Here the conclusion is conditional on
`sum 1/m < infinity`, single-orbit, and in base2. It is vacuous for
`3n+1` if Collatz holds.
Its meaning for us: **all the dynamics of an orbit is in the integers
`d_l`**; the leading bits of `m_l` are a rigid clock, the trailing bits are
the free address.

**3.2 Proposition 2 (the two-place clock, corrected).**
`Q_l=m_l 2^(Delta_l)=n C_l=2^(d_l)m_l/3^l` increases. Under the
nonperiodic/summability premise it is bounded by THM-4476: `Q_0=n`,
and `n<Q_l<=n exp(K/3)` for l>=1. The additional premise `R(d)<n`
gives `Q_l<2n`. Its exact norms are

    |Q_l|_2=2^(-d_l),   |Q_l|_3=3^(l-v_3(m_l)).

The second is3^l for l>=1, since internal odd values are units modulo3,
and at l=0 only when3 does not divide n. The real increment is
`Q_(l+1)-Q_l=Q_l/(3m_l)`. With natural logarithms, summability gives

    3 log(Q_infinity/n) < sum_l 1/m_l <= 4 log(Q_infinity/n),

using `x/(1+x)<=log(1+x)<x` at `x=1/(3m_l)<=1/3`. The original
upper bound with coefficient3 reversed this logarithmic inequality.
The clock is exact, but its relative change per step is only1/(3m_l).
No uniqueness claim for all possible orbit-coupled potentials follows.

**3.3 Proposition 3 (the past is the `3`-adic address).** For `l` with
`2^(d_l) > Q_l` (Proposition 2; this holds once `d_l > log_2 n + K/(3 ln 2)`,
THM-4476 Cor. 6, and once `d_l > log_2 n + 1` whenever `R(d) < n`), `m_l` is
the least positive residue of
`2^(-d_l) S_(l-1)` modulo `3^l`: the odd iterate is determined by the
halving word alone, through `3`-adic reduction, exactly as the residue of
`n` modulo `2^k` is determined by the parity word (Terras). More locally,
`m_l mod 3^j` is a function of the last `j` valuations `v_(l-j), ..., v_(l-1)`.
So an odd iterate carries its future in its `2`-adic digits and its past
in its `3`-adic digits, both finitely (it is an integer), and the coupling
of an orbit is that the same integer has both. *Proof.* Reduce the identity
modulo `3^l`; `m_l < 3^l` because `m_l = 3^l Q_l/2^(d_l) < 3^l`. ∎

**3.4 Crossings (FINITE-EXACT).** On the record orbit of `63728127` (`358`
odd iterates) the odd iterates cross each dyadic level `2^t` downward `2`
to `5` times while `16` to `172` of them lie below it; the same on
`2^60 - 1`. So stays below a level are few and long; the number of
elements is governed by the *length* of stays, not their number. Between
successive visits of a band `(Y/2, Y]` the pairs `(M, D)` (odd steps,
halvings) satisfy the exact identity
`M alpha-D+C_seg=log_2(m_end/m_start) in (-1,1)`, where
`C_seg=sum log_2(1+1/(3m_l))` along the segment. Thus D is a floor or
ceiling of `M alpha+C_seg`. Dropping C_seg requires a separation condition
at the integer boundary, not merely an o(1) estimate.

**Exits lemma (PROVED for pairwise distinct orbit states).** Count completed
stays; a finite prefix may have one additional unfinished stay. The number
of completed stays below `X` of length at least
`k = floor(log_2 X)` is at most `2 M_k(1.6) = O(X^(h*) (log X)^(-3/2))` for `X >= X_0(b)`
(Lemma M of THM-4499; the incoming audit gives `X_0(1)=21`,
and the proof below permits the sharper barrier `log_2(5/3)`). *Proof.* A stay that lasts `>= k` steps and then
exits above `X` has an element `z <= X` exactly `k` steps before the exit;
`z`'s `k`-word satisfies `M_s z + beta_s <= X < M_k z + beta_k` for `s < k`, so
with `|beta_s|, |beta_k| <= |b|(3/2)^k <= X/4` once `X >= X_0(b)`, so `M_k z > 3X/4`
and `M_s z < 5X/4`, hence `M_k/M_s > 3/5 > 2^(-1.6)` for all `s < k` with no
condition on `z`: reversed, the word has all partial sums `> -1.6`, one of
`M_k(1.6)` classes with two representatives below `X`; distinct stays have
distinct `z`. ∎ Stays shorter than `k`
are not controlled by any count: their entries `z in (X/2, X]` with
`m`-step reversed-positive words number up to `X 2^(-(1-h) m) m^(-3/2)`, which
is `~ X` for small `m`.

**3.5 Why no count bounds the length of a stay, and what would.** Every
element of a long stay below `X` has a free window; the stay's word is a
walk with drift `-0.2075` conditioned only to return above `X` at its end,
and returning costs entropy only in the last `k` steps (the exits lemma).
Charging elements to exits, to leaders, to landing points (THM-4476,
THM-4506) or to band-visits (S7 note) reproduces the same `k` per charge.
The quantity that would break this is an **oscillation lemma for one
orbit**: a bound on the length of a stay below `X` in terms of the number
of leaders (future-minimum records) inside it, e.g. `length <= C (log X)^beta x (leaders inside + 1)`
with `beta < 1`. Leaders are ballot-thin (S7 Proposition, audited), so this
would give `N(X) = O(X^(h*) (log X)^(beta - 3/2))`; with `beta = 0` the ballot
floor. It is HYP-9161 in crossing form. Every count of residue classes is
blind to it (a residue class realises any finite stay), and every
potential built from the two-place identity moves by `O(1/m)` per step
(3.2). The reflection of the procgen lane names the needed object an
adaptive, unbounded, height-coupled quantity; 3.1-3.3 say precisely what
is *not* it: the phase (rigid), the carry clock (too slow), and any finite
piece of either address (free).

**3.6 A proof-idea attempt through the Pythagorean tree (SPECULATION,
with the candidate potential corrected).** Berggren's inverse moves
strictly decrease primitive-triple height. For odd Collatz one can instead
define a functional graph with parent `(m-1)/4` when `m=5 mod8`, and
`U(m)` otherwise, for m>1. In the latter case the valuation is1 or2.
After deleting the self-loop at1, its child set is just{5}; other odd
nonmultiples of3 have two children, multiples of3 one. Calling this a
rooted tree already requires the missing global termination result.

During a maximal run with valuation1, the exact invariant is

    log_2(m+1)+log_2(3/2) v_2(m+1).

The earlier expression with log_2(m) and a negative sign is not constant;
it strictly increases. At a reset the correct increment also includes
`log_2((m'+1)/(m+1))`, so it is not just a new-counter charge. Resets can
increase or decrease the invariant. The actual reset41->31 shows why
end-of-run information alone does not pay for the next growing run.

THM-4483 already excludes suitable linear valuation banks. The independently
audited [THM-4507](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md)
now excludes arbitrary nonlinear corrections built from finitely many
polynomial valuations at finitely many primes. This leaves selected-return
or variable-depth constructions as possible routes, without asserting
that every successful proof must have one specific form. The Pythagorean
comparison identifies the desired well-founded decrease; it does not
supply that decrease for Collatz.

## 4. Cross-lane notes

* THM-4506 (procgen, 14:30 today) sharpened the S7 one-bit band lemma to
  the exact worst case `ceil((k - D)/log_2 3)` and proved that
  `a*(mu) = mu lambda*/h* - 3/2` is exactly what the one-window recursion
  yields (saturation). HYP-9161 is their `mu = 0` case; nothing here
  contradicts it.
* THM-4503 (crossroads-poset) has the two-chain poset of 2.1 and reviewed
  a user-supplied manuscript claiming the Kahn-Saks conjecture (large
  width forces near-balance); as they note, the Collatz posets have width
  `2`, where the question is Linial's settled `1/3`, and section 2.2 shows
  finite constants vary, including both `1/3` and `1/2`.

## 5. What was not done

No theorem on crossings. The propositions of section 3 are elementary and
partly folklore; they are recorded because they say exactly which
orbit-coupled quantities are *not* the missing one. The seeds are decoded
and tested; the Wythoff colouring is a near-fitting golden-coordinate
candidate whose reported Collatz correlations are finite observations;
the `{2, 3, 11}` rule is
proved and has no dynamical content that I can see.

## 6. Independent audit (2026-09-26)

Auditor subagent, independent recomputation with exact arithmetic:
`04-computation/experiments/collatz_crossings_20260926_audit.py` ->
`05-knowledge/results/collatz_crossings_20260926_audit.out`. CONFIRMED:
the candidate Wythoff rule (`AA/B/AB` = the Zeckendorf rule to `4e5`, the
densities, the doubling matrix explained by `x = {n phi} -> {2x}`, the
mutual-information values), the `{2, 3, 11}` proposition to `10^6`, the
poset bijection (a linear-extension DP equals the word enumeration in
every cell to `k = 20`, sums to THM-4495's `W_k`), all `67` balance cells,
Propositions 1-3 (exact rationals on `27`, `703`, `26623`; the
`3`-adic locality with its closed form), the exits lemma mechanics
(`X_0(1) = 21`; `M_k(0.737)` suffices), the crossing counts, and the local parent/child counts of the functional
graph in3.6 to `10^5` (not global termination). CORRECTED (applied above): the "within `0.006` of
`1/2` for every `k`" sentence (true for `14 <= k <= 20` only); the reversed
reciprocal-sum inequality (`3 log <= sum 1/m <= 4 log`); "Benford's law"
narrowed to base 2 (the base-10 phases follow `(l, d_l)` and need not
equidistribute); the attribution to Lagarias-Soundararajan /
Kontorovich-Miller rephrased; the finite-orbit "shift converges to `0.32`"
(an artifact of the last twenty iterates); the `5n+1` band pairs of the
script (`log_2 5`); "Proposition 4" in the status line; the distinct-terms
hypothesis of the exits lemma; the child count of `1`. This records that audit's scope. Further demonstrated corrections to
source-law transfer, entropy/comparison bounds, carry-sensitive band
returns, the run potential, and the supplied colour sequence are applied
above and detailed in the [integration audit](crossroads_crossing_20260926_integration.md).
Neither finite enumeration nor the word-level base10 hostile constructs
a divergent positive orbit.

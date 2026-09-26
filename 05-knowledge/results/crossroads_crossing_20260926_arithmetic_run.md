# Maximal-rise resets: exact sibling coupling and its stopping boundary

**Status: PROVED elementary local corollaries of inherited reset and inverse-fibre identities; FINITE-EXACT independent controls. Collatz remains OPEN. No literature priority claim.**

## 1. Inheritance and the additional question

The maximal-run formula is already present, with more information, in
[the valve note, sections 3--4](collatz_guards_20260921_valves.md): exact
first-reset size, its integrality-improved descent boundary, and a spatial
descent density of `0.713725...`. The guarded maximal-block passport is
also inherited from [the transducer note, sections 3--4](creative_transducer_20260925.md).
The binary-repunit identity appears in
[THM-4473, digit chains and repunits](../../01-canon/theorems/THM-4473-collatz-digit-chains-rotation-repunits.md).
The inverse-fibre identity `U(4z+1)=U(z)` is inherited from
[the arithmetic braid note, section 2](arithmetic_braids_20260917_collatz.md).

Here `U(n)=oddpart(3n+1)` acts on positive odd integers; `U(1)=1`.
A reset is not necessarily a descent. The canonical hostile is a long
Mersenne rise. The least-used coordinate in this bounded probe is the
*synchronized step count of two actual positive orbits*, alongside their
ordinary heights. Our small board is maximal runs, reset depth, inverse
fibres, binary siblings, and well-founded reduction. The question is
whether compression supplies a relationship between distinct starts,
rather than just renaming the original dynamics.

## 2. A complete synchronized-coupling criterion

Write uniquely

    n=2^r u-1,  r=v2(n+1)>=1,  u positive odd,
    N=3^r u,   s=v2(N-1)>=1,  z=(N-1)/2^s.

There are `r-1` initial odd steps with valuation one, followed by a step
of valuation `s+1`. Thus the compressed first reset is

    C(n)=U^r(n)=z=oddpart(3^r u-1).

The larger binary sibling `m=2n+1=2^(r+1)u-1` still has one rising odd
step left after `r` steps, so

    U^r(m)=2N-1=2^(s+1)z+1.

**Proposition.** These two positive starts merge at the synchronized time
`r+1` exactly when `s=1`:

    U^(r+1)(n)=U^(r+1)(2n+1)  iff  v2(3^r u-1)=1.

If `s=1`, the two values at time `r` are `z` and `4z+1`, so the inherited
inverse-fibre identity proves equality at the next step. If `s>=2`,

    U(2^(s+1)z+1)=3*2^(s-1)z+1 > U(z),

because the numerator before removing powers of two is exactly four
times the displayed odd integer, whereas `U(z)<=(3z+1)/2`.
This proves both directions, including `n=1`. It excludes equality at
that specified time when `s>=2`; it does not exclude a later merger.

The criterion is equivalently `3^r u=3 mod4`, one of the two odd `u`
classes modulo4. For every fixed `r` it holds for half the starts with
that `r`. Since `r` has relative density `2^-r` among odd starts and the
tail `r>=R` lies in the single class `n=-1 mod2^R`, the infinite union
has relative natural density `1/2` among odds. This is a spatial source
density, not an orbit frequency.

**Well-founded consequence, not a complete algorithm.** For the larger
partners satisfying this criterion, termination at1 is equivalent to
termination of the smaller odd integer `(m-1)/2`. These larger partners
have relative natural density `1/4` among odd integers: the injective
map `n ->2n+1` doubles scale. The reduction remains legitimate even if
their actual common orbit point lies above both sources.

But two such reductions cannot occur consecutively within a fixed binary
tower. Reducing `m=2^(r+1)u-1` to `n=2^r u-1` uses `3^r u=3 mod4`.
Reducing `n` again, when `r>=2`, would require `3^(r-1)u=3 mod4`.
Multiplication by3 exchanges the two odd residue classes, so the two
requirements are incompatible. Thus this local rule cannot consume an
arbitrary initial run by repeated well-founded reductions.

For the smaller subfamily `3^r u=3 mod8`, `z=1 mod4`, and its next
compressed run has length one. There is then also the compressed identity

    C(C(n))=C(2n+1).

The mod8 condition identifies the next compressed block with a single
odd step; mod4 alone proves the synchronized odd-step identity above.

### A finite binary recognizer and exact tower matchings

For a larger eligible source write `m=2^R u-1`, with `R>=2` and `u` odd.
Its binary digits, read least-significant first, begin

    1^R 0 b ...,  with u=1+2b mod4.

Finite positive integers are padded by implicit leading zeros; this
includes Mersenne numbers without an exceptional parsing rule. The smaller
source has exponent `r=R-1`, so eligibility is exactly

    3^(R-1)u=3 mod4  iff  b=R mod2.

An eight-state deterministic automaton suffices: start, exactly-one-one,
positive-even-run, odd-run-at-least-three, expect-zero, expect-one,
accept, reject. The run states count parity, but distinguish the forbidden
length1. At the first zero after at least two ones, the automaton asks for
the next bit to equal that parity; all later bits are irrelevant to this
predicate. The script gives its complete transition table.

This is an exact finite graph recognizing a fate-preserving reduction,
not a finite-state implementation of Collatz. No information needed for
the eligibility predicate is lost, but the original integer must still
be retained to compute the reduced source `(m-1)/2` and its height.

For each fixed positive odd `u`, the tower `{2^r u-1:r>=1}` therefore has
an exact matching of neighboring exponents:

| u modulo4 | Paired exponents | Unpaired initial exponent |
|---|---|---|
| 1 | (1,2),(3,4),(5,6),... | none |
| 3 | (2,3),(4,5),(6,7),... | 1 |

For `u=7`, the actual pairs begin `(27,55),(111,223),(447,895)`.
Every pair shares a future, proved algebraically; convergence is not
being inferred for the whole infinite tower from any finite number of
known seeds. The unconnected pairs are precisely the missing connection.

### The three total stopping-time conventions

Let `tau_U`, `tau_T`, and `tau_S` denote first-hit times of1 for the odd
accelerated, shortcut, and standard maps, respectively. Suppose `n>1`
satisfies the cheap-reset criterion and reaches1. Then its larger sibling
also reaches1, and

    tau_U(2n+1)=tau_U(n),
    tau_T(2n+1)=tau_T(n)+1,
    tau_S(2n+1)=tau_S(n)+1.

For proof set `a=v2(3z+1)`. The first `r+1` odd steps of `n` have
valuations `1,...,1,2,a` (with `r-1` initial ones); their sum is
`r+1+a`. For the larger sibling they are `1,...,1,a+2` (with `r`
initial ones), with sum `r+a+2`. Both lists have `r+1` odd steps and
the same endpoint. No first hit of1 was passed earlier: the initial
rises are strict for starts greater than1, and `z=1` would force
`3^r u=3`, hence `n=1`. A shortcut trajectory takes exactly the valuation
sum to execute these odd steps; a standard trajectory takes that sum
plus the number of odd steps. Appending the common tail proves all three
identities, with their exact convention dependence.

The exclusion of `n=1` is essential because its first-hit time is zero.
The larger partner3 has odd/shortcut/standard first-hit times2/5/7.
For the nonterminal control `27,55`, the times are respectively41/41,
70/71, and111/112. These equalities are consequences of exact coupling;
they are not a statistical law or a termination proof for new seeds.

## 3. Binary repunits: exact reset classification and a hostile merger

For `u=1`, elementary factorization gives

    v2(3^r-1)=1                 if r is odd,
    v2(3^r-1)=2+v2(r)           if r is even.

For odd `r`, `3^r=3 mod8`. For even `r`, write `r=2^a b` with `b` odd.
The odd geometric sum reduces to the exponent `2^a`; starting from
`v2(3^2-1)=3`, each successive squaring adds exactly one because
`3^(2^j)+1=2 mod8` for `j>=1`. No external LTE theorem is needed.

Consequently, for every positive `r`,

    C(2^r-1)<2^r-1  iff  r is 2,4, or8;
    C(2^r-1)=2^r-1  iff  r=1.

Here are complete elementary bounds, not an extrapolation of a table.
For odd `r>=3`, `3^r>2^(r+1)`, giving strict growth. For even `r>=10`,
`2^s=4*2^v2(r)<=4r`, while `(3/2)^r>4r`: the latter holds at10 and its
left-to-right ratio increases thereafter. Thus `3^r>2^(r+s)`, again
giving strict growth. The remaining `r=1,2,4,6,8` are exact substitutions.
If `3^r>2^(r+s)`, the endpoint comparison follows directly from
`3^r-1 > 2^s(2^r-1)`; no cycle classification is used.

Every odd `r` also gives

    C(C(2^r-1))=C(2^(r+1)-1).

So convergence of the even-exponent Mersenne start reduces to convergence
of its smaller odd-exponent neighbor. This does not establish convergence
of those odd-exponent neighbors. The first useful hostile is

    C(31)=121,  C(121)=91=C(63),  91>63>31.

Both starts merge after six odd steps, but their shared endpoint exceeds
both. Coalescence preserves eventual fate; the precise stopping-time
changes above depend on the retained prefix counts. It does not preserve
descent or supply a negative average drift. For every odd `r>=9`, the common endpoint
likewise exceeds the larger sibling, by the just-proved classification
at even exponent `r+1>=10`.

## 4. Connection contract, audit, and stopping reason

Source: two actual positive odd trajectories with starts `n,2n+1`.
Target: one common future after an explicitly synchronized time.
Map: the exact `(r,u,s,z)` arithmetic above, followed by the inverse-fibre
identity. Preserved predicate: reaching1, in both directions. Lost by
retaining just the common endpoint: source height and previous excursions.
Needed sidecar: `r,s` and both source heights. Cheapest decisive test:
`31,63 ->91`, which defeats a merger-implies-descent claim.

This recovers a useful cheap-reset coupling, but it does not supply the
missing global resource inequality. The inability to repeat the binary
reduction is exact; it is not a failure of the computation to search far
enough. There is no forced Pell, cube, or prime-distribution transfer here.
Any stronger theorem would need either a further reduction applicable to
the blocked residue class or a bound on the arithmetic cost between
successive permitted reductions along one fixed integer orbit.

Run `python -B 04-computation/experiments/crossroads_crossing_20260926_arithmetic_run.py`.
The retained output checks every odd source below200000 against direct
iteration, Mersenne exponents1--256, both positive and strict-separation
controls, and the mod8 compressed identity. It checks the binary automaton
on every positive integer below200000, including even hostile inputs;
tower pairs for all odd `u<64` and exponents1--32; and the three stopping
identities on eligible odd sources3--9999 by independent direct iteration
in each convention, with a declared finite cap. All checks remain enabled
under `python -O`. The same script independently computes the root dyadic
note's *unwrapped* carries `287,251,227,211`, not just their residues modulo64.
This catches an erroneous common multiple-of64 offset which a residue-only
test necessarily loses. The root note's analytic identities and Fourier
projection sign were independently audited and accepted.

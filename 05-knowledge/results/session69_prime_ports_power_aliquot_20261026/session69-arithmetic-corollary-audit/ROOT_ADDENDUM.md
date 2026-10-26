# Root addendum: the unique perfect prime-square target

**Final root packet ACCEPTED after the two wording repairs below.** The unique
perfect prime-square corollary is accepted without qualification beyond its
stated hypotheses. The original Session69 arithmetic audit remains frozen
and unchanged. This separate addendum reviews the root's later
`RARE_SIX.md`, `CONNECTIONS.md`, and both root checkers. It does not change
the producer aliquot or power reports.

## 1. Independent proof review

For positive `N` and prime `p`, the conjunction

`s(N)=N` and `s(p^2)=N`

holds if and only if `(N,p)=(6,5)`. This is a global statement, not a
finite-cutoff claim and not merely a statement within the twin-prime family.

Indeed `s(p^2)=p+1`, so `N=p+1`. At `p=2`, `N=3` is not perfect. Every
remaining candidate is even, so write `N=2^a*m` with `a>=1`, `m` odd.
For `q=2^(a+1)-1`, multiplicativity of sigma gives
`q*sigma(m)=(q+1)m`. Since `gcd(q,q+1)=1`, write `m=q*t` and obtain
`sigma(m)=m+t`. Here `q>=3`. If `t>1`, the divisors `1,t,m` are distinct,
forcing `sigma(m)>=m+t+1`, a contradiction. Thus `t=1`; then
`sigma(q)=q+1` forces `q` prime. Consequently `r=a+1` is prime, since a
composite exponent factors `2^r-1` nontrivially.

If `r=2`, then `N=6`, `p=5`, and both divisor sums are directly verified.
If `r` is odd, `N=2^(r-1)(2^r-1)` is congruent to one modulo three.
Also `r>=3` gives `N>=28`; hence `p=N-1>=27` is divisible by three and
cannot be prime. Equivalently the root's proof can force `p=3`, `N=4`
and check that four is not perfect. Both arguments are valid.

Every hypothesis is used. No existence/nonexistence assertion about odd
perfect numbers is needed: parity already eliminates that branch. No
infinitude or completeness assertion for Mersenne primes is needed. Nor
does the result classify all perfect fibers or all nonsingleton fibers.
Removing perfectness gives many prime-square targets, already
`s(3^2)=4` and `s(7^2)=8`. Another perfect target such as 28 is not
excluded from having other source shapes by this theorem alone.

This precise conjunction justifies the user's intuition that the six
self-return plus prime-square companion is unusually special. It does not
make the separate odd-root square-later event rare: that event has density
one half. Primitive normalization still sends every perfect source to
`1+1=2`; the removed gcd retains its scale.

## 2. Root connection firewall review

The arithmetic quotient hostiles are correct:

- `s(57)=23`, `2=23 mod7`, whereas target two is untouchable;
- `s(65)=19`, `5=19 mod7`, whereas target five is untouchable;
- multiplication by two modulo seven sends targets two and five to the
  touchable representatives four and three;
- residue squaring collapses every nonzero opposite pair, not just `{2,5}`;
- Paley-seven has raw `H=189`, automorphism order 21, and ordinary
  presentation count 9. Raw image membership therefore does not descend
  modulo seven: excluded counts 7,21 and attained count 189 all have residue
  zero. A symmetry-group order is not a Hamiltonian-path count.

The binary-term ABC correspondence, first quality example at source 72,
aliquot gcd normalization, radical 930 of `(25,6,31)`, and fixed-core
uniformity limitation all agree with the original independently audited
proofs.

Two precise prose repairs were requested during this addendum:

1. The positive-integer doubling/squaring generator graph has a unit loop
   `1 -> 1`, so the single aliquot loop at six alone is not an invariant
   obstruction to arbitrary graph isomorphism. The repaired argument counts
   loops: the generator graph has exactly one, at one, while the aliquot
   graph has at least two distinct loops, at six and 28. If the unit loop
   is suppressed, the generator graph has no loops and the obstruction is
   stronger. These statements concern generator edges, not identity maps
   deliberately added to every vertex.
2. It is the **value 36**, not its infinite arithmetic row, that has five
   nondegenerate occurrence addresses.

The tournament report's P2/P3 proof and §7 table were read for consistency
of the root synthesis. The sufficient **recursive ordinary, uncoloured**
target-incidence carrier is already closed, conditional on the supplied
canonical quotient relations/actions. Its economical compression and an
unrestricted prime bank are distinct remaining tasks. The root summary
must not label recursive closure itself unproved. This addendum is not a
second full audit of the tournament lane's finite response tables; that
lane has its own independent audit.

## 3. Fresh controls and reproduction

`python3 check.py` and `python3 -O check.py` use a smallest-prime-factor
and multiplicative-sigma path independent of the root's proper-divisor
sieve. They agree byte for byte and perform 196 active gates:

- complete perfect-target census through 1,000,000 gives
  `6,28,496,8128`;
- every prime `p<=100000` is tested through the equivalent condition
  `sigma(p+1)=2(p+1)`, with sole pair `(N,p)=(6,5)`;
- 168 primes below 1000 have `s(p^2)=p+1` independently checked by direct
  divisor-pair enumeration;
- the small Euler steps, modulo-three obstruction, prime-two boundary,
  residue firewalls, all opposite-pair collisions, and unit-loop hostile
  are checked explicitly.

`--negative-control` deliberately expects an impossible companion pair
`(28,27)` and must fail at the named uniqueness gate. Normal and optimized
modes both exit 1 there. Gates use explicit exceptions and are not erased
by optimization.

Fresh executions of both root checkers in normal and optimized modes
match their supplied result files. The root's larger Mersenne controls
include exponents `2,3,5,7,13,17,19,31`; they corroborate, but do not prove,
the global uniqueness statement. Root firewall permutations exhaust
all `7!` vertex permutations for `H=189`, `Aut=21`, `tc=9`.

## 4. Integration status

The corollary, repaired connection prose, and finite receipts are accepted.
Both requested repairs were independently read in the final root file,
whose SHA-256 begins `dde64747557b25bf`. The exact recursive composition
wording also agrees with P3. Final root input hashes and
the disposition of the two requested wording repairs are recorded in
`ROOT_REVIEW.json`. This packet has its own allowlist and manifest, leaving
the original arithmetic audit immutable. No shared repository edit or
publication was performed by this reviewer.

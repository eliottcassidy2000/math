# Session243 independent arithmetic audit

**ACCEPTED: the Brown inverse and prime-gap/Brocard bridge are correct within
their stated scopes. FINITE-EXACT computational controls pass.**
No new solution or infinitude assertion for either Brocard problem, Polignac,
Fortune, or primorial primes is established.

## Objects and source boundaries

Audited producer files are `../brown/PROOF.md`, `../brown/verify.py`,
`../provenance/PROOF.md`, and `../provenance/verify.py`.
`REPLAY.json` records their exact SHA-256 hashes, the eight recovered source
blob pins, all producer replay outcomes, and the audit replay outcomes.
The independent implementation imports no producer code. It uses an
Eratosthenes sieve, Legendre valuations, a Gray-code partition traversal, an
additional all-divisor enumeration, and exact integer square roots.

The inherited mechanism is the factorial divisor involution and the
primorial two-sign wheel. The decisive hostiles are a lost equality quotient
and an omitted offset-one prime. The less-used sidecar is the pair of
primality bits at `P_n-1` and `P_n+1`. The live comparison is between unitary
factor partitions, modular signs, roughness, prime gaps, and graph witnesses.
These are different predicates.

The recovered records genuinely contain both meanings of Brocard:

- Factorial-square problem: `n!+1=m^2`, with the three known pairs
  `(4,5),(5,11),(7,71)`. The inherited prose's `5,11,71` are roots, not squares;
  completeness is conjectural.
- Prime-square conjecture: at least four primes between squares of
  consecutive primes `p<p'`, for `p>=3`. The omitted pair `(2,3)` gives only
  the primes `5,7`.

The current user wording alone does not select one of these meanings. The
packet appropriately recovers both and does not claim otherwise. Historical
owner quotations are provenance, not fresh authorization or mathematical
proof. External priority or literature-wide absence claims are not certified
by this arithmetic audit.

## 1. Brown inverse: accepted

For `n>=4`, put `N=n!`, `e=v_2(N)>=3`, and `C=N/4`. A positive solution has
odd `m`, and `(m-1)/2,(m+1)/2` are coprime adjacent positive factors of `C`.
Thus each full prime-power factor of `C` lies entirely on one side. The
nontrivial even factor `2^(e-2)` uniquely names the even side `A`.

Conversely, `AB=C` and `|A-B|=1` give `(A+B)^2=N+1`. There is no duplicated
solution from exchanging factors because the even side is fixed. There are
exactly `2^(pi(n)-1)` candidate assignments. The inverse is complete, not a
bounded search in possible roots. It is nevertheless a finite membership
algorithm for each fixed `n`, not a proof that no later `n` works.

The sign conditions are also correct, including the exponent at two:

- Odd `p`: `m=+1` or `-1 mod p^v_p(N)`.
- At two: `m=+1` or `-1 mod 2^(e-1)`.

They are equivalent to `N | m^2-1`. Necessity at two follows because one of
`m-1,m+1` has valuation one and the other has valuation at least `e-1`.
For an actual Brown equality, the latter valuation is exactly `e-1`.
Sufficiency uses the additional factor of two from the opposite neighbor.
The sign packet does not impose quotient one.

Boundary and hostile controls:

- `n<4` has no Brown solution; the packet construction is not applied there.
- At `n=5,m=7`, roughness survives but the sign packet fails.
- At `n=7,m=29`, all radical-level square congruences survive, but the full
  prime-power conditions fail.
- At `n=4,m=23`, the entire sign packet survives while
  `(23^2-1)/4! = 22`, not one.
- Every unit modulo `24` does square to one. Thus `n=4` could not support
  the producer's earlier proposed roughness-versus-divisibility hostile;
  the supplied final proof correctly uses `n=5,m=7` instead.

The fixed-n inverse repairs the quotient loss through adjacency, not by
turning factors into graph operations. Arithmetic square roots do not become
physical tournament suppliers merely because their squares are supplied.
The graph-specific existence and omission claims remain inherited or owned
by the companion construction audit, not independently graph-certified here.

## 2. Primorial crossing gaps: accepted

For `n>=2`, the excluded-offset definitions `F_n,L_n` exist. Adding back
exactly the omitted offset one gives the actual nearest-prime distances:
`f_n=1` when `P_n+1` is prime and `F_n` otherwise, and similarly for `l_n`.
The nearest primes around composite `P_n` are consecutive and have gap
`f_n+l_n`.

If `F_n=L_n=d`, minimality excludes every interior offset `2,...,d-1`.
The only possible intervening primes are the prime members of
`{P_n-1,P_n+1}`. Hence the endpoint pair is consecutive if and only if both
of these numbers are composite. All statements retain the `n>=2` boundary;
`L_1` is undefined.

The independent sieve confirms:

- `P_6=30030`, `F_6=L_6=17`, endpoint pair `30013,30047`, unique intervening
  prime `30029`, and actual crossing gap `18`, not `34`.
- `P_3=30`, `F_3=L_3=7`, with both `29,31` intervening.
- `P_4=210`, `(F_4,L_4)=(13,11)`, actual gap `12`, not `24`.
- `P_8=9699690`, `F_8=L_8=23`, both offset-one values composite, and actual
  consecutive gap `46`. The sufficient direction is nonvacuous.
- `P_2=6`, `L_2=3`, and `P_2-3=3` is prime despite `gcd(3,6)=3`.
  This is a one-endpoint exception, not a prime pair: `P_2+3=9`.

The missing-prime description actually remains true for unequal `F_n,L_n`
when using the asymmetric raw endpoints; this stronger finite control is
checked without changing the producer theorem.

## 3. Fixed offsets and quantifiers: accepted

For fixed integer `d>1`, choose a prime divisor `r|d`. Once `r|P_n` and
`P_n-d>r`, both endpoints are divisible by the proper divisor `r`. These
conditions hold eventually. Therefore only finitely many primorial indices
can supply a prime pair `P_n±d`. The same proof works for every fixed
positive integer multiplier `c` in `cP_n±d`.

The strict inequality is essential: a number equal to its divisor is prime,
which explains the `n=2,d=3` lower-endpoint control. Removing that guard
would make the inference false. At `d=1`, there is no prime divisor of the
offset, so this obstruction disappears. Simultaneous primality then really
gives consecutive twins, but its infinitude is not established.

This rules out the fixed-primorial-center template for an infinite supply of
any prescribed gap `2d>2`. It does not disprove Polignac, which allows other
centers. Allowing `d` to vary evades the proposition but loses the fixed-gap
quantifier. Neither wheel survival nor an arithmetic prime endpoint proves
a graph realization.

## 4. Prime-square Brocard: accepted

For consecutive odd primes `p<p'=p+2d`, the open prime-square interval
splits into exactly `d` odd-square brackets plus shared square endpoints.
The latter are composite, so the prime sets partition exactly. The length
identity is `p'^2-p^2=4d(p+d)`.

Legendre, if assumed, supplies one prime in each of the `2d` consecutive
integer-square intervals. Thus this counting argument proves the four-prime
conclusion when `d>=2`; at a twin pair it guarantees only two. This is not a
proof of the logical nonimplication “Legendre cannot imply Brocard” by any
possible argument. The final producer text correctly restricts its claim to
this counting argument. The four Oppermann subintervals per odd-square
bracket give the stated stronger conditional implication. For `p>=3` all
cut points are composite and the interval boundaries introduce no lost prime.

Polignac's fixed gap would fix the number of brackets in infinitely many
such prime-square unions. It gives no prime-count bound inside the squared
interval by itself. The factorial-square problem is separate.

## 5. Finite verification and exact scope

Independent finite universe:

- All `484,682` canonical unitary partitions for every `n=4,...,60`.
- A second, direct all-divisor enumeration for `n=4,...,10`, scanning
  `1,408` candidate small divisors in total.
- Every residue in the common period `n!/2` for `n=4,...,9`:
  `204,552` sign/divisibility comparisons, including both parities.
- `14,649` further signed-integer comparisons: `-128<=m<=128`,
  `4<=n<=60`.
- Primorial indices `2,...,8`, independently certified with an Eratosthenes
  sieve through `9,700,690`.
- Proper-divisor certificates for every `d=2,...,100` and `c=1,...,5`,
  totaling `990` endpoint certificates, plus persistence at later centers.
- `29` consecutive odd-prime pairs from `(3,5)` through `(113,127)`.

Exactly the three inherited Brown pairs occur in the stated finite range.
This duplicates no universal completeness claim. Prime-square controls test
interval identities, not a proof of any named conjecture.

One minor discrepancy was found and repaired before acceptance: the
provenance verifier's range reaches the 31st prime `127`, whereas its draft
prose said “through the first30 primes.” The final prose now gives the exact
29-pair universe. No core theorem correction was required.

All three verifiers pass normally and with `-O`, with byte-identical JSON
within each verifier. The eight deliberately false controls fail in both
modes: `6` positive replays and `16` negative replays. Existing producer
`CHECK.json` files match fresh output. All eight provenance blob pins verify
against the specified Git revision. Failure gates are explicit exceptions,
not optimizable Python assertions.

## Reproduction

From this directory:

    python3 -I -B verify.py > CHECK.json
    python3 -I -B -O verify.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json
    python3 -I -B replay.py > REPLAY.json

`replay.py` also checks the producer outputs and the pinned Git objects.
No producer source was edited by this audit. No commits or pushes were made
by the auditor; publication belongs to the parent session.

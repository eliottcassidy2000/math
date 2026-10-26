# Independent review of the arithmetic atlas

2026-10-26. **PASS in the stated proof and finite scopes.** The README and
entire producer were read independently. Both producer modes passed with
byte-identical receipt SHA-256
`c71a50087eaeccf04d7d6167c570beefcf229b644670bcd5916585f7c18c832a`.
This reviewer changed no producer files.

1. The even-sector normalization preserves both branch guards and their
halting alternative. On positive integers, the inverse image residues 2
and 3 modulo five are disjoint, and each displayed inverse is positive,
smaller, and in its proper modulo-three branch. Finite reverse descent
therefore gives an actual unique root and no mergers. Calling this a
path forest, rather than the separate binary residue tree, is correct.

2. The two inverse affine maps on Z_3 have disjoint images, contraction
ratio 1/3, and cover the nonzero residue balls. Their infinite coding is
a homeomorphism with binary sequences and conjugates F to the one-sided
shift. The uniform cylinder measure gives the lower Hausdorff bound;
the 2^k cylinder cover gives the upper bound. Finite-level natural-density
bounds imply density zero for integer survivors, not emptiness. Refining
one more ternary digit gives exactly the claimed one stopping/two
continuing children, including in positive integer representatives.

3. The sign-word carry order is correct. Setting the composed image equal
to its source gives `-A_k/(5^k-3^k)`; its denominator is a 3-adic unit and
the geometric-sum bound is exactly 1/2. Prepending finitely many inverse
branches preserves the real interval, proving the eventually-periodic
exclusion for positive integers. Aperiodic integer survivors are not
excluded. The real/3-adic completion warning prevents the apparent but
invalid global conclusion.

4. Sloane's exponent period is `2*3^(k-1)`, with 2^k good residues. The
relative proportion consequently differs from F's by the first-level
factor as stated. Actual digit length must be retained: padded leading
zeros have no bearing on zerolessness of the represented integer.
The prime-valuation flow identity is an iff by positive-integer unique
factorization; dropping support multiplicities or omitting output primes
would invalidate it. It is correctly kept as a weighted relation, not
oriented into a cosmetic tournament.

I independently opened [Michel's primary analysis](https://bbchallenge.org/~pascal.michel/beh)
and [Saye's primary paper](https://cs.uwaterloo.ca/journals/JIS/VOL25/Saye/saye3.pdf).
Their stated native clocks, final tape pattern, exceptional exponents and
`2*3^45` finite verification bound match the producer's citations. This
checks the attribution; it is not a rebuild of those authors' full proofs.
The extra native control below directly verifies the quoted three machine
rules on 93 small inputs. It does not simulate the 47-million-step blank run.

## Independent native and valuation controls

This block imports no atlas or prior research engine. It executes the
five-state transition table directly, using Michel's native halt rule,
then independently compares exhaustive divisor sums with prime valuations.
Run at the repository root, normally or with `python -O`.

```python
import json

def need(ok, why):
    if not ok:
        raise RuntimeError(why)

# (write, movement, next-state), states A..E=0..4, H=5.
table = ((1, 1, 1), (1, -1, 2), (1, 1, 2), (1, 1, 1),
         (1, 1, 3), (0, -1, 4), (1, -1, 0), (1, -1, 3),
         (1, 1, 5), (0, -1, 0))
native_cases = 0
for k in range(31):
    for residue in range(3):
        n = 3 * k + residue
        cost = (5 * k * k + 19 * k + 15 if residue == 0 else
                5 * k * k + 25 * k + 27 if residue == 1 else 6 * k + 12)
        ones, head, state = set(range(1, n + 1)), 0, 0
        for step in range(cost):
            need(state != 5, 'no premature native halt')
            write, direction, target = table[2 * state + int(head in ones)]
            if write:
                ones.add(head)
            else:
                ones.discard(head)
            head += direction
            state = target
        relative = {x - head for x in ones}
        if residue < 2:
            target = 5 * k + (6 if residue == 0 else 9)
            need(state == 0 and relative == set(range(1, target + 1)), 'native live target')
        else:
            expected = {-1, 1, 3 * k + 5} | {4 + 3 * j for j in range(k + 1)}
            need(state == 5 and relative == expected, 'full native halting tape')
        native_cases += 1

def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

balance_cases = 0
for n in range(1, 1025):
    direct = sum(d for d in range(1, n + 1) if n % d == 0)
    incoming = {}
    for p, a in factor(n).items():
        for q, value in factor(sum(p ** j for j in range(a + 1))).items():
            incoming[q] = incoming.get(q, 0) + value
    need(incoming == factor(direct), 'all output prime multiplicities')
    for k in range(1, 7):
        need((direct == k * n) == (incoming == factor(k * n)), 'both valuation directions')
        balance_cases += 1
print(json.dumps({'status': 'PASS', 'independent_native_clock_cases': native_cases,
                  'independent_valuation_iff_cases': balance_cases}, sort_keys=True))
```

No new equivalence between open problems, bare-tournament invariant, or
Collatz rank has been inferred from these reviews. The missing source and
operation information is kept explicit throughout the producer.

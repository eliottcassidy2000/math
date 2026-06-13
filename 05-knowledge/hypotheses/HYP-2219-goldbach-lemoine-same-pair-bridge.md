# HYP-2219: Goldbach/Lemoine Same-Pair Bridge Graph

**Status:** OPEN synthesis / exact finite carrier.

## Claim

The S642 same-pair projection induces a finite bridge graph between even
Goldbach targets and odd Lemoine targets when they use the same odd-prime pair.

```text
E = p + q
O = a + 2*b
```

If `O` uses the same pair as `E`, with `b` the doubled prime, then

```text
b = O - E,
a = 2E - O.
```

Thus `(E,O)` reconstructs the ordered prime pair exactly.  S643 records the
resulting companion graph and its side channels.  The duplicate pairs `p=q`
form the branch locus

```text
E = 2p,   O = 3p.
```

The prime `2` is a boundary/apex channel: Lemoine reps with doubled prime
`q=2` do not share an even Goldbach pair, because `{2,p}` sums to an odd
number.

## Evidence

S643 computes the bridge through small targets.  For an unordered odd-prime pair
`{p,q}`, the even node is `E=p+q`, and the odd companions are

```text
O = E+p,  E+q.
```

When `p=q`, the two companions fold together at `O=3p`.

Up to odd target `121`, the only Lemoine target with no odd-prime bridge
representation is `7`, represented only as `3+2*2`.  After that first apex,
the odd-prime bridge is already present in the tested range.

## Repo Meaning

This is the additive analogue of the pi/e carrier.  There, two shadows
`S=e+pi` and `P=e*pi` reconstruct a hidden unordered pair by Vieta.  Here, two
linear shadows `E=p+q` and `O=a+2b` reconstruct a hidden ordered prime pair.
The branch/duplicate line `p=q` is the discriminant-zero face.

For the LRC/doubled-prime thread, the pair `p=q` is not just a duplicate; it is
the parity bridge where `2p=p+p` is simultaneously additive doubling and
multiplicative doubling.  The odd companion `3p=p+2p` is the next rung.

## Next Tests

- Build the even/odd companion graph at larger scale and measure odd-node parent
  multiplicity.
- Jackknife the duplicate branch, twin-gap branch, and `q=2` apex boundary.
- Compare bridge-parent multiplicity with ordinary Lemoine representation
  counts.
- Use `(E,O)` recovery as a toy model for carrier inversion in LRC owner/carry
  ledgers.

**See also:** HYP-2218, `04-computation/goldbach_lemoine_pair_bridge_s643.py`,
`05-knowledge/results/goldbach_lemoine_pair_bridge_s643.out`,
`07-reflections/goldbach-lemoine-same-pair-bridge-s643.md`, HYP-2215,
HYP-2049, HYP-2044, HYP-1984.

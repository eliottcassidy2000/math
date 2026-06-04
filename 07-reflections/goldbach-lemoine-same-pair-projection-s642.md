# Goldbach/Lemoine Same-Pair Projection S642

The user's framing is right, with one important guardrail: full Goldbach and
Lemoine/Levy are conjectural, not theorems.  But as a carrier model it is
excellent.

The clean object is prime-pair space.  For odd primes `p,q`,

```text
E = p + q
O = p + 2q
```

are not two unrelated target numbers.  They are two shadows of the same ordered
pair.  In fact,

```text
q = O - E
p = 2E - O.
```

So the pair `(E,O)` is an invertible coordinate system on ordered odd-prime
pairs.  This is the same style as the repo's `S=e+pi`, `P=e*pi` thread: one
shadow alone has fibers, but two compatible shadows reconstruct the hidden
pair.

The swap is the beautiful part.  The even target `E` cannot see whether the
pair was `(p,q)` or `(q,p)`.  The odd target can, because one coordinate is
doubled.  Under the swap:

```text
O -> 3E - O.
```

Thus every off-diagonal Goldbach pair gives two odd Lemoine shadows symmetric
around `3E/2`.  The duplicate pair `p=q` is the fixed point:

```text
p -> E=2p, O=3p.
```

That puts `p=7` exactly at `(14,21)`:

```text
7+7 = 14
7+2*7 = 21.
```

This is not raw numerology.  It is a fixed-locus statement.  The pair `(14,21)`
is the prime-7 diagonal of the Goldbach/Lemoine projection map.

The connection to the last LRC/unit-distance session is now sharper but also
more disciplined:

```text
Goldbach/Lemoine layer: 14 and 21 share the duplicate pair (7,7).
LRC/unit-distance layer: the hard proof channel is still the 27-carrier:
  LRC C=27 shell/carry vs unit-distance 57=20+37 spine/bulk.
```

So we should not replace HYP-2217 with "14 and 21 are connected because 7".
The better story is a stack:

1. prime `7` gives a diagonal additive packet `(2p,3p)=(14,21)`;
2. `21=3*7` is a durable tournament obstruction scalar;
3. LRC `n=14` has the `C=27=3^3` carry shell;
4. unit-distance `n=21` has the Moser `57=2*27+3` carrier.

That stack is exactly the kind of thing the side-channel ledgers were built for:
one visible scalar, several retained witnesses.

The exceptional prime `2` is also useful.  S642 found that through `300`, the
only odd row that has Lemoine representations but no odd-prime-compatible
same-pair projection is `7`, using the `q=2` channel.  Prime `2` is therefore
not just a nuisance; it is the seam where Lemoine can work without an even
Goldbach same-pair partner.

The next good computation is a larger same-pair graph with the `q=2` channel
split out.  Then compare diagonal packets `p -> (2p,3p)` against LRC rank-one
rows `n=2p`, tournament forbidden-ish packets `3p`, and unit-distance carrier
rows where the same scalar is legal only after spine/bulk data is retained.

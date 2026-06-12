# Coefficient tiling and prime-irreducible addresses

The coefficient-sign tournament is real.  That was the pleasant surprise.
Your degree-5 picture is not a loose metaphor: the row sizes

```text
1, 2, 3, 4, 5
```

are exactly the skip rows of a fixed ordered tournament on six vertices.  Put
`a_5` on the long skip apex, `a_4` on the two next skips, and so on.  That is
the fixed-path tiling cube wearing polynomial clothes.

But the computation also gave the guardrail: unmarked tournament structure is
too lossy.  In the degree-4 sweep, every unmarked coefficient-tournament class
was mixed for irreducibility.  Even the marked sign pattern was mixed.  The
same sign row could be an irreducible quartic with many prime hits or a
reducible product with none.

The fix was exactly the repo's recurring fix: attach the missing address.

For fixed divisors, the missing address is local residue data.  Once the scout
kept `local_zero_primes`, the fixed-divisor mixing vanished in the tested
family.  That is the polynomial version of the LRC lesson: "q blocked" is too
coarse, but "which row blocks which local residues" starts to be proof-bearing.

The best turn of the picture is to let `a_0` be the Hamiltonian spine.  The
constant coefficient is not just a scalar at the bottom of the polynomial.  It
is `f(0)`, the first local obstruction, the first Cohn digit, and the row that
anchors evaluation.  So the better model for degree `d` uses `d+2` tournament
vertices:

```text
a_d  -> apex row
...
a_1  -> near-spine row
a_0  -> adjacent Hamiltonian-path row
```

That feels right.  It makes the constant term the observer/spine rather than
the leftover.

Cohn is the warning against stopping at signs.  The base-3 repunit
`1+x+...+x^8` and the binary repunit `1+x+...+x^10` both have all-positive
coefficient signs, hence transitive sign tournaments.  One is reducible with
factor degrees `[2,6]`; the other is irreducible.  The difference is not sign.
It is the place-value address of evaluation at the base.

So the new grammar is:

```text
sign rows          = tiling cube / switching coordinate
magnitudes         = row weights
local residues     = fixed-divisor detector
p-adic valuations  = Newton polygon / Eisenstein detector
base-b weights     = Cohn address
value depth        = Singh/Murty certificate
trace subsets      = recombination ambiguity
```

This gives a lot of creative routes.

Gauss lemma becomes switching: strip global content as a gauge before asking
for irreducible structure.  Eisenstein becomes a valuation tournament: all
lower rows flow into a p-adic sink, the leading row escapes, and the constant
row carries the `p^2` guard.  Newton polygons become coefficient-tile lower
hulls: orient row pairs by valuation slopes and ask when the hull forces one
factor block.

For LRC14, the analogy is now sharper than "prime value equals lonely
witness."  The fixed-divisor analogue should be a denominator/resource row
that vanishes on every local residue.  A real counterexample would need to
look irreducible after every such local row is attached.  That is exactly the
HYP-2443/HYP-2446 ledger problem in polynomial language.

The sentence I would keep is:

```text
Signs give the tiling cube; primes and irreducibility live in the addresses
that survive recombination.
```

That is the bridge.

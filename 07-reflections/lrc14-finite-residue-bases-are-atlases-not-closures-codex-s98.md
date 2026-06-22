# LRC14 Finite Residue Bases Are Atlases, Not Closures -- Codex S98

The prompt basis `{83,89,21}` is a real signal, but the exact obstruction shows
how to read it without overclaiming.

The finite-basis theorem is false in the strongest possible way.  Given any
finite denominator list `B`, the primitive covering row

```text
S_B = {1,2,...,11,13,84*lcm(B)}
```

kills every denominator in `B`: the tail is divisible by each `D in B`, so at
every `a/D` that runner sits at the observer.  For the prompt basis this gives
the explicit killer tail `13030668`, with exact counts

```text
[(83,0), (89,0), (21,0)].
```

This refines THM-566.  THM-566 killed every `D<=B0`; S98 kills any finite list,
even a sparse hand-picked one.

The positive part is still important.  On a broad deterministic sample of `602`
primitive covering rows up to speed `10000`, `{83,89,21}` certifies `591` rows.
The first miss is not mysterious; it has a small replacement certificate
`D=19` with `4` unit numerators.  So the right object is not a fixed finite
certificate basis.  It is a finite residue atlas that is allowed to move when
the row has killed a modulus by divisibility or resonance.

The character-count scout gives the clean language.  For each `D`,

```text
N(S,D) = #{a mod D : gcd(a,D)=1 and all s*a/D are safe}
```

has a main term plus the resonance packet

```text
sum_s k_s s == 0 mod D.
```

The main term can be healthy while a specific denominator has `N=0`; that zero
is the local resonance/divisor event.  A basis is useful because several
denominators rarely dip together in ordinary samples.  A divisor-loaded tail
manufactures simultaneous dips, so a fixed basis cannot be global.

The apex floor also becomes cleaner.  Covering kills every reduced denominator
`2..14`, not just `14`.  This is the exact boundary between THM-523's easy
non-covering split and the hard covering core.  Once every `q<=14` divides
some speed, all those small rational observers are dead.

Assumption challenge:

- Runners are the wrong tournament vertices for this question.
- Speed-subset deletion is the wrong minor order, because loneliness is not
  deletion-monotone.
- Useful candidate vertices were denominators, numerator residues, character
  modes, divisor-loaded tails, covering obligations, even-graph holes, and
  proof obligations.
- The selected carrier is the scaled residue/character-count atlas.  It
  preserves the witness predicate exactly at each modulus and records when
  information is destroyed by a quotient.

The odd-hole / `{7,21}` / `E_7` analogy should be kept as an address diagnostic.
It may identify the low-resonance shapes that create simultaneous denominator
dips.  But HYP-2872 remains the guardrail: an equinumerous even-graph quotient
is not a proof carrier until it preserves the quantity being proved.

The next proof target should be phrased as:

```text
coherent divisor/resonance packets -> scaled finite atlas,
incoherent packets -> HYP-2875 bandlimited L2 tail.
```

That aligns the residue-basis idea with the current spectrum-sum route instead
of trying to revive a bounded-denominator shortcut that THM-566 already killed.

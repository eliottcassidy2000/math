# Pi Sum/Product/Fraction Trinity (S662)

The useful move was to stop treating `pi` as a scalar and treat it as a carrier
with three complementary representations.

```text
sum      -> additive packet / moment / cancellation order
product  -> local factor / zero / norm / sieve channel
fraction -> recursive boundary / convergent / owner state
```

HYP-2229 already had the first two faces in a strong form: the sine product
keeps elementary factors, and Newton/log-derivative turns them into zeta power
sums.  The missing face is the continued fraction.  Brouncker is numerically
slow in the small audit, but that is not a weakness for the carrier metaphor:
it makes the recursive boundary state explicit.

The S662 script makes the representation faces into a tournament.  The metric
triples are:

```text
sum         (3,1,2)
product     (2,3,1)
fraction    (1,2,3)
raw_decimal (0,0,0)
```

The result is exactly the right shape:

```text
sum -> product -> fraction -> sum,
```

with all three beating `raw_decimal`.  The tournament has one directed
3-cycle, SCCs `{sum, product, fraction}` and `{raw_decimal}`, and three
Hamiltonian paths.  So the representation trinity is not a ranking.  It is a
small nontransitive proof lens: each face repairs what another forgets.

The LRC `n=14` transfer is the best immediate use.  S661 asks for a no-leak
owner-derivative theorem, and S660 says unpaired derivative decks leak.  The pi
trinity says what to keep:

```text
odd-wall / pair-sum packets        (sum face)
C=27 gcd shell / local obstruction (product face)
carry-owner derivative recursion   (fraction face)
```

In that language, a floor atom is not just a scalar maximin value.  It should
be a triune record:

```text
(additive walls, local product shells, recursive owner state).
```

This also gives a way to revisit OCF.  `H(T)` is the scalar evaluation, the
independence polynomial is the retained packet object, and a continued
fraction or continuant encoding might remember boundary state during deletion,
substitution, or Hamiltonian-path gluing.  That is the new tangent worth
testing: not "find another formula", but "which formula face preserves the
side channel the proof needs?"

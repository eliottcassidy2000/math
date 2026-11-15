# Chan–Pak: the carrier checksum is an odd-cycle count, and why literal span-3 cones need rank

Session klein-2026-11-15-ocf. **Chan–Pak remains OPEN.** Owner prompt: build the syllable-adaptive cone certificate
for literal span 3; look for creative, non-branching routes; ask what the novel mathematics fundamentally represents.
(The prompt referred to connections "below"; none arrived with the message, so past themes were used.)

1. **[The unifying frame](UNIFYING.md).**
   - *Known.* THM-002 (OCF): H = I(Ω, 2), where Ω is the odd-cycle conflict graph (Grinberg–Stanley).
   - *New, the apex identity.* A transitive tape with one apex has H = 1 + p(u), where u marks the positions that beat
     the apex. So the carrier checksum p is twice the number of odd cycles through the apex. Exhaustive for tapes of
     length ≤ 9.
   - *Consequences.*
     - p is even: Rédei's theorem.
     - THM-5232's reverse-complement invariance is converse symmetry, so reverse-complement collisions are partial
       converses.
     - Carrier-row surcharges count the odd cycles inside a tape.
   - *Reach.* One apex (equivalently, odd cycles with a common vertex) realizes a stable 50% of the odd numbers
     (misses 7, 13, 21, 29, …).
   - *Program.* Chan–Pak is surjectivity of I(·, 2) on realizable conflict graphs. A closure lemma for two-apex
     overlaps would replace more branches.
2. **[Literal cones](LITERAL_CONES.md).** The syllable-adaptive certificate **cannot work**, which corrects the
   previous plan.
   - *Wrong patterns.* About 99% of the cheap cycles that need long memory are not syllable concatenations.
   - *Uniform memory stalls.* Even six code units of uniform memory leave a cycle of growth 1.38 uncut, below the
     floor 5^(1/3).
   - *What is needed.* The literal language has slow bottlenecks. A certificate must track the rank of the prefix
     product as a monotone potential; that design is specified.

## Replay

| command | checks |
|---|---|
| `python3 oddcycle.py` | the odd-cycle form and the apex identity |
| `python3 apexspec.py` | the one-apex spectrum |
| `python3 memneed.py 1.765`, `python3 sylparse.py`, `python3 rhocrit.py` | the literal-language data; these read `span3_intrep.txt` from `ribbon_spectrum_20261104`, or from `$SPAN3_INTREP` |

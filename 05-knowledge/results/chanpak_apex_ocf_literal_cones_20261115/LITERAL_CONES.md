# Literal span 3: why syllable-adaptive reachable cones cannot certify the growth bound

The task was to build the syllable-adaptive cone certificate proposed in
[SPAN3_LITERAL.md](../chanpak_mod8_synthesis_span3_20261115/SPAN3_LITERAL.md) §3.
- *Setting.* Bottleneck-free literal span-3 words are words over the 21-word code {0^i a : i ≤ 2, a ∈ {1, …, 7}};
  bottleneck-free means rank ≥ 2 in the minimal 62-dimensional transfer.
- *Cones.* Reachable cones (images of the nonnegative 119-state orthant) are invariant under any rule that sends a
  state to a suffix of its extension.
- *What decides it.* So the whole question is which tag set (state memory) makes every cycle of the automaton grow at
  least ρ.

**The measurements show that plan cannot work, and say what can.**

## 1. Language data (modular ranks)

- 18 of the 21 code units are bottleneck-free; the three "00a" units have rank ≤ 1.
- There are 305 bottleneck-free code pairs and 5,186 triples.

## 2. Cheap cycles are mostly not syllable concatenations

Take the 31,648 cyclic code words of up to 4 units that grow more slowly than 1.765. Every one is a bottleneck: its
powers reach rank ≤ 1. Its *bottleneck length* is the number of code units after which that happens, the worst case
over rotations; the distribution runs from 1 to 12.

Few cycles that need memory of 5 or more units parse as cyclic concatenations of the 13 floor syllables
(`memneed.py`, `sylparse.py`):

| bottleneck length | 5 | 6 | 7 | 8 | 9–12 |
|---|---:|---:|---:|---:|---:|
| syllable-composed | 15 | 19 | 11 | 4 | 4 |
| other | 3,546 | 1,905 | 582 | 293 | 44 |

Memory adaptive along the syllables therefore leaves about 99% of the long cheap patterns uncut. **The plan of
SPAN3_LITERAL.md §3 was wrong.**

## 3. Uniform memory stalls far below the floor

The cheapest cycle a uniform memory of m code units leaves uncut, meaning its bottleneck length is above m + 1
(`rhocrit.py`):

| memory m (code units) | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| cheapest uncut growth | 1.000 | 1.000 | 1.170 | 1.246 | 1.308 | 1.380 |
| uncut cheap cycles | 29,060 | 21,843 | 6,423 | 2,862 | 938 | 345 |

Examples of the cheapest uncut cycles: 2|7|01 at m ≤ 2 (growth 1.000), and 4|01|7|6 at m = 6 (growth 1.380).

Even six code units (up to 18 letters) of memory cannot certify more than 1.38. That is far below the floor 5^(1/3).
The literal language has **slow bottlenecks**: cycles of growth as low as 1 that need up to 12 code units to become
rank one. A bottleneck-free word can contain one bounded stretch of such a cycle, and finite-suffix memory cannot
tell when the stretch must end.

## 4. What a certificate needs: rank as a potential

Along a word the rank of the prefix product never increases. It stays between 2 and 62 for bottleneck-free words.
Every slow-bottleneck stretch that is long enough lowers it. So rank is a bounded, monotone potential.
- *The automaton to build.* Its states should be (short tag, current rank r ≤ 62). Transitions may lower r. The cones
  should be the reachable cones restricted to products of rank r. An LP would then charge each rank drop a bounded
  loss, and require growth ρ between drops.
- *Why it should work.* At most 60 drops can occur, so the bound H ≥ C·ρ^n survives, with C absorbing the drops.
- *The open step.* Making the rank-r transition relation exact needs the image subspaces. Over-approximating it
  needs a proof that every cheap stretch lowers the rank, which is a statement about the 62-dimensional semigroup
  not checked here.

This is the precise remaining obstruction for HYP-9742 in literal span 3. The periodic evidence (no cheap periodic
word up to 12 letters is bottleneck-free) still supports the conjecture.

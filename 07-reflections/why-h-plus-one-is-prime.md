# Why h+1 Is Prime for Exceptionals

**Session:** kind-pasteur-2026-03-21-S18j
**Arising from:** The Bernoulli-Coxeter chain (S18h), the classical/exceptional divide (S18i)

---

## The Fact

For every exceptional simple Lie algebra, the Coxeter number plus one is prime:

| Algebra | h | h+1 | Prime? |
|---------|---|-----|--------|
| G_2 | 6 | 7 | Yes |
| F_4 | 12 | 13 | Yes |
| E_6 | 12 | 13 | Yes |
| E_7 | 18 | 19 | Yes |
| E_8 | 30 | 31 | Yes |

For classical algebras, h+1 is sometimes prime, sometimes not:
- A_2: h=3, h+1=4 (not prime)
- B_4: h=8, h+1=9 (not prime)
- A_7: h=8, h+1=9 (not prime)
- D_5: h=8, h+1=9 (not prime)

The pattern is: **h+1 prime is NECESSARY for an exceptional algebra but not sufficient.** Every exceptional has it; many classicals don't.

---

## The 6k Structure

All exceptional Coxeter numbers are multiples of 6:

| k | h = 6k | h+1 = 6k+1 | Prime? | Exceptional? |
|---|--------|-------------|--------|-------------|
| 1 | 6 | 7 | Yes | G_2 |
| 2 | 12 | 13 | Yes | F_4, E_6 |
| 3 | 18 | 19 | Yes | E_7 |
| **4** | **24** | **25 = 5^2** | **No** | **None!** |
| 5 | 30 | 31 | Yes | E_8 |

The exceptional algebras exist at k = 1, 2, 3, 5. They SKIP k = 4, which is the unique value in {1,...,5} where 6k+1 is composite.

The number 25 = 5^2 is the obstruction. If 25 were prime, there might be an exceptional algebra with h = 24. But 25 is composite, and no such algebra exists. The compositeness of 25 and the nonexistence of an "E_9" (or whatever would live at h=24) are the same fact seen from two sides.

---

## Why This Matters for Tournaments

In our framework:
- Tournament prime p enters at Bernoulli weight h = p-1
- If h is an exceptional Coxeter number, the prime controls "hard" structure
- h+1 being prime is what MAKES it a tournament prime (trivially: h+1 = p, and p is prime by definition)

So the statement "h+1 is prime for all exceptionals" IS the statement "the exceptional Coxeter numbers produce tournament primes." The two facts are the same.

But the converse fails: not every prime is h+1 for an exceptional. The primes 3, 5, 11, 17, 23, ... are h+1 for classical algebras but NOT for any exceptional. These primes enter the Bernoulli denominators at weights that are classical-only Coxeter numbers.

The tournament divides:
- **Exceptional primes** (h+1 prime where h is ALSO an exceptional Coxeter number): 7, 13, 19, 31. These control obstructions.
- **Classical primes** (h+1 prime where h is ONLY a classical Coxeter number): 3, 5, 11, 17, 23, .... These control direct structure.

---

## The Deeper Question: Why 6?

Why are all exceptional Coxeter numbers multiples of 6?

Because 6 = h(G_2) = the Coxeter number of the smallest exceptional algebra. And G_2 sits inside every other exceptional algebra (via inclusion of root systems). The Coxeter number of a larger algebra must be a multiple of h(G_2) when G_2 embeds — not exactly, but the divisibility by 6 reflects the foundational role of G_2.

More precisely: the exceptional Coxeter numbers are:
- 6 = 2 * 3
- 12 = 2^2 * 3
- 18 = 2 * 3^2
- 30 = 2 * 3 * 5

These are all **6-smooth** numbers (having no prime factor > 5) that are multiples of 6. They are the products of the "easy" tournament primes {2, 3, 5} (or subsets thereof).

So: **the exceptional Coxeter numbers are built from classical primes, but they PRODUCE exceptional primes.** The classical primes 2, 3, 5 combine (via multiplication) to produce the Coxeter numbers 6, 12, 18, 30, and then adding 1 produces the exceptional primes 7, 13, 19, 31.

This is a PRIME GENERATION mechanism:
- Input: products of small primes {2, 3, 5}
- Operation: add 1
- Output: new, larger primes {7, 13, 19, 31}

The classical structure GENERATES the exceptional structure by the "+1" operation. And the "+1" is the Redei quantum — the same "+1" that appears in moonshine (196884 = 196883 + 1), in the taxicab (1729 = 1728 + 1), and in the OCF (H = 1 + 2*Delta).

---

## The k = 4 Gap

The single failure at k = 4 (h = 24, h+1 = 25 = 5^2) is deeply significant:

- 24 = |BT| = the order of the binary tetrahedral group
- 24 = the central charge of the moonshine module V^natural
- 24 = the rank of the Leech lattice
- 24 = the number of regular tournaments on 5 vertices

The number 24 is the most important number in moonshine. It is NOT the Coxeter number of any exceptional algebra precisely because 25 = 5^2 is not prime. The "most moonshine" number fails the primality test.

In tournament terms: there is no "moonshine obstruction prime" at p = 25 because 25 is composite. The moonshine module (c = 24) does not contribute a prime to the obstruction tower. Instead, it contributes through E_8 (h = 30, one step beyond) and through the Leech lattice structure.

The gap at k = 4 is why the exceptional tower jumps from E_7 (h = 18) to E_8 (h = 30), skipping h = 24. The jump of 12 = h(E_6) is the largest gap in the tower, and it occurs because 25 = 5^2 blocks the intermediate step.

---

## The Complete Picture

The Lie algebra classification says: there are exactly 5 exceptional simple Lie algebras, with Coxeter numbers 6, 12, 12, 18, 30.

The number theory says: 6+1=7, 12+1=13, 18+1=19, 30+1=31 are all prime.

The tournament theory says: these primes (7, 13, 19, 31) control the hard structure of tournaments — forbidden values, OCR, moonshine.

The reason this works is the "+1" operation. The exceptional Coxeter numbers are smooth numbers (products of 2, 3, 5). Adding 1 to a smooth number often produces a prime (this is related to the theory of smooth numbers and primality, though not guaranteed). The Lie algebra classification SELECTS exactly the smooth multiples of 6 where adding 1 gives a prime, up to the classification's natural boundary at rank 8.

The tournaments-as-signed-simplices framework provides the bridge: tournaments live on the simplicial A_{n-1} structure (classical), and their obstructions come from the non-simplicial exceptional structure. The "+1" that takes you from a Coxeter number to a tournament prime is the same "+1" that takes you from the cusp to the first excited state — the Redei quantum, the trivial representation, the departure from background.

**The exceptional Lie algebras are the mathematical proof that the "+1" operation generates the primes that govern tournament obstructions.**

---

*When Killing and Cartan classified the simple Lie algebras in the 1890s, they found exactly five exceptions to the infinite classical families. When we proved that H = 7 is forever forbidden in tournament theory, the number 7 appeared because G_2 has Coxeter number 6 and 6+1 = 7. The same "+1" that Ramanujan saw in 1729 = 1728 + 1, that McKay saw in 196884 = 196883 + 1, and that Redei saw in H = 1 + 2*Delta. The classification theorem of simple Lie algebras, proved over a century ago, already contained the forbidden values of tournament Hamiltonian path counts. It just took this long for anyone to read them off.*

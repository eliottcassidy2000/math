# Krenn pairing and the no-cancellation grammar

**kind-pasteur-2026-07-27-S136.** Provenance note, not truth source.
Seed: Mario Krenn's 3000-euro question (bicolored complex-weighted
graphs; perfect matchings induce inherited vertex colorings; the
question is whether monochromatic sectors can survive with all
mixed sectors cancelling, at given dimension/vertex parameters).

## 1. The dictionary

Krenn's frame is the repo's involution-parity grammar (S131
reflection) wearing quantum optics clothes:

```text
Krenn                          repo
-----                          ----
mixed-color PM cancellation    free involution pairing (sigma on
                               black lines; conjugation on complex
                               fibers; swap on parent pairs)
monochromatic survival         fixed-point locus carries the count
                               (blues odd; all-blue sector; R(k)
                               parity = half-membership)
inherited vertex coloring      per-node colour census of incident
                               structure (pure-blue/mixed/pure-black
                               IS an inherited coloring of nodes by
                               their line fibres)
no-cancellation obligations    LRC graft survival; klein's nonzero
                               C7xC13 tensor cells; HYP-9027 odd
                               Jelonek exponents
```

The hard direction on both sides is identical: proving a
monochromatic/odd sector does NOT cancel. Krenn's partial results
(max-degree <= 3; d > n/sqrt2 exclusions) are inventory arguments
of exactly the THM-2465/THM-1370 type.

## 2. A new mini-theorem from the dictionary (verified)

**Line-multigraph pairing law.** In the merged line multigraph
(nodes = merged classes M_n, edges = complement lines, loops
excluded, coloured blue/black), sigma fixes every node and every
blue line and freely pairs parallel black lines (same endpoints).
Since a perfect matching cannot contain two parallel edges, sigma
acts on PMs preserving the black-count, with fixed points exactly
the all-blue PMs. Hence EVERY black-count stratum with k >= 1 has
even cardinality; and since blue lines touch only SC nodes, any NS
node kills the all-blue sector:

```text
for all n >= 5:  #PM(line multigraph of M_n) is EVEN, stratum by
stratum, with the all-blue (monochromatic) sector empty.        (1)
```

Verified exactly at n = 5 (M = 10 nodes, 32 lines, 8 blue): PMs =
{2 black: 8} -- total 8, every PM uses exactly TWO black lines
(the two NS nodes must be black-matched; the SC rest matches
all-blue), all strata even, all-blue zero. In Krenn language: the
merged metagraph is ANTI-monochromatic -- its mixed sectors pair
off destructively mod 2 and its monochromatic sector is
structurally empty. The lonely survivor grammar is inverted
relative to Krenn's GHZ goal, which is why the repo's hard
obligations are NO-cancellation statements while Krenn's is a
cancellation-design statement: the two problems are the two signs
of one pairing calculus.

## 3. Speculative transfers (flagged, untested)

- The LRC noncirculant-graft no-cancellation could be posed
  Krenn-style: colour contribution words by owner; survival = the
  monochromatic owner sector staying nonzero under the evolution
  pairing. If the graft's pairing involution can be shown to have
  a fixed point in the owner sector, survival follows the same way
  all-blue carries parity here. One session to formalize against
  THM-2445's owner conditioning.
- Krenn's d > n/sqrt2 bound smells like a rank/dimension inventory
  (matchings needed vs available) -- the same bookkeeping as the
  BKK/parameter counts in THM-2446 SS4. Not pursued.

Files: verification inline in this reflection's commit
(krenn_line_pm_pairing test run, n=5 exact); the structural proof
of (1) is the two-line sigma argument above and needs no census.

# CRT fiber products, not modular cusps, are the exact clock carrier

The failed HYP-8880 analogy was still pointing at a real operation, but it
named the wrong object. THM-2057 never needs a modular curve. It chooses a
phase on a core clock and a phase on a tail clock and asks whether those two
residues lift to one integer. That is a generalized-CRT fiber product.

THM-2059 makes the carrier lossless. For a row `aC union {w}` and a proposed
clock `N`, form the safe core packet modulo `N` and the safe tail packet modulo
`h=Na/gcd(w,Na)`. Reduce both packets modulo `d=gcd(N,h)`. Their compatible
pair count is exactly the dot product of the two reduction histograms, and
after the deterministic lift multiplicity it is exactly the number of safe
grid phases `k/(Na)`.

This is more than a repackaging of the small-clock proof. THM-2057 used
`N<=14`, where every nonzero residue is automatically at threshold. The CRT
formula has no such restriction: it computes the surviving subset for any
`N`. The audit found `34854` certificates with `15<=N<=50` in the explicit
larger-clock family tested, while also replaying the small missing-clock
specialization on `44761` rows.

The decisive information is not either packet size. Two nonempty packets can
reduce to disjoint classes modulo `d`. The correct sidecar is the pair of
histograms, or equivalently the bipartite compatibility graph. This also
explains why forcing a tournament here would be destructive: compatibility is
symmetric and bipartite, and the ties are the theorem.

## How this changes the finish target

THM-2058 currently reserves a primitive phase-order and longitudinal-interval
carrier but remains a claimed stub in the checked tree. Its eventual packet
counts should feed directly into THM-2059's `alpha` histogram. The tail
interval supplies `beta`; the CRT dot product becomes the exact acceptance
test. This suggests the following finite pipeline inside every THM-2053 deck:

```text
transverse template
  -> primitive core packets by phase order
  -> reduction histograms modulo gcd(N,h)
  -> tail interval histograms
  -> positive CRT dot product, or an explicit disjoint-support obstruction
  -> on obstruction, split at the next owner/divisibility wall.
```

The obstruction is now concrete: a failed clock is not merely `P_N=0`, but a
pair of disjoint residue supports. Those supports can be transported through
Farey cells and compared under deletion. That is a plausible place for the
signed wall-word/Euler machinery to act, because deletion changes one packet
by an explicit set intersection rather than changing an unnamed modular
object.

## Assumption challenge

The useful quotient changes with the operation. Kelvin/Farey geometry needs
hull owners; the clock layer needs residue histograms; Euler deletion needs
owner-labelled safe components. None of runners, primes, cusps, total packet
sizes, or a forced tournament retains all three. The creative route is a
small typed pipeline of lossless quotients, with each handoff checked by an
exact formula.

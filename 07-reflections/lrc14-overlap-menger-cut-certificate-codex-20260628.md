# LRC14 Overlap-Tax Menger-Cut Certificate Reflection

HYP-3437 is the first place where HYP-3434's overlap tax stops being just a
positive correction term and becomes a finite object.

The useful object is the atomized one-branch arrangement:

```text
E_safe component
  -> endpoints from even gates and branch-0 odd bad intervals
  -> atomic cells with odd-blocker multiplicity
```

Then the algebra is exact:

```text
multiplicity 0  -> branch0 survivor
multiplicity 1  -> ordinary scalar cover mass
multiplicity >=2 -> overlap-tax mass
```

The evidence is stronger than the previous harmonic-sieve statement.  On the
same `150` rows, the `59` negative naive-slack rows all have a strict rescue
cut.  The rank histogram is:

```text
{2:7, 4:2, 5:48, 6:2}
```

The two rank-`6` rows are the canonical AP-with-84 row and its duplicate.  That
matters: the worst rank is exactly the row already controlled by HYP-3431's
canonical corridor-fence proof.  The apparent noncanonical obligation is lower
rank in this bank.

Proof obligations that now look worth isolating:

1. Prove the canonical rank-`6` core `(3,5,7,9,11,13)` is the HYP-3431
   corridor-fence case, not a generic obstruction.
2. Prove noncanonical negative-slack rows have a rank-`<=5` odd-blocker
   overlap core, or name the first exception.
3. Attach endpoint labels to the cut core so HYP-3429's rank-2 survivor spine
   and HYP-3427's wall alphabet can take over after the cut is found.
4. Keep HYP-3436's two-color bad-core extractor separate: HYP-3437 is
   one-branch overlap tax; HYP-3436 should classify simultaneous branch
   failure.

The creative analogies are now disciplined.  Schwarz-Christoffel prevertices
are exactly the interval endpoints.  Bring radical language names branch
monodromy debt.  BDH/Mertens estimates are allowed only as averaged statements
over already labelled cut cores.  None of those routes can replace the atom
sidecar.

The main caveat: this is still a finite stress bank and a brute-force subset
search.  The next proof-facing move is to derive the rank bound from interval
geometry, endpoint order, and the two-adic split, not from enumeration.

# Colorings and Kakeya amplifiers: the 91-stalk closure is an anti-Kakeya theorem

**kind-pasteur-2026-07-27-S136.** Provenance note, not truth source.
Continuation of the Krenn-pairing dictionary; seed directive: see
everything as colorings, merge in Kakeya amplifiers.

## 1. The coloring inventory (one table, cross-thread)

```text
coloring                        where               law it carries
--------                        -----               --------------
blue/black (sigma-fixedness)    line multigraph     T-join/even; PM parity (S136)
inherited node coloring         pure-blue census    alphabet law THM-2454
defect coloring m_Y(h,r)        91-stalk rows       THM-2436 'cannot stay vertical'
owner coloring                  LRC graft words     THM-2445 conditioning; the
                                                    open no-cancellation
residue channels mod 5          twin ancestry       dead-channel law (exact)
chi_7 x chi_13 characters       stalk harmonics     mixed unit character; klein
                                                    tensor cells
quadratic-residue coloring      Paley/T_p           the nonrigid pure-blue atom
```

The repo's strongest moves recolor a set until exactly one color
class carries the obstruction, then pair off the rest (Krenn
grammar) or floor the survivor (amplifier grammar).

## 2. The stalk closure, in Kakeya language

Ordinary words are LINES in Z/91 (13-term APs, unit step);
blockers are vertical lines; the guard is a 26-slab. Verified
descent facts (three-line computation, this session):

```text
every unit-step 13-term AP covers Z/13 EXACTLY ONCE;
the guard covers each F_13 class EXACTLY TWICE;
so guard + five words have mod-13 multiplicity 2+5 = 7 = |F_7|
on every class -- THM-2436's row identity IS a Kakeya descent:
the entire covering defect lives in the F_7 coordinate.        (1)
```

THM-2431/2436's engine is then an **anti-Kakeya theorem**: five
lines in five DISTINCT directions (plus slab and verticals) can
never cover the punctured stalk -- a repeated direction is forced
(distinct_cover = 0, two independent engines). Classical Kakeya
says many-directions forces LARGE sets; the stalk says the
punctured set is too large for few-distinct-directions to cover --
the covering-dual statement. The forced repeated direction is then
fed to the rounding transfer: the amplifier step, where ONE
certified coincidence (the repeated pair's fixed integer lift n)
amplifies into a 91-cell measure cap. That two-step shape --
anti-Kakeya forcing, then amplification of the forced coincidence
-- is the deep-c_3 closure's actual skeleton, and it is the same
skeleton as the bush argument in Kakeya proper (one intersection
point amplified through every needle it meets).

## 3. Amplifier grammar, catalogued

```text
certified nonzero              amplified to
-----------------              ------------
repeated-step lift n != 0      mu(P) <= |C_d|/91 (THM-2436)
one mixed unit character       parent mass >= 1/252
klein tensor cell J*A != 0     owner-resolved mass >= delta/338
odd real collision fiber       d odd for the wild map (HYP-9030)
one blue tiling (Redei)        blues-odd everywhere (THM-643)
```

Open transfer worth one session: a Dvir-style polynomial method
over F_7 for the descended covering problem -- the defect coloring
is a function on 13 rows x 7 columns with row sums fixed; a
degree argument on F_7 might re-prove the repeated-direction
forcing without enumeration, upgrading the atlas's exhaustive
verification to mechanism. Flagged for the LRC lane (codex owns
the graft; this would be a sidecar, not a collision).

## 4. Honest boundary

Z/91 is not a field; Dvir does not apply directly, and the descent
(1) is exactly what routes around the compositeness (13-side
trivializes, 7-side is a field). No claim is made that the
polynomial method WILL close it -- only that the reframe puts a
classical tool within one session's test of a lane that has so far
run on exhaustive atlases.

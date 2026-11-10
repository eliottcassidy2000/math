---
id: THM-5100
title: "Sharp strong elevenfold symmetry gap from the complete two-path-cover census"
status: "PROVED from inherited prime-action reduction and independently audited FINITE-EXACT inputs. Chan-Pak completeness OPEN."
session: Session240, 2026-11-10
source: 05-knowledge/results/session240_seven_eleven_symmetry_20261110/gaps/PROOF.md
audit: 05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_gap/AUDIT.md
---

# THM-5100: the next strong elevenfold carrier costs 1,848,759 paths

[Full proof and census](../../05-knowledge/results/session240_seven_eleven_symmetry_20261110/gaps/PROOF.md)
· [Independent audit](../../05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_gap/AUDIT.md).

Let H(T) count ordered directed Hamiltonian paths. If a finite strong tournament
T admits an actual automorphism of order11, then either

- |T|=11 and H(T) is one of92411,93027,93467,95095; or
- |T|>11 and H(T)>=1848759.

The second threshold is sharp at13 vertices: substitute the eleven-circulant
with connection set{2,3,4,5,10} and two singletons into a directed triangle.
The inherited general prime-action reduction is THM-5099
(`value-refined-scc-recurrence-and-sharp-seven-symmetry-gap`). With two or more
moving eleven-orbits it gives H>=92411^2+22=8539792943; with one, the orbit is
a module and a quotient triangle supplies an induced C3[R,1,1]. The inherited
path-cover kernel gives H(C3[R,1,1])=3H(R)+F2(R), where F2 is the ordered
two-path-cover count. The new complete32-orientation table has rows

    (H,F2,multiplicity) =
    (92411,1571526,10), (93027,1574958,10),
    (93467,1577598,10), (95095,1586310,2).

The core minimum92411 and sevenfold threshold2163 are inherited, not new.
Normal/-O producer checks, independent walk inclusion–exclusion and direct
inflation DP agree. A separate uncached DFS counts the sharp witness.

In particular, below1848759, a strong T with |T|!=11 has no automorphism-group
factor11. Equality is classified within the triangle-inflation family only;
no complete ambient equality classification is claimed. Strongness and the
actual action are essential. No stronger full-group denominator, arbitrary
H realization, or Chan-Pak completeness follows.

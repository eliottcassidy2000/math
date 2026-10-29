# Quasirandomness with arithmetic structure and exact local clocks

2026-10-29. **PROVED scoped deductions / CITED analytic and group inputs /
FINITE-EXACT checks / OPEN extensions.** Same-agent proof review; no claim
of literature novelty or a solution of Collatz.

Concurrent [Session107](../session107_quasirandom_arithmetic_20261029/characters/CHARACTER_STRUCTURE.md)
independently proves and audits the same fixed-colour theorem, with a stronger
unequal-fibre estimate and signed bases. It also independently recovers the
group and Euclid-offset mechanisms. Those shared results are credited jointly.
This packet adds the explicit sigma/Hamming duality, all-odd-term Jacobi
extension, composite-base valuation map and all-level19-adic clock proof.

Quasirandomness here means that large vertex sets have no macroscopic
orientation bias, and every fixed tournament has its random limiting
frequency. It can coexist with exact algebraic relations in a separate
observer. The useful distinction is between **relations to retain** and
**remaining modes whose cancellation must be proved**.

The [proofs](PROOFS.md) connect the requested numbers through concrete maps:

| Object | Mechanism and consequence |
|---|---|
|Prime-prefix Legendre tournament|A quadratic large-sieve norm bound suppresses large coherent cuts. This is inherited P1.|
|Fixed square-class decorations|**Proved here:** their feasible code patterns are equidistributed, jointly with every fixed tournament pattern. Named fixed prime attachments are included.|
|{4,2,5,3,7,13}|The sigma-character observer has3 independent bits,8 patterns and3 exact relations. Direct character tests instead have5 independent bits.|
|{2,5,7,19,29,59,89}|The old square-total Hamming code and the prime-character simplex code are literal orthogonal duals. Pairwise independent bits have exact Fano triple correlations.|
|2,3,7,43,1807,3263443,...|**Proved here:** all odd Sylvester terms form a transitive Jacobi tournament. The future residue1 creates perfect alignment, defeating quasirandomness.|
|13 and139|Their product is Phi_6(43);43 has order6 modulo both. This is an Eisenstein discriminant/order condition.|
|19 and2->3->7->5->2|**Proved here:** the seed's periods modulo19^k are4,4 and12*19^(k-3) for k>=3. A compatible19-adic four-cycle is a different source.|
|1093 and3511|The inherited doubling tower delays one digit of order growth. Paley3511 still has exact square-root cancellation.|
|Fermat numbers and primes|Later terms are2 modulo earlier terms. Fermat prime rigidity and absence of Wieferich delay do not imply mixing.|
|Monster and Baby Monster|Gowers's representation-degree theorem gives a genuine quasirandom-group route. The carriers are Cayley incidence graphs, not tournament automorphism groups.|
|3n±1|FLT supplies local inverses and finite multiplicative orders. Actual Collatz paths still require source, ordered carry and valuation clock.|

The same modulus appears in two retained systems:

    f^4(2)-2 = 1805 = 5*19^2 = 2^11-3^5.

An inherited Collatz word has numerator475 and denominator1805, missing
one19 factor needed for integrality. This shared arithmetic gives a useful
precision test; it does not give a conjugacy of the maps.

The six-label Monster display also has an honest integer interpretation:
extend v_p(R) to e_R(b)=max{j:b^j divides R}. At R=|M|/|B| it gives the
five-cycle4->2->5->3->7->4 and13->2. Passing to the true involution-class
index R/2 changes it to4<->2. The central factor must remain visible.

The next targets are uniformity for a growing decoration bank, prime-factor
selection from Sylvester terms, and joint lifting of clocks with actual
response/carry data. Fixed-rank equidistribution does not prove these.

[AUDIT.md](AUDIT.md) records recovery, primary sources, scope and reproduction.
The standard-library [check.py](check.py) verifies exact CRT fibres, the
Hamming/simplex duality,55 Sylvester pairs, Paley Gram matrices through139,
the19 clock through19^6, and the Wieferich/Fermat controls. Empirical colour
counts use4806 admitted primes below100000; the asymptotic theorem uses
its analytic proof, not this sample. The81/27 asymmetry survives within
every fixed feasible colour class.

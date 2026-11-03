# Surjective seven-adic Hamilton observation in the transported family

**Proved.** For every integer m≥1 and every r modulo 7^m, infinitely many k≥0 satisfy h(q2^k U)≡r (mod 7^m). Consequently each of the two actual strong-tournament families with Hamilton counts 3h(q2^k U)+2 and 3h(q2^k U)+4 separately hits every residue modulo 7^m, and every odd class modulo 2·7^m, infinitely often. This is congruence coverage, not equality to a prescribed integer or primality of a quotient.

Use the notation of [TRANSPORT, sections T3 and T8](https://github.com/eliottcassidy2000/math/blob/f6d5304849e50b29df6e368239e92fffa92e57e3/05-knowledge/results/golden_inverse_defects_20261103/TRANSPORT.md). Set A=[[11,18],[4,7]], z0=(25,34), z_k=z0 A^k=(Y_k,X_k), and h_k=Y_k+6. Modulo 7, A has a single orbit on the 48 nonzero vectors (TRANSPORT section T8, or the included direct check). A short integer matrix calculation gives

    A^48 = I+7C,     C mod 7 = [[6,6],[6,0]].

For each r modulo 7, select a nonzero vector z=(Y,X) with Y+6=r and Y+X≠0. There are at least six choices: when Y=0, choose any nonzero X; when Y≠0, exclude just X=−Y. The 48-cycle gives a corresponding k0. The derivative of Y under the 48-step clock is

    D = (z C)_1 = 6(Y+X) ≠ 0 (mod 7).

For s≥1, binomial expansion gives (I+7C)^(7^(s−1))=I+7^s C (mod 7^(s+1)). Thus, for any current index k≡k0 (mod 48), adding j·48·7^(s−1), 0≤j<7, changes h_k modulo 7^(s+1) by j·7^s D and does not change its earlier residue. As D is a unit, exactly one j supplies the next desired base-7 digit. Induction constructs k modulo 48·7^(m−1) for any target residue mod 7^m. Repeating k by multiples of 48·7^(m−1) leaves h modulo 7^m unchanged; the h_k are unbounded (the positive transfer strictly grows), giving infinitely many distinct counts.

Each target y modulo 7 has exactly six suitable k0 (including y=0, because the missing zero vector already accounts for the unique derivative zero). Thus at least 6 distinct k residues modulo 48·7^(m−1) map to each specified h residue mod 7^m; the count is a lower bound because the six remaining singular base branches were not analyzed.

Both transported count maps have slope 3, a unit modulo every 7^m, so surjectivity transfers. Every actual tournament Hamilton count is odd; CRT gives the odd classes mod 2·7^m. Source-word replacement and even triangle parity guards remain as in TRANSPORT: here the prefix contains no b, so its parity is even. No closure of output states under repeated ternary lifts is asserted.

**Exact valuation corollary.** For either companion and any t≥0, infinitely many family members have exactly 7-adic valuation t: prescribe its residue to be 7^t modulo 7^(t+1). The two companion counts differ by 2, so 7 never divides both for the same k. More generally, the pair of residues is exactly the affine diagonal {(r,r+2): r∈Z/7^m}; it is not an independently selectable pair.

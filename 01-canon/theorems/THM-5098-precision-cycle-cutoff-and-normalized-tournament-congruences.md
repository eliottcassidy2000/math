---
id: THM-5098
title: "Precision-dependent cycle cutoff and mod-four/mod-eight normalized tournament-count congruences"
status: "PROVED for all stated odd orders and precisions, independently audited; exact computations are additional controls."
session: Session239, 2026-11-10
source: 05-knowledge/results/session239_sequence_engines_20261110/precision/PROOF.md
audit: 05-knowledge/results/session239_sequence_engines_20261110/audit/README.md
---

# THM-5098: choose precision before enumerating cycle types

[Full proof](../../05-knowledge/results/session239_sequence_engines_20261110/precision/PROOF.md)
· [Independent audit](../../05-knowledge/results/session239_sequence_engines_20261110/audit/README.md).

For odd n>=1 set U_n=A000568(n)/2^((n-1)/2). Its inherited integrality and
oddness for n>=3 are THM-305 (`two-adic-valuation-tournament-count`); n=1 is direct. For an odd
partition with r cycles and multiplicities m_j, the normalized contribution
has valuation

    delta=sum_(i<j) gcd(l_i,l_j)-(r-1)/2-sum_j v2(m_j!)
         >= (r-1)^2/2-v2(r!).

Hence b bits require only r with the right side<b; r=O(sqrt(b)), independently
of n. This gives a polynomial-in-n algorithm at fixed precision, not an
algorithm polynomial in log n, and not a full-count speed theorem.

For every odd n>=3,

    U_n == -n (mod4).

For every odd n>=5,

    U_n == n^(-1)+2(n-2)^(-1)+4 C_n (mod8),

where C_n counts unordered distinct positive odd pairwise-coprime triples
summing to n. Inverses are taken modulo8. The n=1 mod4 and n=3 mod8 boundaries
are essential; the all-one triple at n=3 has a different centralizer.

The proof classifies the surviving one- and three-cycle strata and discards
all higher strata by a proved valuation inequality. This is a uniform proof
extracted from efficient finite-precision evaluation, not inference from a
finite numerical pattern. No Hamiltonian-path spectrum follows from it.

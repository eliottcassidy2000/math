# PROP-001: Arc-Flip Identity (Equivalent to OCF/Claim A)

**Type:** Theorem (PROVED — follows from OCF)
**Certainty:** 5 — PROVED for all n
**Status:** PROVED (corollary of OCF / Grinberg-Stanley)
**Last reviewed:** kind-pasteur-2026-03-05-S12
**Tags:** #arc-flip #ocf #claim-a #proof-strategy

---

## Statement

For any tournament T on n vertices, any arc i->j in T, let T' be the tournament
obtained by flipping arc i->j to j->i. Then:

```
H(T') - H(T) = I(Omega(T'), 2) - I(Omega(T), 2)
```

Equivalently, E(T) := H(T) - I(Omega(T), 2) is invariant under arc flips.

---

## Why This Implies OCF and Claim A

1. E(transitive tournament) = 1 - 1 = 0.
2. Any tournament is reachable from the transitive tournament by arc flips.
3. If E is invariant under arc flips, then E(T) = 0 for all T.
4. E(T) = 0 means H(T) = I(Omega(T), 2), which is the OCF formula.
5. OCF + Claim B (proved) => Claim A. (See Corollary in paper.)

---

## Algebraic Decomposition of delta_I

When flipping arc i->j to j->i:

**Lost cycles:** odd cycles C of T using arc i->j. All contain {i,j}, forming a clique in Omega(T).
**Gained cycles:** odd cycles C' of T' using arc j->i. All contain {i,j}, forming a clique in Omega(T').
**Common cycles:** odd cycles existing in both T and T' (not using i->j or j->i).

By the A-clique argument (same technique as the Claim B proof):

```
delta_I = 2 * [sum_{gained C'} I(R_{C'}, 2) - sum_{lost C} I(R_C, 2)]
```

where R_C = conflict graph of common cycles vertex-disjoint from C.

**Key observation:** Since V(C) contains {i,j}, the complement T[V\V(C)] is the same in T and T'. So R_C = Omega(T[V\V(C)]).

By strong induction (assuming OCF for tournaments on < n vertices):

```
I(R_C, 2) = H(T[V\V(C)])
```

Therefore:

```
delta_I = 2 * [sum_{C' gained} H(T[V\V(C')]) - sum_{C lost} H(T[V\V(C)])]
```

---

## The Target Identity

To prove PROP-001, it suffices to prove:

```
H(T') - H(T) = 2 * [sum_{C': odd cycle using j->i in T'} H(T[V\V(C')])
                     - sum_{C: odd cycle using i->j in T} H(T[V\V(C)])]
```

combined with strong induction on n. This is a purely combinatorial identity
about Hamiltonian path counts and subtournament path counts.

The LHS can be expressed as:
```
delta_H = #{paths using j->i in T'} - #{paths using i->j in T}
        = sum_{path-pairs (P1,P2) of V\{i,j}} [T[a][j]*T[i][b] - T[a][i]*T[j][b]]
```
where a = end(P1), b = start(P2).

---

## Verification Record

| n | Tests | delta_H = delta_I | Method |
|---|-------|-------------------|--------|
| 3 | 24 | 100% | exhaustive |
| 4 | 384 | 100% | exhaustive |
| 5 | 10240 | 100% | exhaustive |
| 6 | 600+ | 100% | random |

Also verified: the algebraic delta_I formula matches direct computation (100%).

R_C values observed: {1} at n<=5, {1, 3} at n=6.

---

## Relationship to Other Results

- Equivalent to Claim A (CONJ-001) and to OCF (THM-002).
- Uses the same A-clique technique as the Claim B proof (THM-003).
- Supersedes OPEN-Q-009 (arc-reversal invariance): PROP-001 is the E(T) version
  of arc-reversal invariance, which is cleaner because it avoids the vertex v.
- The "simple formula" delta = 2*(#gained - #lost) holds at n<=5 but fails at n>=6
  due to non-trivial I(R_C, 2) values.

---

## Resolution

PROP-001 is now PROVED as a corollary of OCF (THM-002). Since H(T) = I(Omega(T), 2)
for all tournaments (Grinberg-Stanley, arXiv:2412.10572 Corollary 20), E(T) = 0
identically, so it is trivially invariant under arc flips.

The direct combinatorial approaches listed below are no longer needed, but remain
of independent interest:
1. Transfer matrix / path-pair analysis
2. Generating function approach
3. Involution argument

---

## Source

kind-pasteur-2026-03-05-S5 (original); kind-pasteur-2026-03-05-S12 (resolution)

---
id: THM-5120
title: "Lawful rank-one reset transport and prime-word avoidance in the physical source grammar"
status: "PROVED, independently audited; shortest-reset list FINITE-EXACT. Full reset ideal and Chan–Pak coverage OPEN."
session: Session246, 2026-11-10
source: 05-knowledge/results/session246_unifying_principles_20261110/bridges/RESET_PRINCIPLE.md
audit: 05-knowledge/results/session246_unifying_principles_20261110/audit/AUDIT.md
---

# A reset preserves every continuation, with a multiplicative charge

Use the inherited physical source grammar, seed t=(1,0,2,1), and chronological
letters O,Z,J on x=(D,E,U,d):

    O(x)=(D,0,2U,d),
    Z(x)=(D−E+U,d,2U−d,d),
    J(x)=(D−E+U,D−E,D+U,D).

Every word constructs an actual rooted tournament; D is its Hamilton-path
count. The initial graph is an arc and each letter adds a vertex.

**Exact reset theorem.** The two words OJ and JZ have the same full matrix.
For m>=0, k>=1, either choice in

    r=(OJ or JZ) Z^m O^k J

has the rank-one action

    M_r x = (D+2^(m+1)U) v_k,
    v_k=(1+2^k,1,1+2^k,1)=M_(O^(k−1)J)t.

Therefore every physical prefix u and arbitrary suffix w satisfy

    H(u r w)=(D(u)+2^(m+1)U(u)) H(O^(k−1)J w).

The proof keeps a physical representative of v_k, the exact full state,
and the same lawful future transfers. Associativity then proves the identity
for every finite suffix, including every finite prefix of an infinite word.
No infinite Hamilton count or limiting count is defined by this statement.

Both factors exceed one: reachable states have D>=1,U>=2, and D never
decreases under the letters. Hence a prime-count word cannot contain this
reset sublanguage or a word obtained from one of its factors by OJ↔JZ
replacements. This is a necessary constraint within this grammar, not a
construction of primes or an exclusion from the tournament spectrum.

The four shortest rank-one words are OJOJ,OJJZ,JZOJ,JZJZ, all with matrix
(3,1,3,1)^T(1,0,2,0). Exhaustive exact ranks of the 120 nonempty words
through length four prove that no shorter reset exists. This does not
classify the full two-sided rank-one ideal.

The six-vertex word OJOJ has outgoing row masks [54,60,56,17,32,8], is
strongly connected, and has H=15. Thus the factorization need not arise from
an ordinal sum of strongly connected components. It is a continuation-state
reset with an unavoidable retained multiplier.

The elementary common mechanism, **lawful reset transport**, also describes
the precisely indexed modular reset of THM-3182
(`THM-3182-factorial-gauss-manin-rank-one-reset-and-two-transverse-smith-bands.md`)
and the spatial boundary reset of THM-4204
(`THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md`). Their rings,
legal continuations and observers differ and remain explicit in the full
proof. Neither source problem's open claims are transported to tournaments.

[Proof, sampled provenance, literal graph checks and exact verifier](../../05-knowledge/results/session246_unifying_principles_20261110/bridges/RESET_PRINCIPLE.md).

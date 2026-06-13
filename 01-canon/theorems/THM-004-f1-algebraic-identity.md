# THM-004: F1 — Algebraic Identity for inshat

**Type:** Lemma
**Certainty:** 5 — PROVED
**Status:** PROVED
**Last reviewed:** SYSTEM-2026-03-05-S1
**Disputes:** none
**Tags:** #inshat #type-ii #algebraic-identity #signature

---

## Statement

For any path P' of T−v with binary signature s (where s[j] = 1 if v→P'[j], else 0):

```
(inshat(v, P') − 1) / 2  =  #{Type-II positions in P'}
```

where a Type-II position at index j means s[j]=1, s[j+1]=0.

This is a **pure algebraic identity** for any binary signature. It does not depend on the tournament structure beyond the definition of inshat.

---

## Proof / Proof Sketch

Let:
- s = signature (binary string of length n−1)
- b = boundary term = s[0] + (1 − s[n−2])
- #01 = #{j : s[j]=0, s[j+1]=1} (Type-I count)
- #10 = #{j : s[j]=1, s[j+1]=0} (Type-II count)

Then: inshat = b + #01 + #10

We want to show: (inshat − 1)/2 = #10, i.e., b + #01 − 1 = #10.

**Telescoping identity:** #10 − #01 = s[0] − s[n−2] (the number of 1→0 transitions minus 0→1 transitions equals the change from first to last).

So: #10 = #01 + s[0] − s[n−2].

**Substituting:** b + #01 − 1 = s[0] + (1 − s[n−2]) + #01 − 1 = s[0] − s[n−2] + #01 = #10. □

## Notes & History

Discovered during computational analysis (scripts 1-5). The proof is elementary but the identity is foundational — it converts the counting question about inshat into a question about Type-II positions, which then connects to 3-cycles via THM-005.

**Consequence:** The per-path identity (inshat−1)/2 = Σ_C μ(C) over 3-cycle embeddings is equivalent to: #Type-II = Σ_{3-cycles C=(v,a,b), (a,b) consecutive in P'} μ(C).

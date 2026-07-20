---
id: THM-1560
title: "THE HALVING DICTIONARY: {−1,0,1} AND {1,1/2,0} ARE ONE OBJECT IN TWO COORDINATES, AND THE 1/2 IS EXACTLY WHAT DIES MOD 2. (A) The map x ↦ (1+x)/2 carries the sign encoding to the fraction encoding, and for tournaments it IS the standard A = (J−I+S)/2 — the repo's GF(2) machinery (tilings, switching classes, cut⊕cycle, even graphs, two-graphs) and its sign machinery (skew adjacency, Pfaffian, charge, the nullcone) are the same data one halving apart, invertible over ℤ[1/2] and undefined in characteristic 2. (B) THE COLLAPSE PRINCIPLE, a theorem: S ≡ J−I (mod 2) for EVERY tournament, since +1 ≡ −1, so every invariant that is an INTEGER polynomial in S is constant mod 2. Pfaffian-oddness (THM-1475) is the one-line corollary F = Pf, giving (n−1)!!. (C) THE CENSUS confirms the split exactly: Pf(S), det(S), tr(S²), tr(S⁴) are constant mod 2; c₃(T) and tr(A³), whose expression in S carries a 1/2, both vary. (D) THE RÉDEI CLASS IS SMALL: hp needs the 1/2 yet is mod-2 blind — but Hamiltonian cycles, out-arborescences, transitive triples and 4-cycles all vary. So Rédei's oddness is a genuinely special fact, not an instance of anything structural. (E) THE ASYMMETRY IS THE CONTENT: mod-2 blindness is automatic in the sign world and therefore carries no information, and rare in the fraction world and therefore carries a lot — which is precisely why Rédei is a theorem while Pfaffian-oddness is a triviality. (F) On the nullcone side the same triple recurs: Gamma shapes a = ν/2 ∈ {1, 1/2, 0} for ν ∈ {2,1,0} real Gaussians, i.e. n ∈ {4,3,2} and 1/a ∈ {1,2,∞} — degenerating at a = 0, which is exactly where GMC(2) lives."
status: >
  (A) PROVED (an identity) + VERIFIED entrywise at n = 4,5,6.
  (B) PROVED.  One line, and it is the correct general statement behind THM-1475 §1.
  (C) VERIFIED-EXHAUSTIVE at n = 4,5,6 (all tournaments, or 4096 sampled at n=6).
  (D) VERIFIED on five counts at n = 4,5,6.  A bounded census, not a classification: it shows
      the four tested companions of hp are NOT mod-2 blind, it does not prove hp is unique.
  (E) READING of (B)-(D), not a separate theorem.
  (F) STRUCTURAL OBSERVATION, deliberately under-claimed.  I do NOT assert a canonical
      bijection between charge and Gamma shape; what is checkable is that both ladders
      degenerate at the same end.
  Advances no open problem.  It is a dictionary plus one small theorem plus a census.
source: klein-2026-07-20-S349 (owner: still think nullcone; {−1,0,1} and {1,1/2,0} are functionally equivalent but have each shown up many times — this is what I mean by even/odd vs positive/negative)
depends_on:
  - THM-1475  # Pf(S(T)) is odd -- now a corollary of (B)
related:
  - THM-1430  # Seidel switching of graphs lands on E_n (the GF(2) side)
  - THM-1470  # even tournaments / two-graphs
  - THM-1500  # the GMC master theorem and the Gamma ladder
  - THM-1550  # the toral nullcone criterion
script: 04-computation/halving_dictionary_klein_S349.py (+ .out)
---

# THM-1560 — the halving dictionary

## A. The two coordinates

```text
x ↦ (1+x)/2      sends   {−1, 0, 1}  ⟶  {0, 1/2, 1}
```

For tournaments this is exactly the standard conversion between the **skew** adjacency `S`
(`S_ij = ±1`, `S_ii = 0`) and the **0/1** adjacency `A`:

```text
A = (J − I + S)/2,        S = 2A − (J − I).
```

Verified entrywise at `n = 4,5,6`: `A_ij = (1 + S_ij)/2` off the diagonal, and the diagonal
`S_ii = 0` is precisely the `1/2`. So the repo's two halves —

| fraction world `{0, 1/2, 1}` | sign world `{−1, 0, 1}` |
|---|---|
| tilings, switching classes, cut ⊕ cycle | skew adjacency `S` |
| even graphs, two-graphs (THM-1430/1470) | Pfaffian, `det S` (THM-1475) |
| GF(2) parity, even/odd | charge/weight grading, the nullcone |

— are **one object in two coordinate systems**. The dictionary costs a `1/2`, and `1/2` does
not exist in characteristic 2. *That* is the even/odd versus positive/negative split.

## B. The collapse principle

Every off-diagonal entry of `S` is `±1`, and `+1 ≡ −1 (mod 2)`. Hence

> **`S ≡ J − I (mod 2)` for every tournament on `n` vertices** — mod 2 the skew adjacency
> forgets the orientation completely.

**Theorem.** If `I(T) = F(S)` for a polynomial `F` with integer coefficients, then
`I(T) ≡ F(J−I) (mod 2)` is the *same value for every tournament*.

**Corollary (THM-1475 §1).** `Pf(S) ≡ Pf(J−I) = #{perfect matchings of K_n} = (n−1)!!`, odd.
What I proved there as a fact about Pfaffians is the general principle specialised to `F = Pf`.

## C. The census — the split is exactly where the theorem says

| invariant | integer polynomial in `S`? | values mod 2 | constant? |
|---|---|---|---|
| `Pf(S)` | yes | `{1}` | **yes** |
| `det(S)` | yes | `{0}` (odd `n`), `{1}` (even `n`) | **yes** |
| `tr(S²)`, `tr(S⁴)` | yes | `{0}` | **yes** |
| `c₃(T)` | no — needs `1/2` | `{0,1}` | no |
| `tr(A³)` | no — needs `1/2` | `{0,1}` | no |

Exhaustive at `n = 4,5` and over 4096 tournaments at `n = 6`. The invariants that can *see*
the tournament mod 2 are exactly those whose definition passes through the halving.

## D. The Rédei class is small

`hp(T)` needs the `1/2` (it is a count in the `0/1` world) and yet is mod-2 blind — that is
Rédei's theorem. Is that an instance of something? **No.**

| count | values mod 2, `n = 4,5,6` | mod-2 blind? |
|---|---|---|
| `hp` (Hamiltonian paths) | `{1}` | **yes** — Rédei |
| Hamiltonian cycles | `{0,1}` | no |
| out-arborescences at a fixed root | `{0,1}` | no |
| transitive triples | `{0,1}` | no |
| 4-cycles | `{0,1}` | no |

Being a fraction-world count is **not** enough. Rédei's oddness is special.
*(Bounded census — it shows these four companions fail, not that `hp` is unique.)*

## E. The asymmetry is the content

- **Sign world:** mod-2 blindness is *automatic* — every integer polynomial in `S` has it, by
  §B. It therefore carries **no information**, and Pfaffian-oddness is a one-liner.
- **Fraction world:** mod-2 blindness is *rare* — of five natural counts, only `hp`. It
  therefore carries **a great deal**, and Rédei's theorem is a theorem.

This is the useful form of the owner's even/odd versus positive/negative distinction: the two
encodings are equivalent, but *the same property means opposite things in them*. A parity
statement proved on the sign side is usually free; the same statement on the fraction side is
usually not.

## F. The same triple on the nullcone side

The GMC witnesses (THM-1490/1500/1510) use `P = (1+Z)(W − g(Z)U)` with `U ~ Gamma(a)`,
`a = ν/2` for `ν` real Gaussians, and the master equation forces `1 + s·g(s) = (1+s)^{1/a}`,
so `g` is a polynomial iff `1/a ∈ ℤ₊`:

| `a = ν/2` | `1/a` | `G(x) = (1+x)^{−a}` | `n = 2+ν` |
|---|---|---|---|
| `1` | 1 | `1/(1+x)`, rational | 4 |
| `1/2` | 2 | `(1+x)^{−1/2}`, square root | 3 |
| `0` | ∞ | `e^{−cx}`, **degenerate** | 2 |

and the witnesses' charges are `{−1, 0, +1}`. So `{1, 1/2, 0}` counts *half-Gaussians* and
`{−1, 0, 1}` counts *charge*, both indexing `ν ∈ {2,1,0}` under different normalisations.

**Deliberately under-claimed:** I am *not* asserting a canonical bijection between charge and
Gamma shape. What is checkable, and what matters, is that **both ladders degenerate at the
same end** — `a = 0`, `n = 2`, `1/a = ∞` — and that end is exactly where GMC(2) lives and
where THM-1550's `M = 1` argument stops applying.

*Files: `04-computation/halving_dictionary_klein_S349.py` (+ `.out`).*

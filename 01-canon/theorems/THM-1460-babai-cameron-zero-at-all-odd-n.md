---
id: THM-1460
title: "BABAI–CAMERON REMARK 7.4 IS ZERO AT EVERY ODD n, NOT JUST n = 1 MOD 4 — and the proof is that the score-PARITY VECTOR is a COMPLETE INVARIANT within a switching class of tournaments when n is odd. Switching at W changes s_v by |W^c|-2a (v in W) or |W|-2b (v not in W), so for n ODD exactly one of |W|,|W^c| is odd; since W ~ W^c we may take |W| even, and then the parity vector flips exactly on W. W ranges over the 2^(n-1) EVEN subsets while the class has exactly 2^(n-1) members, so the map is a BIJECTION: each reachable parity vector occurs EXACTLY ONCE. Since #odd scores = C(n,2) mod 2, the reachable set is the EVEN-weight coset when n = 1 mod 4 (containing 0, the ALL-EVEN anchor) and the ODD-weight coset when n = 3 mod 4 (containing all-ones, the ALL-ODD anchor, since n odd). BOTH anchors are permutation-invariant, so any automorphism of the switching class fixes that member — hence the Babai–Cameron count is ZERO for EVERY odd n. The owner's n = 1 mod 4 case is the all-even half; n = 3 mod 4 works identically with the all-ODD member. At EVEN n the argument COLLAPSES: |W| and |W^c| share parity, so either no vertex flips or all do, leaving only TWO parity vectors — verified, all-even counts of 0, 4, 64 rather than 1. THIS IS THE 3-vs-1 MOD 4 SPLIT THE OWNER NOTICED: 3, 7, 11 carry an all-ODD anchor and 5, 9, 13 an all-EVEN one — the same split that gives Paley TOURNAMENTS at p = 3 mod 4 and Paley GRAPHS at p = 1 mod 4"
status: PROVED (the bijection and both anchors) + VERIFIED (n=5 exhaustive 64/64; n=3,7,9,11 sampled; even-n control fails as predicted)
author: opus-2026-07-20-S409
depends_on: [THM-474 (mac-mini; cites Babai–Cameron), THM-1445 (skew/symmetric = odd/even), THM-1430, THM-1455]
credits: owner (the n = 1 mod 4 observation and the Remark 7.4 pointer)
---

# THM-1460 — Babai–Cameron Remark 7.4 is zero at every odd `n`

## 1. The question

Babai–Cameron, *EJC* **7** (2000) #R38, Remark 7.4 asks one to enumerate the switching
classes of tournaments in which **some automorphism fixes no member**, and states:
*"We cannot do this."*

The owner observed that at `n ≡ 1 (mod 4)` the answer is **zero**, because every
automorphism must fix the unique **even** member. That is correct, and it **extends to every
odd `n`** — the `n ≡ 3 (mod 4)` case has its own anchor.

## 2. The parity vector is a complete invariant (odd `n`)

Switching at `W` reverses all arcs across the cut, so for `v ∈ W` the score changes by
`|W^c| − 2a_v` and for `v ∉ W` by `|W| − 2b_v`. Modulo 2:

> the score-parity vector flips **on** `W` iff `|W^c|` is odd, and **off** `W` iff `|W|` is odd.

**Let `n` be odd.** Then exactly one of `|W|, |W^c|` is odd. Since `σ_W = σ_{W^c}`, we may
always take `|W|` **even**; then `|W^c| = n − |W|` is odd, and the parity vector flips
**exactly on `W`**. So switching acts on parity vectors by adding `1_W`.

`W` ranges over the `2^{n−1}` even-cardinality subsets, and the switching class has exactly
`2^{n−1}` members. **The map is a bijection: each reachable parity vector is realised by
exactly one member of the class.** The score-parity vector is a *complete invariant* inside
a switching class.

## 3. Which coset, and the two anchors

`Σ_v s_v = C(n,2)`, so `#{odd scores} ≡ C(n,2) (mod 2)`. Hence the reachable parity vectors
form a coset of the even-weight subspace:

| `n mod 4` | `C(n,2)` | reachable coset | permutation-invariant member | anchor |
|---|---|---|---|---|
| `1` | **even** | even-weight | `0` (weight 0, even ✓) | **unique ALL-EVEN tournament** |
| `3` | **odd** | odd-weight | `1…1` (weight `n`, odd ✓) | **unique ALL-ODD tournament** |

The only permutation-invariant vectors in `𝔽₂ⁿ` are `0` and `1…1`; for odd `n` exactly one
of them lies in the reachable coset, and §2 says it is hit exactly once.

**Theorem.** *For every odd `n`, each switching class of tournaments on `[n]` contains a
unique member `A` whose score vector is all-even (`n ≡ 1 mod 4`) or all-odd (`n ≡ 3 mod 4`).
Any automorphism `g` of the switching class permutes its members and preserves the
score-parity type (relabelling does not change the multiset of scores), so `g(A)` is such a
member, whence `g(A) = A`. Therefore `g` fixes a member, and the Babai–Cameron count is
**zero** for every odd `n`.* ∎

## 4. Verification

| `n` | `n mod 4` | all-EVEN members per class | all-ODD members per class |
|---|---|---|---|
| 3 | 3 | `{0: 200}` | **`{1: 200}`** |
| 5 | 1 | **`{1: 200}`** (and `64/64` exhaustively) | `{0: 200}` |
| 7 | 3 | `{0: 200}` | **`{1: 200}`** |
| 9 | 1 | **`{1: 40}`** | `{0: 40}` |
| 11 | 3 | `{0: 40}` | **`{1: 40}`** |

`n = 5` exhaustive: 1024 tournaments, 64 switching classes, **every one** with exactly one
all-even member.

## 5. Even `n`: the argument collapses (control)

For `n` even, `|W|` and `|W^c|` have the **same** parity, so switching either flips **no**
score parity or **all** of them. Only **two** parity vectors occur per class instead of
`2^{n−1}` — the invariant is far too coarse and there is no canonical anchor. Verified:

```
n = 4 : all-even members per class  {0: 51, 4: 9}
n = 6 : {0: 60}
n = 8 : {0: 57, 64: 3}
```

Never uniquely 1. **So oddness of `n` is exactly the hypothesis**, and Remark 7.4 remains
genuinely open at even `n` — which is presumably where Babai–Cameron's difficulty lives.

## 6. The 3-vs-1 mod 4 split the owner noticed

The residue does **not** decide whether an anchor exists — it decides **which** anchor:

- **`p ≡ 3 (mod 4)`: 3, 7, 11** — anchor is the **all-ODD** tournament. `C(p,2)` is odd, so
  an all-even tournament is *arithmetically impossible*. This is the **skew/odd** world, and
  it is exactly where **Paley tournaments** live.
- **`p ≡ 1 (mod 4)`: 5, 13** (and `9`) — anchor is the **all-EVEN** tournament. This is
  where **Paley graphs** live.

So 3 and 7 really do behave alike, and 5 groups with 9 and 13, exactly as the owner said —
and the mechanism is the parity of `C(n,2)`. This is the **same** mod-4 dichotomy that
THM-1445 identified as skew-vs-symmetric: `−1` is a non-residue iff `p ≡ 3 (mod 4)`, which
is what makes the complement an anti-automorphism and puts tournaments in the odd
eigenspace. **One dichotomy, three faces: `C(n,2)` parity, Paley tournament vs Paley graph,
and skew vs symmetric.**

## 7. On the owner's `k = 2` ↔ JC₂ suggestion — recorded, NOT claimed

There is a real rhyme: **THM-1350** proves a Keller counterexample collision can never be a
**double** (minimum 3), and **THM-1455** finds `k = 2` unrealisable. Both are "2 is
excluded". But THM-1350's exclusion is a *proved parity argument* (odd fibre ⟹ a fixed
point), whereas the `k = 2` gap is **still unexplained**. I am not asserting they are the
same mechanism.

**A quick argument I tried and had to discard**, recorded so nobody repeats it: I reasoned
that flipping one arc toggles both triples of a 4-subset that contain it, so `ε(T)` would be
flip-invariant and `k` constant. **False** — flipping an arc of a *transitive* triangle does
not always toggle it (`a→b, a→c, b→c` with `(a,b)` flipped is still transitive). The
premise "flipping an arc of a triangle switches cyclic ↔ transitive" holds for cyclic
triangles but not for all transitive ones. `k = 2` remains open.

## 8. Open

1. **`k = 2`** (THM-1455) — still unexplained; §7 kills the easy route.
2. **Even `n` for Remark 7.4** — the genuinely open half; §5 shows why parity cannot reach it.
3. **Does the anchor give a canonical form?** The all-odd/all-even member is a *canonical
   representative* of each switching class at odd `n`. That is a free normal form for the
   tournament censuses (THM-1415/1455) and may cheapen them considerably.

## Verification

`04-computation/babai_cameron_mod4_opus_S409.py` (mod-4 table; `n=5` exhaustive; `n=9`
sampled; `n=3,7` impossibility), `04-computation/parity_anchor_all_odd_n_opus_S409.py` (both
anchors at `n = 3,5,7,9,11`; even-`n` control). Outputs in `05-knowledge/results/`.

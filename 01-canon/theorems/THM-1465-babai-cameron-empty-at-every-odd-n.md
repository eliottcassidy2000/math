---
id: THM-1465
title: "BABAI–CAMERON REMARK 7.4 IS EMPTY AT EVERY ODD n — not just at n ≡ 1 mod 4, and the odd-cycle indicator is the missing half. (0) Babai–Cameron (EJC 7 (2000) #R38) ask to enumerate the switching classes of tournaments admitting an automorphism fixing NO member, and write 'We cannot do this'. klein-S338 settled n ≡ 1 mod 4 via the unique all-even-score member. (I) THE THEOREM: for EVERY odd n, every automorphism of a switching class fixes a member, so the failure set is empty. Proof: klein's parity map is a σ-EQUIVARIANT bijection from the 2^{n−1} members onto a coset K of the even-weight code, so σ fixes a member iff σ fixes a vector of K; σ-fixed vectors are the unions of cycles; at n ≡ 1 mod 4 take the empty union (0 ∈ K), and at n ≡ 3 mod 4 take any ODD-LENGTH cycle — one exists because n is odd — whose indicator has odd weight and so lies in K. (II) COROLLARY, with Moon's odd-order theorem: at odd n every permutation stabilising a switching class has ODD order. Verified: the stabilising orders are exactly {3}, {3,5}, {3,5,7} at n = 3,5,7, while even orders appear at n = 4 ({2,3,4}) and n = 6 ({2,3,5,6}). (III) SO THE QUESTION IS REALLY ABOUT EVEN n, and the first terms are computed here: the failure set has 2 classes at n = 4 and 4 at n = 6 (up to isomorphism), from 8 and 640 labelled. (IV) I ALSO ACCEPT klein's correction of my THM-1415: I compared switching classes against even GRAPHS (A002854) when the right comparand is even TOURNAMENTS, all scores even, counted by A049313 — my own 1,2,2,6 was A049313 all along"
status: >
  (I) PROVED.  Step 1 is klein-S338's parity law (THM-1440), used as an input and
  credited; steps 2-6 are three lines and are mine.  Independently VERIFIED by exact
  computation: the failure set is empty at n = 3, 5, 7 (all 32768 cosets × 5040
  permutations at n = 7).
  (II) PROVED as a corollary (a stabilising σ fixes a member, hence lies in some Aut(T),
  which has odd order by Moon).  VERIFIED exactly at n = 3,5,7 and shown to fail at
  n = 4,6, which is the control that the statement has content.
  (III) VERIFIED-EXACT by the same reduction.
  HONEST: the argument is elementary and may be folklore — klein flagged the same about
  their half.  Babai–Cameron's "We cannot do this" refers to the general enumeration, and
  what is settled here is only the odd-n half of it.  The even-n enumeration, which is
  where the question actually lives, remains open beyond the two terms below.
source: kind-pasteur-2026-07-20-S128c116 (owner: Babai–Cameron Remark 7.4; at n ≡ 1 mod 4 the answer is zero; 3 and 7 are alike mod 4 while 5 goes with 1 and 9)
depends_on:
  - THM-1440    # klein-S338: the score-parity law and the n = 1 mod 4 case
related: [THM-1415]
external:
  - "Babai–Cameron, Electron. J. Combin. 7 (2000) #R38, Remark 7.4: 'We cannot do this.'"
  - "Moon: the automorphism group of a tournament has odd order."
script: 04-computation/babai_cameron_failure_set_kps_S128c116.py, babai_cameron_n7_kps_S128c116.py (+ .out)
---

# THM-1465 — the odd-n half of Remark 7.4, closed

## 0. The question, and what was already done

Babai–Cameron ask for the number of switching classes of tournaments admitting an
automorphism that fixes **no** member of the class, and say plainly: *"We cannot do
this."* klein-S338 settled `n ≡ 1 mod 4`: such a class contains a unique tournament with
all scores even, any automorphism must fix it, so the failure set is empty there.

That argument is specific to `n ≡ 1 mod 4` — it is where `Σ s_v = C(n,2)` is even. The
owner's observation is the right frame: **`3` and `7` sit together mod 4 and `5` goes
with `1` and `9`**, and klein's theorem covers only the second family.

## I. The theorem

> **For every odd `n`, every automorphism of a switching class of tournaments fixes a
> member. The Babai–Cameron failure set is empty at every odd `n`.**

*Proof.* Let `C` be a switching class and `σ` an automorphism of it.

1. *(klein-S338)* For odd `n` the score-parity map `T ↦ p(T) ∈ 𝔽₂ⁿ`, `p(T)_v = s_v mod 2`,
   is a **bijection** from the `2^{n−1}` members of `C` onto a coset `K` of the
   even-weight code: `K =` the even-weight vectors when `n ≡ 1 (mod 4)`, the odd-weight
   vectors when `n ≡ 3 (mod 4)`.
2. It is **σ-equivariant**: `p(σT) = σ·p(T)`, since relabelling permutes scores.
3. Hence `σ` fixes a member of `C` **iff** `σ` fixes some vector of `K` acting by
   coordinate permutation.
4. The `σ`-fixed vectors of `𝔽₂ⁿ` are exactly the indicators of **unions of cycles** of
   `σ`; the weight of such a vector is the sum of the chosen cycle lengths.
5. `n ≡ 1 (mod 4)`: take the empty union. `0 ∈ K`, and `σ·0 = 0`. ∎
6. `n ≡ 3 (mod 4)`: `n` is odd, so `σ` has **at least one cycle of odd length**. Its
   indicator is `σ`-fixed and has odd weight, hence lies in `K`. ∎

Step 6 is the missing half, and it is where the two residue classes rejoin: at
`1 mod 4` the witness is the *empty* union of cycles, at `3 mod 4` it is a *single
odd cycle*. Both exist precisely because `n` is odd.

**Verified independently**, by a reduction to linear algebra over `𝔽₂`: encoding a
tournament by the edge set it reverses relative to `i → j`, switching classes are the
cosets of the cut space `C = span{δ(U)}`, `σ` acts affinely as `x ↦ Px + c`, and with
`A := P + I`, `v := Ax + c`,

> `σ` stabilises `x + C` ⟺ `v ∈ C`;  `σ` fixes a member ⟺ `v ∈ A(C)`

(and `A(C) ⊆ C`, since `A δ(U) = δ(σU) + δ(U)`), so failure is exactly `v ∈ C ∖ A(C)`.
Exhaustively: **0 failures at `n = 3, 5, 7`** — the last over all `32768` cosets and all
`5040` permutations.

## II. Corollary, via Moon

Moon: the automorphism group of a tournament has **odd order**. A permutation stabilising
a switching class fixes a member by (I), hence lies in that member's automorphism group.
So:

> **At odd `n`, every permutation that stabilises a switching class has odd order.**

Verified, and the control shows it has content:

| n | orders of stabilising permutations | all odd? |
|---|---|---|
| 3 | {3} | yes |
| 5 | {3, 5} | yes |
| 7 | {3, 5, 7} | yes |
| **4** | **{2, 3, 4}** | **no** |
| **6** | **{2, 3, 5, 6}** | **no** |

Read the other way: an even-order permutation lies in no `Aut(T)`, so it fixes no member
*automatically* — and (I) says no such permutation can stabilise a class at all when `n`
is odd.

## III. The question is really about even n

| n | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| switching classes (A049313) | 1 | 2 | 2 | 6 | 12 |
| **failure set, up to iso** | **0** | **2** | **0** | **4** | **0** |
| failure set, labelled | 0 | 8 | 0 | 640 | 0 |

So Remark 7.4 is a question about even `n`, and `2, 4` are its first two terms. At even
`n` klein's parity vector is pinned only up to global complement, there is no canonical
member, and even-order permutations do stabilise classes — all three obstructions to the
argument above appear together.

## IV. Accepting a correction

klein-S338 correctly corrects my THM-1415, which concluded "the graph two-graph theorem
has no tournament analogue" from `1,2,2,6 ≠ A002854`. The comparison was against the
wrong object: the analogue of an even **graph** is an even **tournament** (all scores
even), not an even graph, and the right sequence is **A049313** — which is exactly what
my own `1, 2, 2, 6` was. The counts were right and the conclusion drawn from them was
wrong. THM-1415 §II should be read as superseded.

## Named next

- The even-`n` enumeration. `n = 8` needs `2^{21}` cosets × `40320` permutations, which
  the same vectorised reduction should reach; that would give a third term and make an
  OEIS search worthwhile.
- Whether the even-`n` failure classes have a structural description — at `n = 4` there
  are only two, and they should be identifiable by hand.

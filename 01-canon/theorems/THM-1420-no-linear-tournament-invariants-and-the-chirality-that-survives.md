---
id: THM-1420
title: "THERE ARE NO F_2-LINEAR TOURNAMENT INVARIANTS AT ALL, AND THE REASON IS THAT 'FORWARD' IS NOT CANONICAL — so the search for a descending star-type invariant is provably empty, for EVERY subgroup, base-path-free or not. (D) The isomorphism-difference group Gamma_min (the smallest subgroup whose orbits are unions of iso classes) is ALL of F_2^E: for an adjacent transposition the only pair whose endpoint order reverses is (k,k+1) itself, so the affine shift is c(tau_k) = e_{(k,k+1)}, and S_n is edge-transitive, so Gamma_min contains every basis vector. Verified dim Gamma_min = C(n,2) exactly at n = 3..7. Hence NO proper subgroup descends — not the star group, not switching (THM-1415), not \\cap Gamma_P — and THM-1405's transversality is FORCED, not accidental. (D') The mechanism is exactly the affine twist: with the linear part alone the quotient is 1-dimensional (total parity), and c kills precisely that one dimension; over A_n, where c always has even weight, the dimension survives — quotient dim exactly 1 at n = 3..7. (D'') CHIRALITY: |Aut(T)| is odd so Aut(T) <= A_n, hence EVERY iso class splits into EXACTLY TWO A_n-classes, separated by inversion parity — verified 2,4,12,56 -> 4,8,24,112 at n = 3..6. (A) The tiling fibre over a node is t = H/|Aut|, and since H is odd (Redei) and |Aut| is odd, EVERY fibre is odd — so #iso classes = 2^m mod 2 is EVEN for n >= 3, a derivation rather than an observation. (B) The metagraph edge set computes from Aut-orbits on arcs, never from 2^E: measured speedups 24x, 116x, 698x, 4941x at n = 4..7."
status: >
  (D) PROVED (one line, all n >= 2) and VERIFIED-EXACT (dim Gamma_min = C(n,2) at n = 3..7
  by direct F_2 span computation over all n! relabelings).
  (D') PROVED and verified: A_n quotient dimension exactly 1 at n = 3..7.
  (D'') PROVED (|Aut| odd => every element has odd order => Aut <= A_n) and VERIFIED-EXACT
  by full enumeration at n = 3,4,5,6: A_n-class counts are exactly twice the S_n-class
  counts and inversion parity separates the two halves of every class.
  (A) PROVED from Redei + |Aut| odd; verified n = 3..7 (sum of fibres = 2^m exactly).
  (B) VERIFIED-EXACT n = 4..7, with the involution consistency identity
  (n!/|Aut(C)|)*N(C,C') = (n!/|Aut(C')|)*N(C',C) holding at every n.
source: mac-mini-2026-07-20-S129 (owner: figure out how tiling sets map to iso class nodes
  exactly; compute edges/nodes/tilings from each other efficiently, look for tricks; and
  "a descending star-type invariant has to come from a base-path-independent subgroup —
  the natural candidate being \\cap Gamma over all spanning paths")
depends_on:
  - THM-1405  # the star quotient is the cycle space, transverse to isomorphism
credits:
  - THM-1415  # kind-pasteur-S128c110, CONCURRENT and first to publish: the base-path-free
              # fix is to stop restricting, and the canonical star quotient is SWITCHING.
              # Part C below was derived independently and is CEDED to THM-1415; what is
              # new here is (D), which explains why switching's coarseness is not special.
related:
  - the GF(2) Cut (+) Cycle split (CLAUDE.md)
  - A003141 / min-FAS (the inversion parity here is NOT min-FAS — see Honest scope)
script: 04-computation/tiling_node_edge_dictionary_macmini_S129.py (+ .out)
---

# THM-1420 — no linear invariants, and the chirality that survives

**One line.** The reason no star-type invariant descends to iso classes is not that the base
path is in the way. It is that **`F₂`-linear tournament invariants do not exist at all** — and
they fail to exist for one reason: *which way is "forward" is not canonical.*

## (A) The fibration, and what oddness forces

The tiling fibre over an iso class is `t(C) = H(C)/|Aut(C)|` (orbit–stabilizer on the path
fibration; already canon as "tilings × |Aut| = H"). Two classical oddness facts now bite:

- **Rédei:** `H(T)` is odd.
- **Classical:** `|Aut(T)|` is odd (a tournament has no automorphism of even order).

Hence **every fibre `t(C)` is odd**, and `Σ_C t(C) = 2^{C(n−1,2)}`. A sum of `#classes` odd
numbers is `≡ #classes (mod 2)`, so:

> **`A000568(n)` is EVEN for every `n ≥ 3`** — derived, not observed.

Verified exactly (`Σ t = 2^m`, every `t` odd, every `|Aut|` odd) at `n = 3..7`:
classes `2, 4, 12, 56, 456`, fibre sums `2, 8, 64, 1024, 32768`.

In the **merged** metagraph the fibre over a node is `(2 − [C is SC])·t(C)`, since `H` and
`Aut` are both complement-invariant.

## (B) The dictionary — computing nodes, tilings and edges from each other

| want | from | how |
|---|---|---|
| tilings over a node | `H`, `Aut` | `t(C) = H(C)/|Aut(C)|` |
| merged fibre | `t`, SC-status | `(2 − [SC])·t` |
| **metagraph edges** | `Aut` | flip **one arc per `Aut(C)`-orbit**, weight by orbit size |
| global check | — | `Σ_C H(C)/|Aut(C)| = 2^{C(n−1,2)}` |

The edge trick is the useful one: the arc-flip metagraph never requires touching `2^{C(n,2)}`
labelled tournaments. Measured canonicalization counts vs brute force:

| `n` | Aut-orbit canonicalizations | brute force | speedup |
|---|---|---|---|
| 4 | 16 | 384 | 24× |
| 5 | 88 | 10 240 | 116× |
| 6 | 704 | 491 520 | 698× |
| 7 | 8 912 | 44 040 192 | **4941×** |

with the involution consistency identity `(n!/|Aut(C)|)·N(C,C') = (n!/|Aut(C')|)·N(C',C)`
verified at every `n` — a strong check that the orbit weighting is right.

## (C) The base-path-independent subgroup — ceded to THM-1415

The proposed candidate `⋂_P Γ_P` over all spanning paths is **trivially zero**: `Γ_P` is
supported off `E(P)`, and every edge of `K_n` lies on some Hamiltonian path (verified: all
`E` edges covered at `n = 4..7`), so the intersection is supported on `∅`.

The right move is not to intersect but to **stop restricting**, because

> `Γ_P = Cut(K_n − P)` is exactly the **restriction of `Cut(K_n)`** to `E ∖ E(P)`
> (incidence rows restrict to incidence rows; the restriction is injective for `n ≥ 4`,
> since a cut `δ(S)` inside `E(P)` would need `|S|·|S^c| ≤ 2`),

and `Cut(K_n)` is `S_n`-invariant and base-path-free — it is the **Seidel switching group**.
This identification is **THM-1415 (kind-pasteur-S128c110), which published first**; the
derivation here was independent and concurrent and is ceded. The lemma above — that `Γ_P` is
literally a *restriction* of the switching group — is the bridge showing THM-1405's
base-path-dependent group and THM-1415's canonical one are one object seen twice.

## (D) The theorem: no proper subgroup descends

Let `Γ_min` be the smallest subgroup with `P·x + x ∈ Γ` for all tournaments `x` and all
relabelings `P`. Orbits of `Γ` are unions of iso classes **iff** `Γ ⊇ Γ_min`.

The `S_n`-action on the arc space `F₂^E` is **affine**: `x ↦ P(x) + c(P)`, where `c(P)` marks
the pairs whose endpoint order `P` reverses.

> **Theorem.** `Γ_min = F₂^E` for every `n ≥ 2`.
>
> *Proof.* For the adjacent transposition `τ_k` the only pair whose endpoint order reverses is
> `(k,k+1)` itself, so `c(τ_k) = e_{(k,k+1)}`. `Γ_min` is `S_n`-invariant and `S_n` is
> edge-transitive on `K_n`, so `Γ_min` contains `e_f` for every edge `f`. ∎

Verified: `dim Γ_min = C(n,2)` exactly at `n = 3,4,5,6,7` (`3, 6, 10, 15, 21`), quotient `0`.

> **Corollary.** *No proper subgroup of the arc space has isomorphism-invariant orbits.*
> Not the star group, not switching, not `⋂_P Γ_P`, not any subgroup whatever. THM-1405's
> transversality is **forced**, and THM-1415's "switching is nearly trivial" is not a special
> weakness of switching — **every** candidate is dead for the same reason.

### (D′) The mechanism is exactly the affine twist

Drop `c` and keep only the linear part: `span{e_{P(f)} + e_f}` is the **even-weight
subspace** (edge-transitivity), of codimension **1**, and the surviving functional is total
parity `Σ_e x_e` — the parity of the number of arcs pointing against the reference order
(the *inversion parity*). So there is exactly one candidate linear invariant, and `c` kills
precisely it, because `Σ_e c(P)_e = inv(P) ≡ sgn(P)`.

Over `A_n` the shift `c(P)` always has even weight, so it cannot kill that functional:

| `n` | `C(n,2)` | `dim Γ_min` (`S_n`) | quotient | `dim` (`A_n`) | quotient |
|---|---|---|---|---|---|
| 3 | 3 | 3 | **0** | 2 | **1** |
| 4 | 6 | 6 | **0** | 5 | **1** |
| 5 | 10 | 10 | **0** | 9 | **1** |
| 6 | 15 | 15 | **0** | 14 | **1** |
| 7 | 21 | 21 | **0** | 20 | **1** |

> **The non-existence of linear tournament invariants is the non-existence of a canonical
> orientation.** Inversion parity is a perfectly good invariant of a tournament *plus a
> reference order*; relabeling by an odd permutation flips it. It is a **chirality**, not an
> invariant.

### (D″) Chirality: every class splits in exactly two

`|Aut(T)|` is odd, so every element of `Aut(T)` has odd order, so every element is a product
of odd-length cycles — each an **even** permutation. Hence `Aut(T) ≤ A_n`, and

`|S_n-orbit| = n!/|Aut|`,  `|A_n-orbit| = (n!/2)/|Aut|`,

so **every iso class splits into exactly two `A_n`-classes**, separated by inversion parity.
Verified by full enumeration:

| `n` | `S_n` classes | `A_n` classes | ratio | parity separates |
|---|---|---|---|---|
| 3 | 2 | 4 | 2.0 | yes |
| 4 | 4 | 8 | 2.0 | yes |
| 5 | 12 | 24 | 2.0 | yes |
| 6 | 56 | 112 | 2.0 | yes |

So the one thing that *does* descend is a `ℤ/2`-torsor: every tournament has a **handedness**,
well-defined up to even relabeling, and the two hands are exchanged by any odd relabeling.

## Honest scope

- (D) is proved for all `n ≥ 2`; the table is confirmation, not the proof.
- (D″)'s splitting is proved; the *separating* role of inversion parity is verified at
  `n ≤ 6` and follows immediately (`e(P·T) = e(T) + inv(P)`), so it is safe at all `n`.
- **Inversion parity is not min-FAS.** It is the parity of back-arcs against a *fixed reference
  order*, not the minimum over orders. No claim is made about `A003141` here.
- (A)'s "`A000568` even for `n ≥ 3`" is certainly not new as a fact; what is offered is the
  derivation from Rédei-oddness via the fibration.
- (B) is an engineering speedup with exact verification, not a theorem.
- (C) is **ceded to THM-1415**; only the restriction lemma is claimed here.
- Nothing here bears on LRC or the Jacobian threads.

*Artifacts:* `04-computation/tiling_node_edge_dictionary_macmini_S129.py` (+out).
*Credits:* THM-1415 (kind-pasteur-S128c110) for the switching identification, published first
and independently; THM-1405 for the transversality this explains; CLAUDE.md's Cut ⊕ Cycle
split; Rédei and the classical odd-automorphism fact, which do all the real work in (A)/(D″).

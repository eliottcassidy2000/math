# FINISHED: the D_{2p} irrep decomposition of Paley Hamiltonian paths (THM-131 generalization)

**Author:** opus, 2026-06-29 (dihedral<->tournament thread, prompted to finish proofs).
**Status:** PROVED (character theory) + VERIFIED exactly (p=3,7,11,19).
**Closes:** the explicitly-open "Generalization Conjecture" in
[THM-131](01-canon/theorems/THM-131-d14-irrep-decomposition.md) ("a, b need more careful analysis").
**Artifacts:** `04-computation/paley_dihedral_irrep_general_opus_20260629.py` + `.out`.

## Theorem
Let `p ≡ 3 (mod 4)` be prime, `T_p` the Paley tournament, `H = H(T_p)` its number of
(directed) Hamiltonian paths, and `f` the number of HPs fixed by a reflection of the
dihedral symmetry `D_{2p}` (equivalently the **anti-palindromic** HPs: `v_i = −v_{p−1−i}`,
forcing the middle vertex `0`). Then the HP permutation representation decomposes as
```
   V_HP  ≅  a·ρ_triv  ⊕  b·ρ_sign  ⊕  c·(ρ_1 ⊕ … ⊕ ρ_{(p−1)/2}),
   a = (H + p f)/(2p),   b = (H − p f)/(2p),   c = H/p   (each 2-dim irrep).
```
Consequences: `a + b = H/p` (= #`Z_p`-orbits), `a − b = f`, `a = ` #`D_{2p}`-orbits of HPs
(Burnside), and `b` is the **chirality / sign-isotypic** multiplicity.

## Proof
`G = D_{2p} = ⟨r,s | r^p=s^2=1, srs=r^{-1}⟩`, `p` odd. Conjugacy classes: `{e}`;
`{r^j,r^{-j}}` (`j=1..(p−1)/2`); and ONE class of all `p` reflections (`p` odd). Irreps:
`ρ_triv`, `ρ_sign` (`+1` on rotations, `−1` on reflections), and `(p−1)/2` two-dim `ρ_k`
with `χ_{ρ_k}(reflection)=0`. For the permutation character `χ_V(g)=#{HP : g·HP = HP}`:
- `χ_V(e)=H`.
- `χ_V(r^j)=0` (`j≠0`): a directed HP fixed by the shift `+j` would map its unique source
  `v_0` to `v_0+j=v_0`, forcing `j≡0`. (Rotations are orientation-preserving, so no
  reversal escape.)
- `χ_V(reflection)=f`, the same for every reflection (single conjugacy class).

Then `m(ρ)=\tfrac1{2p}[H·χ_ρ(e) + f·\sum_{p\ reflections}χ_ρ(refl)]` (the rotation terms
vanish since `χ_V(r^j)=0`), giving `a,b,c` above. ∎

## Verification (exact)
| p | H(T_p) | f | a=triv | b=sign | c=2dim | b/(a+b) |
|---|---|---|---|---|---|---|
| 3 | 3 | 1 | 1 | 0 | 1 | 0 |
| 7 | 189 | 9 | 18 | 9 | 27 | 1/3 |
| 11 | 95095 | 185 | 4415 | 4230 | 8645 | 0.4893 |
| 19 | 1172695746915 | 573057 | 30860700921 | 30860127864 | 61720828785 | 0.5000 |

All multiplicities integral; `a+b+2c·(p−1)/2 = H`; `a−b = f`. (`T_7` reproduces THM-131's
`18,9,27`.) The only ingredient left unclosed in general is a closed form for `f` (the
anti-palindromic / self-reverse-complement HP count, sequence `1, 9, 185, 573057, …`),
which appears genuinely hard — like `H(T_p)` itself it likely has no elementary formula.
But the **decomposition is complete given `f`**, which is what THM-131 left open.

## Chirality asymptotic (clean corollary)
`b/(a+b) = (1 − p f / H)/2`. Since `f = o(H/p)` (data: `pf/H = 1, 1/3, 0.021, 9·10^{-6}`),
**`b/(a+b) → 1/2`**: the Paley HP-representation becomes **maximally chiral** as `p→∞`. The
sign isotype is *absent* only at `p=3` (`b=0`, achiral) and present for all `p≥7`. This is
the representation-theoretic form of THM-139's "H-maximization rewards chirality."

## Synthesis: the dihedral analogy tournament ↔ LRC(14)
The same group `D_{14}=D_7` governs both sides, and on both the **sign isotype is the hard
part**:
- **Tournaments (T_7, here):** HPs `= a·triv ⊕ b·sign ⊕ c·(2-dim)`; `b=9` is the chirality,
  forced by the reflection being an **anti-automorphism** (`−1` a non-residue mod 7, THM-127/139).
- **LRC(14) (mac-mini S78 / kps S31av):** the cap `= mean ⊕ sign ⊕ de-Moivre-body`; the
  **sign-isotypic component is the Borsuk–Ulam / odd-degree obstruction** (the imaginary
  Gauss sum `i√7`), precisely because the complement is the *same* heptagon anti-automorphism
  acting freely. The trivial+2-dim (even/real) part is SOS/Brouwer-provable; the sign part is
  the genuine obstruction.

So the tournament-side "chirality multiplicity `b`" and the LRC-side "sign-isotypic
obstruction" are two instances of one phenomenon: the `7≡3 (mod 4)` anti-automorphism puts
weight on the sign representation of `D_7`. The finished tournament theorem makes the
*provable* shadow of the still-open LRC obstruction completely explicit.

---

## How it DIFFERS by p mod 4 (the predictable dichotomy)

Investigating the companion case `p ≡ 1 (mod 4)` (n=10,26) makes the whole picture
predictable from a single bit: **is `−1` a square mod p** (Euler). That bit decides whether
the heptagon/p-gon reflection (negation `v↦−v`) reverses arcs, and it propagates everywhere.

**The dihedral action on HPs is built differently:**
- `p≡3 (mod 4)` — negation is an **anti-automorphism**, so the reflection must be
  **negate + reverse** to land back in `HP(T_p)`. Reversal admits **anti-palindromic** fixed
  paths (middle vertex 0) ⟹ **`f > 0`**.
- `p≡1 (mod 4)` — negation is an **automorphism**, so the reflection is **negate only** (no
  reversal). A pure relabeling `v↦−v+a` fixes a unique *vertex* `a/2`, never a whole path
  ⟹ **`f = 0`** ⟹ **`a = b = H/(2p)` exactly**.

**Verified (exact):**
| p | p%4 | negation | H | f | a | b | b/(a+b) |
|---|---|---|---|---|---|---|---|
| 5 | 1 | aut | 10 | 0 | 1 | 1 | 1/2 |
| 13 | 1 | aut | 1579968 | 0 | 60768 | 60768 | 1/2 |
| 3,7,11,19 | 3 | anti-aut | … | >0 | >b | <a | <1/2 ↗ 1/2 |

**A free/fixed DUALITY between two spaces** (the reflection flips character):
| space the reflection acts on | p≡1 mod 4 | p≡3 mod 4 |
|---|---|---|
| tournament **orientations** | FIXES `T_p` (Brouwer) | FREE, `T_p↔T_p^op` (Borsuk–Ulam) |
| **Hamiltonian paths** | FREE, `f=0` | FIXED, `f>0` (anti-palindromic) |

So the same `p mod 4` that gives the LRC certificate dichotomy is a *theorem* on the
tournament side. The unified table:

| | p≡1 mod 4 (n=10,26) | p≡3 mod 4 (n=14,22,46) |
|---|---|---|
| `−1 mod p` | square | non-square |
| negation | automorphism | anti-automorphism |
| Paley `T_p` | self-complementary | self-converse |
| reflection-fixed HPs `f` | **0** | **>0** |
| HP multiplicities | **a = b** (balanced) | **a > b** (deficit `f`) |
| Gauss sum `Σχ(n)ζ^n` | `√p` (real) | `i√p` (imaginary) |
| quadratic field | `Q(√p)` real | `Q(√−p)` imaginary |
| LRC(2p) certificate | Brouwer / SOS / **even** | Borsuk–Ulam / odd degree / **imaginary** |

**Prediction for the LRC frontier.** Even `n = 2p` splits by `p mod 4`: the `p≡1` cases
(n=10, 26) are the **real / Brouwer / even-degree** side; the `p≡3` cases (n=14, 22, 46) are
the **imaginary / Borsuk–Ulam / odd-degree** side — intrinsically the hard ones — with the
class number `h(Q(√−p))` grading difficulty *within* that column (`h=1` for p=7,11 → n=14,22
gentlest; `h=3` for p=23 → n=46 harder). **n=14 (p=7≡3, h=1) is exactly the simplest hard
case** — which is why it is the first open one. The dichotomy is now a verified tournament
theorem, not just an LRC heuristic.

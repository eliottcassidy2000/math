# The reversal involution and the "half" principle: one Z₂ quotient under palindromic HPs, half-tilings, the half-circle, and palindromic polynomials — with its sign as the master discriminant

*opus-2026-06-29. Prompted to connect palindromic Hamiltonian paths to the half-tiling model
and surmise freely. The thread closes on a single object — the reversal involution `R` — that
the project keeps meeting under different names, and on the fact that its **sign**
`ε=(−1)^{C(n,2)}` is the same `p mod 4` discriminant that makes LRC(14) the hard/imaginary case.*

## One involution, four "half" objects

The **reversal involution** `R: P ↦ P^rev` (THM-084/088/051) sends `fwd ↦ (n−1)−fwd` and
position `j ↦ n−1−j`. It is the anti-diagonal symmetry, and the whole project's "half"
language is its fundamental domain:

| domain | the symmetric object | its **half** (fundamental domain of R) | status |
|---|---|---|---|
| **HPs (tournament)** | anti-palindromic HP `v_i=−v_{p−1−i}` | first `(p−1)/2` vertices + fold | **proved/verified here** |
| **tilings** | grid-symmetric (blue) tiling | half the staircase (one side of the anti-diagonal) | repo (`isGridSym`) |
| **LRC time-circle** | lonely set `L = 1−L` | `L ∩ (0,1/2)` | **verified here** |
| **F-polynomial** | (anti-)palindromic `F`/`SF` | coefficients `0..⌊(n−1)/2⌋` | THM-084/088 |

All four are the **same `Z₂` quotient by `R`**; "half-tiling," "anti-palindromic HP," "mirror
pair `{t,1−t}`," and "palindromic coefficient vector" are one structure in four costumes.

## A genuinely new reframing of `f` (the half-tiling principle, tournament side)

The dihedral reflection-fixed count `f` (which controls `a=(H+pf)/(2p)`, `b=(H−pf)/(2p)`) had
looked formula-less. The half principle re-expresses it: an anti-palindromic HP of `T_p`
(`p≡3 mod 4`) is `v_0…v_{m−1}, 0, −v_{m−1}…−v_0` (`m=(p−1)/2`), **fully determined by its first
half** — the second half and all its arcs are forced by the anti-automorphism. So

> **`f` = the number of Hamiltonian half-paths**: orderings of one representative from each
> antipodal pair `{±a}` with consecutive differences in `QR_p` and the fold condition
> `−v_{m−1}∈QR_p`. (Verified `f_half = f_brute`: 1, 9, 185 at p=3,7,11.)

So `f` is intrinsically a *half-sized* HP count — the tournament-side incarnation of the
half-tiling "the symmetric object is determined by half." This halves the computation and is
the right shape for a possible recursion (the half-system is itself a `QR`-difference digraph
on `(p−1)/2` signed pairs).

## The master discriminant: `ε = (−1)^{C(n,2)} = sgn(R)`

Reversing an `n`-ordering is the longest permutation `w₀`, `sgn(w₀)=(−1)^{C(n,2)}`. For
`n=p` prime this is `+1` (`p≡1 mod 4`) / `−1` (`p≡3 mod 4`). This single sign is what every
earlier dichotomy was measuring:

- **SF** is palindromic (`ε=+1`) vs **anti-palindromic** (`ε=−1`) — *proved*, THM-088.
- the dihedral HP rep has `f=0, a=b` (`ε=+1`) vs `f>0, a>b` (`ε=−1`) — *proved*, this session.
- Gauss sum `√p` real vs `i√p` imaginary; field `Q(√p)` vs `Q(√−p)`.
- LRC(2p) certificate Brouwer/SOS/**even** vs Borsuk–Ulam/odd-degree/**imaginary**.

So the LRC(14) obstruction being "odd/imaginary/sign-isotypic" — the wall every route this
week hit — is `ε(7) = (−1)^{C(7,2)} = (−1)^{21} = −1`, i.e. the **sign of reversing 7 things**.
The half-circle's reversal `t↦−t` carries the same sign obstruction: for a covering set its
fixed points (`t=0`, `t=1/2`) are non-lonely, so the lonely witnesses are forced into mirror
**pairs** with no symmetric representative — exactly the "antipodal pair, certified by odd
degree, not a fixed point to construct" picture, here made elementary.

## Surmise (free, unproven)

If `ε` is genuinely the obstruction class, then the LRC "signed certificate" everyone needs
should be expressible as the `R`-equivariant (sign-isotypic) part of the cap, and the
even/`R`-symmetric part — the "half" data on `(0,1/2)` — should be the SOS/Brouwer-provable
piece. The proof split the geometry already predicts (mac-mini S78: even/real by Euler/Brouwer,
odd/imaginary by Borsuk–Ulam) is then *literally* the `±1` eigenspace split of `R`. The
half-tiling lift→half quotient (HYP-3244) is the same `R`-quotient on the tournament side; the
"missing proof coordinate" it warns about is the `ε=−1` (sign) eigenspace that the half-chart
discards. **The thing to keep across any half/quotient compression is the `R`-odd part.**

## Status
- **Proved/verified:** the `f` half-reframing; `L=1−L` mirror-pairing; the `ε`-unification of
  THM-088 with the dihedral `f`-deficit (both `=(−1)^{C(p,2)}`).
- **Synthesis:** the four-costume `R`-quotient; `ε` as master discriminant across tournament,
  polynomial, number-field, and LRC-certificate dichotomies.
- **Open:** the LRC sign-isotypic certificate itself (= LRC(14)). Unchanged — but now named
  precisely as the `R`-odd / `ε=−1` eigenspace, the one piece the "half" compression drops.

Artifacts: `04-computation/half_principle_reversal_opus_20260629.py`,
`signed_hp_bridge_opus_20260629.py` (+ `.out`s). Related: THM-084, THM-088, THM-051, THM-127,
THM-131 (+ this session's generalization), HYP-3244, mac-mini S78, kps S31av, OPEN-Q-108.

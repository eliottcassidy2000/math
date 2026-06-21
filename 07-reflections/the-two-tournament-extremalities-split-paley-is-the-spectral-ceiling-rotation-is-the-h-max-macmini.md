# The two tournament extremalities split: Paley is the spectral ceiling, the rotation is the H-max, and Doyle-Holt picks the wrong one

*mac-mini-2026-06-21, THREAD B (the user's Doyle-Holt lead, tournament side).*

The user's unifying picture (HYP-2747) said: the LRC consec-max and the tournament
"Paley / H-max" are the **same** genuinely-aggregate theta' ceiling. Thread B set
out to test whether the Doyle-Holt / half-arc-transitive lens *gives* the Paley
extremality. Working it literally on the tournament side sharpens the picture and,
honestly, **corrects one leg of it**: there are **two distinct tournament
extremalities**, Paley realizes only one of them, and the Doyle-Holt
(arc-transitive) lens characterizes the *other* object than the one the user named.

## The two extremalities, decoupled

| extremality | object | symmetry | invariant type |
|---|---|---|---|
| **det(I+S) Hadamard ceiling** (THM-472) | **DRT** (Paley when it exists) | **arc-transitive**, SC, \|Aut\|=n(n-1)/2 | **spectral / weakly-regular** |
| **H-max** (Hamiltonian paths) | **rotation** (consec/carousel) for n≥13 | only cyclic, \|Aut\|=n, **arcO = m = (n-1)/2** | **strictly finer than spectrum** (THM-499) |

At n=7,11 these coincide on Paley (so the user's identification looked exact). They
**split at n≥13**: canon LEM-004 already proves the *circulant* H-max is the
**rotation** at n=13,15,17,19, with Paley only **3rd** at n=19
(H_rot(19)=1,184,212,824,763 > H_2nd > H_paley=1,172,695,746,915). The det(I+S)
ceiling stays on the DRT/Paley side throughout — a genuinely different, *spectral*
extremal (THM-499: H is not a function of the spectrum).

Computed directly this session (lrc14_threadB_rotation_vs_paley_symmetry_macmini):
Paley_n is arc-transitive (arcO=1, \|Aut\|=21,55,171 = the Frobenius groups), a DRT,
and attains det(I+S)=(n+1)^{(n-1)/2}; the rotation is **never** a DRT and **never**
arc-transitive (arcO = m), yet is the H-max. At n=15 (Paley-free) the rotation
(H=198,464,295) even beats the tower DRT (H=198,335,025).

## Why arc-transitivity is the *spectral* end, mechanically

For a circulant tournament, the cyclic group gives one arc-orbit per connection-set
element (per "skip"); the **multiplier group** can fuse them. Arc-transitivity ⟺ the
multiplier group acts **transitively on the connection set** ⟺ the set is a single
multiplicative orbit ⟺ **QR ⟺ Paley** (verified p=7,11,19). The rotation set
{1,…,m} is the *opposite* extreme: trivial multiplier action, m arc-orbits, minimal
symmetry. So "maximal multiplier symmetry" is exactly "spectral flatness =
Ramanujan = DRT" (THM-126), which is the **det** ceiling — **not** the H-max.

## The Doyle-Holt trichotomy (this is the real deliverable)

kps (HYP-2748) showed half-arc-transitive = vertex-transitive **non-self-converse**
(the converse Z_2 unrealized), needing a non-abelian carrier (Holt 27; F_21 at
n=21). Pairing that with the arc-transitivity mechanism above gives a clean
**trichotomy of the DRT's available symmetry across the orders n ≡ 3 (mod 4)**:

- **prime-power n** (7,11,19,23,27,31,43): Paley/QR exists → the DRT is
  **arc-transitive & self-converse**. Maximal symmetry. *The converse Z_2 IS realized.*
- **n = pq, p \| q−1** (21,39,55): no Paley, but a non-abelian Frobenius group of
  order n exists → vertex-transitive **NS** tournament. **Half-arc-transitive
  regime**: the converse Z_2 is **unrealized** (kps's F_21).
- **n = pq, p ∤ q−1** (15,35,51): only Z_n exists, and there are **0 circulant
  DRTs** → the DRT is forced **non-vertex-transitive, NS, fully defective**
  (tower DRT(15): vO=3, arcO=7, \|Aut\|=21).

Proof of the third bucket's obstruction is one line of group theory: a regular-VT
DRT on n=15 vertices would need a group of order 15; the only such group is Z_15
(15=3·5, 3∤4 ⇒ no non-abelian order-15 group); a Z_15-circulant DRT would be one of
the 2^7 antisymmetric sets, and exhaustive search finds **none**. So at 15 the
det-ceiling object *cannot* be arc-transitive or even vertex-transitive.

## The verdict for the user's question

**Doyle-Holt is the obstruction, not the route — for H.** The arc-transitive
(Doyle-Holt "symmetric, converse realized") end is precisely the **det(I+S) / DRT
spectral ceiling**, available only at prime-power orders, and is **not** the H-max
for n≥13. The H-max is the **rotation = consec = the linear/anti-MDS stratum** —
exactly the LRC extremizer (HYP-2602, additive-energy/Freiman, THM-538). So the
user's unification, sharpened, reads:

> The LRC consec-max and the tournament **H-max** are the same **linear/anti-MDS
> (rotation, low-symmetry)** extremality; the tournament **Paley/arc-transitive**
> object is a *different*, **spectral (det-ceiling)** extremality. The shared theta'
> wall — if there is one — runs through the **rotation/consec** side, not Paley.

This is consistent with mac-mini-S14's own collapse finding (theta'=L_y, no PSD
lift helps) and with the apex-prime reflection's "absolute \|E_2\| carrier vs signed
shadow": the H-max / consec object is the *aggregate, non-spectral* one, and the
spectrum-flat Paley is the clean-but-wrong shadow. Three pictures, one demotion of
Paley from "the extremizer" to "the spectral special point."

Open, honest: the **global** (non-circulant) H-max at n≥13 is not known
(search space 2^{C(n,2)}); at n=7 it *is* Paley (189). So the "rotation is H-max"
statement is proven only within the circulant family. But the *decoupling* of the
two extremalities, and the trichotomy of where arc-transitivity is even available,
are exact.

→ HYP-2747 (theta' unification, sharpened), HYP-2748 (Doyle-Holt = converse Z_2),
LEM-004 (rotation is circulant H-max), THM-472 (det ceiling = DRT), THM-499 (H
finer than spectrum), THM-126 (Paley spectral-flat), HYP-2602/THM-538 (LRC
consec-max), the everything-is-the-triangle foundation.

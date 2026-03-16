# Lossy and Lossless: What Geometry Determines

## The divide

**Lossless** = bijective homomorphism = isomorphism. No residue.
**Lossy** = fails bijectivity or homomorphism. Irreducible residue.

## Lossless compressions (no information lost)

- **Multiplication → rapidity**: ln(a·b) = ln(a)+ln(b). Unique factorization guarantees recovery. The rapidity vector determines n exactly.
- **Q: (-1,1) → (0,∞)**: bijective. Every velocity ↔ unique Q-value.
- **arctanh: (-1,1) → R**: bijective. Every velocity ↔ unique rapidity.
- **The formal group**: F_h(x,y) invertible given either input. The homomorphism f(a(+)b) = f(a)·f(b).
- **Ω(T) → H(T)**: one direction of the OCF. Given the conflict graph, H is determined.

## Lossy compressions (information destroyed)

- **H(T) → T**: many tournaments share one H-value. Non-injective.
- **Addition n+m → s**: the sum doesn't determine the summands. Non-homomorphic with multiplication.
- **8D → 3D**: tournament space to transfer matrix eigenspace. Loses 5 = pentagon dimensions.
- **Z[ω] → Z**: hexagonal to integer line. Loses the twist. Residue: Cassini ±1.
- **Cycle topology → H-spectrum**: some H-values (7, 21) are missed. Residue: the forbidden.

## What determines lossy vs lossless

A map is lossless iff it is an **isomorphism** of the source and target structures.

Failure mode determines residue type:
- **Non-injective** (many-to-one): residue = KERNEL = equivalence classes, orbits, Vitali atoms. This is the **symmetric** residue.
- **Non-homomorphic** (structure mismatch): residue = OBSTRUCTION = ±1, forbidden values, commas. This is the **intrinsic** residue.

## Q is the last lossless step

Q is a bijective homomorphism from F_h to multiplication. Q itself loses nothing.

Everything AFTER Q — discretization to integers, projection to H-values, restriction to tournaments — is lossy. The lossy part is the discretization, never Q.

## The geometry: curvature determines loss

- **Zero curvature** (flat, conformal): lossless. The map preserves angles.
- **Nonzero curvature**: lossy. The projection bends, creating distortion = residue.

The curvature quantum = 3. The first bend. The first loss.
The forbidden = 7. The bend too sharp to follow. Full obstruction.
The ±1 = the curvature at the atomic scale. The irreducible distortion.

Below 3: lossless. At 3: first loss. At 7: obstruction. Past 7: permanent ±1 residue.

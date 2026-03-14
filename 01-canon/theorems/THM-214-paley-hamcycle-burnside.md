# THM-214: Fixed Hamiltonian Cycles in Paley Tournaments

**Status:** PROVED
**Found by:** opus-2026-03-14-S89c
**Verified in:** `04-computation/pi_hamcycle_89c.py`

## Statement

For the Paley tournament P_p (p ≡ 3 mod 4, p prime), the cyclic shift σ: i ↦ i+1 mod p acts on directed Hamiltonian cycles. Exactly **(p-1)/2** directed Hamiltonian cycles are fixed by σ.

## Proof

A directed Hamiltonian cycle C = (v₀, v₁, ..., v_{p-1}, v₀) is fixed by σ iff the shifted cycle (v₀+1, v₁+1, ..., v_{p-1}+1) is a rotation of C. Since C visits all vertices, this forces C to be an arithmetic progression: v_i = v₀ + id (mod p) for some step d.

For C to exist in P_p, we need each arc v_i → v_{i+1} to be in the tournament, i.e., d ∈ QR(p). Similarly, the closing arc v_{p-1} → v₀ has difference d mod p.

Since -1 is a NQR when p ≡ 3 mod 4, the reverse direction -d is always a NQR when d is a QR. So each QR d gives exactly one directed fixed cycle (not its reverse).

There are (p-1)/2 quadratic residues mod p, hence (p-1)/2 fixed directed Hamiltonian cycles. □

## Corollary (Burnside)

The number of orbits of directed Hamiltonian cycles under Z/pZ is:

orbits = (hc + (p-1) × (p-1)/2) / p

where hc is the total count of directed Hamiltonian cycles.

## Verified Cases

| p | hc | Fixed | hc mod p | (p-1)/2 |
|---|-----|-------|---------|---------|
| 3 | 1 | 1 | 1 | 1 ✓ |
| 7 | 24 | 3 | 3 | 3 ✓ |
| 11 | 5505 | 5 | 5 | 5 ✓ |
| 19 | 34358763933 | 9 | 9 | 9 ✓ |
| 23 | 374127973957716 | 11 | 11 | 11 ✓ |

---
id: THM-1152
title: THE COMPATIBLE RELATIONS ENUMERATED MOD 4 — and my own "4|e ⟹ full concentration" criterion REFUTED by (1,2,7). (I) THE ENUMERATION, which is the question asked and is exact: the centre is c = (1/4,1/2,3/4), so a relation m = (m₂,m₃,m₄) is compatible iff m·c ∈ ℤ, i.e. **m₂ + 2m₃ + 3m₄ ≡ 0 (mod 4)**. This is a sublattice Λ ⊂ ℤ³ of **index 4** — exactly **16 of the 64 classes mod 4**, listed in full. All observed concentrating relations check out: (−2,1,0) and (−3,0,1) for (1,2,3) give 0 and 0; (1,1,−1) for (5,9,14) gives 0; while (9,−5,0), also a relation of (5,9,14), gives −1 and is NOT compatible — matching its merely partial 6.33×. (II) A FRAMING CORRECTION: I had said (1,2,3) is special for carrying "two independent compatible relations". That is wrong — the relation lattice R(d) = d^⊥ always has rank 2, and R(d) ∩ Λ has index ≤ 4 in it, hence rank 2 as well, so **every** direction carries two. The right criterion is whether **all** of R(d) lies in Λ, since the geodesic closure is the annihilator of R(d): the geodesic passes through c iff R(d) ⊆ Λ, i.e. iff e = 0 or 4 | e, where e generates {m·(1,2,3) : m ∈ R(d)}. (III) BUT THAT CRITERION DOES NOT GIVE THE MAXIMISER, and this refutes my own proposal: **96 non-proportional triples with d ≤ 25 satisfy 4 | e**, and they really do hit — **(1,2,7) has e = 4 and passes through c exactly at u = 3/4** (frac(−3/4), frac(−3/2), frac(−21/4)) = (1/4,1/2,3/4). Yet its sojourn is **7.27×**, not 28×. Other hitting directions: (1,6,7) 10.10×, (2,4,14) 7.28×, (1,6,11) 6.17×, (1,2,11) 4.29×, (1,10,15) 3.77×. Hitting the centre is **necessary but far from sufficient**. (IV) AND MY REPAIR ATTEMPT ALSO FAILS: I proposed sojourn ≈ hits × 2ρ/d_max, which would make sojourn·d_max invariant on the proportional family. It is not — it runs 0.2857, 0.5714, 0.8571, 1.4287, 2.8590 for m = 1,2,3,5,10, growing linearly, while the sojourn itself stays flat at 0.0952
status: (I) PROVED and complete — this is the enumeration requested, and it is exact. (II) PROVED — a genuine correction to my own THM-1151 framing. (III) REFUTES my own criterion, by explicit witness. (IV) refutes my own repair. **The maximiser claim remains unproved**: (1,2,3) is maximal among all directions tested, but the reason is still not identified. Uniform r=5 remains OPEN. Also recorded: my script's summary line hardcoded "NONE" for the 4|e search while its own output listed 96 hits — a reporting bug in my own code, caught on read
source: kind-pasteur-2026-07-18-S128 (cont.80; owner: enumerate the compatible relations mod 4)
depends_on:
  - THM-1151    # whose relation-compatibility rule this enumerates and whose criterion it refutes
related: [THM-1150, THM-1149, THM-1148]
script: 04-computation/compatible_relations_kps_S128c80.py, hit_vs_sojourn_kps_S128c80.py (+ .out)
---

# THM-1152 — the enumeration, and what it does not settle

## (I) The compatible relations, enumerated

The centre is c = (1/4, 1/2, 3/4), so m = (m₂,m₃,m₄) is compatible iff m·c ∈ ℤ:

> **Λ = { m ∈ ℤ³ : m₂ + 2m₃ + 3m₄ ≡ 0 (mod 4) }**, a sublattice of **index 4**.

That is **16 of the 64 classes mod 4**:

```
(0,0,0) (0,1,2) (0,2,0) (0,3,2)   (1,0,1) (1,1,3) (1,2,1) (1,3,3)
(2,0,2) (2,1,0) (2,2,2) (2,3,0)   (3,0,3) (3,1,1) (3,2,3) (3,3,1)
```

Checks against the observed concentrating relations:

| m | m·(1,2,3) | ≡0 mod 4 | origin |
|---|---|---|---|
| (−2,1,0) | 0 | ✓ | d₃−2d₂ = 0, for (1,2,3) |
| (−3,0,1) | 0 | ✓ | d₄−3d₂ = 0, for (1,2,3) |
| (1,1,−1) | 0 | ✓ | d₂+d₃−d₄ = 0, for (5,9,14) |
| (9,−5,0) | −1 | **✗** | also a relation of (5,9,14) |

(5,9,14) has one compatible relation and one incompatible one — matching its partial 6.33×.

## (II) A correction to my own framing

THM-1151 said (1,2,3) is special for carrying "two independent compatible relations". **That
is wrong.** R(d) = d^⊥ always has rank 2, and R(d) ∩ Λ has index ≤ 4 in it, so it is rank 2
too — every direction carries two. The right criterion is whether **all** of R(d) lies in Λ,
because the geodesic closure is the annihilator of R(d):

> geodesic passes through c ⟺ R(d) ⊆ Λ ⟺ e = 0 or 4 | e,

where e generates {m·(1,2,3) : m ∈ R(d)}, and e = 0 ⟺ d ∝ (1,2,3).

## (III) The criterion does not give the maximiser

**96 non-proportional triples with d ≤ 25 satisfy 4 | e.** And they genuinely hit:
**(1,2,7) has e = 4 and passes through c exactly at u = 3/4** —
(frac(−3/4), frac(−3/2), frac(−21/4)) = (1/4, 1/2, 3/4). Yet:

| direction | hits c | d_max | sojourn | ×\|B\| |
|---|---|---|---|---|
| (1,2,3) | ✓ | 3 | 0.09523 | **28.28** |
| (1,6,7) | ✓ | 7 | 0.03401 | 10.10 |
| (2,4,14) | ✓ | 14 | 0.02450 | 7.28 |
| **(1,2,7)** | ✓ | 7 | 0.02449 | **7.27** |
| (1,6,11) | ✓ | 11 | 0.02078 | 6.17 |
| (1,10,15) | ✓ | 15 | 0.01270 | 3.77 |

**Hitting the centre is necessary but far from sufficient.** My criterion identifies a large
family, not the maximiser.

## (IV) My repair attempt also fails

I proposed sojourn ≈ hits × 2ρ/d_max, which would make sojourn·d_max invariant on the
proportional family. It is not:

| m | 1 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|
| sojourn | 0.0952 | 0.0952 | 0.0952 | 0.0953 | 0.0953 |
| sojourn·d_max | 0.2857 | 0.5714 | 0.8571 | 1.4287 | 2.8590 |

The sojourn is flat; the product grows linearly. The scaling story is wrong.

## Honest status

The enumeration asked for is **done and exact** — Λ has index 4, 16 classes, listed. But it
does not finish the maximiser proof, and two of my own claims died in the process: the
"two independent relations" framing (II) and the "4|e ⟹ full concentration" criterion (III),
plus the scaling repair (IV). (1,2,3) is maximal among everything tested, and I still cannot
say why. **Uniform r=5 remains open.**

Also recorded against myself: the first script's summary line printed "NONE — e = 0 is the
only way" while its own output immediately above listed 96 hits. I had hardcoded the
conclusion into the print statement before running it. Caught on read, but it is the same
error I logged in cont.53 and said I would not repeat.

## Named next
- The invariant that actually separates (1,2,3) from the other 96 hitting directions is
  still unidentified. Candidates worth measuring: the number of distinct centres hit (six
  are available), the chord geometry at each hit, and the local density of hits rather than
  d_max alone.
- Given three failed mechanisms and two failed criteria on this one question, the honest
  move may be to stop proposing and instead compute the sojourn exactly as a function of d —
  it is a piecewise-rational quantity and a closed form would settle the maximiser outright.

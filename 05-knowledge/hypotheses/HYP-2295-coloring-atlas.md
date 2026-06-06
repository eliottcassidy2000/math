---
id: HYP-2295
title: The coloring atlas — every problem as a graph coloring (vertex / edge / both); chromatic polynomial = the master partition function
status: OPEN (atlas/synthesis); the LRC χ/α=C_n, the pair-sum edge-coloring, and P(C_n,k)=Potts are VERIFIED
source: claudebox-2026-06-03-S633
related:
  - HYP-2292  # tie-deleted graphs (C_n, Cayley/circulant) — the graphs we color
  - HYP-2245  # partition functions everywhere (chromatic poly = Potts = the master Z)
  - HYP-2288  # trienement: tie=resonance, symmetry controls resonance
  - THM-401   # pair-sum sieve modulus 2n-1 (= the edge-coloring / 1-factorization)
  - HYP-2130  # rigidity = orbit-type
---

# HYP-2295 — everything as a graph coloring

Building on S632 (the tie/conflict graphs are circulants/Cayley graphs) and S626 (partition
functions everywhere): **the chromatic polynomial `P(G,k)` is a partition function (zero-temperature
Potts), so "colorings everywhere" = "partition functions everywhere".** The atlas below reframes our
problems as colorings — **vertex**, **edge**, or **both** — with the controlling chromatic invariant.

## VERTEX colorings (color the nodes)

- **LRC sieve = `χ` of the tie-graph; corrector = `α`.** Tight LRC tie-graph `= C_n` (S632);
  **VERIFIED:** sieve arity `= χ(C_n) = 2` (n even) / `3` (n odd); corrector size `= α(C_n) = ⌊n/2⌋`.
  `P(C_n,2) = 0` for odd n (no 2-sieve) / `2` for even n — the single-vs-multi-sieve dichotomy
  (HYP-2075) as 2-colorability, the **even/odd 2-adic seam as `χ`**.
- **Unit distance = the chromatic number of the plane (Hadwiger–Nelson).** Native vertex coloring;
  `χ(plane) ∈ {5,6,7}` (de Grey ≥5). The grid disproof = maximize edge density = a `χ`-forcing /
  Hadwiger–Nelson construction (HYP-2230/2288).
- **Tournament metagraphs:** `χ(G_n) = n-1`; `χ(E_n)` (even graphs) `= 2,3,5,10,28` (canon). Coloring
  iso-classes.
- **Collatz / shell tower:** color `ℤ` (resp. `(ℤ/m)*`) by 2-adic / doubling-orbit class — the
  multiplier ±-pairs (S625) are a 2-coloring; the "free pair" (dodge) = a free color.

## EDGE colorings (color the pairs)

- **The pair-sum sieve (THM-401) = a proper edge-coloring of `K_n`** by `(v_i+v_j) mod (2n-1)`.
  **VERIFIED** (n=5,6,7): edges sharing a vertex get distinct labels — i.e. the **round-robin
  1-factorization / circle-method schedule**. The sieve modulus `2n-1` is the edge-coloring modulus;
  looseness = a "free" edge-color (a sum-residue off the band).
- **Tournament arc 2-coloring:** the cut/cycle (score/Walsh) split = a 2-coloring of arcs (canon).
- **Round-robin scheduling / Latin squares:** edge-coloring `K_n` = 1-factorization — the same object
  as the pair-sum sieve. (Engineering: scheduling, code design.)

## BOTH (total / mixed colorings)

- **The trienement** (S631/HYP-2288): nodes = points/runners, edges = the `<`/`=`/`>` relation; the
  tie (`=`) is an *edge color* (the resonance), the order is an *arc color*, and loneliness/unit-
  distance is a *node* property — a genuinely mixed (total) coloring.
- **The metagraph G_n/Z_2** (canon spine/ribs/sea): SC-SC / SC-NS / NS-NS is an **edge 3-coloring**
  by endpoint type, layered over a vertex 2-coloring (SC vs NS). Both at once.

## The unification: chromatic polynomial = the master partition function

`P(G,k) = Σ_colorings 1` is the `q=k`, zero-temperature **Potts partition function**; the covering-
depth `Z(z)` (S626) is the same object deformed (fugacity ↔ temperature). So the master bivariate
partition function (HYP-2245) specializes to:
- `z→0` / ground state ↔ proper colorings (loneliness = a free color);
- `χ` ↔ the smallest palette with a proper coloring (sieve arity / Hadwiger–Nelson);
- `α` ↔ a single color class (corrector / lonely set);
- chromatic index `χ'` ↔ edge sieve / round-robin.

**Resonance = an edge (a conflict); the problem is to color around the conflicts; the symmetry of the
conflict graph (Cayley/circulant, the perspective key HYP-2130) sets `χ`, `α`, `χ'`.** Trienement +
symmetry (HYP-2288) + colorings (here) + partition functions (HYP-2245) are four views of one object.

## To do
1. Prove "tight LRC ⟺ tie-graph `= C_n`" and read the sieve/corrector as `χ`/`α` of `C_n`; recover
   the pair-sum sieve arity (THM-401) as the chromatic index of the `K_n` edge-coloring.
2. Make `P(G,k) = Z(z)` precise (which specialization) for the covering-depth graph.
3. Hadwiger–Nelson ↔ LRC: is the chromatic number of the plane the continuous limit of `χ(C_n)` /
   the shell circulant `χ`? (the discrete-circle unit-distance graph → the plane.)

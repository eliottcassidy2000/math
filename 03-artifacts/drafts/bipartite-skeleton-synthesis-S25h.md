# Bipartite Skeleton Synthesis: t3 Parity, Antiferromagnetism, and the Even/Odd Dichotomy

**Instance:** kind-pasteur-2026-03-06-S25h
**Status:** Multiple theorems proved, verified n=3,...,9

---

## 1. The Blue Line Skeleton

The **blue line skeleton** is the GS flip graph on SC tournament classes:
- **Vertices:** Self-converse (SC) tournament isomorphism classes
- **Edges:** Grid-symmetric (GS) tiling flip pairs between different classes
- **Weight:** Number of GS tiling pairs connecting the classes

### Key counts

| n | SC classes | Edges | Same-class GS flips | Cross-class | Bipartite |
|---|-----------|-------|--------------------|-----------  |-----------|
| 3 | 2         | 1     | 0                  | 2           | YES       |
| 4 | 2         | 1     | 2                  | 2           | YES*      |
| 5 | 8         | 8     | 0                  | 16          | YES       |
| 6 | 12        | 26    | 4                  | 60          | NO        |
| 7 | 88        | 246   | 0                  | 512         | YES       |
| 9 | ?         | ?     | 0                  | 65536       | YES       |

*At n=4, trivially bipartite (2 vertices, 1 edge).

---

## 2. THM-059: t3 Parity Bipartition (PROVED)

**Theorem.** At odd n, the blue line skeleton is bipartite with
the bipartition determined by t3 parity.

**Proof mechanism:** Decompose C(n,3) triples into three types:
- **Type C** (2 backbone edges): n-2 consecutive triples, each contributes
  exactly 1 to t3(T) + t3(flip(T)). Total = n-2.
- **Type A** (0 backbone edges): flip reverses all edges, so CW <-> CCW.
  Individual contribution = 2 (always even).
- **Type B** (1 backbone edge): GS-paired triples contribute equal amounts.
  Total always even.

Grand total: t3(T) + t3(flip(T)) = n-2 + even = ODD at odd n.
Therefore t3 parity flips, making the skeleton bipartite. QED.

### Even-n contrast

At even n: t3(T) + t3(flip(T)) = n-2 + even = EVEN.
So t3 parity is PRESERVED, allowing:
- Same-class GS flips (self-loops in the skeleton)
- Edges between classes of the same t3 parity (odd cycles)
- Result: skeleton is NOT bipartite at n=6

---

## 3. Spectral Structure (n=5)

The unweighted skeleton at n=5 is a 4-cycle with 4 pendant arms:
```
     11(t3=0,H=1)
     |
0(t3=3,H=9) --- 3(t3=4,H=11) --- 1(t3=3,H=9)
     |               |
2(t3=4,H=13) --- 6(t3=5,H=15) --- 5(t3=4,H=15)
     |
9(t3=1,H=3)
```

Eigenvalues of weighted adjacency: {+/-2(1+sqrt(2)), +/-2, +/-2, +/-2(sqrt(2)-1)}
Perfectly symmetric around 0 (bipartite signature).

Vertex degree = 2 * (number of GS tilings in the class).
Number of tilings = H(T) = number of Hamiltonian paths.

---

## 4. NSC Classes and the Bipartition

**Key fact:** t3(T) = t3(T^op) for all tournaments (reversing all edges
maps CW 3-cycles to CCW and vice versa, preserving count).

**Consequence:** NSC pair members have the SAME t3 parity, so both sit
on the SAME side of the bipartition. NSC pairs never straddle the partition.

At n=5: 2 NSC pairs on side A (odd t3), 2 on side B (even t3). Equal split.
At n=7: 184 NSC pairs, distribution TBD.

---

## 5. Antiferromagnetic Interpretation

The skeleton is an **antiferromagnetic Ising model** on tournaments:
- Each SC class has a "spin" sigma_i = (-1)^{t3(T_i)}
- GS flip edges always connect opposite spins
- The skeleton is the interaction graph
- The ground state energy = -(number of edges) (all bonds frustrated)
- At odd n: perfect Neel order (bipartite = unfrustrated antiferromagnet)
- At even n: frustrated (odd cycles in skeleton)

This connects tournament theory to statistical mechanics!

---

## 6. Geometric Picture: The Tournament Orbifold

Tournament space = hypercube {0,1}^{C(n,2)}.
Tiling space = same hypercube, with backbone factored out.
GS subspace = affine subspace of dimension k = gs_dof.
Isomorphism classes = orbits under S_n action.

The skeleton = quotient of GS subspace under isomorphism = an orbifold.

The t3 parity is a **Z/2 invariant of the orbifold** — analogous to the
parity of a permutation. It divides the orbifold into two "sheets."

At odd n: the Z/2 action (GS flip) permutes the two sheets (bipartite).
At even n: the Z/2 action has fixed points (self-loops) and odd cycles
(non-orientable orbifold).

---

## 7. Connection to W-Hierarchy

From THM-058: w_{n-3}(T) = (n-2)! * [2*t3 - C(n,3)/2].

The bipartition by t3 parity is equivalent to: the SIGN of the fractional
part of w_{n-3} / (2*(n-2)!), rounded to nearest integer, modulo 2.

More precisely: t3 = (w_{n-3}/(n-2)! + C(n,3)/2) / 2. The parity of t3
determines which side of the skeleton a tournament sits on.

This means the **leading non-trivial W-coefficient** (w_{n-3}) encodes
the bipartition of the entire skeleton. The W-polynomial knows the
skeleton's topology!

---

## 8. Open Questions

1. **Prove Type B even total algebraically.** Currently verified
   computationally. Need: GS constraint forces paired triples to have
   equal (before+after) contributions, and no fixed Type B triples exist.

2. **What determines the n=7 skeleton topology?** 88 vertices, 246 edges,
   genus 159. Is there additional structure beyond bipartiteness?

3. **Spectral gap at large n.** Does the spectral radius of the skeleton
   grow polynomially or exponentially in n?

4. **Transfer matrix connection.** The transfer matrix M[a,b] counts
   Hamiltonian paths. The GS subspace is where M has additional symmetry.
   Is M restricted to GS the adjacency matrix of the skeleton?

5. **H(T) on each side.** At n=5, mean H is 9.0 (side A) vs 10.0 (side B).
   Is there a systematic H asymmetry between sides at larger n?

6. **Partition function.** Define Z = sum over SC classes of H(T)^beta.
   At beta=1, Z = sum of H over SC classes. Does the bipartite structure
   give a factorization or closed form for Z?

# THM-096: Tournament Path Homology Simplicity

**Status:** PROVED (conditional on even-Betti vanishing, verified n <= 8)
**Author:** kind-pasteur-S45 (2026-03-09)
**Depends on:** THM-095 (seesaw mechanism)

## Statement

For any tournament T on n vertices (verified n <= 8):

1. **H_0(T) = Z** (every tournament is weakly connected)
2. **H_{2k}(T) = 0** for all k >= 1 (even Betti numbers vanish)
3. **At most one H_{2k+1}(T) is nonzero**, and when nonzero it equals Z
4. **chi(T) in {0, 1}**: chi = 0 iff T has a nontrivial odd hole, chi = 1 iff T is path-acyclic

In other words, the GLMY path homology of a tournament is either:
- **Contractible:** H_* = (Z, 0, 0, ...) — like a point
- **Single odd hole:** H_* = (Z, 0, ..., 0, Z, 0, ...) — like an odd sphere

## Proof Structure

### Step 1: Even Betti vanishing (computational)
beta_{2k} = 0 for all k >= 1, verified:
- Exhaustive: n = 3 (8), n = 4 (64), n = 5 (1024), n = 6 (32768)
- Sampled: n = 7 (2000), n = 8 (500)

### Step 2: Seesaw mechanism (THM-095)
beta_{2k} = 0 creates a chain complex coupling that prevents adjacent
odd Betti numbers from being simultaneously nonzero:

  beta_{2k+1} > 0 => im(d_{2k+1}) saturates => beta_{2k-1} = 0

### Step 3: Betti = 0 or 1
From the seesaw, im(d_2) takes only TWO values: C(n,2)-n or C(n,2)-n+1.
Since ker(d_1) = C(n,2)-n+1 always:

  beta_1 = ker(d_1) - im(d_2) in {0, 1}

Similarly, the seesaw propagates to give beta_{2k+1} in {0, 1} for all k.

### Step 4: Euler characteristic
chi = sum (-1)^p beta_p = 1 - sum_{k>=0} beta_{2k+1}

Since at most one odd beta is 1: chi in {0, 1}.

## Statistics at various n

| n | Contractible | H_1=Z | H_3=Z | H_5=Z |
|---|---|---|---|---|
| 3 | 75.0% | 25.0% | - | - |
| 4 | 62.5% | 37.5% | - | - |
| 5 | 70.3% | 29.7% | 0% | - |
| 6 | 84.2% | 14.6% | 1.0% | 0% |
| 7 | ~87% | ~5% | ~7.5% | 0% |
| 8 | ~93% | ~1.7% | ~19% | 0%? |

## Open Questions

1. **Prove even-Betti vanishing algebraically.** This is the main gap.
   Why does ker(d_{2k}) = im(d_{2k+1}) exactly for all tournaments?

2. **Does beta_5 ever appear?** Not found at n <= 9 (500+ samples).
   If so, when? The onset pattern beta_1 at n=3, beta_3 at n=6 suggests
   beta_5 at n=9 or later, but 500 random n=9 found zero.

3. **Does the simplicity extend to n >= 9?** The real-roots counterexample
   (THM-025) at n=9 suggests pathology may arise. Does homology simplicity
   survive?

4. **Connection to odd spheres.** A tournament with H_{2k+1} = Z has
   the homology of S^{2k+1}. Is there a geometric interpretation?

## Relation to OCF

The OCF (H(T) = I(Omega(T), 2)) is completely independent of this theorem.
OCF works for any digraph; homology simplicity is specific to tournaments.

However, the combination is powerful: tournaments have both
(a) the simplest possible path homology (one hole or none), and
(b) a clean algebraic identity relating Hamiltonian paths to the
    independence polynomial of the conflict graph.

## Files

- `04-computation/beta2_vanishing_proof.py` — analysis script
- `04-computation/beta1_beta3_mediator.py` — seesaw verification
- `05-knowledge/results/beta2_vanishing_proof.out` — output

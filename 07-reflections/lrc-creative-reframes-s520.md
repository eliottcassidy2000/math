# Reflection: LRC through 10 creative lenses

**Session:** opus-2026-06-01-S520
**Date:** 2026-06-01

## The big picture

Nine radically different mathematical frameworks were applied to LRC, with
computation backing each one. The tournament analysis of the reframes gives
a strict ranking: **topological > spectral > {LLL, relativistic, thermo}**.

## The five most promising new angles

### 1. Rapidity threshold = 0.5 · ln(n-1)

Map each runner's gap g_i = ||v_i t|| to rapidity r_i = atanh(1 - 2g_i).
The lonely condition g_i ≥ 1/n becomes r_i ≤ 0.5·ln(n-1).

The formal group F(x,y) = (x+y)/(1+xy) = tanh(atanh(x) + atanh(y)) is
ADDITION IN RAPIDITY SPACE. So the "total rapidity load" on the observer is
just the sum of individual rapidities. The lonely condition is that all
individual loads are small.

The rapidity budget is (n-1) · 0.5 · ln(n-1), which grows as O(n·log(n)).
This is the project's central object (the formal group) applied to the
project's central open problem (LRC). The connection has not been made
this cleanly before.

### 2. Billiard: good box volume → e^{-2}

The "good box" (all runners in [1/n, 1-1/n]) in the (n-1)-torus has volume
((n-2)/n)^{n-1} → e^{-2} ≈ 13.5%. This is CONSTANT for large n —
a striking fact. The rational-line trajectory massively under-hits this box
(often 0% for initial segments, wall-only witnesses). Why? The trajectory
lies in a low-dimensional rational subtorus that barely grazes the box.

The billiard question: can a rational line in the torus always hit a
box of volume ~13.5%? This is a simultaneous Diophantine approximation
question with a very specific target set.

### 3. Thermodynamic energy bound

Define E(t) = -Σ log(||v_i t||). LRC is equivalent to min_t E(t) < (n-1)·log(n).

At all tested n (3-14), the minimum energy is well below threshold. The
"heat capacity" (Var[E]) grows linearly with n. The temperature parameter
β gives a Boltzmann distribution exp(-βE) that concentrates on the
ground state (most lonely time) as β → ∞.

This reframe suggests: use concentration of measure on the "energy landscape"
to prove the minimum is always below threshold.

### 4. Odd/even parity complementarity

For odd speed v: ||vt|| + ||v(t+1/2)|| = 1/2 exactly.
For even speed v: ||v(t+1/2)|| = ||vt|| (invariant).

So for n ≥ 4: if an odd-speed runner is bad at t, it's automatically good
at t+1/2 (with gap 1/2 - ||vt|| ≥ 1/2 - 1/n ≥ 1/n).

Combined with THM-387: the gap race after a wrap-around at time t_0 can be
analyzed separately for odd and even speeds. Odd speeds that are bad at t_0
are guaranteed good at t_0+1/2. This might reduce the gap-race analysis
to even speeds only.

### 5. Left-right bid anti-correlation

In the "auction" reframe, the correlation between left and right winning
bids is strongly negative: -0.85 (n=3) to -0.48 (n=6). This is the
auction-theory version of THM-387's directed fiber flow. The anti-correlation
decreases with n, suggesting the "joint auction" becomes more competitive
at larger n but never reaches independence.

## What DOESN'T work

- **LLL fails** for initial segments (p=2/n is too large relative to the
  dependency degree). The LLL approach needs sparse speed sets with many
  coprime pairs.
- **Spectral approach is inconclusive** for initial segments (the DC component
  is too small relative to the AC sum). Needs improvement, perhaps via
  a weighted Fourier basis.
- **Chip-firing** doesn't add much beyond the gap-race language.

## The mathematics pointing beyond itself

The rapidity connection is the deepest finding. The formal group
F(x,y) = (x+y)/(1+xy) — which controls Hamiltonian path counts H(T),
the fiber fraction, the Walsh transform, and the Cayley-Dickson tower — also
controls the "loneliness intensity" in the LRC problem. The threshold
rapidity 0.5·ln(n-1) is the Riemannian distance on the hyperbolic line
from the observer to the safe arc boundary.

This is the first time the formal group has been connected to LRC in a
quantitative, computation-backed way. It suggests that the H(T) invariant
(Hamiltonian path count) and the LRC lonely measure may be controlled by
the same underlying structure: the hyperbolic geometry of the formal group.

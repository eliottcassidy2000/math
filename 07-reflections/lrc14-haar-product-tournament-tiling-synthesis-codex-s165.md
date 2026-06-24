# LRC14 Haar Product / Tournament Tiling Synthesis - Codex S165

This pass connected three threads that had been developing in parallel:

```text
discrepancy theory and Haar fronts,
the 2D Haar product rule,
the fixed-path tournament tiling cube.
```

The common structure is not a visual analogy.  It is an algebraic warning:
products are exact before quotienting and dangerous after scalarization.

## Product Rule

For dyadic Haar packets,

```text
h_{I x J} h_{I' x J'} = (h_I h_I') tensor (h_J h_J').
```

In one dimension, `h_I h_J` is zero, an indicator, or a signed finer Haar atom,
depending on whether the intervals are disjoint, equal, or nested.  The 2D
rule is just the coordinatewise product of those cases.

For fixed-base tournament tilings, the corresponding basis is Boolean Walsh:

```text
chi_A chi_B = chi_{A xor B}.
```

The staircase tile address is the support of the character.  Strip counts,
Hamming weights, coefficient profiles, and tournament isomorphism classes are
quotients of that support address.

## Why This Touches LRC14

The recent agents have all been saying the same thing in different languages.
HYP-2984 tracks when kernel deformation preserves a certificate or creates a
boundary defect.  HYP-2985 routes packets by admissible smoothing clocks.
HYP-2987 turns the proof stack into a labelled handoff atlas.  HYP-2986 says
AP/GW are boundary cocircuit packets, not empty scalar mass.

The Haar product rule explains why this convergence is structural.  A
discrepancy coefficient is not just a number.  It has a rectangle address, a
sign, a support boundary, and a parent packet.  If we collapse it to a strip
count or raw tournament class, products no longer know where the residual
coefficient lives.

The tournament tiling model already learned this: the quotient from the
staircase cube to a profile is useful only when the fiber is homogeneous or the
lost support coordinate is restored.  That is exactly the quotient guardrail
from HYP-2978/HYP-2979.

## The Missing Picture

The LRC14 proof object should be a Haar/Walsh packet sheaf:

```text
q/Farey scale
regular-open Haar mass
closed boundary cocircuit debt
endpoint owners
Ramanujan/exact-period phase
Fejer/Toeplitz coefficient packet
kernel/smoothing clock
state-lift or F7 residual
```

The theorem target is not "find the right scalar discrepancy bound."  It is:
prove that every quotient used by the proof is a homomorphism for the relevant
packet product, or else emits the lost coordinate as a named residual.

That makes the tournament tiling model useful again.  Its fixed Hamiltonian
path is the tie path; its triangular support grid is the Boolean Haar address;
its scalar shadows are allowed only after product-rule descent is certified.

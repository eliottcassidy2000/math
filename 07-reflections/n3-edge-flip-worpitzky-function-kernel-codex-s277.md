# n=3 Edge-Flip Worpitzky Function Kernel Reflection

Source: codex-2026-06-28-S277  
Hypothesis: HYP-3147
Script: `04-computation/tournament_n3_edge_flip_worpitzky_codex_s277.py`

The useful small object is the collapsed edge-flip kernel:

```text
P(C -> T)=1
P(T -> C)=1/3
P(T -> T)=2/3
```

with stationary distribution `(C,T)=(1/4,3/4)` and nontrivial eigenvalue
`-1/3`.  This is just the three-coin split after quotienting the cube by the
two n=3 tournament score classes: straight words `HHH,TTT` are cyclic, and
all `2:1` mixes are transitive.

The edge-level asymmetry is the part worth keeping: in a cyclic triangle,
every edge flip breaks the cycle; in a transitive triangle, exactly one
minority edge flip closes the cycle while the other two flips stay transitive.
So the two-node class graph is too coarse unless the minority edge is retained.

Worpitzky enters as the first refinement after entering the transitive fiber.
The six transitive states have descent counts `1,4,1`, exactly Eulerian
`A(3,k)`, verifying `x^3 = C(x+2,3)+4C(x+1,3)+C(x,3)`.  In this reading,
Worpitzky is a function-basis for ordered path data inside the `T` fiber, not
a replacement for the `C/T` transition kernel.

The four functions from the prompt split perfectly:

```text
a+b, a*b   commute and forget orientation;
a^b, b^a   swap under orientation reversal.
```

This gives a concrete quotient rule for LRC packets: symmetric edge functions
are useful shadows, but they cannot certify edge-flip behavior.  A proof-facing
edge packet needs an ordered function channel or an explicit orientation
sidecar.

Pull into current LRC work:

- HYP-3141 edge witnesses should carry `edge_flip_class_kernel`,
  `minority_edge_gate`, `worpitzky_descent_word`, and
  `ordered_function_payload`.
- HYP-3142's antisymmetric shell should be treated as an order-sensitive
  payload rather than residual noise.
- HYP-3143's packet-order audit should recurse through n=3 subtriangles before
  trusting an n=4 quotient.
- HYP-3129's signed SPEC low modes should be checked against the local `-1/3`
  eigenmode.

Assumption challenged: the vertices are not runners, and not merely tournament
isomorphism classes.  The live vertices are edge flips, coin words,
minority/majority roles, descent words, and ordered endpoint functions.

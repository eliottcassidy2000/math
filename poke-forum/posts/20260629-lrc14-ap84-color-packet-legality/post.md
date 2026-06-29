# LRC14 AP84 color-packet legality matrix

HYP-3459 adds a small exact bridge from the older coloring/discrepancy work to
the current AP84 tail.

The moral is precise: a coloring is a quotient, so it must declare what LRC
predicate it preserves and what coordinate it forgets.

For the AP84 family

```text
S_m={1,2,...,11,13,84m}
```

the current colors are not interchangeable:

```text
gate_bucket_hist={'both_outer_inner':24,'outer_only_7_gate':6,'inner_only_5_gate':4,'clean_only_lcm35':1}
correction_hist={0:1,1:22,2:12}
both_outer_inner -> d values [1,2]
d=1 -> both_outer_inner, inner_only_5_gate, outer_only_7_gate
mixed_gate_plus_correction_fibers_count=0
```

So gate color and floor color each forget something, but the pair is clean for
the AP84 local target.  Then HYP-3457/HYP-3454 add the nonperiodic phase
warning:

```text
m=1 and m=36 share residue r=1,
but m=1 is a finite mixed endpoint transient,
while m=36 is rank-one E/E.
```

Suggested next concrete object:

```text
rows    = HYP-3438 survivor gates/components
columns = gate, floor, height, endpoint_phase, branch_mirror,
          incident_C3_Qsqrt, zipper_cocycle
```

Target theorem: the AP84 local-to-global splice into HYP-3439 is a homomorphism
for this seven-color packet product.  If not, the first failed color should
route to HYP-3453/HYP-3455, owner-current, two-adic descent,
exact-period/state-lift, or SPEC debt.

This is the Haar zipper lesson in finite AP84 form: margins and colors are
useful only while their lost cocycles are either reconstructed, killed, or
named.

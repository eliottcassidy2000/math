# Planar-JC concept board — 2026-08-26

Status labels below are local session labels; only cited canon is truth.

| lane | object / representation | predicate | invariant / operation | lost coordinate | cheapest decisive test |
|---|---|---|---|---|---|
| Anchor | hidden lattice `L`, degree-24 shell | twelve attachments do not map to `O` | `mu_6 x <T>` orbits; reciprocal denominator | exact marked ratio after modular reduction | enumerate orbits, bound denominator degree, compute symbolic/modular reciprocal gcd |
| Niche | full lattice `M`, visible degrees `12,16,28` | twelve points cannot share a fibre | visible/hidden projection and fibre degree | two-torsion glue means projection collapse is only necessary | classify low visible maps and attachment fibres integrally |
| Wildcard | cyclic differences `(T^k-1)m` | every difference vanishes on a collapsed orbit | exact integral `C_12` action and degree form | rational projectors can erase the `O/2O` glue | enumerate response shells and find every positive difference degree `<12` or hidden degree `12` |
| Hostile | index-four class `h=(v+omega^2 f+g)/2` | quotient arguments survive gluing | `O/2O` residue line | integral lift | verify the `T` action on `[u,f,g,h]` and invariance of `q` vector-by-vector |
| Entry | exact-`M=12` seam outside `W=0` | reduce or exclude hidden-`E0` locus | Prym variation / wall degeneration | specialization does not classify the whole Hom locus | inspect incoming THM-4248/4249 and test a transverse parameter derivative |

## Current exact signal

On the normalized basis, the lawful action is

```text
T(u)=-omega u,  T(f)=g,  T(g)=-omega f,
T(h)=omega^2 h-omega f.
```

If `m` collapses the attachment orbit, every `(T^k-1)m` vanishes on all
twelve points.  Exact shell enumeration finds:

```text
q(m)=34: 288 vectors have q((T^3-1)m)=12, all at (q(v),q(ell))=(112,24)
q(m)=42: 144 vectors have q((T^3-1)m)=12, all at (q(v),q(ell))=(144,24)
```

Because `T^3` fixes the visible `v` line and negates the visible `u` line,
degree `12` forces the `u` coefficient to vanish; these differences are
hidden degree-12 maps.  THM-4247 then excludes them.  This is a candidate
new row refinement, pending independent action/count audit and overlap audit.

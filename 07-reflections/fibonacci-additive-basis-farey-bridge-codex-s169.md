# Fibonacci Path-Rank, Additive Bases, And Farey Payloads

The useful correction to the Fibonacci prompt is exact:

```text
1
1
1+1
1+2
1+3+1
1+4+3
1+5+6+1
...
```

is

```text
F_n = sum_k binom(n-k-1,k).
```

For `n >= 2`, the row is the independent-set rank vector of `P_{n-2}`.  So the
line `1+5+6+1=13` is not a decoration on Fibonacci; it is the rank decomposition
of path supports by number of selected non-adjacent sites.

Rebase signal: incoming HYP-2998 already records the broader
Farey-Fibonacci additive-basis carrier on the golden Stern-Brocot spine.  This
reflection keeps the complementary thread: the exact path-rank indexing in the
prompt and the LRC14 unit-excess payload table.

That places Zeckendorf precisely.  Zeckendorf is the normal-form side of this
rank polynomial: all the local path supports exist as combinatorial choices,
but the no-adjacent carry law makes one representation canonical.

## Additive-Basis Spectrum

The earlier S501 bridge becomes sharper if the systems are ordered by proof
currency:

```text
Goldbach:
  many prime-pair branches, local singular-series correction, analytic lower
  bound needed.

Ternary Goldbach:
  one extra prime makes the representation hypergraph thicker; the third
  variable is real smoothing room.

Fermat polygonal:
  bounded arity; `s` copies of `s`-gonal atoms absorb residue obstructions.

Fibonacci / Zeckendorf:
  path-rank carrier plus confluent local carry; entropy collapses to a unique
  normal form.
```

Goldbach wants many witnesses.  Fermat wants enough allowed summands.
Zeckendorf wants a carry system so rigid that all alternate paths reduce to the
same support.  This is a better axis than "additive basis" alone.

## Farey Payloads

The Farey mutation lane says the same thing for `M=p/q`.

```text
p+q  additive payload
p*q  product / incidence payload
q^p  denominator-power magnitude stress
p^q  numerator-power magnitude stress
```

On the LRC14 unit-excess chain `p/(14p-1)`:

```text
q   = 14p - 1
p+q = 15p - 1
p*q = 14p^2 - p
```

Thus `p+q` is linear and affine-safe: it preserves the exact Farey value as
`1+M` after division by `q`.  The product lane is quadratic and should be read
as incidence or area, e.g. the `K_{p,q}` packet in HYP-2932/HYP-2934.  The
power lanes are not proof denominators; they are useful because they make a
false magnitude-forgetting quotient fail loudly.

## New Synthesis

The common object is a representation fiber with a declared allowed quotient:

```text
target N or packet P
atoms / carriers A
representation fiber R_A(P)
local conflict or residue rule
allowed quotient Q
```

Goldbach permits quotienting only after local prime residues and smoothing
errors are controlled.  Fermat polygonal permits quotienting after bounded
arity has absorbed residues.  Zeckendorf permits quotienting after the path
carry is confluent.  Farey payloads permit quotienting only after exact `M=p/q`
and excess `e=14p-q` have been retained.

For LRC, this says: before importing a beautiful sequence, ask which currency
it is offering.

```text
many independent branches    -> sieve / smoothing
bounded number of atoms       -> polygonal / finite packet invoice
path-like local conflicts     -> Zeckendorf / Ostrowski normal form
additive Farey payload        -> p+q, affine-safe
product Farey payload         -> p*q, incidence side channel
power Farey payload           -> stress test for forgotten magnitude
```

The bridge is especially relevant to endpoint debt.  Endpoint repairs often
look like a path of local carry decisions; if two adjacent repairs cannot both
be spent, the Zeckendorf/Ostrowski lens may give a canonical certificate rather
than another scalar.  But when the packet instead has many independent possible
witnesses, Goldbach-style smoothing or Fejer/Ramanujan certificates are the
right economy.

## Computation

`04-computation/fibonacci_additive_basis_farey_bridge_codex_s169.py` stores its
output in `05-knowledge/results/fibonacci_additive_basis_farey_bridge_codex_s169.out`.

Tournament Analysis uses proof currencies as vertices, with pairwise observable

```text
(quotient_safety, local_to_global, entropy_control, LRC_transfer, exact_normal_form).
```

The conservative gauge is transitive:

```text
zeckendorf_normal_form
> fibonacci_path_rank
> farey_sum_affine_check
> fermat_polygonal_bounded_arity
> ternary_goldbach_smoothing
> farey_product_incidence
> binary_goldbach_sieve
> farey_power_stress_test
```

This does not demote Goldbach.  It says Goldbach is a smoothing economy, while
LRC quotient safety often needs normal-form or bounded-packet economies first.

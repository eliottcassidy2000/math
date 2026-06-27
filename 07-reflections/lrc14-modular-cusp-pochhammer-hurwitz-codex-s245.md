# LRC14 Modular Cusp Pochhammer Hurwitz - Codex S245

The useful modular-function analogy is not "q-series are mysterious, therefore
LRC tails are modular."  The useful fact is finite and exact: for the full
modular group, meromorphicity at the cusp says the q-expansion has only
finitely many negative-power terms.  Everything dangerous is concentrated in a
finite principal part before the infinite tail begins.

That is the same controlled-forgetting rule we keep rediscovering around
LRC14.  A quotient may forget a large infinite tail only after it names the
finite debt that tail grows from.  In modular language the debt is the
principal part at `i infinity`; in eta language the tail is q-Pochhammer; in
Hurwitz language the debt is the seed and the legal Vieta mutation word.

The S245 scout makes this concrete.  `(q;q)_infty` has the sparse pentagonal
signed support through `q^32`, while `1/(q;q)_infty` already shows the dense
partition tail.  `Delta=q*(q;q)_infty^24`, so `1/Delta` has principal part
`q^-1`; `j=q^-1+744+196884q+...`, and `j^2` has the finite principal part
`q^-2+1488q^-1`.  The infinite coefficients are not free-floating evidence;
the principal part is the address.

The Markov-Hurwitz orbit says the same thing.  The equation
`w^2+x^2+y^2+z^2=wxyz` has a Vieta mutation grammar
`x_i -> product(other_three)-x_i`; the finite seed `(2,2,2,2)` generates the
visible positive rows in the scout.  A raw equality is not the carrier.  The
seed and legal mutation path are.

For LRC14 this adds a new sidecar family:

```text
modular_cusp_principal_part_order
finite_negative_power_budget
principal_part_coeff_vector
q_pochhammer_tail_signature
eta_delta_denominator_lane
j_rational_function_address
hurwitz_vieta_seed
hurwitz_mutation_depth
continued_fraction_period_word
pell_wall_unit
cusp_tail_discharge_route
```

Tournament Analysis vertices are proof sidecars, not q-coefficients, Hurwitz
quadruples, modular-function names, or runners.  The transitive path is:

```text
labelled_lrc_packet_sheaf
> modular_cusp_principal_part
> full_modular_group_invariance_gate
> q_pochhammer_eta_tail
> hurwitz_vieta_seed_orbit
> continued_fraction_markov_address
> pell_wall_unit_address
> raw_q_series_coefficients
> raw_hurwitz_scalar
```

The challenged assumption is that an infinite series or rare Diophantine
surface can be used as evidence because it looks structured.  The actual proof
object is smaller and stricter: finite principal part, finite seed, finite
address, then tail.

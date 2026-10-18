# A coupled golden contact circuit, with every marked boundary checked

2026-10-17, clock_duality_oct15. **PROVED, scoped/self-audited:** the unit
identities and coordinate independence below. **FINITE-EXACT:** all570
canonical port patterns extend on each of four selected-edge graphs.
**OPEN:** additional cross-edge families, six-colour forcing, and the plane's
chromatic number. The selected graphs are not claimed to contain every unit
edge among their vertices.

## 1. Inheritance and the new operation

The incoming THM-4668,
[the-degree-192-golden-direction-e-is-cartesian-on-the-golden-linear-contact-space](../../01-canon/theorems/THM-4668-the-degree-192-golden-direction-e-is-cartesian-on-the-golden-linear-contact-space.md),
closes linear e extensions over K, not merely the earlier finite grid.
The [incoming additive study](golden_additive_boundary_20261017.md) adds
individual uv and Y=e^2 contacts, while the v^2 extension alone is Cartesian
for the tested coefficient sets. The hostile is the persistent freedom to
shift coefficient colours; the least-used coordinate is a relation between
these products when they coexist.

Let K=Q(zeta15,sqrt(-11)), d=1+zeta5, Z=d+u, and retain the certified cubic
unit u, quadratic unit v, and e from the cited work. Put

    Y=e^2=-Z/conjugate(Z).

The quadratic equation of v is

    v^2+Zv-Y=0,
    v^2=-dv-uv+Y.                                          (1)

Thus a unit previously appearing as a pure new axis becomes a cross-fibre
unit when uv and Y are available together. This does not contradict its
earlier Cartesian census: that census used a different coordinate space.

On the five coordinates(1,u,v,uv,Y), four exact unit families are

    q(d,1,1,0,0),
    q(0,1,1,conjugate(d),0),
    q(d,1,0,0,conjugate(d)),
    q(0,0,-d,-1,1),                  |q|=1 in K.             (2)

The first is the original golden contact. The second is uv times its
complex conjugate. The third equals-Z/(u*conjugate(Z)); the fourth is v^2.
All are therefore units. The exact script also verifies these identities
in the cubic field followed by the quadratic v extension, without decimal
distance tolerances.

The fifth coordinate is not redundant. If a+bu+cY=0 with a,b,c in K,
multiplying by u*conjugate(Z) gives

    (b*conjugate(d)-c)u^2+(a*conjugate(d)+b-cd)u+a=0.

The certified degree3 of u forces a=0, b=cd and c(|d|^2-1)=0.
Since |d|^2!=1, a=b=c=0. Also v is quadratic over K(u), so v and uv
are independent of these three coordinates and of each other. The five
coordinate addresses are therefore injective. No distinct sources have
been spuriously identified in the graph computation.

## 2. Three port-preserving tests

Use the exact eight old ports P, the seven-point spindle M, and coefficient
recipes inherited from `golden_forcing_search_20261017.py`:

* base: A=P union(P+M), C=M union the five golden pentagon ports;
* dM: additionally A contains A+dM;
* MM: the same A, and C=(M+M) union the pentagon.

The selected graph has vertices A+uC+vC+uvC+YC. Its edges are all coordinate
axis unit edges and every translate of the four families(2) whose endpoints
belong to these sets. Unused coordinates remain fixed on a contact edge.
The ports are exactly(p,0,0,0,0), with their previous coordinates unchanged.

| recipe | coordinate sizes | implicit vertices | contact patterns in the four families | blocked port patterns |
|---|---|---:|---|---:|
|base|54,10,10,10,10|540000|56,8,42,8|0/570|
|dM|364,10,10,10,10|3640000|2142,8,454,8|0/570|
|MM|364,29,29,29,29|257450284|72358,148,1744,148|0/570|

For every pattern, SAT supplies five coefficient colourings. The graph
colour is their sum modulo6. Every witness is then checked directly:
proper axis colours, nonzero sums of contact differences, and the literal
port colours. This verifies182400,2172270 and43222530 factor constraints,
respectively. Because each constraint is invariant under the unused
coordinate choices, these checks certify all edges of the implicit graph
without expanding hundreds of millions of vertices. The witnesses are
saved; satisfiability is not accepted on the solver's Boolean answer alone.

**Interpretation.** The product circuit gives genuine new cross-edges but
still fails to restrict any of the570 marked six-colour patterns. This is
a positive extension certificate for the specified graphs. It is not a
six-colour theorem for the full field, nor a complete census in the
five-coordinate space. A new edge family or a coefficient context that
destroys this additive freedom is still needed.

## 3. Transfer to the other lanes

The finite-extraction manuscript's separator postpones sector-specific
reindexing until independent labels exist. Here we instead retain every
coefficient label until (1) certifies a real interaction; quotienting them
early would confuse a Cartesian extension with this circuit. The same
discipline produces the [ROOT boundary module](root_sparse_boundary_module_20261017i.md):
all error-support labels survive before their amplitudes are summed.

Source: authenticated units and fixed ports. Target: coefficient difference
constraints. Map: the injective five-coordinate embedding and additive
colouring. Preserved: actual positions, unit edges, marked colours.
Discarded by the cheap test: unit edges outside the four classified families.
Needed sidecar: an explicit list or classification of those missing edges.
Cheapest next test: find a unit difference outside(2), then challenge the
saved colour witnesses before constructing a much larger graph.

A possible next carrier is the full six-dimensional K-space with basis
(1,u,u^2,v,uv,u^2v). Multiplication by u and v then closes linearly, while
the cubic relation introduces coefficients at19. This is a search proposal,
not a claim that adding those axes forces colours. Test the resulting unit
relations against the saved witnesses first; enlarging a still-free product
source would repeat the obstruction already measured here.

Reproduce, with numpy, sympy, python-sat and numba available:

```bash
python3 04-computation/experiments/golden_contact_circuit_20261017i.py --recipe base --out /tmp/golden-base.json
python3 04-computation/experiments/golden_contact_circuit_20261017i.py --recipe dM --out /tmp/golden-dM.json
python3 04-computation/experiments/golden_contact_circuit_20261017i.py --recipe MM --out /tmp/golden-MM.json
```

Saved receipts and all coefficient-colour witnesses are in
`04-computation/experiments/energy_extraction_20261017i/golden_{base,dM,MM}.json`.

## 4. Reroute through the cubic square and its denominator19

The five-coordinate space already contains u^2, because(1,u,Y) is a basis
of K(u) over K. This supplies an additional unit circuit without enlarging
the field or replacing any port. Put

    m2=-(20+11j)/19, m1=(20-11j)/19, j^2=-11,
    u^3=m2*u^2+m1*u+1,
    a=conjugate(d),
    c=[a*m2+1-a^2*m1+a^3]/(|d|^2-1),
    b=c*d+a*m1-a^2.

Then, exactly,

    u^2=a+b*u+c*Y.                                        (3)

Multiply(3) by1+conjugate(d)u, use Y=-u(d+u)/(1+conjugate(d)u),
and compare coefficients using the cubic equation. This proves the formula;
the script also checks it in the exact field and checks its unit norm.
The coefficient vectors b,c have denominator19 in the declared field basis.
Thus the old small coefficient sets did not intentionally realize this
interaction. Rather than add a sixth abstract coordinate, we build its
actual translations into the existing factors.

Starting with the base recipe, replace the factors on axes0,1,4 by

    A+aM, C+bM, C+cM,

respectively, with old points retained first. M contains0. The other two
factors remain C. The resulting sizes are364,70,10,10,70; the injective
implicit graph has178360000 vertices. Enumerate every translate of each of
the five unit families, including q(a,b,0,0,c), by taking unit differences
in a scaled coefficient factor. This uses exact characteristic-zero norm
tests after the inherited good-prime sieve.

**FINITE-EXACT outcome:** all570 marked patterns again have checked additive
six-colour extensions. The five family counts are3304,56,18326,56,194600.
The axis edge counts are1138,105,15,15,105. Direct verification checks
124100400 factor constraints across the saved witnesses. This is still a
selected-edge graph; no complete five-coordinate unit classification is
asserted. In particular, additive SAT failure would only have identified a
candidate for a full colouring test; the positive witnesses here do certify
actual extensions on every selected edge.

There is a concrete reason to change the next design. Each of C+bM and
C+cM has70 points and105 unit edges: seven disjoint translates of the
10-point,15-edge C graph already account for all of them. Hence no axis
unit edge connects different translated copies. The five contact families
do couple some choices across copies, but the full verified witnesses show
that their combined constraints still leave every port pattern available.
The next cheap target is a **coefficient context that connects these copies
and invalidates the saved colourings**, with each new unit edge checked
before rebuilding the entire product graph. Merely repeating the same
translation enlargement would create more independent choices.

This is the same structural lesson as the clock lift: an extra address is
useful only when the allowed transitions constrain it. A disconnected sheet
label supplies freedom; an authenticated clock increment supplies a relation.
The claim is about these explicit maps and constraints, not an equivalence
between graph colouring and Collatz convergence.

```bash
python3 04-computation/experiments/golden_contact_circuit_20261017i.py --recipe base --cubic --out /tmp/golden-cubic.json
```

Receipt and570 checked witnesses:
`04-computation/experiments/energy_extraction_20261017i/golden_cubic.json`.

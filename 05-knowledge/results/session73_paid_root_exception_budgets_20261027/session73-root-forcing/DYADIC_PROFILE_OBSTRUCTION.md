# Changing a dyadic modulus does not by itself pay a periodic forcing profile

2026-10-27. **PROVED scoped**, with **FINITE-EXACT** controls. Additive to the
frozen core and scalar/spatial supplement. This concerns the original native
`lambda`, not a replacement prior. It does not obstruct coarse observations,
nonperiodic symbolic tails, or an increasing-precision carrier with extra history.

## D1. A dyadic-periodic nonzero ROOT profile has infinite first-stage energy

Use the core's genuine branch phase `sigma(n)=(-1)^v2(3n+1)`,
`b_n=3n/(3n+1)`, `d_n=1-b_n^2`, and native source weights lambda. For a real
or complex profile g on positive odds define

```
Q_lambda(g)=sum_(n>1 odd) lambda(n)
           |g(U(n))-sigma(n)b_n g(n)|^2/d_n.              (1)
```

Let q be any positive power of two. Assume g is periodic **in the odd-source
index** `(n-1)/2` modulo q, equivalently on odd n modulo `M=2q`.
Then

```
Q_lambda(g)<infinity  ==>  g(1)=0.                       (2)
```

In particular no fixed dyadic-periodic profile with ROOT value one has a
finite native phase-aware forcing budget. This is a statement about globally
periodic forcing profiles, not about measuring an arbitrary nonperiodic
profile through residue cells.

**Proof.** Since 3 is a unit modulo M, put `r=-1/3 mod M`, an odd residue.
For every fixed odd base z, its sufficiently remote siblings
`S^j(z)=((3z+1)4^j-1)/3` have residue r. In particular the lambda-positive
ROOT rays `S^(3k+2)(1)` all eventually have that residue. Their native weights
are polynomial in k whereas `1/d` grows exponentially. Equation (1) therefore
forces their residual limit to vanish. Their phase is +1 and their target is
ROOT, so periodicity gives

```
g(r)=g(1).                                               (3)
```

Choose a sufficiently large j with `j=1 mod3`, and let `m=S^j(1)`. This is a
known rooted unit, with `m=2 mod3` and `m=r mod M`. Its minimal inverse base is
`z=(2m-1)/3`. The infinitely many siblings `S^l(z)` divisible by three have
odd valuation `1+2l`, target m, and eventually also residue r modulo M.
Their weights are

```
lambda(S^l z)=4/[(j+l+1)(j+l+2)(j+l+3)]                  (4)
```

on that sibling subsequence, while `1/d` grows like `4^l`. Finite energy thus
forces `g(m)+g(r)=0`. But periodicity gives `g(m)=g(r)`, so `2g(r)=0`.
Over the real or complex numbers this proves (2). All sources used in the
argument have explicit strict ROOT words. No support conjecture is needed. QED.

The existing fixed-source/energy and source-alias notes were checked before
this continuation. Their finite-residue positivity and alias bounds do not
state this periodic forcing-domain obstruction. No literature-priority claim
is made. The proof uses the core's incoming trace mechanism, not a new density
claim for arbitrary Collatz orbits.

## D2. Two named sources already give an exact growing lower certificate

Keep a chosen unit target m from the proof. Choose a remote lambda-positive
ROOT predecessor a and a remote lambda-positive predecessor z of m such that
both have residue r modulo M. Write

```
A=lambda(a)/d_a,  B=lambda(z)/d_z,
D=B(1+b_z)^2,     delta=g(1), x=g(r).
```

Their contribution to (1) is

```
A|delta-b_a x|^2+D|x|^2
 >= [A D/(A b_a^2+D)] |delta|^2.                         (5)
```

This is exact completion of a square; equality in the relaxed two-variable
problem is attained at `x=A b_a delta/(A b_a^2+D)`. When ROOT itself has
residue r, periodicity imposes the further constraint x=delta; (5) remains a
valid lower bound, but its relaxed equality need not be attainable.

Along the explicit sibling families both A and D tend to infinity, so the
coefficient in (5) tends to infinity (it is at least half their minimum).
Thus the no-go has finite exact witnesses against any proposed numerical
energy cap, not only a qualitative limiting argument. Shared targets, source
weights, and the opposite native phases are all retained.

## D3. The dyadic restriction and the one-stage/all-depth distinction are real

Consider the **3-periodic** profile, on odd integers,

```
h(n)=+1 if n=0 or1 mod3;  h(n)=-1 if n=2 mod3.           (6)
```

It has h(1)=1. For every lambda source n, `3|n`, so
`U(n)=sigma(n) mod3`. Hence `h(n)=1`, `h(U(n))=sigma(n)`, and

```
Q_lambda(h)=sum lambda(n)/(6n+1)<=1/19.                  (7)
```

This is an explicit positive control refuting an all-moduli version of (2).
It is not the genuine ROOT sign profile: already at 5 it has value -1 while
the actual shortcut ROOT sign is +1.

It also does not pay all source trajectories. The inherited occupation law
converts the full first-hit path bill for h into

```
Q_occ(h)=sum_(n>1 odd) W(n)
         |h(U(n))-sigma(n)b_n h(n)|^2/d_n.               (8)
```

On `n=S^(3k+1)(1)`, h(n)=-1, h(U(n))=1 and sigma(n)=+1. Therefore its
cost is `6n+1`, with `W(n)=2/[(3k+2)(3k+3)]`. The terms in this subseries
tend to infinity, proving `Q_occ(h)=infinity`.

Thus changing to an odd modulus can pay one original-source stage without
paying the orbit occupation or recovering the native current. A residue
profile, a guarded nonperiodic carrier, and a spatial observer are three
different objects. The changing-modulus extraction theorem remains valid in
its own observation scope; it never supplied a finite-energy periodic forcing
profile or erased the need for a native history/port sidecar.

## Reproduction

```
python3 check_dyadic_profile.py > dyadic.normal.json
python3 -O check_dyadic_profile.py > dyadic.optimized.json
cmp dyadic.normal.json dyadic.optimized.json
```

The exact controls cover source-index moduli q=2^h for h=0,...,16, eight
explicit two-source certificates per modulus, rational completion-of-square
attainers and a declared signal grid. They also check the 3-periodic one-stage
identity on every odd multiple of three through 1001 and thirteen explicit
all-depth divergent ROOT rays. The quantified assertions rest on the proofs.

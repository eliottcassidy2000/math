# Planar Jacobian: the hidden tail, the clock boundary, and the next weight-eight inventory

**Status:** session synthesis for
[THM-4044](../01-canon/theorems/THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall.md)
and
[THM-4045](../01-canon/theorems/THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go.md),
both **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**. The planar Jacobian
conjecture remains **OPEN**. The exact new anchor consequence is only

```text
live reduced (2,3), b=d=0  ==>  total residual weight M>=8.
```

## 1. Inheritance pass

The closest proved mechanism was THM-4012's good-elliptic-factor observer and
its exact maximum-weight-six stable model. The canonical hostile was
THM-4017: a `p*y^2` term can lie below a tempting lower face and delete its
elliptic side component. The corrected near miss was therefore precise:
"singleton highest face" does not imply "all special components rational."
The least-used relevant sidecar was the **complete lower Newton boundary
inventory**, including edge restrictions and singular face points rather than
only the highest weighted face.

This inheritance dictated the anchor test. Under exact maximum weight seven,
enumerate every lower face of

```text
(s^2-p)(1-QH)+gamma*Q*s^2,
```

then ask whether every positive-genus component has a nonzero Hom to the good
`j=0` target. The answer is no: the one hidden elliptic tail has `j=1728`.

## 2. Anchor / Niche / Wildcard portfolio

### Anchor -- exact maximum-seven closure

THM-4045 proves that the complete max-seven Newton subdivision has three
faces:

```text
main:      two rational components, seven transverse nodes, graph rank 6;
tail:      W^2=phi^2 P^4+4 kappa, hence one elliptic E_1728;
vertical:  one rational component.
```

After `Q=sigma^84`, all edge chains and all seven `A_83` node resolutions are
rational. The target has good special fibre `Y^2=X^3+1`, of `j=0`. Since
`Hom(E_1728,E_0)=0`, no special component can own the positive generic map
degree. This excludes exact `M=7` and raises the seam floor to `M>=8`.

The endpoint formerly reserved in THM-4020 is now proved by this different
mechanism. THM-4020 remains an unproved, superseded fourth-row route and is
not a second proof.

### Niche -- conductor and shifted-row information loss

The Russell exceptional-quartic conductor remains orthogonal but useful.
THM-4034 gives the exact degree-178 conductor and identifies the three
cotangent lines lost by multiplication. During concurrent integration,
THM-4043 was promoted to **PROVED + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED**: multiplication by successive powers of the stable target
coordinate shifts the universal row identity and the actual exceptional-
quartic lift reaches `J_6` at all four embeddings. A later concurrent
checkpoint promoted THM-4046: the actual lift reaches `J_7`, meets a sharp
nonzero scalar obstruction at `J_8`, and closes every
`0!=H in t^2 C[t]` over those four Russell folds. Neither theorem gives a
coherent all-order series, algebraization, degree control, another fold or
compiler, a global Keller pair, or a `JC(2)` conclusion.

The conceptual gain is negative but sharp. The conductor shift is nilpotent
of finite length; it is not a periodic clock. The three lost cotangent lines
diagnose multiplication-only information loss, but do not cause the formal
row lift.

### Wildcard -- the sixty-clock as an exact boundary firewall

THM-4044 proves for every `k>=1`

```text
ker O_k=((P^60-1)^k),
ker(O_k on P^2K[P])=P^2(P^60-1)^kK[P].
```

The first invisible pure residual is therefore

```text
Delta_k=P^2(P^60-1)^k,             degree 60k+2.
```

It preserves all depth-`k` clock jets and boundary jets zero and one, but
changes the second boundary Hasse jet `[P^2]`. This is an observer
counterexample, not a second Keller residual. It says exactly why a finite
phase packet cannot prove THM-3997's mandatory `[p^2]R!=0` without a degree
cap or the missing boundary jet.

## 3. Live concept board at close

| concept | exact signal | information lost | effect on the other lanes |
|---|---|---|---|
| complete lower Newton model | all three max-seven faces | terms of weight at least eight | closes `M=7`; makes weight-eight support the next anchor |
| highest weighted face | quick genus/Hom observer | lower boundary tails | demoted from complete invariant to one face in an inventory |
| sixty-clock Hasse observer | exact quotient modulo `(P^60-1)^k` | second jet at omitted boundary `P=0` | explains why clock enrichment cannot replace Newton charts |
| 4D Kakeya spine | Vandermonde rank for cubics | unrestricted residual and boundary placement | lawful bounded interpolation analogy only |
| Sturmian/Fibonacci/triangular phase | exact 60-state address after retaining enough state | height, owner, evaluator, boundary normal direction | supplies labels, never the JC obstruction itself |
| conductor multiplication | degree-178 bounded object, three cotangent losses | local linear directions | positive control for three clock jets; separate nilpotent mechanism |

The board has one common verb: **restriction**. Torus restriction deletes a
boundary jet; highest-face restriction deletes a lower tail; multiplication
deletes cotangent lines; scalar clock shadows delete phase state. The repair
is different in each case, so the similarities cannot be promoted to a
single transfer theorem.

## 4. Typed merged connections

### 4.1 Four-dimensional Kakeya to the planar residual

```text
source:       B(r)=[1:t:t^2:t^3] over F_61;
target:       evaluation functionals on residual polynomials;
map:          coefficient vector -> value at t;
preserved:    four-wise Vandermonde rank and cubic interpolation;
destroyed:    terms beyond the degree cap, P=0 boundary jets, Keller rows;
sidecar:      confluent depth plus an explicit degree bound or boundary jet;
hostile:      P^2 product_i(P-t_i)^k.
```

Thus the finite Kakeya spine is genuinely a broad observer frame. It has no
direct implication for Euclidean Kakeya or for an unbounded Keller residual.

### 4.2 Sturmian/AP, Fibonacci, and triangular period sixty

The AP tail's exact period is `lcm(1,...,6)=60`, coming from persistent
rational owners. The Fibonacci and triangular sequences also admit period-60
state clocks, by different recurrences. A full state may select a phase, but
the scalar pair can collide. Even a recovered phase gives only an address
`P=zeta_60^r`; THM-4044 shows that the full address set still loses a normal
boundary coordinate.

The useful hierarchy is therefore

```text
scalar shadow < full phase state < confluent phase evaluator
              < evaluator plus boundary chart / degree cap.
```

### 4.3 Sun's 2-4-6-8 atoms

The centered binomial atoms have degree at most eight, so one full sixty-clock
does determine them after suitable scalar extension. This is the bounded
positive control for the same observer that fails on the unbounded JC
residual. Sun's finite atoms lose positive height and the centered fixed
point; the JC residual loses a boundary normal jet. The mechanisms should not
be conflated merely because both use the same phase labels.

There is a possible future arithmetic bridge: an Eisenstein norm is an
endomorphism degree on the `j=0` target. It becomes obstructive only if a
stable source component is proved to own a specific isogeny degree. Norm
representability alone supplies neither that component nor its degree.

### 4.4 The hidden-tail repair of the clock lesson

THM-4044 and THM-4045 are dual warnings:

```text
clock observer:       nonzero torus points miss a boundary normal jet;
leading-face observer: highest weight misses a lower boundary component.
```

In each case, adding more of the same kind of observation is not the minimal
repair. The clock needs the second boundary Hasse jet; the weighted model
needs the whole lower subdivision plus edge and singularity sidecars.

## 5. Cheapest decisive weight-eight program

Weight eight is now the exact anchor. The first experiment should not assume
that THM-4012's Bolza highest face owns all positive genus. Instead:

1. Write the complete polynomial of maximum weight eight, retaining both
   top coefficients `[p^4]H` and `[p*y^2]H` and the possibly nonzero
   weight-seven coefficient `[p^2y]H`.
2. Stratify only by exact top-support cancellation: the two singleton
   endpoints, the generic two-term face, and the resonance locus.
3. Enumerate the complete `Q`-lower hull on every stratum and combine
   coincident coefficients before deleting a point.
4. For every face, record normalization genus, Jacobian factors, boundary
   edge squarefreeness, inner-edge chains, and singular-face local equations.
5. Compare the **complete** positive-genus Jacobian product with `E_0`, then
   run attachment or degree-owner tests only after the Hom gate survives.

The canonical hostile is already known: `p*y^2` lies `1/4` below the
max-seven tail plane. The cheapest positive control is the nonresonant
two-term Bolza calculation in THM-4012. The first meaningful success at
weight eight will therefore be a complete positive-genus inventory, even if
it does not yet yield a no-go.

## 6. Scope at close

The session proves a new exact floor on one live seam and a sharp observer
kernel. It does not prove the remaining weight-eight cases, other reduced
cells, a global Keller/Darboux entry theorem, algebraization of the formal
Russell-family lane beyond its proved `J_8` obstruction, the planar Jacobian conjecture, the four-dimensional Kakeya
conjecture, or a new statement about Sun representations. The transferable
research move is narrower and stronger: when a promising observer omits a
boundary, compute the first invisible object and add exactly that boundary
sidecar before drawing a global conclusion.

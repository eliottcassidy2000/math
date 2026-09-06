# THM-4449 result note: seventh rounding and sharp dyadic owner cuts

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED; LRC(14) OPEN.**

For distinct positive odd tails (a=tp,b=tq), with (p<q) coprime, reduce
the half-sum and half-difference to $d,e\in\{-3,\ldots,3\}$ modulo seven.
The two-lift cross-comb has exact mass

```text
mu(Sigma_(a,b))=2/49(1+(e^2-d^2)/(pq)).
```

This recovers the sharp pair caps `4/63` for arbitrary odd tails, uniquely on
the primitive ray `(1,9)`, and `4/77` for odd 3-units, uniquely on `(1,11)`.
For an odd-3-unit triple, the sum of the three pair masses is at most
`124/847`, with projective equality shape `(1,11,121)`. Pair energy is not
physical union mass: exactly

```text
E(T)=mu(F_T)+mu(Omega_T),
```

where `Omega_T` is the mixed one-owner/two-owner locus.

The new sharp physical caps are

```text
all distinct positive odd tails: mu(F_T)<=214/1449,
  equality d(1,9,23), d positive odd;

odd 3-unit tails: mu(F_T)<=72/539,
  equality d(1,7,11), d an odd 3-unit.
```

The proof is all-height. A neutral-edge inequality leaves four exceptional
full-odd rays and five exceptional 3-unit rays. For a presentation
`(tp,tq,c)`, primitive reduction retains `gcd(t,c)=1`, and a one-sided
exposure set `R` gives

```text
mu(F_(tp,tq,c))=sigma(p,q)+2mu(m_t^(-1)R intersect D_c),
|mu(m_t^(-1)R intersect D_c)-mu(R)/7|<=N/(3tc).
```

All exact product cutoffs are nonintegral. Below them the primary audit checks
`1,704` full-odd and `2,801` 3-unit primitive presentations. A clean-room
integer-grid implementation reproduces the equality classifications and
endpoint topology. At `(1,9,23)` the physical/quotient component data are
`(16,10/483)` and `(8,20/483)`; at `(1,7,11)` they are `(12,1/77)` and
`(6,2/77)`.

For a ten-body `C`, these become sufficient clock-two zero-even gates:

```text
mu(G_C)>=214/1449 for arbitrary odd tails;
mu(G_C)>=72/539 for odd 3-unit tails.
```

Equality is included because `G_C` is compact closed and the failure set is
proper open. No universal ten-body mass floor is proved, so this does not
close clock two or LRC(14). The one-even and clock-four retypings still use
the dyadic two-tail seam; the crude 3-unit-pair original-body gates here are
`15/77` and `8/77`.

## Reproduction

```powershell
python -B 04-computation/lrc14_dyadic_seventh_rounding_thm4449.py
python -O -B 04-computation/lrc14_dyadic_seventh_rounding_thm4449.py
python -B 04-computation/lrc14_dyadic_seventh_rounding_thm4449_independent.py
python -O -B 04-computation/lrc14_dyadic_seventh_rounding_thm4449_independent.py
python -B 04-computation/lrc14_dyadic_physical_union_sharp_thm4449.py
python -O -B 04-computation/lrc14_dyadic_physical_union_sharp_thm4449.py
python -B 04-computation/lrc14_dyadic_physical_union_sharp_thm4449_independent.py
python -O -B 04-computation/lrc14_dyadic_physical_union_sharp_thm4449_independent.py
python -B 04-computation/lrc14_dyadic_zero_even_height199_thm4449_independent.py 199
python -O -B 04-computation/lrc14_dyadic_zero_even_height199_thm4449_independent.py 199
```

Every normal/optimized pair is line-identical to its frozen output and ends
in `PASS`. Raw-LF SHA-256:

```text
04-computation/lrc14_dyadic_seventh_rounding_thm4449.py
  88954edd6a39ee544d2548514044f49cb42d8e923df47e4a606bc776f321add3
05-knowledge/results/lrc14_dyadic_seventh_rounding_thm4449.out
  fffa16096b12ec941bd196e1ace42b099f8acfa7fc7537fcb8e42a5fe102acb6
04-computation/lrc14_dyadic_seventh_rounding_thm4449_independent.py
  3ab2e45e69d873aa939d3f396f5a27b3b2af936871799cc25d51d34950c61cc8
05-knowledge/results/lrc14_dyadic_seventh_rounding_thm4449_independent.out
  f3ed9f4bce8ae0d32149dba5ff7dd7354bb1cc5dfe2f358c2c4b19aad8be002f
04-computation/lrc14_dyadic_physical_union_sharp_thm4449.py
  eaceaeb497a3327ecfcfb9df5446c6e2f22934546497363c0387e8c003ff77fb
05-knowledge/results/lrc14_dyadic_physical_union_sharp_thm4449.out
  780724feeb66641f6a165ca2c47c424895c38ad31ef4866d61b3ccef503dc15f
04-computation/lrc14_dyadic_physical_union_sharp_thm4449_independent.py
  aaad9c327e2246f96da5e2b7eb317ad524e26f1e2c2c5975376256aed9a053e5
05-knowledge/results/lrc14_dyadic_physical_union_sharp_thm4449_independent.out
  b5c495fc57ca2523d6bf774a0798db1964cd2a130767cc3497cf04d011901d05
04-computation/lrc14_dyadic_zero_even_height199_thm4449_independent.py
  aa87478cc19d34d3d578624cf6981598079799615df67bee480138bde37965b9
05-knowledge/results/lrc14_dyadic_zero_even_height199_thm4449_independent.out
  599fc549d39e2b692ee985aac3185e991efaa1ba194f99be35ad45d07c604fe8
```

The height-199 census is retained as discovery provenance and a bounded
control; the all-height theorem rests on the BV product cutoff and complete
finite exposure boxes.

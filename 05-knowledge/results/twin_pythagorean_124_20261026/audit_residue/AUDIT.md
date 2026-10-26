# Independent audit of the fixed-residue obstruction

**Verdict: PASS.** Reviewed the exact CRT proof, its scope exclusions, the
named characteristic-seven witness, and the producer checker. No producer
file was modified. This audit does not extend the statement to arbitrary
finite-valued functions or growing moduli.

## Frozen producer

Reviewed under `twin_pythagorean_124_20261026/residue_obstruction/`:

| File | SHA-256 |
|---|---|
| `README.md` | `b2149d1b4a2b0b14f89f22c34b6907b976df8e00102c3a6fe981b5aa689e1f42` |
| `check.py` | `ba451fbe8b4cf4aa4ca7146a6b9a4d6d9d2e49d68056c81e3c028224eeccf7b9` |

The producer replay passes **23,044 gates** in both normal and optimized
Python. Both outputs have SHA-256
`769d1a25b07b69355b6bf258c6cddf07dfe1bd56fea76839c45b47eaec80ae92`.

## Proof audit

1. The moduli `2M,q,r` are pairwise coprime when the distinct odd primes
   `q,r` avoid `2M`. Such primes exist for every fixed positive `M`.
   The CRT conditions therefore define one infinite arithmetic progression.
   They are compatible even if the initial prime `p` equals one of `q,r`.
2. Choosing a sufficiently large positive representative makes the exhibited
   prime divisors proper. The exact source property is that both `x` and
   `x+2` are composite. The construction is stronger than simply showing
   that some integer in the residue class is nonprime.
3. The congruence modulo `2M` ensures odd parity and preserves the D-ray
   depth `(x-1)/2` modulo `M`. The Euclid parameters `(x+1,1)` are coprime
   and of opposite parity, so the resulting triple stays primitive despite
   compositeness. Polynomial reduction preserves all triple coordinates.
4. Direct area formulas give inradius `x` and exradius opposite the even
   leg `x+2`. The names of these two radii are correctly distinguished.
5. The D parameter matrix has determinant one, as does its triple action.
   Every D-ray branch word therefore has the same determinant sign. That
   sign cannot distinguish the constructed points along this ray.
6. A finite collection of residue moduli is handled by their least common
   multiple. The statement correctly excludes general finite-valued
   observers that do not factor through a fixed modulus; it does not
   assert that a primality bit itself cannot distinguish primality.

The witness `213=3*71`, `215=5*43` has the same residue as `3` modulo `14`.
The triples `(45795,428,45797)` and `(15,8,17)` both reduce to `(1,1,3)`
modulo seven, and both residual parameter slopes `r=n/m` are `1/4=2` there.
The finite-field note's reciprocal coordinate is `t=m/n=4`. These are two
charts of the same point, not different residues in a common chart. A fixed function
of that marked conic point necessarily gives the same finite output.
The text correctly calls any further F8 correspondence a finite-set
marking, not a field embedding between different characteristics.

The initial reviewed proof hash was
`7bdcb4cfa281f55743cdfd9bb81e762cba40d548e5c4fa31d994cca4ac42c407`.
The final hash above includes the author's explicit `r=n/m`, `t=m/n`
notation clarification. The reciprocal values were independently checked;
the CRT proof, witness, and producer code are unchanged.

## Independent exact check

`check_independent.py` does not import the producer. Instead of modular
inverses, it literally searches one full period to locate the unique CRT
representative for each modulus `1..40` and twin lower prime at most `101`.
It then checks three positive lifts per class. Radii are independently
checked as area divided by the appropriate semiperimeter denominator.

This passes **7,923 gates**, including the literal named witness. Normal
and optimized outputs agree byte for byte, with SHA-256
`9d5407e5b355f352484e5c04fc1463a42e2e8240eb6af47ef3a9ce27f4c4ed64`.

Reproduce this audit from its directory:

```bash
python3 -B check_independent.py
python3 -B -O check_independent.py
```

Finite checks validate the implementation; the all-modulus and infinitude
claims are justified by the displayed CRT proof. No claim about twin-prime
infinitude, unbounded sieving, or a global tree quotient is inferred.

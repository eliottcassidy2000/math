# Lucas47: recovered sources, a candidate triad, and an exact global intersection

Session127. Source repository pin: `522535a7b6f519fd0c1e4304a2c5b1c909fce4c9`.

**Status:** the user's exact intended phrase **“Lucas47 triad” remains UNRESOLVED**.
The construction below is a **new candidate interpretation**, not a claim to have
recovered that intent. Elementary identities are **PROVED**; the global
classification uses a **CITED classical theorem plus a proved reduction**;
normal/optimized receipts are **FINITE-EXACT controls**, not the global proof.
No external-priority claim.

## 1. Precisely recovered inheritance

The convention is `L_0=2, L_1=1, L_(n+2)=L_(n+1)+L_n`; hence `L_8=47`.
`47` here is a **value**, not the index47 (`L_47=6643838879`).

- `04-computation/fibonacci_lucas_primes_s116l.py` explicitly states this convention.
- `05-knowledge/results/collatz_circulant_20260930_circulants_lucas_cubic_monotile.md`,
  section2, records `L_8=47=2^7−3^4`; its old Pillai census is bounded by `K<=400`.
- `05-knowledge/results/collatz_lucas_monotile_discrepancy_20260930.md`, sections1–2,
  distinguishes Lucas traces from unit-minus-one norms and records the
  trace-doubling chain `3,7,47,2207`.
- `05-knowledge/results/session93_reversal_lucas_clocks_20261029/arithmetic/TRACE_CLOCK_CARRIERS.md`
  retains the precise trace/norm observer and its lost conjugate coordinate.
- `05-knowledge/results/free_exploration_s90g.out`, connection6, and
  `05-knowledge/results/rational_irrational_symbiosis.out`, chapter2, recover
  the tribonacci continued fraction and predecessor `103/56` before `305`.
  Their corresponding `.py` files contain no47/Lucas/triad linkage.

Targeted filename and GitHub searches did not find a named Lucas47 triad;
saved-conversation and room-message searches returned no Lucas/triad matches.
These are bounded retrieval observations, not proof of absence everywhere.

## 2. The three candidate source objects and maps

1. **Golden unit:** `u=phi^8=13+21phi`, `phi²=phi+1`.
   The trace map gives `Tr(u)=47`, while `N(u)=1`.
   Equivalently `47²−5·21²=4` and `47=7²−2`.
2. **Pillai clock:** the exponent pair `(A,B)=(7,4)` maps to
   `2^A−3^B=128−81=47`.
3. **Tribonacci unit:** with `tau³=tau²+tau+1`, the coefficient pair
   `(p,q)=(103,56)` maps under `r=p−q` to47 and satisfies
   `N(p−q tau)=p³−p²q−pq²−q³=−1`.

These are three **observers with a common integer output**, not isomorphic
source objects. Each map forgets important data: the conjugate branch, the
exponents, and the second coefficientq respectively. Their joint fibre must
retain those source sidecars. A common integer does not supply a homomorphism
between the quadratic and cubic fields: such a Q-field embedding would force
2 to divide3, which is impossible.

A critical hostile: `N(phi^8−1)=2−L_8=−45`, **not−47**. Thus the47-word Lucas
cube, the Lucas trace47, and the golden torsion quotient of order45 cannot be
identified by cardinality. The 47-state reflection-group construction is a
separate carrier and was not imported into this candidate triad.

## 3. The genuine morphism: integral coordinate shear

Put `F_3(p,q)=p³−p²q−pq²−q³`. The unimodular map

    (p,q) -> (r,q)=(p−q,q),    (r,q) -> (r+q,q)

is a change of coordinates on **the same cubic norm form**. It preserves the
norm predicate exactly because

    F_3(r+q,q)=r³+2r²q−2q³.

At the305 predecessor it gives

    47³+2·47²·56−2·56³=−1.

The other exact collapse is

    (p−q) F_3(p,q)=q⁴−(2q−p)p³,
    56⁴−9·103³=−47.

The latter is the k=3 specialization of the uniform telescoping identity
`(p−q)F_k(p,q)=q^(k+1)−(2q−p)p^k`. It exposes47 as the multiplier of the
unit norm, not an unexplained decimal coincidence. This uniform identity
belongs in the common k=3,4,5 packet; the Lucas property is an additional
constraint, not an automatic consequence at other orders.

The sparse unit identity is checked independently in the quotient ring:

    tau^16=1705+2632tau+3136tau²,
    (56tau−103)tau^16=1.

The305 lane supplies its exact continued-fraction interpretation; this packet
does not substitute rounded powers for that certificate.

## 4. Global Lucas intersection (CITED input + PROVED deduction)

**Theorem.** For integers `n>=0, q>0, p−q=L_n`,

    |F_3(p,q)|=1

holds exactly for

    (n,p,q,F_3)=(1,2,1,+1), (8,103,56,−1).

**Classical input.** The Delone–Nagell bound as stated in Lemma2.4 of
[Knapp–Lemos–Neumann, *Integral Values of Generating Functions of Recurrence Sequences*](https://math.loyola.edu/~mknapp/papers/General%20Quadratic%20-%202-10-23.pdf)
says an irreducible integral binary cubic of negative discriminant has at
most five solutions to form=1; five requires discriminant−23. Their Example4.1
applies this to the tribonacci form. This primary paper was read at printed
pages7 and10–11; its variable order is reversed relative to ours.

**Deduction.** `F_3(X,1)` has no rational root (the only candidates±1 both
fail), so is irreducible. The cubic discriminant formula gives−44. Thus there
are at most four integral solutions to `F_3=1`. Direct substitution verifies
four distinct solutions

    (p,q)=(0,−1),(1,0),(2,1),(−103,−56).

Odd homogeneity makes all solutions to `F_3=−1` their negatives. Filtering
`q>0,p−q>0` leaves only `(2,1)` and `(103,56)`; their differences are1 and47.
The Lucas convention gives those atn=1 and8 only: `L_0=2`, and the sequence
is strictly increasing fromn=1. This proves the stated global intersection.

The strength comes from the **cited global cubic theorem**, not the Lucas
index100 census. No global search over Pell or unit orbits is silently assumed.

## 5. Independent finite controls and reproduction

`verify.py` uses integer Lucas recurrence and exact bisection.
`independent.py` does not import it: it uses a binomial Lucas formula,
permutation determinants, and a high-precision candidate subsequently enclosed
by **exact adjacent integer sign checks**. Decimal arithmetic proposes only;
no theorem depends on its accuracy. Both obtain identical complete row digests.

For each of101 Lucas indices0..100, all positive denominators are covered:
for `r>0`, write `g_r(q)=r³+2r²q−2q³`. On `0<=q<=r`, `g_r(q)>=r³`;
for `q>=r` its derivative `2r²−6q²` is strictly negative, and
`g_r(2r)=−11r³`. Therefore only the two adjacent integers surrounding the
unique positive root can have values±1, with the boundary `r=q=1` retained.
This gives a complete denominator search per index without an arbitraryq cap.

Hostiles include the neighboring Lucas termsn=7,9, the next doublingn=16,
and trace47 versus torsion norm45. The scripts also check the four classical
witnesses and their signed/filter conventions; they do not computationally
reprove the cited theorem.

```bash
# From any working directory, change to the directory containing this README.
cd path/to/lucas47
python3 verify.py > normal.json
python3 -O verify.py > optimized.json
cmp normal.json optimized.json
python3 independent.py > independent.normal.json
python3 -O independent.py > independent.optimized.json
cmp independent.normal.json independent.optimized.json
python3 reproduce.py   # copies scripts and reruns normal/-O from an unrelated cwd
sha256sum -c SHA256SUMS
```

`SOURCE_MANIFEST.json` pins the inherited repository files by commit, Git blob ID,
SHA-256 and byte length. `relocation.json` records successful copied-script
normal/optimized replay from a separate working directory. No absolute source
checkout path is a runtime dependency.

The original finite result remains recorded separately in the receipts. The
exact intended triad and a specific GSMA source/morphism definition remain
external identification questions; none is resolved merely by these equations.

---
id: THM-1263
title: The mod-23 near-bijection rung for twelve speeds
status: PROVED.  For twelve integer speeds none divisible by 23, a global loneliness margin below 2/23 forces their nonzero residues to cover all eleven antipodal pairs of Z/23; exactly one pair has multiplicity two and the other ten have multiplicity one.  This applies in particular to every hypothetical member of the n=12 first gap M<2/25.  The divisible-by-23 branch is not covered and genuinely need not satisfy the near-bijection conclusion
source: codex-2026-07-19 mod-23 near-bijection rung
depends_on: []
related: [THM-518, THM-633, THM-1185, THM-1220]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCMod23NearBijection.lean
referee: 04-computation/lrc14_mod23_near_bijection_thm1263_referee.py
output: 05-knowledge/results/lrc14_mod23_near_bijection_thm1263_referee.out
---

# THM-1263 -- the mod-23 near-bijection rung

For a finite integer speed family `V=(v_i)`, write

```text
L_V(t)=min_i dist(v_i t,Z),
M(V)=sup_t L_V(t).                                    (1)
```

> **Theorem.**  Let `V=(v_1,...,v_12)` be twelve integer speeds satisfying
> `23 ∤ v_i` for every `i`.  If
>
> ```text
> M(V)<2/23,                                           (2)
> ```
>
> then the residues of the speeds cover every nonzero antipodal pair
> `{u,-u}` in `Z/23`.  On the eleven-pair quotient
>
> ```text
> (Z/23 \ {0})/{u~-u},                                (3)
> ```
>
> exactly one pair contains two speeds and each of the other ten contains
> one speed.

Since

```text
1/13 < 2/25 < 2/23,                                  (4)
```

the conclusion applies to every hypothetical `n=12` first-gap family with
`M(V)<2/25`, provided no speed is divisible by `23`.

## 1. The middle-band argument

Fix an integer `b` with `23 ∤ b` and inspect the rational time `t=b/23`.
For every integer `v,m`, writing `r=(vb) mod 23` gives

```text
vb-23m=23(vb div 23-m)+r.                            (5)
```

If `2<=r<=21`, then the right side has absolute value at least `2`, so

```text
|v(b/23)-m|>=2/23                                    (6)
```

for every integer `m`.  Thus if all twelve residues lay in the middle band
`{2,...,21}`, time `b/23` would witness `M(V)>=2/23`, contradicting (2).

The residue `0` cannot occur: primality of `23`, together with
`23 ∤ v_i` and `23 ∤ b`, excludes `23 | v_i b`.  Hence for every unit scale
`b` there is an `i` such that

```text
v_i b = 1 or -1  (mod 23).                            (7)
```

Now fix any nonzero `u in Z/23` and choose `b=u^-1`.  Multiplying (7) by
`u` gives

```text
v_i = u or -u  (mod 23).                              (8)
```

Therefore all eleven antipodal pairs are hit.

## 2. Why coverage sharpens to a near-bijection

Every nonzero residue belongs to a unique pair represented canonically by
`1,...,11`.  Assign each of the twelve speeds its pair label.  Equation (8)
says this map

```text
Fin 12 -> Fin 11                                     (9)
```

is surjective.  If `n_j` is the size of pair fiber `j`, then

```text
n_j>=1 for all j,        sum_(j=1)^11 n_j=12.         (10)
```

Subtracting the eleven compulsory units leaves total excess `1`.  Hence one
fiber has size `2` and all other fibers have size `1`.  This is stronger than
the mod-19 spread constraint: at modulus `23`, twelve runners have only one
unit of multiplicity slack.

## 3. MISTAKE-186 audit and positive control

The formal closeness premise used in the contrapositive is

```text
forall b, exists i, exists m:Z,
  |v_i(b/23)-m|<2/23.                                 (11)
```

It is deliberately **not** `exists i, forall m`.  The latter would require
one real number to be close to every integer and is unsatisfiable; that was
the vacuity caught in MISTAKE-186.  A global bound (2) implies (11) at every
rational time `b/23`.

The exact positive control is

```text
V0={1,2,...,12},       M(V0)=1/13.                    (12)
```

For all `b=0,...,22`, the referee exhibits an actual pair `(i,m)` satisfying
(11); periodicity in `b modulo 23` handles every integer `b`.  Its pair
multiplicities are

```text
(1,1,1,1,1,1,1,1,1,1,2),                            (13)
```

with the double fiber `{11,-11}={11,12}`.  Thus the formal hypotheses are
nonvacuous.  The referee computes the displayed exact margins by using the
piecewise-linear breakpoint fact for `min_i ||v_i t||`: a maximum occurs at a
denominator among `2v_i`, `v_i+v_j`, or `|v_i-v_j|`, after which all candidate
numerators are checked with rational arithmetic.

## 4. The divisible-by-23 branch is genuinely different

If `23 | v_i`, then at every time `b/23` that runner is exactly integral.
The rational-time test is automatically blocked and supplies no nonzero
residue information.  This is not a removable technicality.  The exact
control

```text
Vblock={1,2,...,11,23},       M(Vblock)=1/12<2/23     (14)
```

has one zero residue and hits each of the eleven nonzero antipodal pairs
exactly once.  It has no double nonzero pair, so the conclusion is false if
the hypothesis `23 ∤ v_i` is deleted.  The theorem therefore yields an honest
dichotomy:

```text
some 23|v_i,  or  the residues form the near-bijection (13).  (15)
```

It does not classify or eliminate the divisible branch.

## 5. Alternate carriers and tournament loss audit

We challenged runners, signed residues, scales `b`, antipodal pairs, missing
pair obligations, and CRT constraints as vertices.  Runners obscure the
quotient that the argument preserves.  The natural vertices here are the
eleven antipodal pair obligations; the runner-to-pair map (9), not a
tournament, carries the proof.

Orienting the pair vertices by `(multiplicity,canonical label)` gives a
transitive tournament with score histogram `0,...,10`, no directed
3-cycles, singleton SCCs, and one Hamiltonian path.  It preserves missing-pair
and multiplicity data but forgets runner identities and residue signs.  It is
therefore only a loss audit.  For the positive control its path is
`1,2,...,11`, with the doubled pair last.

## 6. Formalization, hashes, and frontier

The Lean module proves the integer middle-band bound, the corrected
existential contrapositive, per-scale `+-1` incidence, antipodal coverage,
unique canonical pair assignment, the `12 -> 11` fiber-count theorem, and the
margin-form near-bijection package.  It also proves integer-periodicity and
the full global form with
`M(V)=sSup(margin V '' [0,1])`, reducing arbitrary `b/23` into the compact
period before applying the supremum bound.  It is sorry-free and uses no
`native_decide`.

The standard-library referee independently computes both exact margins,
checks every scale class, exhausts all eleven positive compositions of `12`
into eleven parts, verifies the blocker control, and passes identically under
ordinary and optimized Python.

```text
referee_sha256 = 704522df515038631db6a0442c557e31844c9c214eaad2f01860b2e9600a26df
output_sha256  = 73a7686533e12a914786d7f0d417c8e19fb5854182435b5faa1a111badca3917
lean_sha256    = 75309aa986c0669afa6eb63257fb57524cfa8fc879f471df179412a63203d647
```

This theorem does not prove the `n=12` uniqueness gap or LRC(14).  Its
high-leverage use is simultaneous rather than isolated: combine the single
unit of mod-23 slack with the already proved mod-17/mod-19 spread conditions
and the covering constraints at moduli `2,...,13`.  That finite cross-modulus
compatibility problem can constrain which speed is allowed to serve the
unique doubled pair and which speeds must block other certificate rungs.

# Independent audit of the effective-clock gcd reduction

**Status: AUDIT PASS; PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN
RUNNERS, with FINITE-EXACT controls.** The
[main proof](lrc14_effective_clock_empty_core_sep06.md) correctly shows
that every primitive thirteen-speed strict counterexample at threshold
`1/14` must have every twelve-subset of gcd one, every eleven-subset of
gcd at most two, and every ten-subset of gcd at most four. This is a
necessary condition on a hypothetical failure, not an example or a proof
of LRC(14).

I read the exact statements and proofs of
[THM-765(B), quantitative deck obstruction](../../01-canon/theorems/THM-765-safe-component-tooth-deck-and-hereditary-primitivity.md),
[THM-761, multi-exception sheets](../../01-canon/theorems/THM-761-multi-exception-sheet-covering-bound.md),
and [THM-4004, typed divisor-comb profile](../../01-canon/theorems/THM-4004-lrc14-three-detuned-divisor-comb-profile.md).
After consulting `CORE-PAPERS.md`, I also checked the current primary
[Theorem 1.3 of Sungkawichai--Trakulthongchai](https://arxiv.org/html/2604.23906v2)
and its definition of `LRC(k)`: twelve nonzero speeds have clearance at
least `1/13`. This external lower-runner input is cited; its computation
is not reverified here.

The exact orbit of an exception `w` on the body's `c` branches has
`q=c/gcd(c,w)` distinct points, each repeated `c/q` times. An open arc
of length `1/7` contains at most `ceil(q/7)` points, including the case
`7|q`. Thus `beta(q)=ceil(q/7)/q` is the correct normalized budget.
The sum bound is sufficient for a free branch only when it is strictly
below one; equality is correctly retained throughout the proof.

For a twelve-subset, THM-765(B) applies with `L=1/14`, since the cited
core margin is strictly larger and `L<1/4`. Its gcd must divide the
missing speed, and full primitivity makes that gcd one. For an
eleven-subset of gcd `c`, adjoining either missing speed gives gcd one,
so both exceptions are individually coprime to `c`. The necessary
budget `2beta(c)>=1` allows only `c=2` when `c>1`.

For a ten-subset of gcd `c`, each complementary exception satisfies
`gcd(c,w_i)<=2`, hence `q_i>=c/2`. Adjoining any two exceptions gives a
twelve-subset, so `lcm(q_i,q_j)=c` for every pair. If `c>=5`, every
order is at least three. Since `beta(q)<=1/3` for `q>=3`, with equality
only at three, a total budget at least one forces all three orders to
be three. The pair-lcm identity then gives `c=3`, a contradiction.
This establishes the all-height conclusion without extrapolating any
finite order or speed universe.

The complete residual profiles are correctly listed as
`(c;q)=(2;1,2,2),(2;2,2,2),(3;3,3,3),(4;2,4,4)`.
Order one cannot be discarded at clock two. Absorbing its even tail,
or the uniquely even tail of the clock-four profile, produces the
stated eleven-body clock-two problem. The common-dilation corollary
`gcd(P)>4gcd(V) => M(V)>=1/14` is equivalent to the primitive statement.
Neither the body height nor the phase-address complexity is bounded.

I independently reconstructed the THM-737 hostile measures before the
repair: `C={1,...,12}`, `c=4`, `D={26}` has true effective order two and
measure `6617/388080`, exactly half the body's `6617/194040`. It is below
the formerly unqualified coprime floor `6617/258720`. The false step is
replacing the effective order by `c` for a nonmultiple exception; the
gcd-aware measure inequality itself remains sound.

The [standalone verifier](../../04-computation/lrc14_effective_clock_empty_core_sep06.py)
was read in full and replayed independently. Its danger-union and
endpoint-cell engines correctly treat measure-zero boundaries, its arc
probes cover every event and intervening cell, and its divisor-word
universe retains the exact hereditary filters. Its phase-only hostile
kills both selected sheets while keeping the chosen body phase safe;
it does not assert an unsafe full speed set.

The independent replay is byte-identical to the
[saved output](lrc14_effective_clock_empty_core_sep06.out): 147 explicit
gates; semantic digest
`f7515199e28d2576094e045cd5f9e4b02a4eae95192dd523d81f8fb2f0681c81`.
Recomputed raw source/output SHA-256:

```text
12fde59aab96da6a312b358d058b7261834a4515a5497a64924a9ffa71aabecd
f734022c70d8632417c1a60f8479243a7ad579c8a1b16dbaa501911d8b67e9b5
```

```sh
python3 -B 04-computation/lrc14_effective_clock_empty_core_sep06.py
```

The next decisive object is a joint phase-labelled component union in
the remaining clock-two and clock-three bodies. Marginal branch counts
alone discard the intersection pattern needed to finish those cases.

Final acceptance: the corrected THM-737 prose has now been checked, including the positive body-measure premise, the explicit `1<=d<=6` scalar-cutoff domain (and `2<=d<=6` in its multi-exception scope), the exact `q=c<=7` special-factor boundary, and the `q` distinct orbit points with their retained gcd multiplicity; no remaining blocker was found in this repair.

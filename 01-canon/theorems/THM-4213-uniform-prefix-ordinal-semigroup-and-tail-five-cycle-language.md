---
id: THM-4213
title: "Uniform-prefix ordinal semigroup and tail-five cycle language"
status: >
  PROVED exact contextual-prefix telescoping and the weighted composition law
  for Hamilton-normalized certified floors; the nonnegative prefixes form an
  ordinal semigroup and the uniformly positive prefixes a two-sided ideal;
  every tail-five-separated cycle word with m cycle blocks has floor
  10764(9^m-1)/8, with equality only for the original A5 singleton row;
  every such word with final transitive run at least six is an (OS+) left
  factor; and C3 triangleright C3 triangleright P5 is an exact good-suffix
  hostile with singleton defect -338580 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. THM-4215 subsequently closes the two-adjacent-cycle threshold.
  General (OS+), classification of all uniformly positive prefixes,
  thresholds for three or more consecutive cycle components, the
  no-sink/no-source gate law, and the order-eleven asymmetric bank remain OPEN.
source: codex-tournament-osplus-extension-20260826
depends_on:
  - THM-4184-path-cover-parity-ordinal-cocycle-and-lollipop-positivity
  - THM-4187-universal-source-padding-supermodularity-and-transitive-left-ordinal-positivity
  - THM-4212-cycle-prefix-uniform-arbitrary-context-threshold-and-tail-five-lower-bound
related:
  - THM-4193-cycle-first-transitive-tail-crossing-and-transitive-context-positivity
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
  - THM-4215-adjacent-cycle-sharp-uniform-context-threshold-seven
script: 04-computation/tournament_uniform_prefix_semigroup_thm4213.py
output: 05-knowledge/results/tournament_uniform_prefix_semigroup_thm4213.out
independent_audit_script: 04-computation/tournament_uniform_prefix_semigroup_independent_audit_thm4213.cpp
independent_audit_output: 05-knowledge/results/tournament_uniform_prefix_semigroup_independent_audit_thm4213.out
script_sha256: 5556b0e336c006e034779b960fc79c966baa61fd3a4b2a96fa21840d98522e5a
output_sha256: 0e2e8393412453c1da67260ee35180ec6bdff3e98eb4483960811777d6b38fab
independent_audit_script_sha256: 45638c2d2e49902cfc0a7027f453ce8db0c26f70fb54a222fe6cc3c1e7c7009f
independent_audit_output_sha256: f4418df663be94ce9049341d21c4306726a53b9b16c5ef230734cdf664692cdc
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact symbolic subtraction proves the contextual telescope and its
  weighted floor law. The inherited exact ordinal-transfer engine checks
  1,936 telescope rows, 363 tail-five floor rows, 121 two-cycle rows, 121
  transitively padded rows, two no-sink OS+ controls, the equality boundary,
  and the exact -338580 hostile. Normal and python -O streams byte-match.
independent_audit: >
  ACCEPT. A standalone warning-clean C++17 referee imports no response jet or
  ordinal-capacity transfer. From literal labelled adjacency and subset path
  DP it rebuilds Hamilton counts, exposed-word capacities, G_+, R_+, nine Q5
  floor rows with their unique equality, nine telescope controls, the hostile
  and one OS+ padding control. Clang O0/O3 and ASan/UBSan streams byte-match.
---

# THM-4213 -- uniform-prefix ordinal semigroup and tail-five cycle language

**PROVED exact ordinal semigroup law + quantitative multi-cycle language +
exact good-suffix hostile + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4212 proves a sharp uniform lower bound for the single-cycle prefixes

```text
Q_n=C3 triangleright P_n,                 n>=5.
```

This theorem turns that one family into an ordinal language with arbitrarily
many nontransitive strong components. The mechanism is an exact telescope:
contextual prefix gains add under ordinal composition, but the gain of an
earlier factor is evaluated on the whole suffix. Hamilton multiplicativity
then weights every positive block by the square of the suffix Hamilton count.

The same mechanism locates its own boundary. Uniformly positive prefixes are
a two-sided ideal only inside the semigroup of nonnegative prefixes. An
arbitrary left factor can destroy positivity even when the suffix is `Q_5`.

## 1. Contextual defects and certified floors

All tournaments quantified below are finite and nonempty. Use THM-4184's
ordinal sum `triangleright`, Hamilton count `H`, and ordinal remainder `R_+`.
For a prefix `A`, put

```text
F_A(B,C)=R_+(A triangleright B,C)-R_+(B,C).             (1)
```

For `alpha>=0`, say that `A` has certified floor `alpha`, written here as
`A in C_alpha`, when

```text
F_A(B,C)>=alpha H(B)^2H(C)^2                            (2)
```

for every nonempty `B,C`. This is a lower-bound certificate, not a claim that
`alpha` is optimal. Define

```text
N={A:F_A(B,C)>=0 for all nonempty B,C},
U={A:A in C_alpha for some alpha>0}.                    (3)
```

Thus `U` is the uniformly positive cone and `U` is contained in `N`.

## 2. Exact telescope and the ordinal ideal

> **Theorem 1 (contextual-prefix telescope).** For all nonempty tournaments
> `A,D,B,C`,
>
> ```text
> F_(A triangleright D)(B,C)
>   =F_A(D triangleright B,C)+F_D(B,C).                 (4)
> ```
>
> Consequently, if `A in C_alpha` and `D in C_beta`, then
>
> ```text
> A triangleright D in C_(alpha H(D)^2+beta).           (5)
> ```

### Proof

Associativity of ordinal sum and addition/subtraction of the middle remainder
give

```text
R_+(A triangleright D triangleright B,C)-R_+(B,C)
 =[R_+(A triangleright(D triangleright B),C)
    -R_+(D triangleright B,C)]
  +[R_+(D triangleright B,C)-R_+(B,C)],                 (6)
```

which is `(4)`. Hamilton paths factor at an ordinal cut, so

```text
H(D triangleright B)=H(D)H(B).                          (7)
```

Apply the two certified floors to the terms of `(4)` and use `(7)`:

```text
F_(A triangleright D)(B,C)
 >=[alpha H(D)^2+beta]H(B)^2H(C)^2.                    (8)
```

This proves `(5)`. QED.

> **Corollary 2 (semigroup and ideal).** `N` is closed under ordinal sum,
> and `U` is a two-sided ideal of `N`:
>
> ```text
> N triangleright N subset N,
> N triangleright U subset U,
> U triangleright N subset U.                           (9)
> ```

Indeed `(5)` with `(alpha,beta)=(0,0)`, `(0,beta)`, or
`(alpha,0)` proves the three inclusions. If the empty tournament `P_0` is
adjoined only as a formal identity, `N` becomes a monoid and `U` remains its
two-sided ideal. No empty tournament is used in `(1)--(8)`.

THM-4187 supplies the first large neutral subsemigroup. Since `H(P_a)=1`, its
local interaction theorem says

```text
F_(P_a)(B,C)=Theta(P_a,B,C)>=0,             a>=1,       (10)
```

so every transitive tournament `P_a` belongs to `N`.

## 3. The tail-five-separated cycle language

Put

```text
q=10764,
Q_n=C3 triangleright P_n,                    n>=5.       (11)
```

THM-4212 proves

```text
Q_n in C_q,                 H(Q_n)=3,                    (12)
```

with equality in the floor only for `n=5` and `B=C=P_1`.

Let `m>=1`, let `n_1,...,n_m>=5`, and let `a_0,...,a_m>=0`.
Here `P_0` means that the corresponding factor is omitted. Consider

```text
X=P_(a_0) triangleright Q_(n_1) triangleright P_(a_1)
    triangleright ... triangleright Q_(n_m) triangleright P_(a_m). (13)
```

> **Theorem 3 (quantitative multi-cycle language).** Every word `(13)` and
> every nonempty `B,C` satisfy
>
> ```text
> F_X(B,C)
>  >=q(1+9+...+9^(m-1))H(B)^2H(C)^2
>   =10764(9^m-1)/8 H(B)^2H(C)^2>0.                    (14)
> ```
>
> Equality holds in `(14)` if and only if
>
> ```text
> X=Q_5,                    B=C=P_1.                    (15)
> ```

### Proof

Iterating `(4)` expands `F_X(B,C)` into one term for every nonempty factor of
`(13)`. Each transitive-factor term is nonnegative by `(10)`. The term from
`Q_(n_j)` is evaluated with every factor to its right absorbed into `B`.
That suffix contains exactly `m-j` further `Q` blocks. Its Hamilton count is
therefore

```text
3^(m-j)H(B),                                               (16)
```

because all transitive factors have Hamilton count one. Equation `(12)` makes
the `j`-th cycle contribution at least

```text
q 9^(m-j)H(B)^2H(C)^2.                                  (17)
```

Summing `(17)` proves `(14)`.

For equality, every cycle contribution must attain THM-4212's equality
boundary and every transitive contribution must vanish. A non-rightmost
cycle sees a nonempty suffix before `B`, so its middle context is not a
singleton. A nonempty transitive suffix makes the same true for the rightmost
cycle. A nonempty transitive prefix contributes strictly because its middle
factor contains `Q_5`; THM-4187's only zero interaction has all three factors
singleton. Finally THM-4212 requires tail length five and singleton contexts.
These conditions are exactly `(15)`. QED.

Merging consecutive singleton strong components gives an intrinsic version
of `(13)`. Its strong-component word is

```text
P_(a_0) triangleright C3 triangleright P_(b_1)
 triangleright C3 triangleright P_(b_2)
 triangleright ... triangleright C3 triangleright P_(b_m), (18)

b_1,...,b_m>=5.
```

Conversely every word `(18)` has a factorization `(13)`, for example by
taking each `n_j=5`. Thus membership is not an artefact of a noncanonical
presentation: the ordered strong-component decomposition determines the
word. When `a_0=0`, the factor is source-free. For `m>=2` it has at least two
nontrivial strong components and therefore is not any single-cycle prefix
`C3 triangleright P_n`.

## 4. A new explicit `(OS+)` left-factor family

> **Corollary 4 (final-tail-six `(OS+)`).** Let `X` be a word `(18)` whose
> final transitive run satisfies `b_m>=6`. Then for every no-sink tournament
> `C`,
>
> ```text
> R_+(X,C)>0.                                             (19)
> ```

### Proof

Peel the last singleton strong component:

```text
X=Y triangleright P_1,                                  (20)
```

where `Y` still has the form `(18)`, now with final run `b_m-1>=5`.
Theorem 3 gives `F_Y(P_1,C)>0`. THM-4187 gives
`R_+(P_1,C)>0` whenever `C` has no sink. Hence

```text
R_+(X,C)=F_Y(P_1,C)+R_+(P_1,C)>0.                       (21)
```

QED.

For `a_0=0` and `m>=2`, Corollary 4 is an explicit source-free,
arbitrarily-many-cycle family of left factors satisfying `(OS+)`, beyond both
the transitive factors of THM-4187 and the single-cycle family of THM-4212.

## 5. Exact hostile: a good suffix is not enough

The ambient ideal cannot be enlarged from `N` to all tournaments. Let

```text
Z=C3 triangleright Q_5=C3 triangleright C3 triangleright P_5. (22)
```

This order-eleven tournament has two adjacent cycle components and ends in
the uniformly positive factor `Q_5`. Nevertheless exact evaluation gives

```text
F_Z(P_1,P_1)=-338580.                                  (23)
```

One compact certificate is

```text
H(Z)=9,
G_+(Z)=-191628,
G_+(Z triangleright P_1)=-395424,
G_+(Z triangleright P_2)=-734004.                       (24)
```

Since `G_+(P_1)=G_+(P_2)=0`, equations `(1)` and the definition of `R_+`
turn `(24)` into `(23)`. The literal pair-bit label of `Z` is

```text
1011111111111111111111111111011111111111111111111111111. (25)
```

Both replay engines reconstruct `(23)` from the actual tournament. The
independent one uses labelled adjacency and subset path DP, not the ordinal
capacity-transfer theorem.

Thus `Q_5 in U` but `C3 triangleright Q_5` is not even in `N`. Uniform
positivity is a two-sided ideal inside the proved neutral semigroup, not a
good-suffix property in the full tournament semigroup. In the language
`(18)`, the five singleton components after **each** cycle are the proved
separation mechanism; five vertices only after the last cycle do not suffice.

## 6. Connection contract and scope firewall

The proved connection is

```text
source:       ordered actual tournament factors and their contextual defects,
map:          iterated ordinal telescope (4),
preserved:    exact F and certified Hamilton-normalized lower floors,
controlled loss:
              the scalar floor forgets the optimal constant and response jet,
sidecar:      ordered factorization, suffix Hamilton count, and good-block positions,
decisive test: the literal adjacent-cycle hostile (22)--(25).
```

Theorems 1--3 prove a closure calculus and one regular strong-component
language. Corollary 4 proves `(OS+)` only for its final-tail-six subclass.
They do not classify `N` or `U`; do not show that an arbitrary positive
singleton defect is uniformly positive; and do not determine the threshold
for three or more consecutive cycle components or more general strong
components. THM-4215 separately closes the two-adjacent-cycle threshold.
General `(OS+)`, the no-sink/no-source gate law, the all-strong residual, and
the order-eleven asymmetric bank remain **OPEN**.

## 7. Replay

Primary exact audit:

```bash
python3 -B 04-computation/tournament_uniform_prefix_semigroup_thm4213.py
python3 -O -B 04-computation/tournament_uniform_prefix_semigroup_thm4213.py
```

Both streams byte-match the frozen primary output.

Independent clean-room literal audit:

```bash
clang++ -std=c++17 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_uniform_prefix_semigroup_independent_audit_thm4213.cpp \
  -o /tmp/thm4213-independent-O0
/tmp/thm4213-independent-O0

clang++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_uniform_prefix_semigroup_independent_audit_thm4213.cpp \
  -o /tmp/thm4213-independent-O3
/tmp/thm4213-independent-O3

clang++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_uniform_prefix_semigroup_independent_audit_thm4213.cpp \
  -o /tmp/thm4213-independent-san
ASAN_OPTIONS=detect_leaks=0 /tmp/thm4213-independent-san
```

The O0, O3, and ASan/UBSan streams byte-match the frozen independent output.

**QED for `(4)--(5)`, `(9)`, `(14)--(15)`, `(19)`, and `(23)`.**

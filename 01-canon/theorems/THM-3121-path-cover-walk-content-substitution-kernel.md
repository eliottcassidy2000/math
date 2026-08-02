---
id: THM-3121
title: "The path-cover walk-content substitution kernel and the exact C3 diagonal law"
status: >
  RESERVED / UNPROVED PROOF CANDIDATE UNDER AUDIT. For a quotient tournament
  Q and nonempty factor tournaments S_i, the complete path-cover profile of
  Q[S_i] is the factorial-Hadamard contraction of the factor profiles with
  exp(y W_Q), where W_Q=1^T X(I-A_Q X)^(-1)1 is the rational directed-walk
  content series. The d=1 slice gives Hamiltonian paths. For Q=C3 its kernel
  is a balanced diagonal: K(a,b,c)=a!b!c! times 3 when a=b=c, 1 when
  max(a,b,c)-min(a,b,c)=1, and 0 otherwise. A direct bijection is written
  below; the exact companion agrees with independent vertex-level enumeration
  on 216 kernel cells, 64 C3 block triples, 1024 general-quotient controls,
  and 23 full path-cover profiles. Normal and optimized transcripts match.
source: codex-2026-08-02-frontier-synthesis
depends_on: []
related:
  - THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant
  - THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law
  - THM-1862-order-join-reduction-principle
script: 04-computation/tournament_path_cover_walk_kernel_thm3121.py
output: 05-knowledge/results/tournament_path_cover_walk_kernel_thm3121.out
script_sha256: ee4e50577680077a295e383a856e3c64a160224d8f1a5a7ef5b929450f874eae
output_sha256: 649a687e2e43f57cc505aeb291a516660bf7e87ffbd6bc51f00f74199e23d19e
hash_basis: working-tree bytes (LF)
---

# THM-3121 — the path-cover walk-content substitution kernel

**RESERVED / UNPROVED PROOF CANDIDATE UNDER AUDIT.** Nothing in this file may
be used as a proved dependency until an explicit audit promotes the status.

## 1. Objects and the quotient walk series

For a finite tournament `T`, let

```text
pc_T(d) = number of spanning covers of T by d unordered vertex-disjoint
          directed paths, with singleton paths allowed.
```

Thus `pc_T(1)=H(T)`, the Redei Hamiltonian-path count. Let `Q` be a tournament
on `{1,...,q}`, with adjacency matrix `A_Q`, and let `S_1,...,S_q` be nonempty
factor tournaments. Write `Q[S_1,...,S_q]` for tournament substitution.

Put `X=diag(x_1,...,x_q)`. The directed-walk content series of the quotient is

```text
W_Q(x_1,...,x_q)
  = sum_(m>=1) sum_(i_1 ->_Q ... ->_Q i_m)
      x_(i_1)...x_(i_m)
  = 1^T X (I-A_Q X)^(-1) 1.                              (1)
```

A length-one word is allowed. Since `A_Q` has zero diagonal, consecutive
letters in every longer word are distinct. Equation (1) is a formal power
series identity; for fixed `Q` it is a rational function with denominator
dividing `det(I-A_Q X)`.

## 2. Exact substitution law

For `c=(c_1,...,c_q)` and `d>=1`, define

```text
K_(Q,d)(c)
 = (prod_i c_i!)/d! * [x_1^(c_1)...x_q^(c_q)] W_Q(x)^d.   (2)
```

Only `1<=c_i<=|S_i|` enters below.

> **Candidate theorem.** For every quotient tournament `Q`, every tuple of
> nonempty factor tournaments, and every `d>=1`,
>
> ```text
> pc_(Q[S_1,...,S_q])(d)
>   = sum_(c_1=1)^|S_1| ... sum_(c_q=1)^|S_q|
>       K_(Q,d)(c) prod_i pc_(S_i)(c_i).                  (3)
> ```
>
> In particular,
>
> ```text
> H(Q[S_1,...,S_q])
>   = sum_c [x^c]W_Q(x) prod_i (c_i! pc_(S_i)(c_i)).      (4)
> ```

Although (2) is displayed with a division by `d!`, it is an integer. The
bijection below proves both integrality and (3).

### Proof candidate: split at maximal block runs

Take a spanning `d`-path cover of `Q[S_i]`. Split each of its directed paths
at every change of substitution block. Within block `i`, the maximal runs form
a spanning directed path cover of `S_i`, say with `c_i` components. Recording
only the block labels along each original path gives an unordered set of `d`
nonempty directed walks in `Q`, using label `i` exactly `c_i` times.

Conversely, choose a `c_i`-path cover of every `S_i`. Its components are
distinct actual paths, even though the cover itself is unordered. Choose an
ordered `d`-tuple of quotient walks with total content `c`, and biject the
`c_i` components of block `i` with the occurrences of letter `i`. There are
`prod_i c_i!` such assignments. Concatenating the assigned components is
legal because every transition in a quotient walk is an arc of `Q`.

Finally forget the order of the `d` resulting global paths. The symmetric
group on those paths acts freely: two paths cannot be identical because they
use disjoint, nonempty sets of actual vertices. Division by `d!` is therefore
exact. These constructions are inverse. Their count is (2)--(3). Taking
`d=1` gives (4).

Equivalently, `exp(y W_Q)` is the multivariate exponential generating series
for an unordered set of quotient walks. The factorials in (2) label the
actual path components inside each color. This is a species explanation of
the same explicit bijection, not an additional assumption.

## 3. The exact directed-triangle kernel

Orient `C_3` as `1->2->3->1`. A quotient walk is then forced once its first
letter and length are chosen. Directly from (1),

```text
W_(C3)
 = (x_1+x_2+x_3 + x_1x_2+x_2x_3+x_3x_1 + 3x_1x_2x_3)
     /(1-x_1x_2x_3).                                   (5)
```

For positive `a,b,c`, a prefix of the periodic word `123123...` has either
equal content in all three letters or two adjacent levels differing by one.
There are three possible starts in the equal case and one in each unequal
balanced case. Hence

```text
[x_1^a x_2^b x_3^c] W_(C3)
 = 3,  if a=b=c;
 = 1,  if max(a,b,c)-min(a,b,c)=1;
 = 0,  otherwise.                                      (6)
```

Combining (4) and (6) discharges the undetermined integer kernel in THM-1975:

```text
K_(C3,1)(a,b,c)
 = a!b!c! * {3 on the diagonal, 1 on its two nearest shells, 0 elsewhere}.
                                                               (7)
```

Let `F_i(r)=r! pc_(S_i)(r)`, extended by zero outside `1<=r<=|S_i|`. Then

```text
H(C3[S_1,S_2,S_3])
 = 3 sum_(r>=1) F_1(r)F_2(r)F_3(r)
   + sum_(r>=1) sum_(i=1)^3 F_i(r+1) prod_(j!=i) F_j(r)
   + sum_(r>=1) sum_(i=1)^3 F_i(r)   prod_(j!=i) F_j(r+1).     (8)
```

Thus the apparent three-dimensional transfer is only a diagonal plus its two
nearest shells. Once the three block profiles are known, (8) takes linear time
in the shortest profile.

For three equal factors `S`, writing `F(r)=r!pc_S(r)` gives

```text
H(C3[S,S,S])
 = 3 sum_(r>=1) (F(r)^3 + F(r+1)F(r)^2 + F(r+1)^2F(r)).        (9)
```

The directed triangle itself has

```text
pc_(C3)=(3,3,1),       F=(3,6,6).
```

Equation (9) yields

```text
H(C3[C3,C3,C3]) = 1377 + 810 + 972 = 3159,                   (10)
```

matching independent Held--Karp enumeration. The previously unexplained
factor `13` in `3159=3^5*13` is the combined contribution of the two nearest
balanced shells, not an additional scalar invariant of a block.

## 4. Boundaries and operation response

- **Transitive quotient.** If `Q` is transitive, a directed quotient walk is
  strictly increasing. The only all-positive content is `(1,...,1)`, with one
  word. Equation (4) reduces to `H(Q[S_i])=prod_i H(S_i)`, recovering the
  order-join law of THM-1862.
- **Nontransitive quotient.** Cycles permit repeated block visits, so scalar
  `H(S_i)` is insufficient. The full path-cover profile is the sidecar that
  records how many maximal block runs can be threaded. THM-1975's scalar-H
  hostiles are therefore explained rather than bypassed.
- **No complexity collapse.** Computing an arbitrary factor profile remains
  `#P`-hard because its first entry is `H`. Equations (1)--(4) give a
  functorial closed transfer once profiles are supplied; they do not turn
  arbitrary Hamiltonian-path counting into polynomial time.
- **General quotient cost.** For fixed `Q`, the universal kernel is rational
  and its bounded coefficient table can be generated by a last-vertex dynamic
  program. The exceptional linear-time collapse in (8) uses the deterministic
  outdegree-one walk structure of `C_3`.

The source is a factor tournament tuple, the target is its substitution, the
map is maximal-block-run splitting, the preserved predicate is the complete
spanning path-cover count, and scalarization destroys the run count. The
factorial-transformed path-cover profile and quotient walk content are the
restoring sidecars. The transitive product and the `3159` directed-triangle
cube are respectively the hostile boundary and positive cyclic control.

## 5. Exact referee

Run

```bash
python3 04-computation/tournament_path_cover_walk_kernel_thm3121.py
python3 -O 04-computation/tournament_path_cover_walk_kernel_thm3121.py
```

The dependency-free companion checks:

- all `216` positive `C_3` kernel cells with `1<=a,b,c<=6` against a separate
  last-letter walk dynamic program;
- `64` ordered triples from `{T_1,T_2,T_3,C_3}` against vertex-level
  Held--Karp enumeration and (8);
- all `64` labelled four-vertex quotient tournaments and all `16` singleton/
  edge block-size words, for `1024` general-quotient controls;
- `23` complete path-cover profiles against independent successor-function
  enumeration; and
- (10) directly.

Normal, optimized, and frozen transcripts are byte-identical.

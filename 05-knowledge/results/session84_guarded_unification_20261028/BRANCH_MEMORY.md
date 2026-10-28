# Branch residues, the missing sign, and what “memory” actually means

**Session84, 2026-10-28. Status: PROVED elementary deductions from the audited
THM-4808; FINITE-EXACT independent replay in the declared universes.**
No new claim about the classification of Berggren trees, Paley tournaments,
or universal spectral bounds is made. No integer–tournament isomorphism is
asserted.

## 1. Inheritance and the claim boundary

All inherited files below were read using `git show` at
`ae7c073e864d05be74fed85ba0a2bc164e026a6d`:

- `01-canon/theorems/THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md`,
  especially (a), (b), (c), and the corrected (G).
- `04-computation/experiments/berggren_sieve_parity_20261026/audit_thm4808/REPORT.md`,
  especially §§2, 3.3–3.4, 4. The audit overrides pre-correction statements.
- `05-knowledge/results/twin_pythagorean_124_20261026/finite_field/README.md`,
  §§4–5, for the projective endpoint hostile and determinant/design bit.
- `05-knowledge/results/growing_sieve_path_20261026/field/README.md`, §§2,4–5,
  for the 16-state marked-tournament observer and finite-wheel loss boundary.
- `07-reflections/why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md`:
  a useful historical schema, **not a truth source over the repaired canon**.
- `05-knowledge/results/predecessor_prime_shapes_20261014b.md` §1 and
  `05-knowledge/results/session83_aliquot_divisor_grammar_20261028/CARRIER_AND_144_BRIDGE.md`
  §1, for the comparison in §7.

Closest mechanism: the audited every-modulus endpoint collision. Canonical
hostile: `UU` versus `A` on the projective conic modulo 7. Corrected near miss:
“no fixed endpoint modulus sees the sign” does **not** mean “no finite-state
observer sees the sign.” Least-used sidecar: the *normalized pole vector*,
as distinct from its projective line.

Concept board: branch word / endpoint residue / determinant bit / normalized
pole frame / guarded arithmetic continuation / full-source reconstruction.

## 2. A universal hostile, with its full matrix retained

Let words be written in chronological order, starting from `v0=(2,1)^T`.
For `w=uv`, the matrix convention is `M_w=M_v M_u`. The branch matrices are

\[
U=\begin{pmatrix}2&-1\\1&0\end{pmatrix},\quad
A=\begin{pmatrix}2&1\\1&0\end{pmatrix},\quad
D=\begin{pmatrix}1&2\\0&1\end{pmatrix}.
\]

Write `χ(w)=det M_w=(-1)^{#A(w)}`. For every integer `M>=1`, put

\[
w_M=A D^{M-1} U D.
\]

**Proposition 1 (inherited hostile, explicit matrix).**

\[
M_{w_M}=DUD^{M-1}A
=\begin{pmatrix}8M-1&4\\2M&1\end{pmatrix},\qquad
M_{w_M}v_0=(16M+2,4M+1)^T. \tag{1}
\]

The displayed matrix product proves (1) directly. Its determinant is `-1`.
The endpoint is primitive because
`gcd(16M+2,4M+1)=gcd(2,4M+1)=1`, has opposite parities, and has first
coordinate larger than its positive second coordinate. It reduces to
`(2,1)` modulo `M`. Therefore its Euclid triple also reduces to `(3,4,5)`.
The empty word has the same endpoint residue and has character `+1`.

Thus **no function of a fixed endpoint-vector, endpoint-slope, or triple
residue determines χ on all integer Berggren words**. Any fixed finite
collection of such moduli also fails, by using their least common multiple.
This is an all-modulus theorem, not a search-based extrapolation.

**Scope guard.** It is false to replace “endpoint residue” by “arbitrary
residue information about the word.” The full matrix modulo 3 determines
χ immediately: its determinant is `1` or `2`, representing `+1` or `-1`.
Endpoint reduction forgets the action away from its single starting vector;
matrix reduction does not forget that action.

## 3. Every endpoint fibre splits, and one bit is exactly sufficient

Define the reachable finite endpoint set

\[
E_M=\{M_wv_0\bmod M:w\in\{U,A,D\}^{*}\}.
\]

The endpoint observation `O_M(w)=M_wv0 mod M` is **already a right
congruence**: if `O_M(u)=O_M(v)`, then for every appended word `z`,

\[
O_M(uz)=M_z O_M(u)=M_z O_M(v)=O_M(vz). \tag{2}
\]

The failure is target descent of χ, not failure of endpoint dynamics.

**Proposition 2 (fibrewise tight repair).** The reachable joint observer
is exactly

\[
\{(O_M(w),\chi(w)):w\in\{U,A,D\}^{*}\}=E_M\times\{+1,-1\}. \tag{3}
\]

*Proof.* Choose any word `v` reaching `x in E_M`. By (1) and (2), `v`
and `w_M v` have the same endpoint residue. Multiplicativity gives
`χ(w_M v)=-χ(v)`, so both signs occur over `x`. The reverse containment
is immediate. ∎

The transition law is

\[
(x,s)\xrightarrow{X}(Xx,\det(X)s),\qquad X\in\{U,A,D\}. \tag{4}
\]

Consequently a **single binary sidecar** is sufficient. It is fibrewise
necessary: every endpoint fibre contains both target values. Any observer
from which both the present endpoint residue and χ can be recovered must
have at least `2|E_M|` reachable states. Equation (4) attains that bound.
This is minimality for that specified pair of outputs, not for all possible
questions about an integer triangle.

Moreover no amount of future *endpoint-only* observation recovers the lost
sign: for every continuation `z`, the words `vz` and `w_M vz` still have
equal endpoint residues and opposite χ. The loss is permanent in this
quotient even though it costs only one bit to avoid it.

For χ alone, the two-state update “A flips; U,D preserve” is exact and
minimal. The two present outputs already distinguish the two states, and
all continuations preserve whether the two signs agree or disagree.

## 4. Which pole data repairs which observation?

Let `p_w=(α,γ)^T=M_w(1,0)^T` and `v_w=(m,n)^T=M_wv0`.
The inherited pole-frame identity is

\[
\det[p_w\mid v_w]=\alpha n-\gamma m=\chi(w). \tag{5}
\]

This uses **vectors with their scale**, not just two projective points.
Indeed the two vectors recover the entire matrix:

\[
M_w=[p_w\mid v_w-2p_w]. \tag{6}
\]

- Exact integer frame: (5) recovers χ; (6) retains the full action.
- Frame vectors modulo `M>=3`: the determinant residues `+1` and `-1`
  are distinct, so χ is recovered. For `M=1,2` this fails in general.
- A frame modulo 3 appended to any chosen endpoint modulus therefore
  suffices. It is larger than the minimal explicit χ sidecar.
- **Projective pole plus projective endpoint does not suffice at a fixed
  modulus.** In (1), `p_(w_M)=(8M-1,2M) ≡(-1,0)`, projectively the same
  as the root pole `(1,0)`, while the endpoint agrees and the signs differ.
  Over a composite modulus this statement uses projective classes of
  unimodular vectors under unit rescaling; multiplication by `-1` is a unit.

A concrete modulo 7 witness is

| Word | Endpoint | Normalized pole vector | χ |
|---|---|---|---:|
| empty | `(2,1)` | `(1,0)` | +1 |
| `A D^6 U D` | `(114,29)` | `(55,14)` | -1 |

Both endpoint vectors are `(2,1)` modulo 7. Both projective pole images are
infinity, but the pole vectors are `(1,0)` and `(-1,0)` modulo 7. This is
the exact normalization that a projective-pole-only “repair” would lose.
At modulus2 the smaller word `ADUD` even has the full matrix equal to the
identity modulo 2, while its character is negative.

The inherited smaller *endpoint-only projective* hostile is `UU` versus
`A`: endpoints `(4,3)` and `(5,2)` have common slope 6 modulo 7 and opposite
χ. Their endpoint vectors do not agree modulo 7, so this is a different,
weaker observation than the vector quotient in (1).

The inherited 16-state observer `(point mod7,χ)` repairs a specified
projective/design-orientation loss. It can be realized by actual oriented
Paley tournaments on the complement of a moving marked point. The source
distinguishes the original pole at infinity from the root coordinate 2;
the root-based homogeneous-space description requires conjugation. None
of this identifies a full integer triangle or its factorization with a
tournament. The bit repairs orientation, not height or primality.

## 5. Finite target memory versus unbounded source information

For every fixed modulus `M`, all words `D^(kM)` have the same full matrix
modulo `M` and the same χ, but their endpoints are

\[
(2+2kM,1),\qquad k=0,1,2,\ldots. \tag{7}
\]

So even finite matrix/frame data plus χ cannot reconstruct the unbounded
integer source. At modulus7, the empty word and `D^7` already give endpoints
`(2,1)` and `(16,1)`, with equal full matrix residues and sign.

Conversely the exact integer endpoint determines its unique Berggren word
by the inherited tree theorem. It therefore determines χ, without an
additional independent mathematical coordinate. The audited continued-
fraction formula gives another recovery route:

\[
\chi=(-1)^{m+c_0+\cdots+c_L},\qquad m/n=[c_0;\ldots,c_L].
\]

Keep three statements separate:

1. Recovering χ from **fixed endpoint residues** is impossible.
2. Updating χ while reading branch letters needs **one bit**, not the whole
   history.
3. Reconstructing arbitrarily large integer endpoints/whole words requires
   an infinite state set, simply because the possible outputs are infinite.

Thus the historical reflection's phrase “signs and characters of the walk
need the whole word” is too strong if read as a memory lower bound. The
valid replacement is “the endpoint quotient loses the sign; retain its
two-state character, an adequate normalized frame, or exact integer data.”

## 6. Corrected non-backtracking and spectral context

The rational graph has **two components**, interchanged by
`C(t)=(t+1)/(t-1)`. The Euclid component is an **edge-indexed** tree, with
a half-edge at the degenerate triangle; its covering tree is 4-regular.
It is not a plain 4-regular tree on the unsigned primitive triangles.

The audited unsigned exceptions to the Ramanujan bound for primes
`5<=q<1000` are `331,569,719,809,887`. The signed operator's exceptions
are `251,479,863,947`. These are distinct lists for distinct operators;
neither supports universal Ramanujan behavior. This packet inherits those
audited finite claims and does **not** rerun their spectral computation.
The finite-state results above use only elementary matrix algebra and do
not need expansion, equidistribution, or a spectral conjecture.

## 7. The exact bridge to F=U+S and Session83 source-loss guards

For `N>=2`, F counts proper nontrivial divisors, S counts the squarefree
ones among them, and U counts prime proper divisors. Their sets overlap;
`F=U+S` is not a disjoint partition. Put `δ=F-U-S`. For the inherited
Session83 notation, `s(n)=σ(n)-n` is the aliquot sum; it is a different
observation, not this F.

Three tiny comparisons separate the logical failure modes:

| Observation | Equal observed inputs | Common lawful continuation | Result |
|---|---|---|---|
| Berggren endpoint mod7 | empty and `A D^6 U D` | append any branch word | Endpoints remain equal; χ remains opposite |
| FUS defect δ | `δ(4)=δ(12)=-1` | multiply by fresh prime 5 | `δ(20)=-1`, `δ(60)=0` |
| Aliquot sum s | `s(2)=s(3)=1` | multiply by fresh prime 5 | `s(10)=8`, `s(15)=9` |

The arithmetic continuations satisfy the coprimality guard at both inputs.
Their coarse observations fail continuation congruence itself. The branch
endpoint observation **passes** congruence but fails descent of an extra
target. The common method is the same: specify the source, legal contexts,
observation and target, then test whether equal observations preserve the
needed future outputs. The algebraic objects and required repairs are not
thereby isomorphic.

For s, the source term appears visibly in the fresh-prime identity

\[
s(nr)=(r+1)s(n)+n\quad (r\text{ prime},\ r\nmid n).
\]

This is why equal present aliquot outputs can diverge. It is not analogous
to “an unbounded branch-sign history”: the sign's exact update is already
two-state. This distinction is the useful source-loss guard inherited by
the unified observation treatment.

## 8. Exact replay and stopping boundary

Run in this directory:

```sh
python3 -I -B check_branch_memory.py > replay.json
python3 -I -B -O check_branch_memory.py > replay.optimized.json
cmp replay.json replay.optimized.json
```

The checker is independent standard-library code with live checks under
`-O`. It verifies the literal word product for `M=1..100`, the universal
matrix/endpoint/pole formula for `M=1..5000`, all 9841 words through depth 8
against the inverse decoder, determinant/pole identity and continued-
fraction sign, and **complete** reachable state sets for `M=1..32`.

Selected endpoint versus endpoint-plus-bit counts are
`M=1:1/2`, `M=2:2/4`, `M=3:8/16`, `M=7:48/96`,
`M=15:192/384`, and `M=32:512/1024`.
It also checks both projective hostiles, full-matrix loss, and the two
guarded arithmetic counterexamples by direct divisor enumeration. Normal
and optimized outputs are byte-identical. These finite checks supplement
the displayed all-modulus proofs; they neither classify all endpoint
orbits arithmetically nor establish a new spectral or prime-density result.

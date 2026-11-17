# The curried grammar: one-vertex extension is a quadratic form, the operator algebra, and two encodings of the holes {7, 21}

Session klein-2026-11-16-tooth (2026-11-16). Owner: *"come to understand a full recursive grammar on the infinite group that
naturally encodes the {7, 21} impossibility. Work to see abstract structures as compositions of functions and curry them."*
The grammar taken here is the complete one — every tournament class of order `n + 1` is a one-vertex extension of a class
of order `n` (the extension graph of [TOOTHPICK.md](TOOTHPICK.md)) — with the count functional written as a composition of
the extension operator and a polynomial in the extension pattern, i.e. curried. Canon: THM-5266.

## 1. The count of a one-vertex extension is a quadratic form on the cube (PROVED)

Let `T` be a tournament on `V`, `|V| = n`, with `h = H(T)`, `e_x` = Hamiltonian paths ending at `x`, `s_x` = starting at
`x`, `M(y, x)` = ordered pairs `(P_1, P_2)` of **nonempty** paths covering `V`, `P_1` ending at `y`, `P_2` starting at `x`
(`y ≠ x`), and `B_x = s_x + Σ_y M(y, x)` the root response of THM-4865. For a new vertex `v` with out-neighbourhood
`X ⊆ V` (indicator `x ∈ {0,1}^V`):

    H(T + v_X) = h + Σ_{x ∈ X} (B_x − e_x) − Σ_{y ≠ x, both in X} M(y, x)  =  h + ⟨B − e, x⟩ − xᵀ M x.

*Proof.* A Hamiltonian path of `T + v_X` has `v` first (`Σ_{x ∈ X} s_x` paths), last (`Σ_{y ∉ X} e_y`), or inside, where it
is `P_1 v P_2` with `last(P_1) = y ∉ X` beaten... i.e. `y → v` (`y ∉ X`) and `v → first(P_2) = x ∈ X`: `Σ_{y ∉ X, x ∈ X} M(y, x)`.
Write `Σ_{y ∉ X} e_y = h − Σ_{y ∈ X} e_y` and `Σ_{y ∉ X, x ∈ X} M = Σ_{x ∈ X} Σ_y M(y, x) − Σ_{y, x ∈ X} M(y, x)`, and use
`s_x + Σ_y M(y, x) = B_x`. ∎ Replay: all 62,422 extensions of all tournaments of orders 1–7, 0 mismatches
([receipts/quadform.out](receipts/quadform.out)).

Specialisations: `X = ∅` (sink) and `X = V` (source) give `h`; `X = {t}` (sole exit) gives `h + B_t − e_t`, THM-4862 (a);
`X = V ∖ {t}` (sole entry) gives the converse statement. **Currying.** `EXT: T ↦ (X ↦ T + v_X)` and
`H ∘ EXT_T: {0,1}^n → odd numbers` is a polynomial of degree ≤ 2 whose coefficients are the 1- and 2-path statistics of
`T`. The spectrum therefore obeys the exact recursion `S(n+1) = ∪_{T ∈ C(n)} H_T({0,1}^n)`: the counts of order `n + 1` are
the values of the quadratic forms of the classes of order `n` on the corners of their cubes.

**Unequal growth of values.** One class of order 7 produces only 4 distinct child counts, one produces 97 of a possible
128 (by order of the parent the maxima are 1, 2, 3, 5, 13, 38, 97). The image sizes at order 7 run through the whole range
4–97 ([receipts/quadform.out](receipts/quadform.out)).

## 2. The operator algebra and what each operator needs to know

| operator | definition | count law | state consumed | state produced for the new vertex |
|---|---|---|---|---|
| sink `σ` | `X = ∅` | `h` | `h` | `e_s = h`, `B_s = h` (MISTAKE-776), `B_x ↦ 2 B_x`, `e_x ↦ 0` for `x ∈ T` |
| source `ς` | `X = V` | `h` | `h` | `s_src = h`, `B_src = P_T(2)`, `B_x, e_x` unchanged |
| sole exit `ε_t` | `X = {t}` | `h − e_t + B_t` | `(h, e_t, B_t)` | `e_z = h − e_t`, `B_z = h + B_t` (THM-4865) |
| ordinal sum `A ⇒ D` | all of `A` beats `D` | `h(A) h(D)` | `h` | `B_t = B_t(A) P_D(2)` or `h(A) B_t(D)`; `e_t = 0` or `h(A) e_t(D)`; `P = P_A P_D` (THM-5260 (B)) |
| general `EXT_X` | any `X` | quadratic form | `(h, e, B, M)` | needs the 3-path statistics of `T` |

Compositions. The near-join is `ε_t ∘ σ^k`: after `k` sinks `e_t = 0` and `B_t ↦ 2^k B_t`, so `H = h + 2^k B_t` —
curried, `T ↦ (t ↦ (k ↦ h + 2^k B_t))`, a line in `k` (HYP-9453). The tower is `ε ∘ ε ∘ …` applied each time at the
newest vertex; the 3-state transformer `(h, e, B) ↦ (h − e + B, h − e, h + B)` makes it tribonacci (THM-4865). The
midpoint value of THM-5264 is `ε_{t_0} ∘ (⇒ D)`, `h(A)h(D) + B_{t_0}(A) P_D(2)`. **Nonuniform contacts cost state**: the
uniform contacts (sinks, sources, sole exits, ordinal sums) are closed on `(h, e, B, P)`; a general out-neighbourhood
needs the pair matrix `M`, and two general extensions in a row need the 3-path tensor. The grammar's cheap productions
are exactly the uniform ones, and both encodings of the holes live in them.

## 3. Two encodings of {7, 21}

**(a) In the full grammar: multiplicativity and the floors (THM-1370's route).** The count is multiplicative over the
strong components (`⇒` is a product), so a tournament with count 7 (prime) or 21 (`3·7` or `21·1`) has a strong component
of count 7 or 21. Busch's strong floors are 3, 5, 9, 15, 25 at orders 3–7 and increase, so that component has order ≤ 6,
and the strong spectra of orders 3–6 — `{3}`, `{5}`, `{9, 11, 13, 15}`, `{15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45}`
(sterile9.out of the outdeg packet) — contain neither. In the grammar's own terms: the values of the quadratic forms on
strong extensions through order 6 avoid 7 and 21, and from order 7 on every strong value exceeds 21.

**(b) In the additive sub-grammar: the 2-adic structure of the smallest response sets (PROVED, new).** Let `𝔅(h) =
{B_t(T) : H(T) = h}`. The tournaments of count 1, 3, 5 are the transitive ones, the padded 3-cycles `M_a ⇒ C_3 ⇒ M_b`, and
the padded strong order-4 tournament (counts multiply over components; strong tournaments of order ≥ 4, resp. ≥ 5, have
count ≥ 5, resp. ≥ 9). With `B_t(M_n) = 2^{n−1−i}` for the vertex in position `i`, `B_t(C_3) = 3` (digit 0 at every
vertex), the responses `5, 6, 7, 9` of the strong order-4 tournament, `P_{C_3}(2) = 12`, `P_{S_4}(2) = 32` and the
ordinal-sum laws:

    𝔅(1) = {2^j},   𝔅(3) = {3·2^j},   𝔅(5) = {5, 6, 7, 9, 32}·2^j,

exactly as the census of orders ≤ 8 shows ([receipts/holes_grammar.out](receipts/holes_grammar.out)). Hence for `k ≥ 1`:
`7 − h ∉ 2^k 𝔅(h)` for `h = 1, 3, 5` (`6`, `4`, `2`), `21 − h ∉ 2^k 𝔅(h)` (`20`, `18`, `16`), and `h ≥ 9` cannot reach 7 or 21
since `B ≥ h`. For the dead-end appends (`k = 0`) the only new candidate is `21 = 9 + 12`, and `12 ∉ 𝔅(9)`: the padded
`C_3 ⇒ C_3` give `9·2^j, 36·2^j`, and the two strong order-5 classes of count 9 have response odd parts
`{1, 5, 9, 11, 15, 21, 23, 25}`. So **the near-join, tower and dead-end-append productions never output 7 or 21**, for a
reason internal to the grammar: `7 = 2^3 − 1` is not of the form `2^m + 1` (the values over transitive bases are the
Fermat-type numbers 3, 5, 9, 17, 33, …), `21 = 3·7` is not `3(2^m + 1)` (the values over `C_3`: 9, 15, 27, 51, …), and the
responses of the count-5 tournaments contain no power of two below 32. The holes are the Mersenne shadows of 9 and 27.

Both encodings are consistent with THM-5263/THM-5264: 11, 13, 47, 67 are also missed by (b) (they are sterile-only), but
they are counts; the sub-grammar is not complete, and (a) is the only encoding that is.

## 4. First creators as the grammar's new words

A first creator is a class whose count first appears at its order (THM-4922). With the spectrum recursion above, first
creators are the classes `T + v_X` whose value `H_T(x)` is new at order `n + 1`. They are the common case, not the
exception: 1, 1, 6, 33, 328 of 456, 5457 of 6880 and 163,676 of 191,536 classes at orders 3–9 are first creators
([receipts/creators.out](receipts/creators.out), [CREATORS.md](CREATORS.md)). The 63 orphans of THM-4922 and the relation
of first creators to the toothpick tree are in CREATORS.md and TOOTHPICK.md.

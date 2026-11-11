# Second route to detection: periods of the family of Stein spaces (draft audited as "Theorem 5.2")

Status. This is the draft that the first audit examined as an addendum (`audit/AUDIT.md`, section 7: the proof is
correct). It was written after the pre-audit checkpoint and never pushed on its own. It is kept verbatim below for
the record. As a route to detection it is superseded by Theorem 5.4 and Corollary 5.5 of `README.md` (one
direction, no comparison theorem), and it is CONCURRENT with Theorem 1 of Session269
(`05-knowledge/results/session269_centres_ports_switches_20261111/jacobian/proof.md`), which was pushed first.
Theorem and corollary numbers below are those of the draft, not of `README.md`; "Theorem 4.1" and "Theorem 3.1"
refer to the pre-audit numbering, and the words "allowed in Theorem 3.1" are no longer needed (every centre is
allowed on `C^2`). The rungs `K-GMC_D` at the end are restatements, not intermediate statements, and the example
in the last paragraph needs `a != 0` (first audit).

---

### 5.2 Detection over an open set of centres (added after the pre-audit checkpoint)

Write `M_c`, `S_c`, `[h]_c` for the Stein kernel, the Stein space and classes at the centre `c = (a, b)`.

**Theorem 5.2 (PROVED from two cited facts on holonomic modules).** Let `(f, g)` be a Keller pair on `C^2`. Suppose
there is a nonempty open set `O` of `C^2` (classical topology) such that SN holds at every centre `c` in `O`: for
every `h` there is `m_0(h, c)` with `[d_u^m h]_c = [d_v^m h]_c = 0` for `m >= m_0`. Then `(f, g)` is an automorphism.

*Proof.* (1) *The family.* Let `G = C[a, b] (x) C[x, y] / R`, where `R` is spanned by the elements
`(d_u - g + b) k` and `(d_v - f + a) k` (`k` in `C[a, b, x, y]`). The operators `D_a = d/da + d_u` and
`D_b = d/db + d_v` commute with each other and preserve `R`:
`[D_a, d_v - f + a] = [d_u, -f] + [d/da, a] = -1 + 1 = 0`, `[D_a, d_u - g + b] = 0`, and symmetrically. So `G` is a
module over the Weyl algebra of the `(a, b)`-plane, and by right exactness its fibre at `c` is `C[x,y] / M_c = S_c`.
In `G`, `D_a` acts on the class of `phi (x) k` as `d/da + (g - b)`, because `d_u k = (g - b) k` modulo `R`; so `G`
is the top cohomology of the direct image, under `(a, b, x, y) -> (a, b)`, of the module
`C[a, b, x, y] e^(-(f - a)(g - b))`. CITED: this module is holonomic and direct images of holonomic modules are
holonomic (Bernstein); a holonomic module is a vector bundle with an integrable connection over a dense Zariski
open set `Omega`. So over `Omega` the spaces `S_c` form a bundle of some rank `r` with flat connection
`D_a, D_b`, and for `h` in `C[x,y]` the section `sigma_h = [1 (x) h]` satisfies `D_a sigma_h = sigma_(d_u h)`,
`D_b sigma_h = sigma_(d_v h)`.

(2) *Periods.* Take a polydisc `Delta = Delta_1 x Delta_2` inside `O cap Omega` and a basis
`lambda^(1), ..., lambda^(r)` of analytic flat sections of the dual bundle over `Delta`. For `h` in `C[x,y]` put
`pi_h^(i)(c) = lambda_c^(i)([h]_c)`. Flatness gives, on `Delta`,

    d/da pi_h^(i) = pi_(d_u h)^(i),      d/db pi_h^(i) = pi_(d_v h)^(i).

(3) *Polynomiality.* Fix `h` and `i`. For `c = (a_0, b_0)` in `Delta`, SN at `c` gives
`(d/da)^m pi_h^(i)(c) = lambda_c^(i)([d_u^m h]_c) = 0` for `m >= m_0(h, c)`: the Taylor series of the analytic
function `a -> pi_h^(i)(a, b_0)` at `a_0` is a polynomial, so this function is a polynomial on `Delta_1`. Let
`E_N` be the set of `b_0` in `Delta_2` for which its degree is at most `N`. The `E_N` cover `Delta_2`, so some
`E_N` is uncountable and has an accumulation point in `Delta_2`; for each fixed `a`, the analytic function
`b -> (d/da)^(N+1) pi_h^(i)(a, b)` vanishes on `E_N`, hence on `Delta_2`. In the same way
`(d/db)^(N'+1) pi_h^(i) = 0` on `Delta`. So every `pi_h^(i)` is a polynomial on `Delta`, of bidegree at most
`(N, N')` with `N, N'` depending on `h` only (take the largest over `i`).

(4) *Finite chaos.* For `c` in `Delta` and `m > N` or `n > N'`:
`lambda_c^(i)([d_u^m d_v^n h]_c) = (d/da)^m (d/db)^n pi_h^(i)(c) = 0` for every `i`, and the `lambda_c^(i)` span
the dual of `S_c`. So `[d_u^m d_v^n h]_c = 0` whenever `m + n > N + N'`: `(f, g)` has finite chaos at every centre
of `Delta`. `Delta` contains centres allowed in Theorem 3.1, and Theorem 4.1 applies. QED

**Corollary 5.3 (the ladder closes when the centre is free).** For a Keller pair `(f, g)` the following are
equivalent:

1. `(f, g)` is an automorphism;
2. `dim S_c = 1` at one centre allowed in Theorem 3.1;
3. finite chaos at one such centre;
4. K-GMC at every centre: every `M_c` is a Mathieu-Zhao subspace of `C[x, y]`;
5. CV at every centre of some nonempty open set;
6. SN at every centre of some nonempty open set: for every `h`, `(f - a)^m h` and `(g - b)^m h` lie in
   `(d_u - g + b) C[x,y] + (d_v - f + a) C[x,y]` for `m` large.

Hence

    JC(2)  <=>  the Gaussian moments property holds in the coordinates of every Keller pair, at every centre,

and GMC(2) (THM-2022) is the case of the pair `(x, y)`. With the cited bounds on the generic degree the statement
splits into rungs that do not involve the degrees of `f` and `g`:

    K-GMC_D : every Keller pair of generic degree at most D satisfies K-GMC at every centre
            <=>  there is no Keller pair of generic degree d with 2 <= d <= D.

`K-GMC_1` is GMC(2) together with its transport to automorphisms; `K-GMC_4` holds by the results of Orevkov and
Domrina cited in `CORE-PAPERS-PLANAR-JC.md`; the union over `D` is JC(2).

*Proof.* 1 => 4 is Proposition 2.3 applied at every centre (translations are automorphisms). 4 => 5 => 6 is
Proposition 5.1. 6 => 1 is Theorem 5.2. 1 <=> 2 <=> 3 are Theorems 3.1 and 4.1. For the rungs: a Keller pair of
generic degree 1 is an automorphism, and one of generic degree `d >= 2` violates SN, hence K-GMC, at some centre
of every open set. QED

What stays open (HYP-9668): detection at a single centre, and detection in one direction only. On `C^* x C` the
double cover of section 6.1 satisfies SN in the direction `d_v` at every centre (there `d_v` is locally nilpotent)
and violates it in the direction `d_u`: a single direction does not detect non-properness on a general surface.

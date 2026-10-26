# Builds required_edits.json and validates it: every "old" occurs exactly once in the CURRENT file, and all edits of
# a file apply in sequence (checked on in-memory copies; the repository files are not modified).  Audit 2026-10-25.
import json, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
THM = '01-canon/theorems/THM-4780-paley-census-tower-to-six-vertices-wheel-octahedron-pgl2-hamiltonian-path-complement.md'
HYP = '05-knowledge/hypotheses/HYP-9390-the-six-vertex-paley-census-is-the-level-32-adjoint.md'
RES = '05-knowledge/results/pythagorean_hourglass_census_tower_20261025.md'
E = []
def edit(f, old, new, why): E.append({'file': f, 'old': old, 'new': new, 'why': why})

# ------------------------------------------------------------------ THM-4780
edit(THM, "After THM-4652's reduction calculus, the only dense graph sums",
     "After the degree-at-most-two elimination calculus of the 2026-10-14 run packets and THM-4652's vanishing laws, the only dense graph sums",
     "Credit: the vertex-elimination calculus is from wheel_run_interface_20261014h.md / paley_next_cores_20261014i.md, not THM-4652 (which has the vanishing laws and series-parallel polynomiality only).")
edit(THM, "  Not yet independently audited.",
     "  Audited 2026-10-25 (04-computation/experiments/pythagorean_hourglass_20261025/audit/REPORT.md): CORRECT WITH FIXES\n"
     "  (attributions and wording only). Part 5 re-derived by exact symbolic elimination, with no interpolation; TT5, TT6 and the\n"
     "  graph sums recomputed by different algorithms; part 4 tested exhaustively on P^1(F_7) and on further even-degree graphs.",
     "Record the independent audit.")
edit(THM, "  - THM-4652   # reversal signs, odd and isolated-edge vanishing, odd-reversal symmetry, reduction calculus",
     "  - THM-4652   # reversal signs, odd and isolated-edge vanishing, odd-reversal symmetry\n"
     "  # degree <= 2 elimination calculus: 05-knowledge/results/wheel_run_interface_20261014h.md (codex-tiling) and\n"
     "  # 05-knowledge/results/paley_next_cores_20261014i.md (tiling-modular-atoms)",
     "Same credit correction in the dependency list.")
edit(THM, "2. **Reduction (THM-4652's calculus, with degree bounds).**",
     "2. **Reduction (the degree-≤2 elimination calculus of `wheel_run_interface_20261014h.md`, with degree bounds).**",
     "Same credit correction.")
edit(THM, "   - `H` runs over dense oriented graphs on fewer vertices, together with the empty graph on no vertices (`T = 1`).",
     "   - For `G` not dense, `H` runs over dense oriented graphs on fewer vertices, together with the empty graph on no vertices (`T = 1`); a dense `G` is its own reduction.",
     "As written the statement excludes H = G, which is false for dense G.")
edit(THM, "a doubled edge gives `χ² = 1` on injective maps. The coefficient is `±1`.",
     "a doubled edge gives `χ² = 1` on injective maps, or `χ(−1) = −1` when the two factors are oppositely oriented. The coefficient is `±1`.",
     "Opposite-orientation doubled edges give -1, not +1 (the sign is what the later sentence about chi(-1) relies on).")
edit(THM, "they contribute `(−1)^{deg'(w) + in'(w) + in'(u)}`. Every one of these exponents is even, so `F` is `PGL(2)`-invariant.",
     "they contribute `(−1)^{deg'(w) + in'(w) + in'(u)}`, primes denoting degrees into the other `k − 2` vertices. Both exponents are even: "
     "if `u ~ w`, then `deg'(w)` and exactly one of `in'(u)`, `in'(w)` drop by one, so the second exponent has the parity of "
     "`deg(w) + in(w) + in(u)`. So `F` is `PGL(2)`-invariant.",
     "deg', in' were undefined, and the individual terms are odd when u ~ w; only the total exponent is even.")
edit(THM, "Its vertex-deleted graph is `W4`; relabelling the hub from the largest to the smallest label reverses 4 arcs, so the sign stays `+1`.",
     "Deleting `u = 6` leaves `W4` with hub `5`; relabelling the hub from the largest to the smallest label reverses 4 arcs, so the sign stays `+1`.",
     "Say which vertex is deleted (the relabelling described is the one for u = 6).")
edit(THM, "- HYP-9274 (B) states `U(W4) = −q(q−1)(q + a_q)`. Together with",
     "- HYP-9274 (B) states the injective identity `T_W4 = −q(q−1)(a_q + 4)`; with the line above it is `U(W4) = −q(q−1)(q + a_q)` "
     "(as recorded in HYP-9275 (A); `U(W4)` is `q(q−1)` times Ahlgren–Ono's moment). Together with",
     "Misattribution: HYP-9274 (B) states the injective sum, not the all-maps U(W4).")
edit(THM, "- HYP-9274 (B) gives `U(W4) = −q(q−1)(q + a_q)` at odd primes. It is Ahlgren–Ono's moment `Σ_b χ(b)(Σ_c χ(c(c−1)(c−b)))² = −a_p − p` plus a collision correction.\n- By the decomposition above, `S_W4 = −4 − a_p`.",
     "- HYP-9274 (B) states `S_W4 = −4 − a_p` at odd primes. Its proof: `U(W4) = q(q−1)·Σ_b χ(b)(Σ_c χ(c(c−1)(c−b)))²` exactly, "
     "Ahlgren–Ono evaluate the moment as `−a_p − p`, and the collision correction `U(W4) = T_W4 − q(q−1)(q−4)` converts this to the injective sum.",
     "Same misattribution; U(W4) is exactly q(q-1) times the moment, the collision correction is the U -> T step.")
edit(THM, "There the wheel identity itself is FINITE-EXACT, not cited.",
     "There the wheel identity itself is FINITE-EXACT, not cited; in this Frobenius form it was already FINITE-EXACT, by direct evaluation, "
     "at 9, 25, 27, 49, 81, 121, 125, 169 (`wheel_apery_extension_20261014h.md` §4).",
     "Credit the earlier prime-power tests of the Frobenius-trace wheel identity (q = 27 was already checked there).")
edit(THM, "They count orders 3 and 4 in `k`-th power Paley digraphs. We found no published order-5 or order-6 count.",
     "They count orders 3 and 4 in `k`-th power Paley digraphs. Both values are older facts for the classes containing `P_q`: "
     "`TT3 = n(n−1)(n−3)/8` holds for every regular tournament, and the `TT4` value for every doubly regular one (part 7). "
     "We found no published order-5 or order-6 count.",
     "TT3/TT4 of P_q are special cases of classical regular / doubly regular counts.")
edit(THM, "- **Earlier repo work (2026-10-14).** The vanishing and reduction calculus is THM-4652 (opus-2026-10-14-S34).",
     "- **Earlier repo work (2026-10-14).** The vanishing laws are THM-4652 (opus-2026-10-14-S34). The degree-≤2 injective elimination "
     "calculus is codex-tiling's `wheel_run_interface_20261014h.md`; the ten-edge terminal classification, including the survival of "
     "`K6 − P6` and the vanishing of the other ten-edge six-vertex terminals, is tiling-modular-atoms' `paley_next_cores_20261014i.md`.",
     "Credit correction (reduction calculus and the ten-edge part of the k = 6 survivor census).")
edit(THM, "  - The PGL(2) lemma and the resulting proof of HYP-9275's octahedron identity.",
     "  - The PGL(2) counting lemma and the resulting proof of HYP-9275's octahedron identity. The invariance behind it is classical: "
     "with an even in-degree orientation, `F(x) = χ(∏_{(i→j)} det(x_i, x_j))` in homogeneous coordinates, which is well defined on "
     "`P¹(F_q)` because every degree is even and `GL(2)`-invariant because `|E|` is even (the triangle case is the oriented Paley "
     "two-graph on `P¹(F_q)`; see e.g. Gunderson–Semeraro, arXiv:2204.10775, and Babai–Cameron on switching classes of tournaments).",
     "Novelty overclaim: the Moebius invariance of even-degree chi-products is classical; the new part is the counting consequence and its use.")
edit(THM, "- **Part 7.** It is elementary. The class-count relations are standard for regular tournaments.",
     "- **Part 7.** It is elementary and classical in substance: for regular `T`, `TT4 = n·C(k,3) − Σ_v c3(N⁺(v))` and "
     "`c4 = n(n²−1)(n−3)/48 − Σ_v c3(N⁺(v))`, so the bound is the known lower bound `c4(T) ≥ c4(DR_n)` on 4-cycles of regular "
     "tournaments, with the same equality case (see Savchenko, arXiv:2403.07629, §1, and the references there). The class-count "
     "relations are standard for regular tournaments.",
     "The TT4 lemma is equivalent to a classical 4-cycle bound; cite it.")

# ------------------------------------------------------------------ HYP-9390
edit(HYP, "with |b(q)|^2 the Sym^2 Frobenius trace minus q^2. Given THM-4780 and HYP-9274 (B), this is equivalent at primes to HYP-9276.",
     "with |b(q)|^2 = q^2 minus the Sym^2 Frobenius trace (at primes, |b(q)|^2 = q^2 - b(q^2)). Given THM-4780 and HYP-9274 (B), at "
     "primes p = 3 (mod 4) this is equivalent to HYP-9276 (1) at those primes.",
     "Sign error: |b(q)|^2 = q^2 - S_q (e.g. q = 27: S_27 = -1207, |b(27)|^2 = 1936); and HYP-9390 says nothing at q = 1 mod 4.")
edit(HYP, "a_q^2 and CM forms of Q(sqrt-2), Q(sqrt-3) and",
     "a_q^2 and CM forms of Q(sqrt-2) and",
     "The Q(sqrt-3) columns of tt6_search.py ignore the unit ambiguity and are not the CM newforms, so that negative result is not established.")
edit(HYP, "  (b) FINITE-EXACT at q = 27, 243, 343, 1331, 2187 through the census; PROVED for primes (HYP-9274 (B)).",
     "  (b) FINITE-EXACT at q = 27, 243, 343, 1331, 2187 through the census (and earlier, by direct evaluation, at 9, 25, 27, 49, 81,\n"
     "  121, 125, 169 in wheel_apery_extension_20261014h.md section 4); PROVED for primes (HYP-9274 (B)).\n"
     "  Audited 2026-10-25 (pythagorean_hourglass_20261025/audit/REPORT.md): all FINITE-EXACT values reproduced with independent code\n"
     "  and PARI's mf eigenforms; the Sym^2 reading of |b(q)|^2 corrected (sign).",
     "Credit the earlier prime-power wheel tests; record the audit.")
edit(HYP, "These are the Frobenius traces; `−b(q²)` here is HYP-9276's `Sym²` trace minus `q²`.",
     "These are the Frobenius traces: `|b(q)|² = q² − S_q`, where `S_q = α_b^{2r} + (α_bβ_b)^r + β_b^{2r}` is HYP-9276's `Sym²` trace "
     "(`S_p = b(p²)`, so `|b(p)|² − p² = −b(p²)` at primes). For example `S_27 = −1207` and `|b(27)|² = 1936 = 44²`.",
     "The sentence as written is false: -b(q^2) = -S_q, while S_q - q^2 = -|b(q)|^2.")
edit(HYP, "  - CM forms of `Q(√−2)`, `Q(√−3)` and `Q(√−7)` of weights 3 and 5 give no identity (`tt6_search.out`).",
     "  - CM forms of `Q(√−2)` and `Q(√−7)` of weights 3 and 5 give no identity (`tt6_search.out`). The `Q(√−3)` columns there take "
     "the trace from the first representation `4p = a² + 3b²` and so ignore the unit ambiguity: the weight-3 column gives `−13`, "
     "`−37` at `p = 7, 19`, against `2`, `26` for the newform `η(2z)³η(6z)³`, and the weight-5 column has the same defect. So that "
     "case is untested.",
     "Unsupported negative result (the Q(sqrt-3) columns are not the CM newforms).")
edit(HYP, "At `q = 27, 243, 343, 1331, 2187` it holds in Frobenius-trace form: the census values satisfy THM-4780 (6) with `a_q = α^r + β^r` (`paley_prime_powers.out`).",
     "At `q = 27, 243, 343, 1331, 2187` it holds in Frobenius-trace form: the census values satisfy THM-4780 (6) with `a_q = α^r + β^r` "
     "(`paley_prime_powers.out`). The same Frobenius form was already FINITE-EXACT, by direct evaluation of the wheel, at 9, 25, 27, 49, "
     "81, 121, 125, 169 (`wheel_apery_extension_20261014h.md` §4); the census adds 243, 343, 1331, 2187. Here `a_q` is the Frobenius "
     "trace, not the Fourier coefficient `a(p^r)` (they differ already at `q = 9`, and `a(27) = 152 ≠ 260`).",
     "Credit, and make explicit that a_q is not the Fourier coefficient at prime powers.")
edit(HYP, "- One of them, the 4-regular graph with 14 automorphisms, reduces by THM-4780 (4) to six-vertex sums.",
     "- One of them, the 4-regular graph `K7 − C7` with 14 automorphisms, reduces by THM-4780 (4) to `±(q − 6)·T_{K6−P6}`.",
     "Precision: the reduction is to the K6 - P6 sum itself.")
edit(HYP, "- The other 37 are genuine. Their signed multiplicities `m_H` range over `{±1, ±3, ±4, ±6, ±8, ±12, ±16, ±24, −32, −40, −48}` (`census_level_forecast.out`).",
     "- The other 37 are genuine. Their multiplicities satisfy `|m_H| ∈ {1, 3, 4, 6, 8, 12, 16, 24, 32, 40, 48}`. The signs depend on "
     "the labelling of each `H`; in geng's labellings (`census_level_forecast.out`) the signed values are "
     "`{+1, +3, ±4, +6, ±8, ±12, ±16, ±24, −32, −40, −48}`.",
     "-1, -3, -6 do not occur; signs are labelling conventions.")
edit(HYP, "- If the data lie in that span, the 37 sums cancel in `TT7`. If not, seven vertices carry new arithmetic.",
     "- If the data lie in that span, so does the genuine combination `Σ m_H T_H` (it is not zero: `2^21·TT7(P_7) = 0` while the known "
     "part is `−25536`). If not, seven vertices carry new arithmetic.",
     "The 37 sums do not cancel; the test is membership of their combination in the span.")
edit(HYP, "1. **Whole censuses.** Every 5-vertex class count of `P_q` lies in `Z[1/2][q] + Z[1/2][q]·S_W4`. Every 6-vertex class count lies in `Z[1/2][q] + Z[1/2][q]·S_W4 + Z[1/2]·S_{K6−P6}`, by the same reduction.",
     "1. **Whole censuses.** Every 5-vertex class count of `P_q` lies in `Q[q] + Q[q]·S_W4`, and every 6-vertex class count in "
     "`Q[q] + Q[q]·S_W4 + Q·q(q−1)·S_{K6−P6}`, by the same reduction. Denominators divide `2^{C(k,2)}·|Aut τ|`, so `Z[1/2]` is right "
     "only for labelled counts: already the 3-cycle count is `q(q²−1)/24`, and the regular 5-tournament count has polynomial part "
     "with denominator `5·2^10`.",
     "False as stated: class counts have odd denominators |Aut tau|; the K6-P6 coefficient is c q(q-1), not a constant.")
edit(HYP, "THM-4780 makes this the same event as a failure of HYP-9276 (or, at prime powers, of (b)).",
     "THM-4780 makes this the same event as a failure of HYP-9276 (1) at that prime (or, at prime powers, of (b) or of HYP-9276 (2)).",
     "At prime powers both conjectural inputs can fail.")

# ------------------------------------------------------------------ results note
edit(RES, "Every `.py` there has its `.out` beside it.",
     "Every `.py` there has an output file beside it except the TT5 search-history scripts `tt5_exact_poly.py`, "
     "`tt5_residual_fit.py`, `tt5_scaling.py` and `tt5_modular_search.py` (the last one's output is the first line of "
     "`tt5_exact_f8.out`, which imports it); `paley_tt5_arith.py` and `paley_tt6.py` write `paley_tt5_data.out` and `paley_tt6_data.out`.",
     "Inaccurate reproduction claim.")
edit(RES, "`TT6` = the level-32 adjoint; **THM-4780**, **HYP-9390** | §4 |",
     "`TT6` = the level-32 adjoint (conditional on HYP-9276); **THM-4780**, **HYP-9390** | §4 |",
     "Mark the conditional statement.")
edit(RES, "   - The potential `c` is the Lorentz-invariant height: the descent is a walk in `O(2,1; Z)` with a positive potential.",
     "   - The potential `c` is the time coordinate for the Lorentz form `a² + b² − c²` that the Berggren matrices preserve (the form is "
     "invariant, `c` is not): the descent is a walk in `O(2,1; Z)` with a positive potential.",
     "c is not Lorentz-invariant.")
edit(RES, "   - Collatz has no such potential: no quadratic form is preserved, which is the separation of `ternary_berggren_20260925.md`.",
     "   - Collatz has no such potential. The proved separation in `ternary_berggren_20260925.md` §6 is narrower: `J = tr²/det` is an "
     "integer on Berggren words and never on accelerated Collatz words, so no fixed projective conjugacy identifies them; that note "
     "says explicitly that it does not distinguish convergence behaviour.",
     "Misquotes the cited result (a single affine Collatz word, rescaled, does preserve an indefinite binary form).")
edit(RES, "The mean tracks `(6/π²)(log √c)²`, the Yao–Knuth law for sums of partial quotients.",
     "The mean is close to `(6/π²)(log √c)²` (the Yao–Knuth leading term for sums of partial quotients) only at small `c`: the ratio "
     "is 0.98, 0.92, 0.88 in the three windows and about 0.86, 0.82, 0.80 near `c = 1.5·10⁶, 1.5·10⁸, 1.5·10¹⁰` (audit, "
     "`audit_berggren.out`). The leading constant is not identified.",
     "Not supported by the data; the ratio keeps falling.")
edit(RES, "   - In a regular tournament, `TT4 ≥ n(n−1)(n−3)(n−7)/64`, with equality iff the tournament is doubly regular.",
     "   - In a regular tournament, `TT4 ≥ n(n−1)(n−3)(n−7)/64`, with equality iff the tournament is doubly regular. This is the "
     "classical lower bound on 4-cycles of regular tournaments in another form (`c4 = n(n²−1)(n−3)/48 − Σ_v c3(N⁺(v))`; Savchenko, "
     "arXiv:2403.07629, §1).",
     "Cite the classical equivalent.")
edit(RES, "   - Put `p_i → p_j` iff `(p_j/p_i) = 1`. This is a tournament by reciprocity.",
     "   - Put `p_i → p_j` iff `(p_j/p_i) = 1`. This is a tournament by reciprocity. It is the converse of the orientation in "
     "`bsd_selmer_growth_20261013a.md` §3 (there `i → j` means `[p_j/p_i] = 1`); the values below are invariant under reversal, but "
     "the labels source + C3 and sink + C3 of the two middle witnesses swap.",
     "Orientation convention differs from the cited note.")
edit(RES, "- THM-4652's calculus eliminates vertices of degree `≤ 2`. Its vanishing laws then kill",
     "- The 2026-10-14 elimination calculus (`wheel_run_interface_20261014h.md`, `paley_next_cores_20261014i.md`) removes vertices of "
     "degree `≤ 2`. THM-4652's vanishing laws then kill",
     "Credit correction.")
edit(RES, "- A new PGL(2) lemma reduces the dense graphs whose degrees are all even.",
     "- A PGL(2) counting lemma reduces the dense graphs whose degrees are all even. Its invariance is the classical Möbius invariance "
     "of `χ(∏ det(x_i, x_j))` over even-degree graphs on `P¹(F_q)`.",
     "Novelty overclaim.")
edit(RES, "the adjoint of `32.3.c.a` (THM-4653 + HYP-9276) |",
     "the adjoint of `32.3.c.a` (THM-4653 + HYP-9276; conditional, HYP-9276 is OPEN) |",
     "Mark the conditional statement.")
edit(RES, "every six-vertex class count is affine in `(a_q, |b(q)|²)`. So the four-loop term of HYP-9371 is a rational function of `q` plus a rational function times `a_q`, and the five-loop term adds `|b(q)|²`.",
     "every six-vertex class count is affine in `(a_q, S_{K6−P6})`, hence in `(a_q, |b(q)|²)` if HYP-9276 holds. So the four-loop "
     "term of HYP-9371 is a rational function of `q` plus a rational function times `a_q`, and the five-loop term adds `S_{K6−P6}` "
     "(conditionally `|b(q)|²`).",
     "Mark the conditional statement.")
edit(RES, "the CM forms of `Q(√−2)`, `Q(√−3)`, `Q(√−7)` (weights 3 and 5),",
     "the CM forms of `Q(√−2)`, `Q(√−7)` (weights 3 and 5; the `Q(√−3)` columns of `tt6_search.py` ignore the unit ambiguity of "
     "`4p = a² + 3b²` and are not the CM newforms, e.g. `−37` instead of `26` at `p = 19` in weight 3, so that case is untested),",
     "Unsupported negative result.")
edit(RES, "- The `ℓ¹` norm is the hostile example: triangle inequality without CS.",
     "- The `ℓ¹` norm is the hostile example: it satisfies the triangle inequality but comes from no inner product (it fails the "
     "parallelogram law), so the equivalence has no content there. With the dot product, `|x·y| ≤ ‖x‖₁‖y‖₁` does hold.",
     "As written ('without CS') it is false for the dot product.")
edit(RES, "Around a sign change of `L`, the two branches of `√L` are exchanged, and the double cover is a Möbius band over a loop.",
     "Around a simple zero of `L` (a loop in the complex `t`-plane; on the real line the two real branches meet at the zero instead), "
     "the two branches of `√L` are exchanged, and the associated real line bundle over the loop is a Möbius band.",
     "Real/complex conflation: on the real line nothing is exchanged at a sign change.")
edit(RES, "- Sanchez-Flores (Graphs Combin. 1996) and Neiman–Mackey–Heule (arXiv:2011.00683), ST27.",
     "- Sanchez-Flores, Graphs Combin. 10 (1994) 367–376 (`v(6) = 28`, ST27) and 14 (1998) 181–200 (uniqueness), and "
     "Neiman–Mackey–Heule (arXiv:2011.00683), ST27.",
     "Wrong year; the two Sanchez-Flores papers are 1994 and 1998.")

# ------------------------------------------------------------------ validate
ok = True
by_file = {}
for e in E:
    by_file.setdefault(e['file'], []).append(e)
for f, edits in by_file.items():
    text = open(os.path.join(ROOT, f), encoding='utf-8').read()
    for e in edits:
        c = text.count(e['old'])
        if c != 1:
            ok = False; print(f'COUNT {c} != 1 in {f}: {e["old"][:80]!r}')
    cur = text
    for e in edits:
        if cur.count(e['old']) != 1:
            ok = False; print(f'SEQUENTIAL FAIL in {f}: {e["old"][:80]!r}')
        else:
            cur = cur.replace(e['old'], e['new'])
    print(f'{f}: {len(edits)} edits, sequential application {"OK" if ok else "FAILED"}; length {len(text)} -> {len(cur)}')
json.dump(E, open(os.path.join(os.path.dirname(__file__), 'required_edits.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('total edits', len(E), 'valid' if ok else 'INVALID')
sys.exit(0 if ok else 1)

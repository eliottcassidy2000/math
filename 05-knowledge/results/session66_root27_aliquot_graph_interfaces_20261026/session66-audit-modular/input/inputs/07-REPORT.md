# Session 65: the distinct fifteens, a Frobenius defect, and the ternary kernel

**Date:** 2026-10-26. **Repository base:** `6a8ad5d3bf01a232da4de41bcbff68390a9337ee`.
**Scope/status:** recovered audited repository results; primary-source CITED facts; elementary proofs below; FINITE-EXACT controls. The new synthesis is **pending independent hostile audit**, not canon. Sloane and Erdős remain OPEN. No machine correctness theorem is inferred from testing or a wiki listing. No repository files were edited and no Git publication was performed by this worker.

## Findings first

1. The fifteen Ogg primes form an exact *prime set*, not just a count. Their supersingular rationality condition equals a zero-genus condition and their set equals the prime support of the Monster's order.
2. Sloane's **15 is an exponent cutoff**, not a cardinality. His original 1973 motivation was ternary multiplicative persistence. The conjecture forbids digit **0** in no further powers after exponent 15. Erdős forbids the absence of digit **2** after exponent 8. The historical **15-state** Erdős machine is a third, encoding-dependent fifteen.
3. A quantitative tournament connection survives adversarial typing: if Frobenius has `b` transposed pairs on the supersingular `j` set, its minimum orientation-preservation defect over **all** tournaments on that set is `b = g(X_0^+(p))`. This gives an exact obstruction interpretation, **not a canonical arithmetic tournament** or a new explanation of why the list has size fifteen.
4. The natural common ternary carrier is the tower `(Z/3^k)^×`, with cyclic kernel of order three. Every finite depth of both digit sieves survives: explicitly `n_k=3^(k−1)−1` produces an all-ones tail of length `k`. The compatible inverse-limit exponent is **(even parity, −1 in Z_3)**, not a nonnegative integer exponent and not the ordinary integer `−1`.

## 1. Inheritance and corrections

Read `AGENTS.md`, `START-HERE.md`, the topic-relevant frontier, `RESEARCH-PROTOCOL.md`, `ACTIVE-GUARDRAILS.md`, and `CORE-PAPERS.md`. Files copied into `inputs/` are pinned to the base above; manifest hashes identify the actual bytes.

- `05-knowledge/results/ogg_triangular_chowla_20260927.md`, §2, already derives/rechecks the class-number and Paley half-row reading. Its audit corrects the difference between arithmetic infinity statements and finite checks.
- `05-knowledge/results/supersingular_divisor_shapes_qseries_20261022.md`, §2: total versus prime-field counts and three differently bounded computational routes. **MISTAKE-666** prevents saying all three routes ran through 3000, conflating prime with prime-power level, or declaring an asymptotic split true at every prime.
- `05-knowledge/results/session32_square_monster_20261023/paley/PALEY_MONSTER_BRIDGE.md`: authentic shared quotient `N_M(<7A>)/He ≅ AGL(1,7)=Aut±(P7)`, ordinary subgroup `7:3`. Conjugation on `C7` alone does not supply translations. This is much stronger than equal cardinalities, but discards the Held kernel and is not an action of the entire Monster on seven points.
- The adjacent `monster/MONSTER_BABY_PRIME_PARTITION.md` retains level structure: at `p=29,41,59,71`, all supersingular `j` are rational but there is one conjugate pair of coarse level-two points. Forgetting that level erases the obstruction. **MISTAKE-668** distinguishes order ratio from subgroup index.
- `05-knowledge/results/session52_saturation_residues_metagraph_20261024/fifteen-kernels/README.md`: fifteen matchings ↔ fifteen binary symplectic Lagrangians; six-pentad outer transport; Bolza weighted kernel orbits `6,4,4,1`. Kernel count, edge classes, and codomain classes differ.
- The adjacent `mrb-supersingular/RESEARCH.md`, §3, supplies explicit Ogg and non-Ogg primes in **every unit residue modulo 15**. Thus no nontrivial sound mod-15 residue rejection keeps every Ogg prime. It also distinguishes modular **level** 15 from reduction **modulo** 15.
- **THM-4754**, `three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md`, is audited PROVED and already supplies the signed 3-adic tournament used below. **THM-4580**, `zeroless-powers-of-two-the-5-adic-tree-and-verification-to-1-1e11.md`, is the decimal analogue and explicitly credits Saye's ternary lift.
- Session63 `cancellation/RESULT.md` keeps odd-action orbit polynomials, coefficient characteristic, and legal toggles separate. A root multiplicity is not a symmetry-kernel dimension.
- **THM-4770** and **THM-4771** are `RESERVED` empty cryptid/odd-orbit stubs at this base; neither is a proved dependency. No lemma below relies on their titles.

Closest mechanism: involution pair-orbits, with the ordinary tournament automorphism odd-order obstruction. Canonical hostile: a single swapped pair. Corrected near miss: a quotient forgetting an entire level-two Frobenius pair. Least-used sidecars: actual Frobenius permutation, carry digit, and exponent's parity coordinate.

## 2. Ogg's fifteen, exactly

Write

`O={2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}`.

For a prime `p`, let `S_p` be the **set of geometric supersingular j-invariants**, `s_p=|S_p|`, and `r_p=|S_p∩F_p|`. This is not a set of all elliptic curves, a stack-weighted mass, or the primes of supersingular reduction for one fixed characteristic-zero curve. Frobenius `F:j↦j^p` acts on `S_p` as an involution.

**CITED primary theorem.** Aricheta, *Supersingular Elliptic Curves and Moonshine*, SIGMA 15 (2019), 007, [Theorem 1.1 / 3.1](https://arxiv.org/html/1809.07421v2), states `Q_p(X)/2=g(X^p)−g(X)` at primes not dividing the level. Taking `X=X_0(1)` gives

`b_p := (s_p−r_p)/2 = g(X_0^+(p))`.

Here `Q` counts nonrational points, not pairs. At level one the theorem has no excluded prime. Ogg's classification gives `b_p=0 iff p∈O`, and this is exactly the prime support of the Monster's order. The source explains the moonshine implication from Monster prime divisors to genus zero while distinguishing it from a full conceptual explanation of the absence of extra primes. At `p=2,3`, `S_p={0}` and both sides are zero. The classical class-number formulas used in our check are restricted to `p>3`.

For `p>3`, the inherited classical computational formulas are

`s_p=floor(p/12)+epsilon_p`, with `epsilon=(0,1,1,2)` for `p=(1,5,7,11) mod12`,

`r_p=(h(−4p)+[p=3 mod4]h(−p))/2`.

Our reduced-form census checks every prime through 211 and finds precisely `O` as the zero-defect set in that range. This is not a new proof of the unbounded classification.

## 3. Elementary theorem: the full involution-defect polynomial

Let `σ` be an involution of a finite **labelled** set `V`. Suppose it has `f` fixed points and `b` transpositions, so `n=f+2b`. For a tournament `T` on `V`, define

`D_σ(T)=#{ {x,y}⊂V : the orientation of (σx,σy) differs from that of (x,y) }`.

Each unordered pair is counted **once**, with an ordered representative used only to compare orientations. Reversing that representative does not change the test. Define `m=b(f+b−1)`.

### Theorem (elementary; independent audit pending)

`Σ_T z^{D_σ(T)} = 2^[C(f,2)+b+m] z^b(1+z²)^m`.

Consequently:

- `min_T D_σ(T)=b`;
- every value is congruent to `b mod2`, and all `b,b+2,…,b+2m` occur;
- exactly `2^[C(f,2)+b+m]` orientations attain the minimum;
- a σ-invariant tournament exists iff σ is the identity.

**Proof.** The action on unordered pairs has three types. A pair of fixed vertices is fixed pointwise and permits two orientations, both contributing zero. Each of the `b` transposed vertex-pairs is fixed setwise but reversed; both orientations contribute one. Every remaining pair-orbit has two unordered pairs. Choosing orientations on its two edges gives four possibilities: two are transported consistently (defect zero), and two inconsistently (defect two). There are

`m=[C(n,2)−C(f,2)−b]/2=b(f+b−1)`

such orbits. Multiplying their independent factors proves the formula. There are no loops, ties, or assumed transitivity. QED.

**Automorphism boundary.** If a finite tournament had an ordinary automorphism group of even order, Cauchy's theorem would supply a nonidentity involution, contradicting the last conclusion. Thus ordinary automorphism groups are odd; anti-automorphisms are a different type and can be involutions. This is the parity obstruction used by the repository's signed-action and odd-kernel lanes, not a novelty claim.

### Genus specialization and hostile

Apply the theorem to the canonical `F` on `S_p`. Together with the cited geometric identity,

`min_T D_F(T)=g(X_0^+(p))`.

For `p=37`, explicit arithmetic in `F_37[t]/(t²−2)` gives

`S_37={8, 3+10t, 3−10t}`.

Frobenius fixes `8`, exchanges the other two, and the supersingular polynomial is `(j−8)(j²−6j+31)`. Thus `f=b=m=1`, genus and minimum defect are one, and the polynomial is `4z+4z³`. **No** orientation makes Frobenius an ordinary automorphism. The checker derives these points from the Legendre Hasse polynomial, checks Frobenius by field exponentiation, and independently counts trace-zero curves for every `j∈F_37`, obtaining only `j=8`.

For `p=11`, both supersingular j-invariants are rational, so both tournaments on that two-point set are Frobenius-invariant. This illustrates how little the orientation itself says: Frobenius is simply the identity.

### Transport contract

Source: the arithmetic `S_p` with its canonical Frobenius action. Target: the complete unoriented pair set with that action; then the **space of all** orientation assignments. Map: forget arithmetic labels but retain the permutation, and measure orientation failures. Preserved: fixed-point/pair decomposition and exact minimum defect, hence the cited genus. Lost: isogenies, moduli geometry, prime, and all arithmetic beyond the permutation. Sidecar: field/prime, full labelled Frobenius action, and theorem identifying its pair count with genus. Hostile: `p=37`; forgetting the action or using ordered pairs changes the invariant. This is an optimized obstruction recoding, **not** an intrinsic preferred binary comparison between supersingular j-values, a preferred tournament, or a new moonshine explanation.

## 4. Sloane, Erdős, and the encoding-dependent fifteen

### Original sources and native arithmetic

Sloane's own [OEIS A102483](https://oeis.org/A102483) lists exponents `0,1,2,3,4,15`, credits the conjecture to 1973, and still marks it open. His author-hosted [1973 article, page 98](https://neilsloane.com/doc/persistence.2.jpg) was downloaded and visually inspected. It frames the zero-digit assertion through multiplicative persistence in base three. We rely on independent exact arithmetic for the intermediate example:

`(222222222222222)_3 → 2^15=(1122221122)_3 → 2^6=(2101)_3 → 0`.

The six digits equal to two in `1122221122` give `64`, not `32`. The source's displayed intermediate is not relied on. Sloane's conjecture would imply that the maximum ternary multiplicative persistence is three, and the chain attains three. This is a **one-way consequence**, not an asserted equivalence.

Erdős's different conjecture is that the only exponents whose ternary powers omit **2** are `0,2,8`; indeed `2^8=(100111)_3`. Saye's [2022 paper](https://arxiv.org/abs/2202.13256v2) verifies both conjectures through exponent `2·3^45`, with a proof-backed recursive search; that is a cited computation, **not rerun here**. The direct current check through 1000 recovers both exceptional lists. Omitting digit one instead has exponents `1,3,9` in this range: for `n≥1`, doubling a `{0,1}`-digit word produces a `{0,2}`-digit word without carries, giving the exact shift equivalence between omission of two at `n−1` and omission of one at `n`.

### Explicit machines and freshness (retrieved 2026-10-26)

Stérin–Woods, [*Hardness of busy beaver value BB(15)*, v2](https://arxiv.org/abs/2107.12475v2), Theorems 1 and 2, prove Erdős-counterexample machines with **15 states/2 symbols** and **5 states/4 symbols**. The explicit [15-state transition file](https://docs.bbchallenge.org/other/powers_of_two_15_2.txt) is available. Their historical smallest-known statements are not current universal minimality claims.

The [Cryptids page](https://wiki.bbchallenge.org/wiki/Cryptids) presently lists a **14-state** Erdős machine and the following **6-state/4-symbol Sloane** machine, both dated October 2026 and **not independently verified**:

```
2RB2RB1RA3RA_0LC1RB3RB2RA_0LE1LD2LC3LC_0RB1LD2LD3LD_1RD2RD3RD0LF_1RE2RZ------
```

Its syntax has 21 nonhalting instructions, one halting instruction, and two undefined entries. Thus `BBi(21)` is not a 21-state statement. No independent correctness proof for this Sloane encoding was located in the targeted searches; the wiki claim is recorded, not promoted. A bounded run would not establish the claimed iff. This encoding can be checked formally later by a macro invariant; it cannot be identified with the supersingular prime set.

## 5. Exact ternary kernel, signed tournament, and realization failure

Set `u_k=2·3^(k−1)`. The following elementary lift is also Saye's Lemma 3.1 in the [v2 PDF](https://arxiv.org/pdf/2202.13256v2), inherited by THM-4580.

### Primitive-root and kernel lemma

For every `k≥1`,

`2^{u_k} ≡ 1+3^k (mod3^{k+1})`, and `ord_(3^k)(2)=u_k`.

**Proof.** At `k=1`, `4=1+3`. Cubing `1+3^k mod3^{k+1}` gives `1+3^{k+1} mod3^{k+2}`, including `k=1`. Since any order modulo three is even, and each such nontrivial first-order lift multiplies the order by three, induction gives the exact order. QED.

Consequently exponent reduction is an isomorphism

`Z/u_k Z → (Z/3^k)^×, n↦2^n`,

and the kernel from level `k+1` to `k` is the `C3` generated by `2^{u_k}`. If a lift is `n+a u_k`, its next digit is

`d_{k+1}(2^{n+a u_k}) ≡ d_{k+1}(2^n)+a(2^n mod3) (mod3)`.

The three possibilities are all three digits, so removing one digit retains exactly two lifts. At depth `k`, Sloane's zero-free tails number `2^k`; Erdős's two-free tails number `2^(k−1)`, because the units digit must be one in the latter. Leading zero-padding counts here: these are **fixed-length trailing words**, not complete integer expansions.

### Native signed tournament, and the lost predicate

THM-4754 defines `S_k` on `Z/3^k`: `x→y` when the first nonzero ternary digit of `y−x` is one. Translation preserves orientation, multiplication by two reverses every arc, and multiplication by four preserves every arc. Thus the full power-of-two clock acts **signed**, while its odd-order square subgroup acts ordinarily. The restriction to units is legitimate. Identifying each lift fibre with `F3` gives a directed triangle, with gauge retained by the exponent increment `+u_k`.

But the digit sieve is not invariant under that kernel: modulo nine the kernel cycle

`1→4→7→1`, with padded digits `01→11→21→01`,

contains two zero-free tails and one forbidden tail. A kernel rotation moves safe to unsafe. The matching of ternary carriers does **not** transport the all-digit avoidance predicate as a group-invariant tournament property. Retain the marked forbidden digit, carry, least exponent, and total digit length.

### Explicit surviving branch at every finite depth

The element `2^{u_k/2}` is the unique order-two unit, hence `−1 mod3^k`. Therefore

`2^{3^(k−1)−1} ≡ −1/2 ≡ (3^k−1)/2 (mod3^k)`.

The right side is the length-`k` ternary all-ones word. Moreover,

`n_(k+1)−n_k=u_k`.

So these exponents are compatible across every depth; adding any multiple of `u_k` produces arbitrarily large ordinary exponents with the same safe tail at that fixed depth. No theorem saying one **fixed** finite number of trailing digits must reveal a forbidden digit can prove either conjecture.

The inverse-limit exponent carrier is `lim Z/u_k Z ≅ Z/2 × Z_3`. Our branch has parity zero and 3-adic coordinate `−1`. Any ordinary integer with the latter coordinate would equal `−1`, which has parity one. Thus **no ordinary integer exponent realizes this entire branch**, even though every finite truncation is realized infinitely often. In particular it is not a counterexample. The obstruction is exactly loss of the ordinary-integer/complete-word realization sidecar.

This is structurally comparable to the MB5 safe-cylinder lane's binary survivor tree and vanishing density; it is not an equality of their Cantor sets, and neither zero density nor nonempty inverse limits settle the sparse integer-orbit problem.

## 6. Other actual fifteens: which maps exist?

| Fifteen | Object and status | Genuine bridge / boundary |
|---|---|---|
| Ogg primes | all-prime arithmetic classification | same actual prime set as Monster prime support; no map to ternary exceptions supplied |
| Sloane cutoff | largest known exceptional exponent, conjecturally final | persistence mechanism; six known exceptions, **not fifteen exceptions** |
| Erdős BB machine | historical proved 15-state encoding | compiles a counterexample search; state count can improve |
| `B4=15` | partitions of a labelled four-set | direct block counts `1+4+3+6+1`; integer partitions of four number five |
| `15` matchings / Richelot kernels | matchings of six labels, Lagrangians of `F2^4` | actual equivariant matching→span map; [Florit–Smith](https://arxiv.org/html/2101.00917v2) supplies geometry with characteristic hypotheses |
| `15` nonzero vectors / Pauli classes | points of the symplectic space | commutation pairing retained; product phases need a sidecar |
| `15` A5 involutions | icosahedral half-turn axes | mark a sixth label, map matching to double transposition fixing its partner |
| Mazur's `15` | full rational elliptic torsion-group types | eleven cyclic and four noncyclic types; ground field is essential; recovered cited classification, no new source audit |
| `15` pentagon families | families admitting monohedral tiling | [Rao's exhaustive-search statement](https://arxiv.org/abs/1708.00274); not fifteen individual pentagons |
| `15` sporadic Apéry-like sequences | Straub 2023 known normalized catalogue | [Straub's stated 15 known cases](https://arxiv.org/html/2301.12248v1); not all conceivable recurrence sequences |

Most instructive hostile: a transposition fixes seven duads but only three matchings. Therefore there is **no ordinary S6-equivariant bijection** between these equal-size point/line sets. A six-pentad outer action repairs the relation only by changing the action sidecar. None of these counts supplies a canonical ordering or a reduction modulo fifteen for the Ogg primes.

## 7. Tests, provenance, and exact stopping boundary

Run:

```
python3 check.py > receipt.json
python3 -O check.py > receipt-optimized.json
cmp receipt.json receipt-optimized.json
```

All checks use explicit runtime conditions, not removable assertions. Executed successfully in normal and optimized modes with byte-identical output.

Finite universes:

- every tournament on `n≤6` for one representative of **every** involution cycle type; full defect distributions, not just minima;
- classical formula census for all primes through 211, with `p=2,3` separately typed;
- all `37²` Legendre parameters and all 37 prime-field j-values via independent point counts;
- powers `2^n` for `0≤n≤1000` with all three digit-omission lists;
- every exponent residue, every three-way lift for depths `1≤k≤8`;
- every vertex-pair in `S_k`, `1≤k≤5`, testing `×2` anti and `×4` ordinary;
- parser-only shape/instruction census of the current Sloane wiki code.

No full search through `2·3^45`, no proof of Sloane/Erdős, and no correctness verification of the recent Sloane or 14-state Erdős machine was performed. The quotient/action connections are completely specified above. The Frobenius-defect and inverse-limit conclusions are ready for independent hostile audit; any canon promotion must wait for that audit.

The root-proposed arbitrary-permutation extension is proved and independently replayed computationally in `GENERAL_PERMUTATION.md`; its minimum is the number of even vertex-cycles, not merely the permutation sign.

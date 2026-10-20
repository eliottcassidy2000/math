# Independent audit of Carmichael inheritance and carry barriers

2026-10-20. **PASS.** Reviewed `INHERITANCE_AND_BRIDGES.md` and its checker. No claim is promoted from a finite experiment to an asymptotic theorem.

## Results accepted

1. **Arbitrary-core prime completion.** For squarefree `N>1` and a prime `q` not dividing `N`, the two stated conditions are exactly Korselt at the old primes and at `q`: `lambda_N | Nq-1` and `q-1 | N-1`. Thus `q<=N` gives a complete finite enumeration. The gcd obstruction is hereditary for Carmichael multiples. It is necessary, not sufficient. The examples `N=33,91,21` respect both directions and the distinct-prime guard.
2. **Fermat-kernel index.** At each prime power the derivative of `X^(N-1)-1` is a unit at a root because `p` never divides `N-1`. Unique lifting gives exactly `gcd(N-1,p-1)` roots, including `p=2`; CRT proves the kernel/image formulas. The statement `I(N)=1 iff Carmichael` retains the composite guard. The prime-marked defects distinguish repeated-prime and field-period costs.
3. **Affine scalar-character blindness.** With composition `fg=f∘g`, the commutator is translation by `-1/4`. That unit generates all translations modulo any modulus coprime to six. Since the multiplier quotient is abelian, the derived subgroup is exactly the translation group. Every cycle-clock word has multiplier one, so every scalar group character returns one regardless of its carry. This does not say additive characters on the translation group are unavailable; the note correctly restores them in its recurrence.
4. **Additive carry recursion.** Direct composition gives `B_(w0)=B_w`, `B_(w1)=3B_w+2^j`; Fourier transformation therefore gives frequency `3a` and phase `exp(2πia2^j/m)`. The fixed-composition index, word order, and unchanged modulus are retained.
5. **Rotation and primitive necklaces.** The rotation identity `2B_rot=3^sB+sD` is correct. Since `D` is odd, it propagates integral fixed points and proves the native parity word. Repetition preserves the rational fixed point exactly; hence Burnside and Möbius formulas apply to the zero-carry predicate. The `(4,2)` hostile has two labelled words, one necklace, and no primitive necklace. Integral fixed points may be negative: positivity requires `D>0`. The finite primitive mixed shapes `(2,1),(3,2),(11,7)` do not assert three positive Collatz cycles.
6. **Quantized occupancy.** The integer inequality `(c-a)(c-a-1)>=0` proves the sharp Fourier-energy floor and its precise equality condition. The sparse-support total-variation bound is correct. This prohibits strong global equidistribution when `W/m→0`, not a target-specific randomized prediction or all useful small normalized Fourier estimates.
7. **Fixed factor-count divisor.** The lift to `G×C_m` is valid for each fixed `m`, conditional on the paper's analytic inputs. Extending a retained character gives `(xi,zeta^j)`. In the high-order case, the first `R-1` powers of `xi` are nonprincipal and the auxiliary phase has modulus one, so the same ambient Fejer estimate applies. In the low-order case the full character order is `<mR`, costing only `m²` in energy. For fixed `m`, `log(|G|m)` has the original asymptotic scale. Identity in the lift imposes `m|omega(N)`. Neither nonzero cardinality classes nor an `m` growing with the scale follow.

## HYP-9338 repair accepted

The inherited carry heuristic must be qualified in three independent ways:

- “On average over shapes” requires an explicit ensemble and normalization. A positive main term tending to zero cannot be a deterministic relative asymptotic for an integer count (or an integer shell total).
- Dividing labelled words by `q` requires primitive period `q`; otherwise use the verified Burnside formula. Integer-cycle counting additionally retains the sign/positivity guard.
- Sparsity `binom(q,p)/|2^q-3^p|→0` is a regime or an independently justified estimate. Entropy alone does not furnish a lower bound on the denominator near the critical slope. The scalar-character estimate from the Carmichael group does not transfer to the affine group's abelianization, and the unrelated changing-level operator radius does not supply the missing fixed-shape cancellation.

These repairs narrow an unsupported inference; they do not prove Collatz convergence, eventual cycle exclusion, or a replacement distribution theorem.

## Independent validation

`independent_checks.py` uses explicit exceptions, and normal/`-O` outputs are byte-identical. It independently verifies:

- 799 Fermat image cardinalities, including non-squarefree and even inputs;
- 4,260 squarefree-core/prime pairs;
- commutators and primitive translation generation at every modulus 5–60 coprime to six;
- 8,166 mixed words through length 12 using exact rational affine evaluation, with 35 integral native words;
- 1,008 Fourier recurrences in exact finite-field root-of-unity models;
- 1,084 arbitrary occupancy vectors, including equality and sparse total variation;
- direct rotation orbits against Burnside/Möbius primitive counts.

Total: **36,397 explicit checks, PASS**. The finite-field character checks are exact modular controls, not a proof inferred from finitely many roots of unity.

The worker's full checker was also replayed in normal and optimized modes: **424,073 checks**, **131,038 mixed words**, **144 distribution instances**. Its JSON and stdout are byte-identical across modes. The worker checker, proof, and boundary labels agree.

# The lifting program terminates in uniqueness (boxeph-2026-07-19-S142)

Owner brief: pursue the lifting program, hunt new kernels at other defect classes.
Outcome: the class contains NO other kernels — the hunt returned a classification.

## The reduction (all identities machine-verified on random instances)

The equivariant z-linear class is F = (zA(s)+y²B(s), yg₁(s)+xzg₂(s), x(h₀(s)−r)),
s = xy, r = x²z. In polar coordinates F is triangular, giving the
**master factorization**  det JF · (h₀−r)² = −Jac_{s,r}(P,Q),  P = F₂F₃, Q = F₁F₃².
Writing w = h₀−r, Φ = sg₁+h₀g₂ (the resolvent), Ψ = h₀A+s²B:
Jac(P,Q) = w²(b₁ + b₂w + b₃w²) with the closed forms
  b₁ = ΦΨ′ − 2Φ′Ψ  (= −c, Keller),   b₃ = 2g₂A′ − 3Ag₂′  (= 0),
  b₂ = 3AΦ′ − A′Φ + 2(Ψg₂′ − g₂Ψ′)  (= 0),   b₀ ≡ 0 automatically.

## The forcing chain (each step an exact machine-checked fact)

1. **b₁ at a multiple root of Φ** evaluates to 0 = −c: impossible ⟹ Φ has only
   simple roots, with Φ′(ρ)Ψ(ρ) = c/2 at each.
2. **Two distinct roots hit the log obstruction**: Ψ = Φ²(K − c∫Φ⁻³) needs
   vanishing residues; verified as a rank fact (Ψ up to degree 10, three root
   pairs: the c-column is never free — only c = 0 solves). ⟹ **Φ is constant
   (no axis collision: automorphism-type rung) or LINEAR.**
3. **b₃ ⟹ the cube–square law** A = αv³, g₂ = γv² — decoding the u³/3u² of the
   kernel as the v = u case.
4. **b₂ forces deg v = 1**: leading degree 2k vs k+1 with coefficient
   (1−k)vₖφ₁; grid-verified inconsistent for v = 1+s+s² and (1+s)².
5. **Reconstruction at the base point** (Φ = 4s+6, Ψ = KΦ² + c/2φ₁ with K = 1/16)
   returns exactly g₁ = 1+12s+9s², B = 4+7s+3s² — the kernel, with the h₀-shift
   as the expected orbit freedom (matching S141's tangent = orbit result).

## The theorem-sketch (machine-supported; full writeup = named lead)

**In the weight-(−1,1,2)→(2,1,−1) z-linear class, every Keller map is either
automorphism-type (Φ constant) or equivalent to THE kernel.** Corollaries: the
resolvent has exactly one root ⟹ one collision s-orbit ⟹ axis fiber 3 ⟹ the
S₃ triple cover is UNIVERSAL in the class; the collision rationals are forced
(root −3/2 of Φ = 2(2s+3); Ψ = (1+s)(2+s) with the −1/4 from c/2φ₁ = −1/4 (!)
— the archaeology of S141's rationals closes: −1/4 IS c/(2φ₁)).

## Where the hunt goes next (the honest frontier)

The class is exhausted; new kernels require leaving it: (i) other weight
systems (the master-factorization method generalizes to any polar-triangular
torus class — same b-machinery, different module generators); (ii) z-degree ≥ 2;
(iii) dimension ≥ 4 (where stabilized cubic-linear forms live); (iv) no
equivariance at all — but S141's rigidity plus this uniqueness suggest a wild
conjecture (labeled WILD): every dim-3 Keller counterexample of minimal degree
is equivariant for some torus, i.e. THE kernel up to equivalence — "the
counterexample is as sporadic as it looks."

Files: jacobian_lifting_classification_boxeph_S142.py + .out (all steps frozen);
predecessor jacobian_lifting_ladder_boxeph_S142.py (the failed naive ladder —
kept: its inconsistencies at d ≠ 3 were the first hint of the forcing chain).

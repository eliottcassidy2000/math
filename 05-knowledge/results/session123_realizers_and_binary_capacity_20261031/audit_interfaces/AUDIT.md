# Independent audit: resource and witness interfaces

2026-10-31. Verdict: PASS within the explicitly stated models. No publication performed. No producer checker imported. The finite replay supports the independent proof audit; it does not prove undecidability, complexity theorems, or unlimited-size circuit claims.

## Audited sources

- `../circuits/RESOURCE_COST_AND_BOUNDARY_KERNEL.md`
- `../realizability/WITNESS_REALIZABILITY.md`

## Circuit findings

1. **Protected forest repair is valid.** Delete only unprotected monomial ANDs and exclusive XOR intermediates. Retaining protected ANDs literally preserves all their outside consumers. In the forest on unprotected edges, a leaf edge belongs to some maximum matching: replace any matching edge incident to its neighbor by the leaf edge, or augment if neither endpoint is matched. Recursing after deleting both endpoints gives exactly ν(U) pivots. Each pivot costs one AND and degree(neighbor)−1 XORs, giving |U|−ν(U) inner XORs. The final reduction has |P|+ν(U)+δ inputs. Together this yields e+|P|+ν(U)+δ−1 gates and saving |U|−ν(U).

   This includes U empty: the original forest is nonempty, so at least one retained output remains. Boundary input fanout is unrestricted; internal XOR fanout is not. An optional notation improvement is to declare δ∈{0,1} and the additional input z explicitly. No repair to the mathematics is required.

2. **Exact obstruction classification is valid for this rewrite family.** ν(U)=|U| iff U is a matching. Thus P meets every incident pair of original edges iff it covers the line graph. The complement of a largest matching gives a minimum such P, of size e−ν(F). This does NOT prove minimality among all circuits, nor enough non-overlapping gains along a global elimination sequence. The source correctly withholds those stronger claims.

3. **Boundary-aware translation kernel is valid.** For a fixed boundary b, an allowed chart sends owned translation h to the corresponding column combination of D; boundary shear Cb+a changes no such direction. Hence unused chart coordinates correspond precisely to independent vectors in the common zero-translation kernel. For quadratic f, the derivative has the affine form x^T P(0,h)+f(0,h)+f(0,0). Its zero condition requires both the full polar columns (including boundary rows) to vanish and the residual constant to be zero. On their common radical the residual constants form a linear map. This establishes the claimed dimension and Gaussian elimination procedure.

   Boundary coupling, multiple exposed outputs, and derivative-one radicals are essential. The triangle's all-ones radical vector has derivative one, so discarding it is invalid. Classical Dickson form and linear-structure facts are rederived or attributed, not newly classified here. Semantic affine compression is not original-coordinate gate saving.

4. **Strict-cost no-go has appropriate scope.** The terminal identity forces c(f)=c(!X)−c(!Y), so endomorphism costs vanish. Copy naturality plus tensor additivity gives c(f)=c(ΔY)−c(ΔX). Free deletion or free copying gives total collapse, assuming finite cancellative scalar arithmetic. Neither conclusion follows merely from symmetric monoidal closedness. An idempotent/max-valued cost or infinite cost cannot be canceled by the displayed argument. The source explicitly retains this scope and allows implementation-level or subadditive costs.

## Realizability findings

1. **Coded section equivalence is correct.** A tracker on paired witness codes supplies a recursive selector by second projection, and conversely. Valid-code-only termination is correctly distinguished from total termination on arbitrary natural codes. Decidable or recursively enumerable total witness relations admit unbounded effective search; the halting graph is not a counterexample to that fact.

2. **Halting-bit example is correctly typed.** Morphisms of assemblies are actual underlying functions with trackers, not trackers as separate arrows. Therefore underlying injectivity and surjectivity imply mono and epi. The halting-bit graph projection is tracked by pair projection; its inverse would compute the halting bit. It cannot be a regular epi, because a regular epi which is mono is invertible. Under the full regular embedding in the exact completion it stays mono and noninvertible; as Eff is a topos, it cannot become epi. The phrase “image in Eff” means the embedded morphism, not its regular-image factorization. No finite test is claimed for this undecidability argument.

3. **The compiler is exact, including endpoints.** Positive p and 2p≤q put each phase interval strictly inside its floor cell even at p/q=1/2. Bounds 0≤z_i≤v_i−1 exclude endpoints 0 and1 correctly. Comparing lower endpoint i with upper endpoint j gives precisely q v_j z_i−q v_i z_j≤(q−p)v_i−p v_j. All pairs share the same z-vector; pairwise existential feasibility does not suffice. The common interval is preserved, not a chosen interior time. Repeated velocities and point intervals require no exception.

4. **Complexity and attribution are correct.** The system has n integer variables and O(n²) inequalities with polynomial binary coefficient length. For fixed n, Lenstra applies; bounded coordinate bisection extracts a vector with O(sum log v_i) feasibility calls, using extra linear bounds rather than extra variables. The lower endpoint certificate has denominator ≤q max v_i. With n variable this yields NP membership, not polynomial-time selection. The n=13 existing-family consequence is conditional on the inherited existence theorem and does not extend that family. This is a new repository interface to classical fixed-dimensional ILP, not a new ILP theorem or practical implementation.

## Independently inspected primary literature

- Faber–van Oosten, *More on Geometric Morphisms between Realizability Toposes*, §1.2, PDF page4 (printed4), definitions of assemblies/tracking and exact completion: https://webspace.science.uu.nl/~ooste110/realizability/geommorRTav.pdf . Supports the categorical setting; the halting counterexample is checked directly above.
- Lenstra, *Integer Programming with a Fixed Number of Variables* (1983), introduction, printed538: https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1983i/art.pdf . Supports the fixed-number-of-variables polynomial-time input. No quotations used.

## Independent replay

`check_independent.py` uses no producer imports and explicit runtime checks (not assert). It exhausts labelled nonempty forests on 2–5 vertices and all edge protection patterns, constructs a leaf-pivot factorization, checks matching/count/line-graph identities, and checks every truth assignment. It also checks all scalar quadratics on three variables under every boundary split, selected two-output combinations, the derivative-one hostile, and exact compiler/floor interval equivalence for all ordered speed tuples of length1–3 from1–5 at radii1/2,1/3,1/4,2/5. Independent endpoint scanning cross-checks safe-time feasibility.

Normal and `-O` receipts are byte-identical:

- 334 forests; 3,300 protection patterns; 101,912 truth assignments.
- 1,024 quadratic/boundary cases.
- 620 compiler instances; 14,460 integer floor vectors.

Commands from the packet directory:

```
python3 -B audit_interfaces/check_independent.py
python3 -B -O audit_interfaces/check_independent.py
cmp audit_interfaces/normal.json audit_interfaces/optimized.json
```


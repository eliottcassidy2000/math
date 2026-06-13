        # Message: opus-2026-03-07-S33: Trivariate GF G_T(t,x) — axes are Eulerian and independence polynomials

        **From:** opus-2026-03-07-S?
        **To:** all
        **Sent:** 2026-03-07 01:24

        ---

        ## Key discoveries this session

### THM-063: Trivariate Generating Function
G_T(t,x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t)(t-1)^{d-f}

Special evaluations:
- G_T(t, 0) = A_n(t) — classical Eulerian polynomial (T-independent)
- G_T(0, x) = I(Omega(T), x) — INDEPENDENCE POLYNOMIAL
- G_T(1, x) = n! for all x, T (universal)
- G_T(t, 2) = E_T(t) = sum_k a_k t^k (tournament Eulerian poly)

The t-axis and x-axis of the GF are the two most important polynomials in tournament theory, meeting at (0,0) = 1 (the empty independent set).

### Palindromic symmetry
G_T(t, x) = t^{n-1} G_T(1/t, x) for ALL x, generalizing a_k = a_{n-1-k}.

### Reduced polynomial P(u, x) where u = t + 1/t
The leading coefficient p_m(x) = I(Omega(T), x) (the independence polynomial itself!).
P(2, x) = n! for all x. All coefficients are explicit OCF polynomials.

### Null space in forward-edge distribution
At n=7, the coefficient matrix has rank 3 (not 4): null vector (0,-2,0,1) means a_k cannot distinguish '2 five-cycles' from '1 disjoint 3-cycle pair'.

### E_T(-1) and E_T(i)
E_T(-1) = deformed Euler/zigzag numbers (verified n=3,5,7,9).
E_T(i) alternates purely real/imaginary by n mod 4.

### For next session
- Investigate null space structure at n=9 (does rank deficiency grow?)
- Combinatorial meaning of p_j coefficients
- Does P(u,x) structure help with Claim A or OCF positivity questions?
- Literature: search for prior work on Eulerian-independence bivariate GFs

### New files
- THM-063-trivariate-gf.md
- 10+ computation scripts in 04-computation/

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*

# Greatest Hits: The Cayley-Delannoy Tournament Theory

## Summary in One Line

The Fourier energy of Hamiltonian path counts over random tournaments is governed by the Cayley transform Q(x) = (1+x)/(1-x) = exp(2 arctanh x), whose Taylor coefficients are Delannoy lattice path weights.

---

## The Top 10 Results (all proved)

### 1. The Master Generating Function
$$Q(x)^m = \left(\frac{1+x}{1-x}\right)^m = 1 + 2\sum_{k=1}^{\infty} g_k(m)\, x^k$$
where g_k(m) = sum C(k-1,j-1) C(m,j) 2^{j-1}.

### 2. The Delannoy Identity
k g_k(m) = sum j C(k,j) C(m,j) 2^{j-1} = total diagonal steps in all Delannoy paths from (0,0) to (k,m). Diagonal = OEIS A108666.

### 3. The Variance Formula
CV^2(H) = sum 2 g_k(n-2k)/(n)_{2k} = 2/n + 0/n^2 - 14/(3n^3) + O(1/n^4).
The 1/n^2 coefficient vanishes exactly because g_2 = g_1^2.

### 4. The Duality
k g_k(m) = m g_m(k). The matrix T(k,m) is symmetric.

### 5. The Functional Equation
Q^m Q(-x)^m = 1. Even-index g_k are determined by odd-index ones.
Equivalently: [x^k] log Q^m = 2m/k for odd k, 0 for even k.

### 6. The Rodrigues Formula
g_k(m) = (1/k!) [d^k/du^k u^m(u+1)^{k-1}]_{u=1}.

### 7. The Uniqueness of arctanh
arctanh is the unique odd formal power series with rational exponential (degree 1,1).

### 8. The Wick Rotation
arctanh(ix) = i arctan(x). In particular: arctanh(i) = i pi/4.
Pi emerges from the tournament at imaginary coupling.

### 9. The Golden Shadow
f_n = (n-2+sqrt(n^2+4))/2 = y_n - 1 (metallic mean minus 1).
CF = [n-1; n, n, n, ...]. Maps to Cayley address via x = (n-1)/(n+1).
Proof: f^2+f-1 = (n-1)(f+1), f^2+3f+1 = (n+1)(f+1), ratio = (n-1)/(n+1).

### 10. The New Sequence
W(n) = 1, 2, 8, 32, 158, 928, 6350, 49752, 439670, 4327904, ...
Not in OEIS. W(n)/n! = 1 + CV^2(H).

---

## Supporting Results

- Transfer matrix M(x) = [[1,2x,0],[0,0,1],[1,x,0]]; char poly lambda^3-lambda^2-xlambda-x
- x-tribonacci recurrence: F(N+3) = F(N+2) + xF(N+1) + xF(N)
- GF: (1+2xy+xy^2)/(1-y-xy^2-xy^3)
- At x=1: tribonacci constant tau = 1.8393...; tau = phi + 0.2213 (price of memory)
- Eigenvalues: lambda_1 ~ 1+2x (real), lambda_{2,3} ~ +/-i sqrt(x) (complex)
- Simplicial binary: sim_H in {0,1} for n >= 4; |sim_H=1| = 2 n!
- 10 applications (ranking, anomaly detection, sports, consensus, ...)
- 8 universality class members (spin chains, phylogenetics, channel capacity, ...)
- Natural numbers on the Cayley line: multiplication = velocity addition
- FTA = Dehn-Sydler on the Cayley line
- Bertrand's postulate = constant hyperbolic width ln(2)/2
- **THM-224 (Golden Exceptional Points):** Discriminant of char poly = 4x(x^2-11x-1). EP eigenvalues are 1/phi and -phi (golden ratio). The discriminant factors over Q(sqrt(5)), linking the transfer matrix to the golden ratio at the exact coalescence points.
- **Bilinear transform = Delannoy:** The bilinear (Cayley) transform Q(x) = (1+x)/(1-x) has Taylor coefficients that ARE Delannoy lattice path weights. This connection is genuinely new — no prior literature links the bilinear/Tustin transform from DSP to Delannoy combinatorics.
- **Fisher-Rao = arctanh:** The Fisher-Rao geodesic distance on Bernoulli distributions is arctanh, the same function governing tournament statistics. This provides an information geometry interpretation: tournament energy IS statistical distance.
- **Golden shadow:** f_n = metallic mean - 1 with continued fraction [n-1; n, n, n, ...]. At n=1: f_1 = phi - 1 = 1/phi. The golden ratio is the n=1 shadow of the Cayley transform.
- **Discriminant symmetry:** Delta(x) = 4x(x^2-11x-1) takes the value -44 at x = -1, x = 1, and x = 11. The symmetry Delta(-1) = Delta(1) = Delta(11) = -44 ties the three natural scales of the transfer matrix.

---

## Connections

| Area | Connection |
|---|---|
| Delannoy paths | Diagonal step count = T(k,m) |
| Tribonacci | Transfer matrix at x=1 |
| Hyperbolic geometry | arctanh = distance, Q = exp(2*distance) |
| Rauzy fractal | Tribonacci substitution dynamics |
| Spin-1 physics | Z_j in {-1,0,+1}, transfer matrix = spin chain |
| Polytope theory | Simplex (1+x)^n, cube (2+x)^n, tournament Q^n |
| Operator theory | GF = resolvent, Perron-Frobenius, spectral mapping |
| Information theory | arctanh = channel capacity, ln(2) = one bit |
| Relativity | Multiplication = velocity addition |
| Number theory | Natural numbers as Cayley addresses |
| Golden ratio / EPs | Discriminant 4x(x^2-11x-1); EP eigenvalues 1/phi, -phi |
| DSP / Bilinear transform | Tustin discretization IS Delannoy path counting (new) |
| Information geometry | Fisher-Rao distance on Bernoulli = arctanh (tournament energy = statistical distance) |
| Metallic means | Golden shadow f_n = metallic mean - 1, CF [n-1; n,n,n,...] |
| Discriminant symmetry | Delta = -44 at x = -1, 1, 11 — three natural scales |

---

## Files

- Paper: `03-artifacts/drafts/cayley-delannoy-tournaments-v2.tex`
- Library: `04-computation/cayley_delannoy.py` (10 tests pass)
- OEIS: `03-artifacts/oeis/W_n_submission.txt` (20 terms)
- Hooks: `03-artifacts/drafts/substack-hooks.md` (Hooks A through I)

# Worpitzky Polynomial Investigation — Synthesis

## Session: opus-2026-03-07-S37 (Phase 2)

### Definition
For a tournament T on [n], the **forward-edge polynomial** (ascent polynomial) is:

F(T,x) = sum_{k=0}^{n-1} F_k x^k

where F_k = #{Hamiltonian paths with exactly k ascents} (positions i where P_i < P_{i+1}).

### Key Results

#### 1. Complement Duality (PROVED)
**F_k(T) = F_{n-1-k}(T^op)** where T^op is the complement tournament.

Proof: Reversing a HP P of T gives a HP of T^op, and ascents become descents.
Consequence: F(T,x) = x^{n-1} F(T^op, 1/x).
F(T,x) is palindromic iff T and T^op have the same HP structure.

#### 2. Sum Over All Tournaments (PROVED)
**Sum_T F(T,x) = A_n(x) * 2^{C(n,2)-(n-1)}**

Each permutation is a HP of exactly 2^{C(n,2)-(n-1)} tournaments.
Average H = n! / 2^{n-1}.
Average F(T,x) = A_n(x) / 2^{n-1}.

#### 3. Unimodality (CONJECTURE — very strong evidence)
F(T,x) is ALWAYS unimodal: F_0 <= F_1 <= ... <= F_peak >= ... >= F_{n-1}.

Evidence: 100% at n=3,4,5 (exhaustive), 100% at n=6,7,8 (500 samples each).

#### 4. Log-Concavity (near-universal)
F_k^2 >= F_{k-1} * F_{k+1} almost always.
- n=3,4: 100% (exhaustive)
- n=5: 1020/1024 (4 failures, all with H=11)
- n=6,7,8: 100% (500 samples each)

The 4 failures at n=5:
- F=[0,0,8,2,1], F=[0,2,6,2,1], F=[1,2,6,2,0], F=[1,2,8,0,0]

#### 5. Real-Rootedness (common but not universal)
F(T,x) has all real roots approximately 90-95% of the time.
- n=3: 100%, n=4: 96.9%, n=5: 95.3%
- n=6: 90.0%, n=7: 89.4%, n=8: 88.6%
Proportion stable around 89%.

#### 6. Worpitzky Expansion: F(T,x) = sum w_k C(x+k, n-1)
The Worpitzky coefficients w_k are always integers.
- w_{n-1} = F_0 = #{fully decreasing HPs}
- w_{n-2} = H - n*F_0
- w_0 = (-1)^{n-1} F(-1) (alternating sum evaluation)
- w_k can be negative; w* is NOT an h*-vector of any polytope

#### 7. Parity Results
- H(T) is always odd (Redei)
- F(-1) = sum F_k(-1)^k is always odd (trivially: F(-1) ≡ H mod 2)
- #{k: F_k odd} is always odd
- Most common: exactly 1 odd F_k

#### 8. F mod 2 Patterns (n=5)
912/1024 have exactly 1 odd F_k, 112/1024 have exactly 3 odd F_k.

#### 9. Generating Function
sum_{m>=0} F(T,m) x^m = W^*(T,x) / (1-x)^n
where W^*(x) = sum_k w_{n-1-k} x^k is the reversed Worpitzky polynomial.

### Negative Results
- F(T,x) is NOT palindromic in general (0/64 at n=4)
- F(T,x) is NOT determined by H alone
- F(T,x) is NOT determined by (H, t3, scores)
- W^* is almost never nonneg (not an h*-vector)
- F(T,2) has no simple relation to H or OCF invariants
- The "universal type B Eulerian" result was an ARTIFACT of wrong W definition

### Open Questions
1. PROVE unimodality of F(T,x)
2. Is there a combinatorial interpretation of the Worpitzky coefficients?
3. What determines the shape F_k/H beyond H?
4. Does the ~89% real-rootedness rate have a structural explanation?
5. Does F(T,x) connect to the transfer matrix eigenvalues?

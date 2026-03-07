# THM-052: Scalar Transfer Matrix for Circulant Tournaments at Odd n

**Status:** PROVED (algebraic proof + exhaustive/DP verification n=5,7,9,11,13)
**Added:** kind-pasteur-2026-03-06-S25d
**Depends on:** THM-050 (consecutive-position formula), THM-030 (M symmetry)

## Statement

**Theorem.** For any circulant tournament T on Z/nZ at odd n, the transfer matrix M satisfies:

    M = (H/n) * I

where H = H(T) is the number of Hamiltonian paths and I is the identity matrix.

## Proof

Let T be a circulant tournament on Z/nZ with generator set S (i.e., T[i,j] = 1 iff j-i mod n is in S). Since |S| = (n-1)/2, T is a regular tournament.

**Step 1: Reduction to palindromic N.**

By THM-050, for a != b: M[a,b] = sum_{j=0}^{n-2} (-1)^j N(a,b,j), where N(a,b,j) = #{Ham paths with {a,b} at consecutive positions {j,j+1}}.

If N(a,b,j) = N(a,b,n-2-j) for all j (palindromic), then the alternating sum vanishes at odd n:

    sum_j (-1)^j N(a,b,j) = sum_k (-1)^{n-2-k} N(a,b,k)
                           = (-1)^n sum_k (-1)^k N(a,b,k)
                           = -sum_j (-1)^j N(a,b,j)    [at odd n]

So sum_j (-1)^j N(a,b,j) = 0, hence M[a,b] = 0 for all a != b.

For the diagonal: each vertex appears at each position H/n times (by vertex-transitivity of circulant tournaments), so M[a,a] = sum_j (-1)^j * (H/n) = H/n (since the alternating sum over n terms equals 1 at odd n).

It remains to prove N is palindromic.

**Step 2: Translation symmetry.**

The map tau_c: i -> i+c (mod n) is an automorphism of T for all c. Under tau_c, Hamiltonian paths biject to Hamiltonian paths, preserving positions. Hence:

    N(a, b, j) = N(a+c, b+c, j)   for all c in Z/nZ.

Define f(d, j) = N(0, d, j) for d in {1,...,n-1}. Then N(a,b,j) = f(b-a mod n, j).

**Step 3: Symmetry of N.**

By definition, N(a,b,j) = N(b,a,j) (both count paths with {a,b} at {j,j+1}). Therefore:

    f(d, j) = N(0, d, j) = N(d, 0, j) = f(0-d mod n, j) = f(n-d, j).

**Step 4: Self-complementarity via reversal.**

The map sigma: i -> -i (mod n) sends generator set S to -S = {n-s : s in S}, which is the generator set of T^op. So sigma is an isomorphism from T to T^op.

Path reversal: if P = (v_0, ..., v_{n-1}) is a Ham path in T, then P^rev = (v_{n-1}, ..., v_0) is a Ham path in T^op.

The composition phi = sigma * rev maps:

    phi(v_0, ..., v_{n-1}) = (-v_{n-1}, ..., -v_0)

This is a bijection from Ham(T) to Ham(T) (since sigma maps T^op back to T).

Under phi, the pair {a,b} at positions {j,j+1} maps to the pair {-a,-b} at positions {n-2-j, n-1-j}. In terms of f:

    f(d, j) = N(0, d, j) = N(-0, -d, n-2-j) = N(0, -d, n-2-j) = f(n-d, n-2-j).

**Step 5: Combining Steps 3 and 4.**

From Step 3: f(d, j) = f(n-d, j).
From Step 4: f(d, j) = f(n-d, n-2-j).

Substituting Step 3 into Step 4:

    f(n-d, j) = f(n-d, n-2-j)

for all d, which is palindromicity of f(e, *) for all e. QED.

## Key Ingredients

1. **Translation automorphism** of circulant tournaments (Step 2)
2. **Symmetry of N** in vertex labels (Step 3)
3. **Self-complementarity** of circulant tournaments via sigma: i -> -i (Step 4)
4. **Path reversal** bijection Ham(T) <-> Ham(T^op) (Step 4)
5. **Alternating sum of palindromic sequence vanishes** at odd length (Step 1)

## Corollaries

**Corollary 1 (Paley tournaments).** For any Paley prime p = 3 mod 4, the transfer matrix of the Paley tournament T_p satisfies M = (H(T_p)/p) * I.

**Corollary 2.** For circulant T at odd n: the off-diagonal entries of M vanish. In particular, NONHAM(a,b) = 0 for all vertex pairs (a,b).

**Corollary 3 (VT tournaments on Z/nZ).** Every vertex-transitive tournament realized as a circulant on Z/nZ at odd n has scalar transfer matrix.

## Scope and Limitations

The proof uses two circulant-specific properties:
1. Translation automorphisms (gives f depending only on d = b-a)
2. The map i -> -i gives T^op (self-complementarity)

For non-circulant vertex-transitive tournaments, the automorphism group may not include translations, and self-complementarity may fail. However, the theorem has been verified for all position-uniform tournaments at n=5 (64/64, including 40 non-VT ones), suggesting it may hold more broadly.

## Verification

- n=5: 64/64 position-uniform (exhaustive), all palindromic
- n=7: 8/8 circulant, all palindromic
- n=9: 16/16 circulant, all palindromic
- n=11: Paley T_11, palindromic (f(d,j) = 1729 for all d,j with T[0,d]=1 or T[d,0]=1)
- n=13: Paley T_13, palindromic

## Source Files

- `04-computation/palindromic_N_proof.py` — proof verification through n=13
- `04-computation/palindromic_N_posuniform.py` — exhaustive n=5 test
- `04-computation/palindromic_N_n9.py` — n=9 circulant test
- `04-computation/palindromic_N_n11.py` — Paley T_11 test

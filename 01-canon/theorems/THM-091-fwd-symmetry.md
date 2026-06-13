# THM-091: Forward-Edge Distribution is Symmetric

**Status:** PROVED (algebraic, 3 lines)
**Proved by:** opus-2026-03-07-S46c
**Scope:** All tournaments, all n

---

## Statement

For any tournament T on n vertices, the forward-edge count fwd(sigma) = #{i : sigma(i) -> sigma(i+1) in T} has a distribution symmetric about (n-1)/2 over uniformly random sigma in S_n.

Equivalently: for all k, the number of permutations with fwd = k equals the number with fwd = n-1-k.

Equivalently: ALL odd central moments of fwd are zero.

Note: This is STRONGER than saying E[fwd] = (n-1)/2 (which is the first moment). It says ALL odd cumulants vanish: kappa_3 = kappa_5 = ... = 0.

Note: This does NOT mean F(T,x) is palindromic. F_k counts permutations with exactly k forward steps among CONSECUTIVE pairs in sigma, while the symmetry is about fwd(sigma) + fwd(sigma^rev) = n-1. The F-polynomial is palindromic only when T is self-complementary.

---

## Proof

For sigma in S_n, define sigma^rev(i) = sigma(n-1-i). Then:

fwd(sigma) + fwd(sigma^rev) = n-1

Proof: sigma has n-1 consecutive pairs (sigma(i), sigma(i+1)) for i=0,...,n-2. Reversing sigma maps pair (sigma(i), sigma(i+1)) to (sigma(n-1-i), sigma(n-2-i)). These are the same pair of vertices in opposite order. Since T is a tournament, exactly one direction is an arc, so exactly one of the two permutations gets a "forward" edge from this pair.

Since sigma -> sigma^rev is a bijection on S_n, fwd and (n-1)-fwd have the same distribution. QED.

---

## Corollaries

1. **Zero skewness (THM-090):** kappa_3 = 0 for all tournaments. This immediately gives E[fwd^3] = 3*mu*E[fwd^2] - 2*mu^3 = A(n) + 6*t3/n.

2. **Zero odd cumulants:** kappa_{2k+1} = 0 for all k >= 0 and all tournaments. The distribution is fully characterized by its even cumulants kappa_2, kappa_4, kappa_6, ...

3. **Cumulant hierarchy:**
   - kappa_2 depends on t3 only (THM-089)
   - kappa_4 depends on (t3, t5) at n=5, on (t3, t5, alpha_2) at n=6
   - Each higher even cumulant adds one new cycle complexity level

4. **Moment-Worpitzky bridge:** Since the Worpitzky coefficient c_j is a linear combination of moments M_0, ..., M_{n-1-j}, and each moment M_r is determined by the cumulants kappa_2, kappa_4, ..., the Worpitzky coefficients inherit the graded cycle structure from the cumulant hierarchy.

---

## Verified

- kappa_3 = 0: n=4 (3 classes), n=5 (8 classes), n=6 (24 classes) exhaustive
- Scripts: fwd_cumulants.py

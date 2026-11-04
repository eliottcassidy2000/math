# Complete modular certificates for prescribed odd-perfect supports

## Result and limits

This is a proved **fixed-support decision certificate**, not an odd-perfect-number nonexistence proof. It combines the S179 primitive-divisor domain theorem with exact one-row-per-prime selection. Its added content is an explicit finite modular formulation that is complete, preserves integer choices, and has a directly replayable recurrence. General generating-function and modular dynamic-programming methods are standard; historical novelty and computational superiority are not claimed.

The extension from a bounded prime support to all prime supports remains unproved. In particular, neither checking many supports nor finding stronger certificates on each prescribed support closes that quantifier gap.

## 1. Exact finite rows inherited from S179

Let S be a finite set of s >= 2 odd primes and q in S with q = 1 mod 4. For p in S let D_p = {ord_r(p): r in S, r != p}. Retain lengths n > 1 such that every divisor d > 1 of n belongs to D_p, and n is odd for p != q or n = 2 mod 4 for p=q. There are at most s-1 such lengths per p.

Why this is complete: for each d dividing n, Bang–Zsigmondy supplies a primitive prime divisor r of p^d-1, hence ord_r(p)=d. The exception d=2, p+1 a power of two cannot occur: non-Euler lengths are odd; Euler p=q is 1 mod 4, and q+1 >= 6 is 2 mod 4. The base-two exception is unavailable. Every such r divides sigma(p^(n-1)), so perfectness forces r in S. Thus each genuine exponent length occurs in this finite list.

For each retained n, compute sigma(p^(n-1)) by its geometric sum, strip 2 and all primes of S, and reject the row unless the remaining cofactor is 1 and its 2-adic valuation equals [p=q]. No general-purpose factorization is needed. Define its signed odd-valuation vector

    c(p,n)_r = v_r(sigma(p^(n-1))) - (n-1)[r=p].

Write R_p for the finite list of these vectors, retaining their length labels. If any R_p is empty, reject immediately. Otherwise an OPN with precisely S and Euler prime q exists **if and only if** choosing one row from each R_p gives sum_p c(p,n_p)=0. Necessity follows by valuations; sufficiency follows from exact closure and the unique factor 2. This is the exact selection equation already present in S179, written in signed coordinates.

## 2. No-wrap modular completeness theorem

For each coordinate r define the integer bounds

    L_r = sum_p min_{c in R_p} c_r,
    U_r = sum_p max_{c in R_p} c_r,
    m_r = 1 + max(|L_r|, |U_r|).

Let G be the finite abelian group product_r Z/m_r Z. Modulus 1 is allowed and has one element. The zero-coordinate count in the group algebra product

    product_{p in S} (sum_{c in R_p} [c mod m])

is **exactly** the number of exponent assignments producing an OPN with support S and Euler prime q.

Proof. Expansion chooses exactly one row per source, and its coefficient is a nonnegative integer counting those choices. Every resulting integer coordinate z_r lies in [L_r,U_r], a subset of (-m_r,m_r). The only multiple of m_r in that open interval is 0. Thus a selection sums to the identity of G if and only if every integer coordinate is zero. The exact-selection equivalence above finishes the proof. No cancellation of coefficients, floating-point tests, or unverified integer-program solver is involved.

An interval with L_r>0 or U_r<0 already rejects the support and can be checked first. Smaller arbitrary moduli are also sound **rejection** filters: a genuine integer zero maps to a modular zero. Their acceptance is inconclusive unless a no-wrap bound or a subsequent exact check is supplied.

## 3. Replayable recurrence

Choose any order p_1,...,p_s. Start F_0(0)=1 and F_0(g)=0 otherwise. For i=1,...,s, set

    F_i(g) = sum_{c in R_{p_i}} F_{i-1}(g-c mod m).

Store only states with positive counts. An independent verifier recomputes all rows, all bounds, and this recurrence. With the displayed moduli, F_s(0)=0 is a complete negative certificate; F_s(0)>0 counts actual solutions. Keeping a predecessor pointer also recovers a witness exponent tuple. Counts must be exact nonnegative integers. Reducing the **counts** modulo a prime could falsely erase a positive solution count and is not permitted.

If T=sum_p |R_p|, the dense recurrence uses O(T product_r m_r) additions and O(product_r m_r) stored counts. Sparse storage can do better but has no guaranteed advantage over the S179 bound (s-1)^s on row tuples. This theorem makes no polynomial-time claim. Small-modulus rejection first is optional, not a completeness assumption.

## 4. Why retaining integrality and the no-wrap bound matters

The one-source abstract catalogue R={(-1),(1)} has a fractional convex combination equal to zero, but no integer choice equals zero. Modulo 2 both rows are 1, so the modular test rejects it. This is an abstract integrality control, **not** a claimed arithmetic support that survives every S179 filter.

Conversely, the singleton catalogue R={(2)} maps to zero modulo 2 despite containing no integer zero. The theorem chooses m=3 and rejects it. This prevents promoting a small-modulus acceptance to an exact arithmetic conclusion.

The existing S179 valuation/Farkas certificate remains useful. Its linear separation and the modular selection test reject different relaxations. Neither is a replacement for proving coverage of all supports.

## 5. Verified scope and provenance

Run `python verify.py` and `python -O verify.py`. The deterministic checker compares the modular recurrence with independent exhaustive integer selection for 1,000 signed catalogues, including 174 positive cases. For 189 arithmetic support/Euler-prime cases, it constructs 124 exact closed rows and checks 116 nonempty arithmetic row prefixes (up to three sources), as well as empty-domain handling. Genuine valuation balances for perfect numbers 6 and 28 are positive controls for the generic compiler. These even-perfect controls do not pass through the odd-support theorem.

The small arithmetic supports all have an empty full row domain, so this experiment does **not** demonstrate an arithmetic improvement over S179. The nontrivial modular-vs-convex distinction is presently the explicitly labelled abstract control.

Primary external sources inspected 2026-11-04:

- Gyory and Smyth, *The divisibility of a^n-b^n by powers of n*, Integers 10 (2010), Theorem 13, https://math.colgate.edu/~integers/k27/k27.pdf , for the classical primitive-divisor theorem used by S179.
- Pace P. Nielsen, *An upper bound for odd perfect numbers*, Integers 3 (2003), https://math.colgate.edu/~integers/d14/d14.pdf , establishes an explicit upper bound conditional on a fixed number of distinct prime factors. Fixed-prime-count finiteness is already classical; the present work claims neither that discovery nor a uniform bound on the number of prime factors.

Internal prerequisites are S177 odd-order support, S178 primeflow, S179 oddperfect/PROOF.md, and S179 crossfield/REPORT.md. No forbidden ledger or canon file was opened. The unresolved major theorem is still: exclude or uniformly control every support, not merely decide one given support.

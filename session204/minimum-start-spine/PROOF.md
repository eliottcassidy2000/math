# Minimum-start ancestry has logarithmically many exceptional junctions

Status: PROVED, independently reviewed. This is an all-order structural consequence, not a coverage or finite-state theorem.

## Inheritance and precise gain

The dependency is the independently audited Session203 bounded-entry theorem, `session203/bounded-entry-grammar/PROOF.md`, Sections 2--3. The local reviewed copy is `session203-bounded-entry-release/session203/bounded-entry-grammar/PROOF.md`. Its exact decomposition and its bounds of 103 eligible vertices and 10 actual entry vertices at start budget 26 are inherited, not new here. The gain is to choose the minimum-start entry rather than an arbitrary entry, producing a multiplicative charge and at most four exceptional stages at that budget.

For a nonempty strong rooted tournament (T,v), let k=s_v(T). If T is not a singleton, order the strong components of T-v as C,D_1,...,D_t. Put A=N_T^+(v) intersect V(C), m=|A|, and q=product_i H(D_i), with empty product 1. The inherited exact identity is

    k=q sum_{a in A} s_a(C).

Here A is nonempty and every summand is positive. All H(D_i) are odd positive integers; a nontrivial strong D_i has H(D_i)>=3. Thus q=1 means all tail components are singletons. The inherited case q=m=1 is exactly a physical J/O/Z macro: add a root beating only the single selected a in C, then add the tail singleton vertices in their ordered positions with their recorded incidences to the new root. The intermediate singleton boundary is handled as a literal graph, not by a transfer formula outside its domain.

## The selected-spine theorem

At each inverse step choose a in A minimizing s_a(C), breaking ties by any fixed ordering of the original vertex labels. Replace (T,v) by (C,a), and repeat until a singleton is reached. This produces an induced-core ancestry spine, since C has strictly fewer vertices than T. Let its successive start counts be k_0,...,k_N=1 and its junction data be (q_i,m_i), 0<=i<N. Then

    k_{i+1} <= floor(k_i/(q_i m_i)),
    product_{i=0}^{N-1} q_i m_i <= k_0.              (1)

Call a stage exceptional when q_i m_i>1. Every other stage is the inherited same-budget J/O/Z macro and satisfies k_{i+1}=k_i. Consequently the number E of exceptional stages satisfies

    E <= floor(log_2 k_0).                           (2)

For k_0<=26, E<=4. The selected spine has at most five stretches of same-budget macros, separated by at most four exceptional junctions.

**Proof.** The least of m_i positive integers is at most their average, so the exact identity gives k_i >= q_i m_i k_{i+1}. Multiplying these inequalities and using k_N=1 proves (1). Each exceptional multiplier is an integer at least two, proving (2). If q_i m_i=1, both positive integer factors are one; the exact identity gives equality of the budgets, and the inherited classification identifies the macro. Conversely every exceptional step strictly decreases the selected budget. Termination follows from strict core-order decrease, not from start-budget decrease. QED.

The bound can be refined when tail complexity is recorded. If E_t stages have q_i>1 and E_m stages have q_i=1,m_i>1, then

    3^{E_t} 2^{E_m} <= k_0.                          (3)

This does not double count a stage with both q_i>1 and m_i>1. In particular at budget 26 there are at most two tail-bearing exceptional stages. More generally the total number B of nontrivial strong tail blocks across the entire selected spine obeys

    3^B <= product_i q_i <= k_0,
    B <= floor(log_3 k_0).                           (4)

The proof multiplies H(D)>=3 over those blocks. They are distinct deleted vertex sets along the nested ancestry, so the total really counts graph pieces, not reused copies.

## What this does and does not compile

The theorem is existential and graph-theoretic for every strong rooted tournament, of arbitrary order. It is not a census assertion. It improves the inherited strict-budget descent to bounded multiplicative junction complexity. It does not bound the lengths of same-budget stretches, the intermediate numerical coordinates, or the positions of all other marked entries in the retained core. The unchosen vertices of A remain present in C and must still be marked for reconstructing the exceptional parent; choosing one minimum-start root does not erase their incidence or joint path-cover data.

In particular a spine is not a rooted branching tree with independently realizable children. The remaining entries are vertices of the same physical core. Replacing their joint responses by arbitrary independent profiles would be invalid.

Applied to the strong source part of the proposed F(8,137) strategy, this bounds exceptional ancestry complexity for all required start budgets at most26. Weak sources still require the inherited ordered-component normalization. There remains an infinite arithmetic image/occupancy problem along the arbitrarily long same-budget words and their marked junctions. No eventual threshold, no new certified target, and no proof of F(8,137) follows.

## Finite controls

`verify.py` independently counts Hamilton starts by subset recursion and SCCs by reachability on every labelled tournament of orders 1 through 6. It checks all 136,743 strong rooted spines, including 340,900 inverse steps, 152,504 exceptional steps and 20,320 tail-bearing steps. Both normal and optimized Python runs agree. The finite census is a diagnostic for the identities and implementation; the all-order theorem rests on the preceding proof. Checks use explicit exceptions rather than Python assertions.

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The reviewed dependency is copied into `dependencies/BOUNDED_ENTRY.md`; the independent mathematical audit is `audit/REVIEW.md`.

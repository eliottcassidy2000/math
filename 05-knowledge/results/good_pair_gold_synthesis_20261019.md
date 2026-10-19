# Gold partitions from good pairs: an attributed synthesis

Date: 2026-10-19. Status: arbitrary finite-poset theorem, with a current-measure weighted version. This combines existing injection and Gold-certificate arguments; no historical priority claim is made. The [primary-source comparison](chain_surplus_primary_source_attribution_20261019.md) records the exact matching hypotheses. Packet88 is unchanged.

## 1. Attribution and convention

Zaguia's Definition1 specifies a good pair by the hypotheses below (possibly in the dual), and Theorem2 states its balanced-pair consequence. His Section2 proves the decreasing insertion-bin law by swaps. [Author-deposited paper, arXiv:1610.00809v3](https://arxiv.org/pdf/1610.00809v3).

Peczarski's 2017 Section2 chooses a half-probability crossing on a chain and uses a swap to obtain a balanced triple; Fact1 turns that certificate into Gold. His structural hypotheses concern full principal upsets, rather than just their difference. [Publisher paper, pp.107–109](https://cdm.ucalgary.ca/article/download/62349/46785/178291).

We give their combination explicitly. Neither accessed theorem states this combined good-pair Gold conclusion. That observation is limited to the checked texts and does not assert literature novelty. The earlier monotone insertion argument also occurs in [Zaguia 2012, pp.4–5](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v19i2p29/pdf/).

Use the weak Gold convention: on every realized path, E≥t1+t2, where E is the initial number of extensions and ti the surviving number after i comparisons. Determined comparisons, including repetition, are allowed. For weights, replace all counts by total weights. A zero-weight path requires no condition.

## 2. Theorem and exact hypotheses

Let P be a finite poset and let a,b be distinct incomparable vertices. Assume

    D(a)⊆D(b),
    T=U(b)\U(a) is a chain, possibly empty,
    P(a<b)≤1/2.

Then P has a weak Gold strategy. If T is nonempty and the first inequality in probability is strict, two comparisons between a and consecutive vertices of the chain {b}∪T suffice. The order of those two comparisons is chosen from the extension counts before either answer is received.

If P or its dual has such a pair, the same conclusion holds. No efficient counting algorithm is asserted.

**Very-good-pair corollary.** If D(a)=D(b) and both U(a)\U(b) and U(b)\U(a) are chains, orient the pair so its forward probability is at most one half and apply the theorem. The dual condition also suffices. Thus every poset containing a very good pair has Gold.

## 3. Full insertion proof

Write c0=b<c1<...<ck for {b}∪T, with k≥0. Every ci is incomparable with a. For i≥1, a<ci would contradict ci∉U(a), while ci<a would imply b<a.

Let qi count the extensions having exactly i members of this chain before a, for 0≤i≤k+1. These cells exhaust all extensions. We show

    q0≥q1≥...≥q_(k+1)≥0, q0≤E/2.

For each i=0,...,k, take an extension in bin i+1 and swap a with ci. All predecessors of a precede b because D(a)⊆D(b), so moving a earlier is valid. Moving ci later preserves its predecessor constraints. Every successor z of ci is either in U(a), and hence after the old a, or in T; in the latter case it is a later marked chain vertex and also after the old a. Thus the swap is valid and injects bin i+1 into bin i. The inequality on q0 is precisely the assumed probability bound.

These are Zaguia's insertion maps, with the common-successor check spelled out. No assumption on independence of exterior positions is used.

## 4. Constructive numerical Gold step

More generally take nonnegative q0≥...≥q_(m−1), E=sum qi>0, m≥2, with q0≤E/2, and assume each nontrivial prefix is an available comparison.

If q0=E/2, compare at the first prefix and then make any comparison, including repeating the first. Both branches have mass E/2, so every final mass is at most E/2 and Gold follows.

Otherwise q0<E/2. Let ℓ be the first index whose inclusive prefix reaches E/2 and set

    a0=sum_(i<ℓ) qi, x=qℓ, b0=sum_(i>ℓ) qi.

Then 1≤ℓ≤m−2, a0<E/2, b0≤E/2, and x≤q0≤a0. The upper index bound follows because the prefix before the last bin is at least E−q0>E/2.

If b0≤a0, use prefix a0 first and prefix a0+x second. The three nonzero-path Gold deficits are

    E−2a0, a0−x, a0−b0.

If a0≤b0, reverse those comparisons. The deficits are

    b0−a0, b0−x, E−2b0.

All are nonnegative; equalities and zero cells are harmless. For the poset, prefix through qj is exactly the comparison a<cj. Thus the selected endpoints are c_(ℓ−1),cℓ. This explicitly realizes the adjacent-threshold/half-crossing Gold mechanism.

If T is empty there are only two bins. Monotonicity and q0≤E/2 force q0=q1=E/2, so the half-case applies. No nonexistent second chain endpoint is needed.

## 5. Weighted and conditional version

Fix S={a,b}∪T. Replace every S-label in an extension by a blank, retaining the outside labels and all occupied positions. Give each extension a nonnegative weight depending only on this outside skeleton, with positive total E. The swaps above preserve that skeleton, so the decreasing-bin law remains true.

The theorem holds under the additional **current-measure** hypothesis that the total weight of a<b is at most E/2. That hypothesis cannot simply be carried over from the original uniform measure after conditioning. One may condition on any positive-weight event determined by the skeleton and then test it anew. More generally, any nonnegative reweighting constant on each skeleton is allowed, with the same current-measure test.

For a very good pair, one can use a common skeleton erasing {a,b} together with both upset differences. Both orientations' swaps preserve it. Choose the orientation whose current probability is at most one half; equal downsets make both orientations structurally eligible. Thus its weighted Gold conclusion holds for every such positive-total weighting. Conditional outside comparisons disjoint from this marked set preserve the permitted form. Arbitrary internal reweighting is not covered.

For the narrower chain-surplus class of packet88, the first two bins are equal under its skeleton weights, so the current-measure probability bound is automatic. See [chain-surplus theorem](chain_surplus_gold_theorem_20261019.md).

## 6. Distinguishing the hypotheses

Take five vertices u,v,s,t,w and precisely the generated relations

    u<s<t, u<w, v<w.

For (a,b)=(v,u), the upset difference is the chain s<t, and the downsets are equal. Yet U(u)={s,t,w} is not a chain: w is incomparable with s,t. Thus the designated pair does not meet Peczarski's full-upset-chain requirement. In fact u,v are the only two vertices with equal lower-cover sets; U(v) is a chain but P(v<u)<1/2. Therefore neither orientation meets his p.109 probability-based relaxation. Strictness follows from the critical-pair injection and the extension u,s,t,v,w, which is outside its image.

This example distinguishes that structural statement only. It is not asserted to escape other established Gold results (indeed it is a five-point poset).

Conversely, two disjoint two-element chains have a very good pair formed by their minima, with two nonempty disjoint upset differences. Those upsets are not nested, so this marked pair falls outside packet88 while satisfying the new synthesis and Peczarski's original structural statement.

## 7. Exact finite controls

The accompanying independent stdlib checker enumerates naturally labelled posets through six vertices, verifies the insertion swaps and skeleton monotonicity, tests every structurally eligible orientation whose current count meets the half condition, and checks the very-good-pair orientation rule. It also checks integer monotone laws without first-two equality and explicit weighted skeleton mixtures. Normal and optimized runs have separate receipts. The all-size proof is Sections3–5; finite controls are regression checks only.

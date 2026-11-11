# Exact finite-width boundary compilers for strict product growth

Date: 2026-11-11. **PROVED:** the all-round height bound, complete late rows
for two independent seeds, and eventual affine interval compression near a
singleton face. **FINITE-EXACT:** the separate small controls. No assertion
that these strict numerical gates are the constructors of the full
tournament Hamiltonian-path spectrum is made.

## 1. Complete rows for two independent seeds

Let a,b≥2 be multiplicatively independent. Identify a^u b^k with (u,k),
and start strict synchronous vector addition from (1,0),(0,1). Retain old
values at every step. Define c_1=0 and, for k≥2,

    c_k=k+min_(1≤i≤k−i)[c_i+c_(k−i)+[2i=k]].          (1)

The beginning is

    k:    1  2  3   4   5   6   7   8   9  10  11  12
    c_k:  0  3  6  10  14  19  23  28  33  39  44  49.

**Theorem.** At every generation t, every live (u,k) with k≥1 satisfies

    u≤kt−c_k.                                        (2)

For every t≥k the entire row is exactly

    {(u,k):1≤u≤kt−c_k},                              (3)

plus the original seed (0,1) when k=1. There are no other pure-axis values.
Thus every fixed strip 1≤k≤K is completely represented after time K by K
integers c_k and the clock t. The statement concerns all later t, not a
finite enumeration or a fit to sampled maxima.

### Universal upper bound

The only zero-grade value is (1,0). Appending it increases u by one and
costs a generation, which respects (2). Retention respects (2) as well.
For two positive grades i+j=k, the inductive bounds at generation t−1
give a total at most `k(t−1)−c_i−c_j`. If i=j, the parents must have
different integer first coordinates, so their sum is at most one less
than twice the common bound. Taking the best split gives exactly (1)–(2).
The seed (0,1) satisfies its bound, including at t=0. Induction proves (2).

### Constructing every point of the upper bound

For k=1, the row is 0≤u≤t. For k=2 and t≥2, two different entries of
the grade-one row at t−1 have every possible sum from 1 through 2t−3.
This proves the base cases, with their actual legal parents.

Now use induction on k≥3 and a split i+j=k attaining (1). At time
t−1≥k−1, both lower-grade rows have their asserted complete intervals.

If i≠j, all cross-row pairs are distinct regardless of their first
coordinate. Their sums fill a complete interval through kt−c_k, starting
at most at 2. The point (1,k) is obtained by the serial b-spine in exactly
k generations, so the whole interval (3) is present.

If i=j≥2, the two copies of the row interval give every strict sum from
3 through kt−c_k. The points (1,k) and (2,k) are present by time k:
the first uses the b-spine; the second is

    (1,k−1)+(1,1)=(2,k),

whose parents have different grades because k≥3. The row interval has at
least two entries: the chain choice in (1) gives c_i≤i(i+1)/2 for i≥2,
and t−1≥2i−1 gives `i(t−1)−c_i≥3`. This proves (3). QED.

### The correction is a real equality obstruction

Dropping the term [2i=k] would give c_2=2, predicting the nonexistent
point (2,2) at generation 2. Its only balanced split would repeat (1,1).
In fact (2,2) first appears at generation 3. Different grades do not pay
this tax; it records equality of the full parent vector, not the parity
of either exponent or a generic loss factor.

The constants also have an exact binary-tree interpretation. Take a full
binary tree with k leaves. At each internal vertex charge its number of
descendant leaves, and add one when its two children have equal leaf
counts. Then c_k is the minimum total charge. This follows by splitting
at the root and applying (1). The height-budget proof shows why this
combinatorial minimum measures a physical construction frontier.

For u>k(k−1)−c_k, the exact first generation of a^u b^k is

    ceil((u+c_k)/k).                                  (4)

The inequality ensures that the right side is at least k, where (3)
constructs the target. Equation (2) excludes every earlier generation.

The exact number of mixed values in a fixed strip 1≤k≤K at t≥K is

    K(K+1)t/2−sum_(k=1)^K c_k.                        (5)

Its size grows linearly even while the unrestricted mixed interior grows
on exponentially larger exponent scales. Formula (5) is a boundary count,
not a density estimate for the entire cone or numeric prefix.

## 2. A general singleton-face compiler

Let S be a finite set of nonzero vectors in N_0^r, with strict vector-sum
growth. Suppose an integer linear functional ell is nonnegative on S and
exactly one seed A has grade zero. It may have negative coefficients in
the ambient coordinates; only nonnegativity on the generated cone matters.
All other seed grades are positive integers. Fix K≥1.

**Theorem.** The portion of the full closure with 0<ell≤K lies in finitely
many cosets of ZA. Choose one representative v_q for each such coset.
Then for each coset there are integers m_q,b_q, a positive integer
1≤s_q≤ell(v_q), and a finite T_q such that for every t≥T_q,

    S_t intersection (v_q+ZA)
      ={v_q+nA:m_q≤n≤s_q t+b_q}.                     (6)

All these quantities admit a finite exact algorithm. Thus any fixed-width
strip has a finite collection of affine moving endpoints, even though the
whole cone need not saturate and the number of strips is unbounded.

### Finite transverse classes and exact transitions

The grade-zero closure stays {A}. Every positive-grade construction uses
at most K leaves from S\{A}; all other leaves are A. Hence its coset is
represented by one of the finitely many sums

    v_mu=sum_(s≠A) mu_s s,
    0<sum mu_s ell(s)≤K.                              (7)

Conversely A+v_mu is attained by the mixed-seed theorem, so every such
coset occurs. Relations between the vectors can identify these classes;
the algorithm computes cosets, rather than treating the seed alphabet as
free. Representatives differing by an integer multiple of A are merged.

Write R_q(t)={n:v_q+nA∈S_t}. Exact updates retain R_q, adjoin R_q+1, and
take pairwise sums from positive-grade cosets. If q+r=s and

    v_q+v_r=v_s+kappa(q,r)A,

the new longitudinal coordinate is u+v+kappa(q,r). When q≠r every pair
is legal. When q=r the pair u=v is deleted. Grades add, so positive
grade k depends only on smaller positive grades and its own +1 update.
These are exact finite-layer transitions with integer carry kappa; no
unavailable higher-grade state can feed a lower one.

Each full-closure row has a least coordinate m_q and is exactly the ray
`m_q+N_0`: append A to its least attainable point. The minimum is finite
and computable directly from (7). For each mu, take n=0 if v_mu is a seed
or mu uses two distinct off-face seed types, and take n=1 otherwise.
The candidate is v_mu+nA; minimize its coordinate in each coset. These
are all possible minima by the mixed-seed theorem. Each candidate uses
at most K+1 leaves, so every m_q is already present by generation K.

### Eventual intervals and affine endpoints

Proceed by increasing positive grade. Once all lower-grade rows are
intervals with affine upper endpoints, each positive-pair birth source
for q is an interval with constant lower endpoint and an affine upper
endpoint of integer slope at least 2. The diagonal case trims both ends
by one; choose a finite starting time when the input rows have length≥2.

Choose that starting time t0≥K and retain the finite exact row R_q(t0). Its +A
continuations have maximum `max R_q(t0)+t−t0`, of slope 1. A pair-source
upper endpoint of slope a≥2, born at time r and then continued by A to
time t, has upper endpoint `ar+b+t−r`. This is maximized by the latest
available r. Thus all future maxima are the maximum of finitely many
explicit affine integer functions, with slopes at most the target grade.
Such a maximum eventually equals the line of largest slope and then
largest intercept among ties.

There are no persistent gaps between these sources. The minimum point
m_q supplies the interval `[m_q,m_q+t−K]`. Each pair-source interval
starts at one of finitely many fixed coordinates, and each retained
transient point supplies its own interval by repeated +A. Eventually
the interval from m_q reaches every such fixed starting coordinate.
Their union is then exactly the full interval to the maximum. This
proves (6), with a computable threshold obtained from the finitely many
line crossings, interval lengths, and starting-coordinate comparisons.

The induction is effective: all transient strip states are finite and
can be generated by the exact coset transitions above. It never requires
the unrestricted cone or its rapidly growing interior.

## 3. Transfer boundary

Source→target: graded vector-construction trees→finite transverse cosets
with one longitudinal integer coordinate. Preserved: complete fixed-strip
membership, strict equality of parents, and synchronous time. Lost: higher
grades and numeric order between unrelated cosets. Required sidecars:
the actual quotient by ZA, representative carries, and the retained clock.

For integer multiplication use prime-exponent vectors and retain their
numeric weight sum_p v_p log p when numeric prefixes are requested. For
the singleton extreme rays of the [classification theorem](../classification/PROOF.md),
an integer exposing functional supplies exactly the hypotheses above.
That theorem needs only the first positive layer; (6) supplies every fixed
finite collection of layers and literal recursive witnesses.

Ordinal joins of tournaments multiply Hamiltonian-path counts, but they
also allow equal-count factors. Consequently they do not preserve the
deleted numerical diagonal used in (1). The correction c_2=3 is therefore
not a lower bound for tournament construction depth, and this compiler
does not prove an omitted path total. A transfer to marked tournament
ports must specify a native grade and actual equality/availability rules;
scalar count equality alone is insufficient.

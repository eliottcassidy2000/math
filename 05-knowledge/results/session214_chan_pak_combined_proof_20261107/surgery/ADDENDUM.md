# Two-component structure and a second all-order omission, 27,371

**Status: PROVED structural reduction plus FINITE-EXACT certificate, pending independent audit.** This extends only the class F in `PROOF.md`. Neither 7,477 nor 27,371 is an ambient forbidden Hamiltonian-path count; explicit two-exit tournaments realize both.

## 1. Two nontrivial components force a nearly transitive factor

A **nontrivial** strong component means one with at least three vertices. A two-vertex tournament cannot be strong. Singleton components contribute H=1. Keep the definition of F: a transitive chain, at most two additional vertices a,b, and at most one chain vertex beaten by a.

**Structural theorem.** A tournament in F has at most two nontrivial strong components. If it has exactly two, they have the form

    B => R_m,

apart from singleton components, where B is a strong one-apex tournament and R_m, m>=3, is the transitive m-vertex order with only its extremal arc reversed. Therefore

    H(T) = H(B) (2^(m-2)+1).

Conversely, the dominance product of any one-apex tournament B and R_m belongs to F. It realizes precisely the displayed product. This statement does **not** say every nonstrong tournament has two nontrivial components: a single strong two-apex core with singleton padding is also possible.

**Proof.** Every directed cycle meets at least one of a,b because deleting both leaves a transitive tournament. Thus every nontrivial strong component contains an apex, and there are at most two. If there are two, one contains a and the other contains b. Each contains at least two chain vertices.

If a's component came first in the dominance order, a would beat every chain vertex in b's component, including at least two such vertices. This contradicts the chain-outdegree restriction. Hence a's component is later. Its chain is transitive and has at least two vertices. Its first chain vertex has no incoming chain arc, so strongness forces an incoming arc from a. This uses a's only allowed chain exit. All other chain vertices beat a. Ordering the chain first and a last exhibits exactly R_m.

Every Hamiltonian word in R_m either avoids its only backward arc, giving the original increasing order, or uses the block (last,first). Each of the m-2 interior vertices can lie before or after that block, increasingly within either run. All 2^(m-2) choices work and are distinct. Thus H(R_m)=1+2^(m-2). The component product law proves the formula.

Conversely, choose B's chain before R_m's chain; the combined chain is transitive. Use B's apex as b and R_m's exceptional last vertex as a. All vertices of B beat a. Within R_m, a beats only its first chain vertex. Thus the dominance product is a presentation in F. QED.

## 2. A proved cutoff for the one-apex boundary

**Lemma.** A strong one-apex tournament with L>=3 chain vertices has

    H >= 1 + 2^(L-2).

**Proof.** Write the apex as b. Strongness forces b->x_0, since otherwise x_0 is a global source. It also forces x_(L-1)->b, since otherwise x_(L-1) is a global sink. For every even-cardinality subset S of the L-2 internal chain vertices, the word

    b, x_0, S in increasing order, x_(L-1), b

is a directed odd cycle of length 3+|S|. The cycles are distinct, because their vertex sets recover S. There are 2^(L-3) even subsets when L>=3. The inherited odd-cycle collection formula

    H(T) = sum over disjoint odd-cycle collections F of 2^|F|

counts the empty collection with weight one and each displayed singleton collection with weight two. All omitted contributions are nonnegative. Thus H>=1+2*2^(L-3), as claimed. The formula is inherited from Session213 `path-gadgets/PROOF.md`, ultimately THM-4876; no new identity is claimed. QED.

Consequently, a strong one-apex core with H<=27,371 has L<=16, since L>=17 would give H>=32,769. Adding two successive sink vertices a,z and placing z after its chain embeds it into the exact-one-exit census with chain length L+1<=17, preserving H. This is the same padding mechanism as the initial packet.

## 3. The second exclusion and the exact initial spectrum

**Theorem.** Through 27,371 the values of H attained by F are exactly the positive odd integers other than

    7, 21, 7,477, 27,371.

Here the nonattainment of 7 and 21 is inherited. The exclusion of 7,477 was proved in the initial packet. The other positive values are finite-exact witnesses. The exclusion of 27,371 is all-order in F.

**Proof of the new exclusion.** Suppose H(T)=27,371 for some T in F. A tournament with no nontrivial strong component is transitive and has H=1. If T has two nontrivial strong components, section 1 forces a divisor 2^k+1 with k>=1. None divides 27,371. For a transparent arithmetic certificate, 27,371=101*271, both factors are prime, and none of 101,271,27,371 is one more than a power of two. Equivalently one can directly test k=1,...,14; larger 2^k+1 exceed 27,371.

There is therefore a single nontrivial strong component S with H(S)=27,371. It is in F by hereditary closure. The **CITED** Busch bound from the initial packet gives

    f(20)=9*5^5=28,125 > 27,371,

so S has order at most 19.

If a presentation of S uses both apexes and a has exactly one chain exit, its chain has length at most 17 and is directly included in the census. If a has no chain exit, strongness requires a->b. Deleting b then leaves the transitive order (chain,a), so S is a one-apex tournament. Presentations with just one apex are also one-apex. Section 2 covers all these one-apex cases with the same chain-length-17 census after count-preserving padding.

The exact census in `addendum_verify.py` scans every presentation with chain length 1 through 17, exactly one a-chain exit, either apex arc and every b-cut. It contains no value 27,371. Thus no possible S exists. The census also supplies a witness for every positive odd integer below 27,371 except 7,21,7,477, proving the positive part. QED.

**Dominance-closure distinction.** Allowing arbitrary finite dominance products of members of F still cannot realize 7,477: primality would require a factor with H=7,477. However, this larger closure **does** realize 27,371, since the positive census provides counts 101 and 271 in F and their dominance product has count 101*271. Thus the second exclusion belongs to F itself, not its full multiplicative closure.

## 4. Physical control with two exits

A chain of length eleven, a beating x_0,x_3, b beating a and precisely x_2,x_5,x_8,x_9, has outgoing bit rows (chain first, then a,b)

    [6142,8188,4088,6128,8160,4032,8064,7936,3584,3072,6144,9,2852].

The direct subset DP and run transfer each give H=27,371. Thus both new class omissions are repaired by allowing just two chain exits for a. No claim is made that this larger class is complete.

## 5. Verification and scope

`addendum_verify.py` reuses the inherited counting mechanisms through the frozen `counting.py`, but has a separate receipt. It scans 8,388,612 parameter presentations, reports the four omissions through 27,371, replays one witness for every other odd value, tests the one-apex bound on all strong one-apex cuts of chain lengths 3 through 12, checks the two-component classification on all parameter presentations through chain length six, and verifies the two-exit control by subset DP.

    python3 addendum_verify.py > ADDENDUM_RESULT.replay.json
    python3 -O addendum_verify.py > ADDENDUM_RESULT.optimized.json
    cmp ADDENDUM_RESULT.json ADDENDUM_RESULT.replay.json
    cmp ADDENDUM_RESULT.json ADDENDUM_RESULT.optimized.json

The proof explicitly separates the two-component case from single-core padding. It does not assume that a composite count must have two nontrivial components, nor that a nonstrong presentation necessarily factors nontrivially. This distinction is essential to the reduction.

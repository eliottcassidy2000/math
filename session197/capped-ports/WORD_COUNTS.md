# Fibonacci word counts and the exact first departure from H=1

Session 197 corollary of the moving-port language. These are counts of construction words, not unlabeled tournament classes. The graph-count and random-time claims are exact for one exposed port. No literature-priority claim is made.

## 1. Three transitive states

A nonempty tournament has H=1 exactly when it is transitive. A nontransitive tournament contains a directed triangle, whose three Hamilton paths inject into paths of every induced extension.

Use the full one-port language from `PROOF.md`: right or left birth, exceptional subset empty or the sole port, and forget the old port or newborn. These are eight equally distinguished letters, even if two letters happen to construct the same graph.

A transitive tournament has a unique order. When its order is at least two, call its port a boundary port if it is first or last, and an interior port otherwise. Let B and I denote those two categories. A singleton with its sole port is a separate initial category U.

The exact numbers of one-letter transitions that remain transitive are

    from B: 4 to B, 2 to I, 2 to rejection;
    from I: 2 to B, 2 to I, 4 to rejection;
    from U: 8 to B.

To check this, an empty exceptional subset adds a universal source or sink. Forgetting the old port exposes that new boundary vertex; retaining the old port moves a boundary port to the interior exactly when the birth occurs on its own boundary side. A nonempty exceptional subset is safe only when a right birth selects the current sink, inserting immediately before it, or a left birth selects the current source, inserting immediately after it. Each safe choice has one boundary and one interior forgetting outcome. In every other nonempty-subset case a directed triangle is formed. The two boundary directions have the same aggregated transition counts. U is exceptional because both ends coincide before its first birth.

## 2. Exact word enumeration

Let b_m and i_m count length-m words that keep H=1 when starting from B or I. The empty word counts once. The transition matrix is

    M = [[4,2],[2,2]],
    (b_m,i_m)^t = M^m (1,1)^t.

Thus both sequences satisfy a_m=6a_(m-1)-4a_(m-2), with

    b_0=i_0=1, b_1=6, i_1=4.

For the Fibonacci convention F_0=0,F_1=1, this gives

    b_m=2^m F_(2m+2),
    i_m=2^m F_(2m+1).

Their ordinary generating functions are

    sum b_m z^m = 1/(1-6z+4z^2),
    sum i_m z^m = (1-2z)/(1-6z+4z^2).

For a singleton start, u_0=1 and u_m=8b_(m-1)=2^(m+2)F_(2m) for m>=1, so its generating function is 1+8z/(1-6z+4z^2).

The boundary/interior distinction is essential: transitive graphs with different exposed positions have different word languages even though their H values are all one. Their unmarked isomorphism class alone is insufficient.

## 3. Exact random first-exit times

Choose each of the eight letters independently and uniformly. Let tau be the first positive number of births at which H>1. Then P(tau>m) is b_m/8^m, i_m/8^m, or u_m/8^m for the corresponding initial category.

After leaving U there is at least a 1/4 chance of rejection at every next letter, so tau is finite almost surely and the survival sum converges. Evaluation of the preceding generating functions at z=1/8 gives

    E[tau | B]=16/5,
    E[tau | I]=12/5,
    E[tau | U]=21/5.

These expectations do not depend on the initial transitive order beyond its category. No assertion about a uniformly random tournament is intended: this is the specified uniform random construction-word process.

## 4. The first new Hamilton count is exactly a binary clone value

Write the old transitive order as v_0,...,v_(n-1), with port v_i. A right birth selecting that port, when it is not the sink, creates exactly one nontrivial SCC:

    v_i -> Tr_d -> z -> v_i,   d=n-1-i>=1.

The preceding vertices remain singleton components. A left birth selecting a nonsource port is the converse configuration, with d=i>=1. No other one-port operation can be the first departure from transitivity.

The strong component is the actual triangle wrapper of a transitive d-vertex module. Its count is

    H=2^d+1.

For a direct proof using odd-cycle collections, every directed cycle uses both outside vertices, so no two cycles are disjoint. For every nonempty odd subset of the d ordered vertices there is exactly one odd directed cycle: the subset is traversed in its increasing order between the two outside vertices. There are 2^(d-1) such subsets. Hence H=1+2*2^(d-1)=2^d+1. All exterior singleton components contribute one.

Conversely, every d>=1 occurs as the first departure from H=1, starting from any nonempty one-port transitive seed. First add a universal sink and expose it. Next make a right selected-port birth and expose the newborn; since the old port was a sink, this safe move makes the new port penultimate. Add d-1 universal sinks while retaining the port. Its distance from the last vertex is now d. A final right selected-port birth creates exactly the stated component, and this is the first nontransitive step.

Thus the complete first-exit count set is exactly {2^d+1:d>=1}, independently of the initial transitive seed.

Two immediate spectrum consequences are sharply scoped:

- No first-exit count is divisible by seven, because powers of two modulo seven lie in {1,2,4}, never -1. In particular no first departure produces a 7p value.
- If a first-exit count is prime, d must be a power of two. Indeed d=2^s t with odd t>1 would make 2^d+1 divisible by 2^(2^s)+1. This is only a necessary condition, not a claim that every Fermat number is prime.

These are first-exit obstructions. Later words can have further nontrivial components or more complicated strong components; no corresponding restriction on their full Hamilton-count spectrum is asserted.

## 5. Sources and verification

The operation conventions and the general finite-machine theorem are in the adjacent `PROOF.md`, frozen SHA256 `77be8657a64155c261aebb00e06b3e3db7114ad764d6b417a51c119ff85e4ad1`. The triangle-cloning value 2^d+1 is inherited from the named S196 cloning family and is rederived here to make the first-exit proof self-contained. The new content of this corollary is the exact construction-word enumeration, its random-time consequence and the complete first-exit characterization for this language.

`verify_word_counts.py` checks all physical one-step transitions from transitive orders through nine, direct Hamilton counts for every first-exit graph in that universe, complete accepted word prefixes through depth five for all three initial categories, the Fibonacci recurrence, the exact rational expectations and literal witnesses for the converse construction. It uses explicit exception checks and no external packages. The all-order results come from the state classification and direct component proof, not extrapolation from the finite controls.

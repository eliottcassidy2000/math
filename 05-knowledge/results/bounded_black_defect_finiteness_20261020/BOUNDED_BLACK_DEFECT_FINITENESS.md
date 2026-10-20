# Bounded black defect forces bounded strong order

Status: PROVED from the CITED Hamilton-cycle lower bound of Thomassen and automorphism-group bound of Dixon. The combinatorial lemmas below are proved directly. Independent finite tests are validation, not a replacement for the all-order proof. Research date: 2026-10-20.

## 1. Statement and conventions

For a tournament T, write A=Aut(T), H(T) for the number of directed Hamiltonian vertex words, and h(T)=H(T)/|A|. The A-action on words is free. A word P=(p_0,...,p_(n-1)) is **blue** when its positional reversal r_P(p_i)=p_(n-1-i) is an anti-automorphism of T. Let b(T) count the A-orbits of blue words and set

    Delta(T)=h(T)-b(T).

Blue is invariant under A because r_(aP)=a r_P a^(-1). Thus there are exactly Delta(T)|A| black Hamiltonian words. A tournament that is not self-converse has no anti-automorphism, so b(T)=0. In the ordinary tiling fibre these are exactly the inherited blue-orbit and black-defect counts.

**Theorem.** For each fixed integer d>=0, there are only finitely many isomorphism classes of strong tournaments T with Delta(T)<=d.

In particular the theorem holds for strong self-converse tournaments. Equivalently, the minimum black defect among strong tournaments of order n tends to infinity as n tends to infinity. The singleton, if included as a strong tournament, contributes just one exceptional class with Delta=0. All cycle arguments below concern n>=3; there is no strong tournament of order two.

Here **strong** always means strongly connected. Hamiltonian cycles are counted up to cyclic rotation, with their directed orientation fixed. Such a cycle is equivalently its successor permutation R, a single n-cycle. No division by reversal is made.

## 2. Exact cycle-cut quotient

Fix a directed Hamiltonian cycle C=(v_0,...,v_(n-1)), with indices modulo n and R(v_j)=v_(j+1). Put

    G=A intersect <R>,    m=|G|,    q=n/m.

A tournament automorphism cannot have even order: its involutory power would transpose a pair and reverse that pair's arc. Hence m is odd. As a subgroup of <R>, G=<R^q>.

Cut C just before v_i to obtain P_i=(v_i,v_(i+1),...,v_(i-1)). An automorphism carrying P_i coordinatewise to P_j is necessarily R^(j-i), because the words specify its value on every vertex. Therefore two cuts have the same A-orbit exactly when their indices are congruent modulo q. The cuts occupy exactly q global path orbits; no automorphism outside <R> creates further identifications.

Their positional reversals are

    r_i(v_j)=v_(2i-1-j),    r_i r_a=R^(2(i-a)).

If no r_i is anti-automorphic, all q cut-orbits are black. Otherwise fix one anti-automorphic r_a. Then r_i is anti-automorphic if and only if R^(2(i-a)) belongs to A, equivalently q divides 2(i-a). Consequently the number beta(C) of blue cut-orbits is exactly gcd(2,q) when nonzero, and zero otherwise. In all cases beta(C)<=2.

The q-beta(C) black cut-orbits are distinct global black path orbits. Thus

    q-beta(C) <= Delta(T),    q <= Delta(T)+2.                  (1)

This argument does not require self-converseness. For a non-self-converse T it gives the stronger q<=Delta(T).

## 3. The missing degree coordinate is restored by rotation orbits

The permutation R^q has exactly q vertex orbits, each of odd length m=n/q. On each orbit O it acts transitively by automorphisms of T[O]. Therefore every vertex of T[O] has the same outdegree within O. Summing these m equal outdegrees gives m(m-1)/2, so each is (m-1)/2.

All outneighbors within O are also outneighbors in T. Consequently

    delta^+(T) >= (m-1)/2 = n/(2q)-1/2.

If Delta(T)<=d and s=d+2, equation (1) yields

    delta^+(T) >= n/(2s)-1/2.                                (2)

Strongness supplies at least one Hamiltonian cycle by Camion's theorem. Applying the preceding argument to any such cycle proves (2) for every nontrivial strong T under consideration. No claim that T itself is regular is needed: only its q induced rotation-orbits are regular. This distinction is essential once the defect is positive.

## 4. Few black words and few square roots bound all cycles

Partition the Hamiltonian cycles into those with at least one black cut and those whose every cut is blue.

### 4.1. Cycles with a black cut

A closing Hamiltonian word determines its directed cycle uniquely by adding the last-to-first edge. Distinct cycles therefore have disjoint sets of cut words. Selecting one black cut from each black-containing cycle is injective into the set of black Hamiltonian words. Hence

    HC_black(T) <= Delta(T)|A|.                              (3)

This is a bound on cycles, not cycle-orbits. It does not count a closing word as several cycles, and it does not assume that every black word closes.

### 4.2. All-blue cycles

If every cut of C is blue, then both r_i and r_(i+1) are anti-automorphisms. Their composition is an automorphism, and

    R^2 = r_(i+1) r_i belongs to A.                          (4)

For a fixed permutation g of an n-element set, the number of full n-cycles R with R^2=g is as follows:

- If n is odd, it is one when g is a full n-cycle, and zero otherwise.
- If n=2m is even, it is m when g consists of two m-cycles, and zero otherwise.

**Proof of the root count.** Squaring an n-cycle produces gcd(2,n) cycles of equal length, so these are the only possible types. Any root R commutes with g=R^2. For odd n, the centralizer of the n-cycle g is <g>; writing R=g^a, the congruence 2a=1 modulo n has exactly one solution, which is coprime to n and hence gives a full cycle.

For even n=2m, let g's cycles be X and Y and fix x in X. A full-cycle root must swap X and Y: these are the even and odd positions in its full orbit. Choose y=R(x) in Y. Commutation and R^2=g force, for every j modulo m,

    R(g^j x)=g^j y,    R(g^j y)=g^(j+1)x.

Conversely these formulas define a full 2m-cycle root for each of the m choices of y. Thus there are exactly m roots. This proof includes n=2 as a permutation statement, although order two is irrelevant for strong tournaments. QED.

In particular there are at most n roots over any g. Since a cycle has one well-defined successor R, (4) and the root count give

    HC_blue(T) <= n|A|.                                    (5)

Not every square root counted here need be a directed Hamiltonian cycle of T or have all-blue cuts; the use is only an upper bound. One may replace n by 1 for odd n and by n/2 for even n, but this sharpening is unnecessary. In a tournament the even case is empty unless n/2 is odd, because R^2 has order n/2.

Combining (3) and (5), for Delta(T)<=d we obtain

    HC(T) <= (n+d)|Aut(T)|.                                (6)

## 5. Contradiction between factorial and exponential growth

We use two published inputs with their precise hypotheses:

- **CITED, Thomassen (1985):** a strong tournament whose minimum outdegree is at least 3k+3 has at least 4^k k! Hamiltonian circuits. Carsten Thomassen, *Hamilton Circuits in Regular Tournaments*, in *Cycles in Graphs*, North-Holland Mathematics Studies 115 / Annals of Discrete Mathematics 27, pp. 159-162. DOI https://doi.org/10.1016/S0304-0208(08)73005-6 ; primary publisher record https://www.sciencedirect.com/science/article/pii/S0304020808730056 . Despite the paper's title, the used statement assumes strongness and a minimum-outdegree bound, not regularity.
- **CITED, Dixon (1967), Theorem 1:** every tournament of order n satisfies |Aut(T)|<=3^((n-1)/2). John D. Dixon, *The Maximum Order of the Group of a Tournament*, Canadian Mathematical Bulletin 10(4), pp. 503-505. DOI https://doi.org/10.4153/CMB-1967-048-9 . The primary PDF and its theorem were rechecked; the argument uses the odd-order solvability theorem of Feit-Thompson.

Fix d and put s=d+2. For n sufficiently large define

    k = floor((n-7s)/(6s)).

Then k>=1 and

    3k+3 <= n/(2s)-1/2 <= delta^+(T).

Thomassen, (6), and Dixon therefore imply

    4^k k! <= HC(T) <= (n+d) 3^((n-1)/2).                   (7)

Here k=n/(6s)+O(1), with s fixed. To see directly that (7) eventually fails, observe that k!>=(k/2)^(k/2) for every integer k>=2. The logarithm of the left side divided by n is at least

    (k/n) log 4 + (k/(2n)) log(k/2),

which tends to infinity because k/n tends to 1/(6s)>0. The logarithm of the right side divided by n is

    log(n+d)/n + ((n-1)/(2n)) log 3,

which tends to (log 3)/2. Thus (7) is impossible for all sufficiently large n, uniformly over T. There are finitely many tournaments of each bounded order, proving the theorem. QED.

No numerical cutoff is claimed from this limit argument. The theorem is effective in principle using the explicit factorial inequality (7), but this packet does not spend a large census or compute a purported sharp cutoff.

## 6. Exact consequences for paired-black leaves and the strong core

### 6.1. The deficit-two center alphabet is finite

Let d_j count strong self-converse isomorphism classes Z of order j with Delta(Z)=2, and write D(z)=sum d_j z^j. The theorem with d=2 proves that **D(z) is a finite polynomial**. It does not determine all its coefficients.

The inherited exact leaf reduction says that every mixed paired-black leaf has the unique form A -> Z -> A^op, where A is a word of rigid strong atoms {1,C3} and Z is one of these deficit-two centers. Therefore its series is now an all-order rational function,

    MixedLeaf(z)=D(z)/(1-z^2-z^6).

The previously proved classification of non-self-converse leaves then gives the complete leaf series

    Leaf(z) = (1/2)[1/(1-z-z^3) - (1+z+z^3)/(1-z^2-z^6)]
              + D(z)/(1-z^2-z^6).                          (8)

It follows that mixed-leaf counts eventually satisfy m_n=m_(n-2)+m_(n-6). Total leaf counts eventually satisfy the recurrence obtained from

    (1-z-z^3)(1-z^2-z^6)=1-z-z^2+z^5-z^6+z^7+z^9,

namely

    l_n=l_(n-1)+l_(n-2)-l_(n-5)+l_(n-6)-l_(n-7)-l_(n-9).

These eventual statements do not specify a recurrence onset because deg D remains uncomputed. Inherited exact checks give d_4=1 and d_j=0 for j<=9, j!=4, hence D(z)=z^4+O(z^10). The claim D(z)=z^4 at all orders is NOT proved.

### 6.2. A center-independent leaf asymptotic

Let a_n count all words in atoms of sizes 1 and 3, so a_n=[z^n]1/(1-z-z^3). Since D is a finite polynomial, (8) has the form

    Leaf(z)=1/[2(1-z-z^3)] + F(z)/(1-z^2-z^6)

with a fixed polynomial F(z)=D(z)-(1+z+z^3)/2. The coefficient of z^(2r) in 1/(1-z^2-z^6) is a_r; odd coefficients vanish. The sequence a_r is nondecreasing by appending a singleton. Therefore, for C equal to the sum of the absolute values of F's coefficients,

    |l_n-a_n/2| <= C a_floor(n/2).

Concatenation gives an injection from a word of weight floor(n/2) and a word of weight ceil(n/2) to a word of weight n. The split is recovered at that fixed weight boundary, so

    a_n >= a_floor(n/2) a_ceil(n/2).

Finally a_(3r)>=2^r by independently choosing each weight-three block to be either the atom 3 or three atoms 1. It follows that a_ceil(n/2) tends to infinity, and hence

    l_n/a_n -> 1/2.

The unknown finite center alphabet affects only a lower-order correction. This is an asymptotic for leaf vertices, not an asymptotic equality for the number of bridges.

### 6.3. Bounded paired-black degree also forces bounded strong order

For a self-converse node represented by strong T, the inherited quotient-degree identity is deg_Q(T)=Delta(T)/2. For a non-self-converse pair {T,T^op}, it is deg_Q=h(T)=Delta(T). Consequently

    Delta(T) <= 2 deg_Q(T)

in both cases. For each fixed integer K>=0, only finitely many strong-core nodes have paired-black degree at most K. Equivalently, the minimum paired-black degree among strong-core nodes tends to infinity with order.

This is a degree constraint, not a theorem about the connectedness, neighbor identities, or bridge topology of the strong core. It does not imply that every bridge is a leaf edge, or classify bridges joining nonleaf regions.

## 7. Inheritance, hostile boundaries, and connection contract

Closest proved mechanism: the cycle-cut quotient in the preceding near-purity obstruction packet. New step: its large cyclic automorphism orbit gives an internal degree lower bound, which activates Thomassen even when the tournament is not regular.

The prior pure-blue finiteness argument used full regularity and a direct injection into automorphisms. The present argument tolerates a fixed number of black path-orbits by splitting cycles into an exceptional word-controlled part and a square-root-controlled part. It does not import the earlier asymptotic Cuckler theorem.

Boundary controls:

- C3: H=3, |Aut|=3, h=b=1, Delta=0.
- Regular cyclic C5: H=15, |Aut|=5, h=b=3, Delta=0. Pure does not mean rigid.
- The strong S4 with outneighbor masks (2,4,9,3): H=5, |Aut|=1, b=3, Delta=2. Positive defect does not force global regularity; its q=4 rotation blocks are singletons.
- Removing strongness destroys the conclusion. Every transitive tournament has one Hamiltonian path and one blue orbit, so arbitrarily large nonstrong tournaments have Delta=0. Rigid wrappers around S4 likewise produce unbounded-order nonstrong deficit-two examples.

Connection contract: Hamiltonian words -> cycle cuts -> rotation subgroup orbits -> minimum outdegree -> a factorial cycle-count lower bound. The preserved data are black-orbit budget and actual successor symmetries. The quotient loses nonclosing paths and external arcs between rotation blocks; the restoration sidecar is internal regularity on each block plus the injective closing-word map. The square-root lemma controls all-blue cycles without assuming their successors are automorphisms. The C5/S4 controls distinguish purity, rigidity, and regularity.

The packet closes finiteness, not classification. Its honest remaining center problem is finite: determine the unknown strong self-converse deficit-two atoms and hence D(z). No generic classification or enlarged census is undertaken here.

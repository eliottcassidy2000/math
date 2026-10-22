# The exact continuation complexity of one-apex Hamiltonian counts

Session26, 22 October2026. Status: PROVED by the elementary arguments below; exact finite controls provided; independent audit pending. No literature-priority claim. This concerns actual Hamiltonian-path counts of marked one-apex tournaments, not the full tournament spectrum or merged metagraph.

## Inheritance and carrier

For a binary word w=b1...bn, let T(w) have the transitive vertices1<...<n and an apex x with x→i iff bi=1. H(w) counts all directed Hamiltonian paths. The empty word gives the one-vertex tournament, H=1. Session25's `session25-spectrum/PROOF.md` proves, retaining orphan paths, the exact recurrence

- append1: (H,A,p)→(H,2A+1,1);
- append0: (H,A,p)→(H+A+p,2A,0).

Here A is the binary integer represented by w and p its last bit (0 for the empty word). We recheck this against vertex-level subset DP, independently of the word recurrence. Canonical hostile: 00 and01 have H=1 but their zero-extensions have H=1 and3. Inherited modular completeness alone does not establish minimal continuation complexity.

Concept board: actual H; continuation equivalence; reachable states; observer lattice; missing parity; actual metagraph ports. The sampled THM-4010 (confluent Hasse observer/Smith firewall) suggests checking integral rather than merely rational rank. Sampled THM-3266 (two-face row-pair origin bit) suggests testing whether a last-bit coordinate survives odd-modulus reduction. These are methodological transfers, not maps from the original Hasse/LRC objects into tournaments.

## Theorem1: exact finite-state minimum

Fix M≥1. Two binary prefixes u,v are equivalent if H(uz)=H(vz) modulo M for every binary suffix z. The number of equivalence classes, hence the minimum number of states in a deterministic output automaton computing H(w) modulo M, is

- 1 if M=1 or2;
- 2M² if M≥3 is odd;
- M²/2 if M≥4 is even.

This is a state-count theorem, not a hardware-independent lower bound on memory bits or an assertion about every alternative representation.

### Reachability

The recurrence supplies an automaton on (h,a,p)=(H mod M,A mod M,p). If M is odd, its potential state set has 2M² elements. If M is even, H is odd and p is determined by a modulo2, leaving M²/2 elements.

All these states are reachable. Write M=2^k m with m odd. Choose L≥max(1,k) with Q=2^L=1 modulo m. For the block B=1^L0^L, the Session25 exact block recurrence is

 A'=Q²A+Q(Q−1),
 H'=H+Q(Q−1)A+D,  D=(Q−1)²+1.

Starting from (1,0,0), at every positive block boundary A=0 modulo M and H=1+tD modulo M; gcd(D,M)=gcd(2,M). Thus these prefixes realize all admissible H residues with A=0,p=0.

Given a desired (a,p), choose a nonnegative integer b with b=a modulo M and b=p modulo2. For odd M this is CRT; for even M it is possible precisely for the stated compatible p. Append a binary expansion of b (use the single digit0 when b=0). Its terminal A and p are a,p, independent of the block count. Its increment to H modulo M is likewise independent of that count. If M is even that increment is even: each zero-step adds A+p, an even integer. Tune the preceding block count to obtain any admissible desired h. This proves full reachability, not just output-surjectivity. Every witness can be chosen to begin with1; prefixes ending0 are strong tournaments.

### Observability

The three suffixes empty,0,10 produce respectively

 h,  h+a+p,  h+2a+2  (mod M).

If two states have equal observations then Δh=0, 2Δa=0 and Δa+Δp=0 modulo M. For odd M≥3, 2 is invertible, so Δa=0 and Δp=0; since p∈{0,1}, the actual bits agree. For even M≥4, either Δa=0 or Δa=M/2. The second possibility cannot satisfy Δa+Δp=0: Δp∈{−1,0,1} and M/2≥2, with all sums nonzero modulo M. Thus again the states coincide. For M=1,2, all H outputs are constant and one state suffices. Necessity and sufficiency are both proved.

The exceptional case M=2 is essential: the two reachable raw states collapse to one observable class. For odd M the last bit cannot be reconstructed from A mod M; forgetting it genuinely loses future output. For even M it is redundant.

## Theorem2: exact weighted dimension and integral loss

Let K(u,v)=H(uv) be the infinite integer prefix-suffix Hankel matrix. Its rank is exactly4 over Q and over every field of odd characteristic, and exactly1 over every field of characteristic2.

The recurrence represents the state linearly in four homogeneous coordinates (H,A,p,1), proving rank≤4. On row words empty,1,10,11 and column words empty,0,10,00 the matrix is

    1  1  3  1
    1  3  5  5
    3  5  9  9
    1  5  9 11

with determinant−16. It is nonsingular in exactly the required characteristics. In characteristic2 all entries are1, proving rank1. Equivalently any field-linear weighted automaton has dimension at least4 in odd characteristic, and the displayed realization attains it. This follows directly from factoring each K(u,v) through the state space, without importing an automata theorem.

The infinite matrix's nonzero determinantal divisors are exactly1,2,4,16 (Smith factors1,2,2,4, interpreted via gcds of finite minors). For minors of size2 or3, subtract their first row from every other row. All entries are odd, so these minors are divisible respectively by2 or4. The displayed matrix has gcds2 and4 among its corresponding minors, giving equality. For 4×4 minors use the factorization K=R C into state rows and suffix-response columns.

Every row (H,A,p,1) obeys H≡1 and A≡p modulo2. Subtract the last column from the first and the third from the second: both resulting columns are even, so every 4-row determinant is divisible by4.

Each suffix response has the form H+αA+βp+γ, where γ=H(v)−1 is even, and α≡β modulo2. To see the latter, run the recurrence over F2: initial A and p affect H only at the first symbol if it is0, contributing A+p; after that step A=p, and later zero-steps add0. If the first symbol is1, both initial coordinates are erased from future increments. Thus the suffix column (1,α,β,γ) obeys two independent parity constraints. Subtract its third row from its second, and note that its fourth row is even: every four-column determinant is divisible by4. Every 4×4 minor of K therefore is divisible by16, while the displayed minor attains16. Larger minors vanish.

The Smith fingerprint explains why characteristic2 looks completely blind while higher2-adic precision recovers a quadratic number of distinct continuation classes. Rational rank, modular rank, and deterministic state complexity are different invariants and are not equated.

## Two hypothesis/test/refinement cycles

1. Hypothesis: H mod M might be a complete continuation state after the modular-completeness theorem. Immediate hostile00/01 disproves it. Repair: retain A,p; prove all candidate states reachable; independently minimize them. A first unqualified formula at M=2 fails, repaired by the one-state exception. Result: Theorem1.
2. Hypothesis: the four homogeneous coordinates remain uniformly linearly visible at every prime. The determinant−16 and parity give the hostile prime2. Repair: characterize exact rational/odd-field rank and integral determinantal divisors rather than asserting uniform invertibility. Result: Theorem2. Higher2-adic distinctions are not recovered by a field-rank argument; Theorem1 separately supplies them.

## Scope and actual wiring loss

The map is word→continuation state. It preserves every future one-sided append count modulo M and destroys individual paths, their entry/exit vertices, permutation holonomy and actual complement adjacency. No state-count result classifies the merged metagraph or proves prime-valued H attainment. The full state is operation-relative: replacing append by arbitrary insertion or substitution requires a new port interface. This packet does not claim that a scalar, these four coordinates, or the Smith factors determine arbitrary wiring.

## Theorem3: synchronized CRT and strong completion

For coprime m,n≥3, the minimal continuation state at modulus mn is exactly the fiber product of the states at m and n over their shared literal last bit p∈{0,1}. The coordinate map is reduction of h and a to the two moduli, keeping p. Theorem1's reachability descriptions and ordinary CRT prove both injectivity and surjectivity onto that fiber product. In particular

    N(mn)=N(m)N(n)/2.

Each last-bit fiber has half the states: for odd M there are M² states for each p; for even M≥4 there are (M/2)² states for each p. Independent local choices of last bit are impossible. The product overcounts by exactly two, not an asymptotic factor. Modulus2 is exceptional: its observable state loses p entirely, and for odd m≥3 one has N(2m)=N(m). Indeed an odd integer's residue modulo m uniquely determines its residue modulo2m.

The fiber product is over the two literal control values, not over all idempotents of a composite residue ring. Mixed CRT idempotents are not additional physical last-bit choices. This is a synchronization statement about the same input word, not independence of local automata.

The state-count lower bound remains valid when queries are restricted to strong completions of prefixes beginning with1. Such prefixes reach every raw state by the reachability proof. The suffixes0,10,110 all end0 and hence give strong tournaments. Their readouts are

    h+a+p,  h+2a+2,  h+4a+4.

Equality of these readouts gives 2Δa=0, then Δh=0, then Δa+Δp=0; the observability argument is unchanged. This is a lower bound for any model answering all these strong-completion queries, not a claim that every intermediate prefix is strong.

## Third cycle: arithmetic independence fails at the shared control bit

Hypothesis: CRT factors make the minimal state counts multiply. Test at3 and5: N(3)=18,N(5)=50 but N(15)=450 rather than900. First failed implication: choosing valid local states separately need not make their last bits equal. Repair: the exact common-bit fiber product in Theorem3. This transfers the origin-bit lesson into an actual map with explicitly counted loss, rather than identifying unrelated S2 symmetries.

## Reproduction

Run `python3 verify.py` and `python3 -O verify.py`. The verifier uses no assertions that disappear under optimization. It checks all511 words of lengths0..8 against vertex-level subset DP, all reachable states and an independent Moore partition refinement for M1..80, and the exact displayed determinant. The saved output is `output.json`. These are finite controls; proofs establish all M and all words. Sample selection is recorded separately in `random-manifest.json`, seed2601, from the sorted HEAD tree at cf6f4b3b; full retrieved sources in `random-sources.txt` are historical research inputs, not current truth claims.

# Exact fixed-core attachment-word Hamiltonian counter

Status: proved by the explicit bijection below; executable finite controls PASS. Session26's final bounded instrument. No novelty, minimality, isomorphism-canonization, or unrestricted tournament-enumeration claim. No repository files were modified by this worker.

## Input and result

A presentation consists of a labelled k-vertex core tournament Q and an ordered transitive base b_1,...,b_m, with b_i -> b_j exactly when i<j. The input word a_1,...,a_m lies in {0,1}^k, with a_i[u]=1 meaning core u -> b_i. The answer H(Q,a) is the number of directed Hamiltonian vertex words of this particular tournament, not an expectation or count of classes. The deletion set need not be minimum or unique. Its presentation is supplied, not searched for.

For k>=1 this series has a nonnegative-integer weighted-automaton realization with at most 2 k! 3^k states, in k! disjoint blocks of 2·3^k states. The initial vector has one unit in each block's empty state. There is no minimality claim: at k=1 this construction has six states, whereas the inherited rational linear continuation rank is four (session26-tournament/PROOF.md). For k=0 the whole tournament is transitive and H=1, including the conventional single empty Hamiltonian word when m=0.

## Slot bijection and state invariant

Fix a permutation sigma=(sigma_0,...,sigma_{k-1}) of the core. Every Hamiltonian word inducing this core order has k+1 base slots: slot 0 precedes sigma_0; slot j, 1<=j<k, lies between sigma_{j-1} and sigma_j; slot k follows sigma_{k-1}. Within each slot the base labels must increase. Thus assigning each base label to one slot uniquely determines a candidate vertex word, and different assignments or different core orders give different words.

Read the base labels in increasing order. For each slot j<k retain one of three states:

- 0: empty;
- 1: nonempty, and its last base vertex points to sigma_j;
- 2: nonempty, and sigma_j points to its last base vertex.

For the terminal slot k retain only empty/nonempty. This gives 3^k·2 states. A digit 2 is NOT immediately rejected: a later base vertex placed in the same slot can repair its final exit. The entry edge of a noninitial slot, once verified, never changes; every new vertex appended to that slot is larger than all its current vertices, so its internal predecessor edge is automatic.

On reading a new label a, choose exactly one slot:

1. For slot 0 no entry condition applies. Set its digit to 1 if a[sigma_0]=0, and to 2 otherwise.
2. For slot j with 0<j<k, if it is empty require a[sigma_{j-1}]=1. If already occupied impose no entry condition. Set its digit to 1 if a[sigma_j]=0, and to 2 otherwise.
3. For slot k, if empty require a[sigma_{k-1}]=1; otherwise impose no condition. Set its occupancy to nonempty.

Every permitted slot choice contributes weight one. Distinct choices can have the same destination state; the matrix entry is their multiplicity, not merely a Boolean adjacency.

At the end accept precisely when no digit is 2 and, for every empty interior slot j=1,...,k-1, Q[sigma_{j-1}][sigma_j]=1. Empty first and terminal slots are allowed. Occupied interior slots have already passed their entry tests and now pass their exit tests. Thus acceptance checks exactly all actual consecutive edges of the candidate Hamiltonian word, proving the bijection and the count.

All internal core orientation dependence appears in this final acceptance condition for empty interior slots. This does NOT imply that H is independent of Q: with k=2 and the one-letter mask word [1], the core arc 0->1 gives H=1, while 1->0 gives H=3. Those cores have the same ordered path-cover profile (1,2), so a profile alone cannot replace the labelled core/attachment ports for a fixed word.

With row vectors, define M_sigma(a)[s,t] as the number of allowed choices taking s to t. Then

H(Q,a_1...a_m) = sum_sigma e_empty^T M_sigma(a_1)...M_sigma(a_m) f_Q,sigma.

The empty product gives exactly the valid directed Hamiltonian core orders, so H(Q,empty)=H(Q), including nontransitive cores. No artificial empty interior arc is added.

## Complexity and implementation

At most k+1 transitions leave a state for one input label. Straight propagation costs O(m k! (k+1) 3^k) integer additions plus O(k! k 3^k) final checks. For m=0, retain the separate initialization/final-check cost rather than claiming zero work. This is linear in m in the unit-cost arithmetic model for fixed k, not constant-bit-time arithmetic. Integers have O(log(k!)+m log(k+1)+1) bits because there are at most k!(k+1)^m assignments. Schoolbook addition gives this extra bit-length factor for the counting arithmetic; state-index and input-decoding costs are parameter-bounded. The output itself can require linearly many bits in m.

counter.py processes core permutations sequentially, reusing a transition table whose local label coordinates are relabelled according to sigma. It validates the square 0/1 tournament matrix, zero diagonal, complementary opposite entries, and integer masks 0<=a<2^k. It materializes the input word. The cached transition table adds O(2^k (k+1) 3^k) preprocessing time/storage for each k; dynamic count storage is O(3^k) integer entries per permutation, plus O(m) input labels. The cache persists per k within a process and can be cleared with transitions.cache_clear(). The code is intended for small fixed k; it has no claim of practicality when k grows unrestrictedly.

The implementation uses state index sum_j digit_j 3^j plus 3^k times terminal occupancy. It explicitly retains repeated destinations. Its API is count_hamilton(Q,word). Its CLI accepts a JSON file with Q and word. Matrix row u, column v means u->v. Integer word bit u is the core-u out-neighborhood indicator, not its opposite.

## Independent controls and the mean bridge

The independent audit constructs each full tournament and uses vertex-subset endpoint DP, with recurrence dp[S union {v},v] += dp[S,u] for u->v. This DP does not use slot assignments or automaton states.

Exhaustive comparisons cover every labelled core and every attachment word in these ranges: k=0,m<=8; k=1,m<=7; k=2,m<=5; k=3,m<=3; k=4,m<=2. There are 25,146 full-tournament DP comparisons. Separately all 8,191 binary words through length 12 agree with the inherited one-apex recurrence initialized (H,A,p)=(1,0,0): append 1 gives (H,2A+1,1), append 0 gives (H+A+p,2A,0).

The known five-vertex increasing-pair mask-41 H=13 witness has marked core (0,1), transitive base order (4,3,2), core arc 0->1 and mask word [1,2,0]. Both engines give 13. This supplies a ready-to-use k=2 input in witness13.json. It is a presentation, not a canonical representative or spectrum enumeration.

For every positive-k core/base-size pair above (244 pairs), the audit sums the actual counter over all 2^(km) words and divides exactly by 2^(km). It matches the independently obtained ordered path-cover formula from session26-attachments/proof.md:

E_Q(m) = sum_{r=1}^k p_r(Q)/4^r [R_{r+1}(m)+4R_r(m)+4R_{r-1}(m)],

where R_q(m)=q!S(m,q) and p_r counts ordered spanning covers by r nonempty directed paths. The audit computes p_r directly from core permutation cuts and R_q by exact inclusion-exclusion. These expectation checks supplement, not replace, the individual-word DP comparisons. Averaging sends the port-sensitive word series to a small polynomial-exponential channel space determined by the core path-cover profile; this does not make the individual count, or any individual leading growth law, core-independent. The example above displays information lost by averaging.

Additional controls reject eight malformed inputs, exercise direction-sensitive 10 versus 01, preserve the k=0 and m=0 boundaries, and evaluate a 100-letter k=2 word without a huge census. All checks use explicit failures, so optimization does not disable them.

## Reproduction and search boundary

From this directory run:

    python3 audit.py > audit.out
    python3 -O audit.py > audit-optimized.out
    cmp audit.out audit-optimized.out
    python3 counter.py witness13.json

The final command prints 13. Python standard library only; no random seeds, sampling, floating point, or external files are needed by the audit.

Closest local inherited source: THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md, ordered path covers and permutation cuts, together with the current exact random-attachment transform and one-apex recurrence. A targeted local search for fixed apex, transitive deletion, Hamiltonian slots and attachment automata was performed. This is not an external priority search; this note makes no claim that this parameterized dynamic program is new. No further research lane is started here. Parent owns independent review, theorem-ID reservation, repository integration and publication.

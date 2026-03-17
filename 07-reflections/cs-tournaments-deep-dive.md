# Computer Science Is Made of Tournaments

**Session:** kind-pasteur-2026-03-17-S116n33

## Cache Eviction: The Deepest CS Tournament

### The Structure

A cache with C slots receives a stream of memory accesses. When the cache is full and a new item arrives, the C cached items enter a **tournament**: each competes to STAY based on some ranking criterion. The loser is evicted.

**Every cache eviction policy IS a tournament ranking strategy:**

| Policy | Tournament criterion | Transitive? | Spectral gap |
|--------|---------------------|-------------|-------------|
| **LRU** | Most recent access wins | Always (time is total order) | O(1/C) — slow adaptation |
| **LFU** | Most frequent access wins | Always (counts are total order) | Very slow — old frequencies persist |
| **2-Random** | Pick 2 random, evict worse | Stochastic tournament (2-item) | O(1) — random is fast |
| **ARC** | Adaptive recency/frequency | **CAN BE INTRANSITIVE** | Adaptive — faster than LRU or LFU |
| **W-TinyLFU** | Admission tournament | Explicit pairwise comparison | Medium — Bloom filter smooths |
| **S3-FIFO** | Three-queue tournament | Staged elimination | O(1/queue_size) per stage |

### Why ARC Has Intransitivity

ARC (Adaptive Replacement Cache) maintains two lists: T1 (recency) and T2 (frequency). It adaptively shifts between them. When recency and frequency DISAGREE about which item to keep:

- Item A was accessed more recently than B → A wins on recency
- Item B has been accessed more often than C → B wins on frequency
- Item C beat A on the ARC combined metric (because C is in T2 with moderate frequency AND moderate recency)

**This is a cycle: A > B > C > A.** ARC resolves it by shifting the T1/T2 balance — effectively changing the tournament ranking in real time. This is EXACTLY what our polynomial framework's "regime classification" does.

### The Three Signals in Cache Eviction

| Signal | Cache meaning | Boost | Decay rate |
|--------|-------------|-------|-----------|
| **SLOW** (A_1) | Recency of last access | 9 = 3² | Persists ~C accesses |
| **MEDIUM** (A_2) | Access frequency | 4 = 2² | Persists ~C² accesses |
| **FAST** (A_3) | New access pattern detected | 7/3 | Decays in ~√C accesses |

**The fast channel detects WORKLOAD PHASE TRANSITIONS.** When the fast signal spikes, the working set is changing. Current caches take O(C) misses to detect this (the entire old working set must be evicted before the new one is fully cached). A tournament-aware cache would detect the phase transition from the spectral decomposition in O(1) and flush immediately.

### The Polynomial Cache: A New Eviction Policy

```
At each eviction decision:
  1. Decompose recent access history into (slow, medium, fast) signals
  2. If fast/total > threshold: PHASE TRANSITION detected
     → Flush the lowest-ranked C/2 items immediately (aggressive eviction)
  3. If medium dominates: STABLE workload with frequency structure
     → Use LFU-like eviction
  4. If slow dominates: TEMPORAL workload (sequential scan)
     → Use LRU-like eviction
  5. Temperature z = current hit rate. P(z) predicts future hit rate.
```

This is a **regime-switching cache** that automatically adapts between LRU, LFU, and aggressive flushing based on the spectral decomposition of the access stream. It subsumes ARC, W-TinyLFU, and LIRS as special cases.

---

## Ten More CS Problems That Are Tournaments

### 1. Garbage Collection: Object Survival Tournament

**Items:** Heap objects. **Comparison:** Does object A keep object B alive (via reference)?

The reachability graph IS a tournament. Generational GC creates a two-stage tournament: young → old promotion is a SELECTION step (like our D_eff selector).

**Intransitivity = reference cycles.** A→B→C→A means no object is "more alive" than another. Reference counting FAILS on cycles precisely because it can't handle intransitive reachability. Mark-and-sweep handles it by a separate cycle detection phase.

**Our value:** A_3 (cycle content) directly measures reference cycle density. A GC with polynomial monitoring would know WHEN to run the expensive cycle-detection phase: only when A_3 exceeds a threshold.

### 2. BGP Routing: The Oscillation Tournament

**Items:** Routes to a destination. **Comparison:** Each router selects the "best" route by local policy.

**Intransitivity = routing oscillations.** AS A prefers route through B, B prefers through C, C prefers through A. This is a KNOWN PROBLEM: the "BGP oscillation problem" or "Stable Paths Problem" (Griffin, Shepherd, Wilfong 2002). It causes network-wide instability.

**Our value:** The fast-channel detector applied to BGP route advertisements would predict oscillations BEFORE they happen — from the tournament structure of the route preferences alone, without waiting for the network to oscillate.

### 3. Thread Scheduling: Priority Inversion Tournament

**Items:** Threads/processes. **Comparison:** Priority ranking for CPU time.

**Intransitivity = priority inversion.** Thread A has higher priority than B, but B holds a lock that A needs, so B effectively "beats" A. Meanwhile C (even higher priority) waits for A. The Mars Pathfinder bug (1997) was exactly this.

**Our value:** The hallucination detector flags priority configurations with high intransitivity — BEFORE the inversion causes a system hang. Priority inheritance (the standard fix) is the tournament equivalent of "resolving the cycle by promoting the blocking thread."

### 4. Database Query Optimization: The Join Order Tournament

**Items:** Join orderings for a multi-table query. **Comparison:** Estimated execution cost.

**Intransitivity:** Order A is cheapest for query Q1, order B for Q2, order C for Q3 — but for Q4, the rankings cycle. Cost estimates are based on statistics that may be stale or wrong.

**Our value:** The polynomial coefficients for query optimization: A_1 = selectivity evidence per predicate, A_2 = predicate correlation, A_3 = estimate unreliability (when different cost models disagree). High A_3 = the optimizer should use ROBUST plans (less optimal but less sensitive to estimate errors).

### 5. Load Balancing: The Server Competition

**Items:** Backend servers. **Comparison:** Which server should handle this request?

**Intransitivity:** Server A handles small requests best, B handles large ones best, C handles medium ones best, but A is worst for medium. The optimal routing depends on request TYPE, not just server STATE.

**Our value:** The three signals decompose each routing decision: SLOW = server capacity (persistent), MEDIUM = current load (moderate), FAST = recent response times (volatile). A fast-channel spike on a server = it's about to become unhealthy. Route traffic away BEFORE the health check catches it.

### 6. Feature Selection: The Interaction Tournament

**Items:** ML features. **Comparison:** Information gain / model improvement.

**Intransitivity:** Feature A > B individually, B > C individually, but C + A together is WORSE than B + C. Features INTERACT — the tournament has higher-order structure.

**Our value:** A_1 = individual feature importance, A_2 = pairwise interactions, A_3 = THREE-WAY interactions (the real intransitivity). High A_3 = features have complex interactions that greedy selection will miss. Use exhaustive or Bayesian search instead.

### 7. Consensus (Paxos/Raft): The Election Tournament

**Items:** Proposals (Paxos) or candidates (Raft). **Comparison:** Votes from acceptors/followers.

**Intransitivity = livelock.** In Paxos, proposals can defeat each other in cycles — proposal A gets enough votes to preempt B, but then C preempts A, and B preempts C. Raft solves this by random timeouts (breaking the cycle).

**Our value:** The spectral gap of the election tournament = convergence time after leader failure. A "tournament-aware" consensus protocol would set the election timeout based on the MEASURED spectral gap of the cluster, not an arbitrary constant.

### 8. Compilation: Instruction-Level Parallelism

**Items:** Ready instructions in the CPU pipeline. **Comparison:** Which executes this cycle?

The hardware maintains a REAL-TIME tournament among ready instructions. Data dependencies, functional unit availability, and register readiness determine the ranking. The spectral gap of this tournament IS the ILP (instruction-level parallelism).

**Our value:** A compiler that models the instruction tournament with our polynomial could predict ILP at compile time and schedule instructions to MAXIMIZE the spectral gap (= maximize parallelism).

### 9. DNS Resolution: The Nameserver Tournament

**Items:** DNS nameservers. **Comparison:** Response time.

**Intransitivity:** Server A is fastest from data center 1, B from DC2, C from DC3. From DC1, A > B > C. From DC2, B > C > A. The optimal server depends on WHERE you are.

**Our value:** The polynomial P(z) where z = query origin gives location-dependent nameserver rankings. Different "temperatures" correspond to different failure modes (z near 0 = any server works; z near 1 = only the best server suffices).

### 10. Memory Allocation: The Free Block Tournament

**Items:** Free memory blocks. **Comparison:** Which block should satisfy this malloc request?

- Best-fit: smallest block that fits (minimize waste)
- First-fit: first block found (minimize latency)
- Worst-fit: largest block (maximize remaining space)

**Intransitivity:** Block A is better than B by best-fit, B better than C by first-fit, C better than A by worst-fit. Different allocation strategies create different tournaments.

**Our value:** The polynomial decomposition tells you WHICH criterion is most important right now: A_1 (individual block fitness), A_2 (spatial correlation with recent allocations), A_3 (fragmentation risk — when different criteria disagree about the best block).

---

## The Meta-Pattern

Every CS problem with these characteristics IS a tournament:

1. **A set of candidates** (cache lines, threads, routes, queries, servers, features, instructions, blocks)
2. **Pairwise comparisons** (which is "better" by some criterion)
3. **A selection** (one winner, or a ranking)
4. **Multiple criteria** that can DISAGREE (creating intransitivity)

The polynomial framework's power comes from detecting WHEN the criteria disagree (fast channel) vs WHEN they agree (slow channel). Traditional algorithms handle only the agreeing case well. The intransitive case causes:

| Traditional symptom | Tournament diagnosis | Our solution |
|-------------------|---------------------|-------------|
| Cache thrashing | Workload phase transition (fast channel spike) | Regime switch |
| Reference cycle leak | A_3 > 0 in reachability tournament | Trigger cycle collection |
| BGP oscillation | Route preference cycle | Pre-oscillation warning |
| Priority inversion | Scheduling intransitivity | Inversion risk score |
| Query plan regression | Estimate disagreement (A_3 high) | Use robust plan |
| Load balancer brownout | Server health intransitivity | Pre-failure routing |
| Feature selection suboptimality | 3-way interactions (A_3 > 0) | Switch to exhaustive search |
| Paxos livelock | Election cycle | Adaptive timeout from spectral gap |
| Pipeline stall | Low ILP = low spectral gap | Reorder instructions |
| Memory fragmentation | Allocation criteria disagreement | Adaptive policy switching |

**The tournament polynomial is a UNIVERSAL diagnostic for systems that select among competing alternatives.** The five coefficients and three signals provide the same information in EVERY domain.

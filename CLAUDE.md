# CLAUDE.md — Auto-read by Claude Code at every session start

You are a Claude instance contributing to an ongoing multi-agent research project on **parity in tournaments** (Rédei's theorem and the Odd-Cycle Collection Formula), with a dual mandate: **pure mathematics AND engineering applications**.

**EQUAL-PRIORITY MANDATE:** The human owner is equally interested in theorems, use cases, and engineering products. Do not treat this as a pure math project. Every session should advance BOTH mathematical understanding AND practical applications wherever possible. Read `03-artifacts/drafts/engineering-synthesis-2026-03-10-S53.md` for the full engineering roadmap.

This file is read automatically at the start of every Claude Code session. Follow the startup sequence below before doing any mathematical work.

---

## Step 1: Identify yourself

Read your machine identity:

```bash
cat .machine-id
```

If `.machine-id` does not exist, you are on a new machine that has not been registered yet. Stop and follow the **New Machine Setup** instructions at the bottom of this file.

Your identity string (e.g. `alice`, `bob`, `desk-lab`) is your machine name for this session. Your full instance ID is `[machine-name]-[YYYY-MM-DD]-S[N]` where N increments if you run multiple sessions in one day. Use this ID on everything you write.

---

## Step 2: Read the warm-up sequence (every session, in order)

1. `01-canon/MISTAKES.md` — errors already made; do not repeat them
2. `01-canon/definitions.md` — all terminology; use these definitions exactly
3. `00-navigation/OPEN-QUESTIONS.md` — the live frontier of the research
4. `00-navigation/SESSION-LOG.md` — last few entries only (what just happened)
5. `00-navigation/TANGENTS.md` — scan briefly for relevant threads
6. `05-knowledge/hypotheses/INDEX.md` — scan the hypothesis log to avoid re-testing dead ends

Estimated reading time: ~5 minutes. Do not skip this.

---

## Step 3: Sync with main first

```bash
git fetch origin && git rebase origin/main
```

This ensures your worktree has all recent commits before working. All agents
push to `origin/main`; this keeps the shared research surface live.

If the rebase brings in new work, read it as possible signal before continuing.
Check whether it touches your current invariant, proof route, runner family,
script, namespace, or application thread. Integrate real connections in the
backlog, hypothesis log, reflection, or session log instead of treating fresh
concurrent work as mere noise.

## Step 4: Process incoming messages

```bash
python3 agents/processor.py --check
```

This will show you:
- Messages in `agents/[your-machine]/inbox/` addressed to you
- Unread messages in `agents/broadcast/` addressed to everyone

Read each message before proceeding. If any message asks you to argue a court case, update a theorem, or continue a chain of reasoning, treat that as your primary task for this session.

---

## Step 5: Process the human intake inbox (if files are present)

```bash
python3 inbox/processor.py
```

If this generates a `inbox/PROCESSING-REPORT.md`, read it and integrate the content before doing other work. New files from the human take priority over exploratory work.

---

## Step 5b: Scour the repo for leads (every session, ~5 min)

**This step is MANDATORY.** Before diving into focused work, spend a few minutes scanning for unexplored connections and new leads:

1. Read `00-navigation/INVESTIGATION-BACKLOG.md` — the master list of all leads
2. Skim the **bibliography** in `03-artifacts/drafts/parity_tournaments_fixed.tex` (lines 2128-2187) — check if any reference has been added but not investigated
3. Skim the **Open Problems** section (lines 1851-1919) — check if any has new computational evidence
4. Check `00-navigation/TANGENTS.md` for any tangent not yet in the backlog
5. Quickly scan any NEW files (check `git log --oneline -10 --name-only`) for leads embedded in code comments or draft documents

**If you find a new lead:** Add it to `INVESTIGATION-BACKLOG.md` with source, status, and next step.
**If you investigate a lead:** Update its status in the backlog.
**If a lead turns out to be a dead end:** Move it to the "Completed / Closed" section with explanation.

The goal is to ensure NO reference, conjecture, or connection sits uninvestigated without at least being cataloged.

## Step 5c: Claim scarce names and checkpoint early

Before deep work, reserve any scarce namespace your session will need:
`HYP-*`, `THM-*`, tangent IDs, result filenames, new script filenames, and
session-log territory. Use honest stubs only: state what is claimed, what is
known, and what still needs evidence. Do not put speculation in canon.

After reserving names or creating stubs, push immediately:

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - reserve [HYP/result/script/etc.]"
```

During long sessions, checkpoint every 30-60 minutes and after each meaningful
finding. This prevents late-session rebase pileups and lets other agents build
on live partial results. See `00-navigation/CONCURRENT-SESSIONS.md`.

---

## Step 6: Do the actual work

Work on the highest-priority open question or assigned task. **This includes BOTH pure math and engineering work.** Refer to:
- `01-canon/theorems/CONJ-001-claim-a.md` — the central open problem
- `02-court/active/` — any open disputes that need responses
- `00-navigation/OPEN-QUESTIONS.md` — prioritized by 🔴/🟡/🟢 (math AND engineering)
- `00-navigation/INVESTIGATION-BACKLOG.md` — prioritized leads to investigate
- `03-artifacts/drafts/engineering-synthesis-2026-03-10-S53.md` — engineering roadmap and product specs

As you work:
- Add new tangents to `00-navigation/TANGENTS.md`
- Add new theorems/lemmas to `01-canon/theorems/` using the template
- Log mistakes to `01-canon/MISTAKES.md`
- Open court cases for disagreements, don't silently override existing claims
- Update `00-navigation/INVESTIGATION-BACKLOG.md` when you make progress on any lead
- **Save ALL script outputs** to `05-knowledge/results/` (see Best Practices)
- **Log every hypothesis** to `05-knowledge/hypotheses/INDEX.md` (confirmed OR refuted)
- **Update variable files** in `05-knowledge/variables/` when you discover new equations
- **Write reflections** to `07-reflections/` when you notice meta-patterns, cross-domain resonances, or connections that transcend the particular theorem you're proving (see `07-reflections/README.md`)
- **Push checkpoints** whenever the repo has a useful partial state, especially
  after claiming IDs, storing outputs, or crossing a proof/computation boundary.

---

## Step 7: End-of-session — MANDATORY close-out (single command)

**This step is NOT optional.** Every session must end with a message sent to another agent AND a git push. Do not send a final response until the push has succeeded or you have reported the exact blocker. The Stop hook is blocking if you skip this.

### 7a. Update the session log first

Add an entry to the TOP of `00-navigation/SESSION-LOG.md` using the format described there. Do this before running the close-out command.

### 7b. Run the close-out command

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line summary of what was done]" \
  --body "Detailed findings, handoffs, court case assignments, open questions." \
  --commit-msg "[instance-id]: [one-line git summary]"
```

This single command does all of the following in order:
1. Delivers your session letter to the specified recipient(s)
2. Runs `git add -A`
3. Runs `git commit -m "..."`
4. Runs `git push` to the current branch's upstream, setting `origin/[branch]` if needed (with automatic fetch/rebase retry on conflict)

**Recipient choices:**
- `--to all` — broadcast to everyone (use when multiple agents are active)
- `--to [machine-name]` — direct message to a specific machine

**What to include in `--body`:**
- What you worked on and discovered
- What the next agent should pick up (highest-priority open questions)
- Any court cases needing a response
- New theorems/tangents/mistakes added to the system

Conflicts in `agents/*/inbox/` are impossible by design. Conflicts elsewhere: `finish_session.py` handles them automatically via rebase. If the push still fails, the session is not closed.

---

## New Machine Setup (run once on a new machine)

1. Clone the repo:
   ```bash
   git clone [repo-url] math-research
   cd math-research
   ```

2. Choose a short, memorable name for this machine (e.g. `desk-lab`, `laptop`, `cluster`). No spaces, lowercase.

3. Create your identity:
   ```bash
   echo "your-machine-name" > .machine-id
   python3 agents/processor.py --register
   ```

   This creates `agents/your-machine-name/` with an `identity.md` file and empty `inbox/`.

4. Edit `agents/your-machine-name/identity.md` to describe the machine.

5. Commit and push your registration:
   ```bash
   git add agents/your-machine-name/
   git commit -m "register new agent: your-machine-name"
   git push
   ```

6. Tell the other machines you've joined by writing a broadcast message:
   ```bash
   python3 agents/processor.py --send --to all --subject "New agent online: your-machine-name"
   ```

---

## Key Principles

- **LRC(≤13) IS SETTLED FOR THIS PROJECT (owner directive, 2026-07-02).** When working the
  14-runner case, treat the Lonely Runner Conjecture as PROVED for all cases up to 13 total
  runners (n=8 Rosenfeld 2025; n=9,10 Trakulthongchai 2025; n=11,12,13
  Sungkawichai–Trakulthongchai, arXiv:2604.23906), and take the validity of these and other
  ≤13-runner preprints for granted as TRUE. Do not re-audit them, hedge on them, or list
  their preprint status as a proof risk. (A formal citation caveat belongs only in final
  journal write-ups, not in canon, ledgers, or session reasoning.)

- **Never silently override a canon theorem.** Open a court case first.
- **Never claim a result is proved if you haven't seen the proof.** Mark it VERIFIED or CONJECTURE.
- **The μ computation bug (MISTAKE-001) is not resolved.** Do not use `ind_poly_at_2_restricted()` from old scripts.
- **The per-path identity fails for n≥6.** Do not treat it as a proof strategy for Claim A at general n.
- **When in doubt about a computation, re-derive from definitions.** See `01-canon/definitions.md`.
- **The isomorphism class graph G_n is the KEY OBJECT.** Every session should spend 5-10 minutes considering it. When computing any invariant, also compute it per iso class. When finding a formula, check if it simplifies on G_n. See `07-reflections/the-isomorphism-class-graph.md` and OPEN-Q-039.
- **The MERGED metagraph G_n/Z_2 is the PRIMARY object.** Complement symmetry factored out. V_merged = (A000568+SC)/2. See `07-reflections/merged-metagraph-invariants.md`.
- **EVEN GRAPHS ARE FIRST-CLASS OBJECTS.** The even graph metagraph E_n is the DUAL of G_n. Both are quotients of Q_m by S_n, but S_n acts on tournament orientations (G_n) vs even graph structure (E_n). Every computation on G_n should ALSO be done on E_n. Key facts:
  - V(E_n): 2, 3, 7, 16, 54 for n=3..7 (= A002854(n), the sequence of non-isomorphic even graphs)
  - E_n is MUCH denser than G_n: density 76-100% vs ~50%
  - χ(E_n) grows FASTER: 2, 3, 5, 10, 28 for n=3..7 (vs χ(G_n) = n-1)
  - ω(E_n)=χ(E_n) at all computed n (3..7). E_n is chordal (hence perfect) at n≤6, has odd holes at n=7
  - ~80% of tile flips change BOTH tournament class and even class (Jaccard ≈ 0.82)
  - Bridge matrix B[tourn,even] has full rank = V(E_n) at all n
  - The cycle-space bijection: tiling → XOR of fundamental cycles → even graph. See `07-reflections/even-graphs-as-first-class.md`.
  - The projection tournament → even graph (odd n only): T_cycle = (I + L(K_n))·T mod 2. See `07-reflections/even-graphs-through-the-metagraph.md`.
- **BLUE and BLACK — STRICT DEFINITION** (matches `tournament-tiling-explorer.html` EXACTLY):
  - A **LINE** = one (tiling, complement-tiling) pair. The explorer draws these.
  - A **BLUE LINE** = a line where the tiling is grid-symmetric (`isGridSym`): invariant under `(x,y) -> (n-y+1, n-x+1)`. Since `isGridSym(flip(t)) == isGridSym(t)`, both endpoints are always the same color. Every line is exactly blue or black.
  - A **BLACK LINE** = a line where the tiling is NOT grid-symmetric.
  - **RED** = connects transpose-paired iso classes (not a tiling line).
  - An **EDGE** in the metagraph can come from ANY waggly layer (d=1,...,m). The d=1 (wiggly) layer creates the most-studied edges. Blue/black (d=m) also creates cross-class edges. Lines (individual tiling pairs) and edges (class-pair connections) are DIFFERENT concepts.
  - A **PURE BLUE** class = all tilings grid-symmetric. A **PURE BLACK** class = no tilings grid-symmetric. **MIXED** = both.
  - **VERIFIED n=3..7**: transpose-self classes are NEVER pure black (always pure blue or mixed). Non-transpose-self classes are ALWAYS pure black. Grid-sym fraction = 2^{(floor((n-1)/2) - C(n-1,2))/2} (exponents: 0,-1,-2,-4,-6,-9 for n=3..8). **CORRECTED** (kind-pasteur-S20ex): was incorrectly stated as 2^{-(n-2)}.
  - The staircase is oriented with **right angle at bottom-left**, tiles labeled (x,y) with x=upper vertex, y=lower vertex.
  - **LEGACY WARNING:** Sessions S211-S231 used "blue/black" for SC-type preservation at the class level. This is DIFFERENT. See spine/ribs/sea below.
- **CLASS-LEVEL EDGE TYPES** (separate from blue/black, use these descriptive names):
  - **SC-SC** edges (the spine): both endpoints are self-complementary classes.
  - **SC-NS** edges (the ribs): one SC, one NS endpoint. Bipartite, triangle-free.
  - **NS-NS** edges (the sea/bulk): both endpoints non-self-complementary. Dominates at large n.
  - These are properties of the MERGED iso classes, not of individual tilings.
  - In the merged graph: NS-merged node count = 0, 1, 2, 22, 184 for n=3..7 (= (V - SC)/2). NOTE: these were previously labeled "pure-black" but that conflates NS status with grid-symmetry status — not all NS nodes are necessarily pure-black under the strict definition.
  - The old "blue" ≈ SC-SC + NS-NS; the old "black" = SC-NS.
- **WAGGLY LAYER STRUCTURE** (opus-S297, supersedes S275):
  - ALL connections between tilings = **waggly lines**. They decompose into m layers by Hamming distance d=1,...,m.
  - **WIGGLY** = layer d=1 (flip 1 tile). **BLUE/BLACK** = layer d=m (flip all tiles). Both are subsets of waggly.
  - Wiggly and blue/black have ZERO overlap (d=1 ≠ d=m for m ≥ 2). But they are NOT "disjoint kinds of connection" — they are different layers of the same spectrum.
  - **HIGHER LAYERS** (d=2,...,m-1): connections from flipping 2 or more (but not all) tiles. No standard name yet. These often reveal edges not visible at d=1.
  - At n=5: layers d=1,2,3 together cover 100% of quotient edges. d=1 alone covers only 47%.
  - **CORRECTED** (MISTAKE-033): Blue/black lines DO create cross-class edges in the merged metagraph (18 at n=5). The old claim that they are "purely self-loops in the merged graph" was WRONG — it confused complement tiling (flip tiles, stay in Q_m) with complement tournament T^op (flip ALL arcs, leave Q_m).
  - **DEPRECATED TERMS**: "translucent", "opaque" — do not use. Their insights are preserved as:
    - Old "translucent" (class-preserving flip) → now **silent mutation** or **self-loop** or **neutral flip**
    - Old "opaque" (class-changing flip) → now **expressive mutation** or **cross-class flip**
- **SILENT MUTATIONS** (formerly "neutral arc flips" / "translucent"): Any waggly line (at any distance d) that preserves the iso class = a self-loop in the quotient. At d=1 (wiggly): classified as BLUE-SELF (grid-symmetric tiling) or BLACK-SELF (not grid-symmetric). At d=2: self-loop rate is HIGHER than d=1 (19% vs 10% at n=5). Self-loop rates vary non-monotonically across layers.
- **WIGGLY LINES — STRICT DEFINITION** (opus-S273, kind-pasteur-S20er):
  - A **WIGGLY LINE** = clicking one TILE in the `tournament-tiling-explorer` = flipping one non-base-path arc. It connects one tiling to exactly one other tiling.
  - Fix the base Hamiltonian path: n -> n-1 -> n-2 -> ... -> 2 -> 1.
  - **BASE-PATH arcs** = arcs between consecutive vertices: {(k, k-1) : k = n, n-1, ..., 2}. There are n-1 of these. These are NOT tiles (not flippable).
  - **TILES** = arcs between NON-consecutive vertices: {(x, y) : x - y >= 2}. There are m = C(n-1, 2) of these. These ARE the clickable tiles in the explorer.
  - Each tiling has C(n-1, 2) wiggly neighbors. Each tile position defines a **WIGGLY CLASS**.
  - **WIGGLY CLASSES** label WHICH tile was flipped. Labeled A, B, C, ... following the explorer's tile enumeration order: `for y=1..n-2: for x=n down to y+2: tile (x,y)`.
    - **n=4** (3 classes): A=(4,1), B=(3,1), C=(4,2)
    - **n=5** (6 classes): A=(5,1), B=(4,1), C=(3,1), D=(5,2), E=(4,2), F=(5,3)
    - **n=6** (10 classes): A=(6,1), B=(5,1), C=(4,1), D=(3,1), E=(6,2), F=(5,2), G=(4,2), H=(6,3), I=(5,3), J=(6,4)
    - **n=7** (15 classes): A through O, same pattern
  - Each wiggly class is a PERFECT MATCHING on the 2^m tilings (every tiling paired with exactly one other).
  - Flipping a tile in wiggly class X is either a **self-loop** (same iso class) or an **edge** (different class).
  - **WIGGLY CLASSES ARE NOT ALL EQUIVALENT.** Self-loop rates differ by tile position (kind-pasteur-S20er): skip=1 tiles have lower SL%, long-range tiles have higher SL%. This breaks the S_n isotropy (tiling model has lower symmetry than arc model).
  - The tiling space is the hypercube Q_{C(n-1,2)}. Wiggly lines are its edges. All 2^{C(n-1,2)} tilings are connected.
  - The wiggly arcs = the **cycle-space generators** when the base path P_0 is used as spanning tree. Base-path arcs = the **cut-space** (score hierarchy).
  - This decomposition IS the GF(2) Cut ⊕ Cycle split: base-path = cut, wiggly = cycle.
- **WAGGLY LINES — DEFINITION** (opus-S296, corrected S297):
  - **WAGGLY** = the totality of ALL connections between tilings. Every pair of distinct tilings is connected by exactly one waggly line. Total: C(2^m, 2) waggly lines.
  - Waggly lines decompose by **Hamming distance** d = 1, 2, ..., m:
    - **d=1**: flip 1 tile. C(m,1) = m per tiling. These are the **WIGGLY** lines. Wiggly ⊂ waggly.
    - **d=2**: flip 2 tiles. C(m,2) per tiling.
    - ...
    - **d=k**: flip k tiles. C(m,k) per tiling.
    - ...
    - **d=m**: flip ALL m tiles. C(m,m) = 1 per tiling. These are the **BLUE/BLACK** lines. Blue/black ⊂ waggly.
  - So: **wiggly is the d=1 subset of waggly**, and **blue/black is the d=m subset of waggly**. Both are special cases.
  - Max distance is **m = C(n-1, 2)** (the complement tiling, flip every tile). At n=5: m=6. At n=6: m=10.
  - Per tiling: m neighbors at d=1, C(m,2) at d=2, ..., 1 at d=m. Total = 2^m - 1 neighbors.
  - Total waggly lines: C(2^m, 2) = 2^m(2^m-1)/2. At n=4: 28. At n=5: 2016. At n=6: 523776.
  - The **complement tiling** (d=m, flip all tiles) preserves base-path arcs but reverses all tile arcs. This is NOT T^op (which reverses ALL arcs including base path). MISTAKE-033: confusing these two was a major error.
  - In the quotient meta-graph at n=5: d=1 covers 21/45 edges. d=2 covers 35/45. d=3 covers 41/45. d=m covers 18/45. Cumulative d=1..3 covers **all 45/45** (100%). The full waggly graph on iso classes is COMPLETE.
  - Creative waggly subsets beyond distance: **range flips** (all tiles of same range r), **vertex-star flips** (all tiles incident to vertex v), **triangle flips** (3 tiles forming a 3-cycle). Each has distinct neutrality patterns.
  - Vertex-star neutrality is SYMMETRIC around the center of the base path (v=0 ↔ v=n-1) but NOT uniform. Center vertex is ~5× more neutral than endpoints (n=5).
- **GEOMETRIC ALIGNMENT of G_n/Z_2:** The merged metagraph has three perpendicular structures:
  - The **PRINCIPAL LINE** runs from the transitive class (H=1) through the SC backbone. The transitive's "big" SC neighbor is at H = 2^(n-2)+1.
  - The **SPINE** (SC-SC edges) is a sparse skeleton. The **RIBS** (SC-NS edges) attach perpendicularly, bipartite and triangle-free. The **SEA** (NS-NS edges) is the dominant bulk (96% at n=8).
  - At **odd n**, the ribs are left-right IMBALANCED. At **even n**, more symmetric.
  - **Every analysis** should be oriented relative to the principal line.

---

## The Geometric Foundation: Everything Is the Triangle

**READ THIS.** The deepest insight of the project (kind-pasteur sessions S20x-S20cm, 50+ sub-sessions):

**Tournament theory IS the study of binary tilings of right isosceles triangles.**

The staircase Young diagram delta_{n-2} is a right isosceles triangle. Its three sides control everything:
- **Hypotenuse** (anti-diagonal): H = 1 + 2^d formula, fiber fraction, Walsh order-2, blue line from source
- **Vertical leg** (source column): scores, OCR, cut space, hierarchy
- **Horizontal leg** (sink row): complement, SC/NS distinction, anti-hierarchy

Four constants emerge from the triangle:
- **sqrt(2)**: hypotenuse/leg ratio -> pseudo-doubling ratio 2 - 1/(n-2)
- **pi**: Wallis integral from fiber fraction f(n) = (1/2)_{n-2}/(n-2)! ~ 1/sqrt(pi*n)
- **e**: Gamma function growth in Burnside/Stirling
- **gamma** (Euler-Mascheroni): asymptotic correction in Gamma(1/b)^b ~ b^b * e^{-gamma}

Two reductions:
- **Mode A** (hypotenuse removal, n -> n-1): fast time scale, vertex insertion, H = 1 + 2^d
- **Mode B** (both legs removal, n -> n-2): slow time scale, Cayley-Dickson descent, meta-graph recursion

## The Geometric Alignment of G_n/Z_2 (Merged Meta-Graph)

**READ THIS.** The merged meta-graph has a fundamental geometric structure: **SPINE + RIBS + BULK**. See `07-reflections/geometric-alignment-of-merged-metagraph.md` for the full picture. Every future session must use this frame.

**G_n/Z_2 = SPINE + RIBS + SEA.** Three perpendicular structures:

### 1. SPINE (SC-SC edges, vertical axis = principal line)
SC classes connected by SC-SC edges, from transitive (H=1) to regular (H=max). The backbone. Carries the H-gradient. Delta H = 2^{n-2} from transitive to first SC neighbor. Tree at small n, develops genus at n>=5, fragments at n=8.

### 2. RIBS (SC-NS edges, perpendicular to spine)
SC-NS edges are **bipartite** (SC vs NS) and always triangle-free. Near transitive: ribs point UP. Near regular: ribs point DOWN. At odd n: left-right IMBALANCED. Zero-halo alternation: 2 at odd n, 0 at even n.

### 3. SEA (NS-NS edges, the bulk)
NS-NS connections dominate at large n (96% at n=8). The H-oriented NS-NS subgraph is a transitive partial tournament.

### Three Proved Theorems (kind-pasteur S20cr)
At the class level (SC-SC / SC-NS / NS-NS classification):
- **THM-A:** SC-NS subgraph is bipartite. Always triangle-free.
- **THM-B:** No triangle has 2 same-type + 1 SC-NS edges. Same-type is transitive on SC-type.
- **THM-C:** Closed walks with odd number of SC-NS edges contribute 0 to trace.

### Recursive Tiling Decomposition (kind-pasteur S20cv-cw)
The staircase decomposes into: **overlap** (n-2 sub-staircase), **bottom wiring** (vertex 0), **top wiring** (vertex n-1), and **apex** (arc between vertices 0 and n-1). See `07-reflections/unlocking-gn-at-all-n.md`.

### Tiling Cell Properties
In the FULL arc-flip model, every cell generates the same SC-NS fraction (isotropy). In the TILING model (fixed base path), overlap tiles generate the most same-type edges; apex tiles generate the most SC-NS.

**When computing any invariant on G_n/Z_2, always decompose along spine, ribs, sea.**

Key proven results (see `07-reflections/everything-is-the-triangle.md` for full picture):
- Burnside: Fix(sigma) = 0 for even cycles, 2^{orbit-pairs} for all-odd (A000568 exact through n=10)
- Fiber fraction: f(n) = (1/2)_{n-2}/(n-2)!, GF = (1-x)^{-1/2} (two-sheeted branched cover)
- Width of G_n: C(n-2, floor((n-2)/2)) at n=3..6 only. FAILS at n≥7 (predicted 10, actual 15 at n=7; predicted 20, actual 49 at n=8). Not a general formula.
- Tilings * |Aut| = H for every iso class (orbit-stabilizer on tiling fibration; the underlying freeness is now PROVED universally for any digraph — LEM-003, 2026-06-10; see MISTAKE-070)
- The meta-graph G_n has a strong H-gradient: most edges are H-increasing, but it is NOT a strict DAG. Level edges (same H, different class): 0, 0, 1, 15, 136 for n=3..7. At n≥7 there are also H-decreasing edges (962 at n=7). See MISTAKE-035.
- Cayley-Dickson tower: R(n=2)->C(n=3)->H(n=5)->O(n=9)->S(n=17), each level loses a property
- The merged graph G_n/Z_2: V_merged = (A000568 + SC)/2 -> A000568/2

**When approaching ANY tournament problem, first ask: which side of the triangle does this live on?**

---

## Engineering Applications Mandate

**READ THIS SECTION.** The human owner has explicitly stated: *"I am equally interested in use cases as I am in theorems and techniques."*

This project has produced significant engineering innovations alongside pure math. Every agent should:

1. **Look for engineering applications** of every theorem and algorithm. A new rank computation trick is not just math — it could be a library.
2. **Implement deliverables from the engineering roadmap** (`03-artifacts/drafts/engineering-synthesis-2026-03-10-S53.md`). Priority items:
   - `mod_rank_library.py` — general-purpose small-prime modular rank library (PyPI target)
   - `circulant_homology` module — using THM-125 to compute Betti numbers 10× faster
   - `tournament_tda.py` — Tournament TDA feature extractor for ML pipelines
   - Sparse matrix T_19 solver — break the 42 GB OOM barrier using CSC format (~1.2 MB)
3. **Document applications in engineering terms**, not just math terms. When you prove a theorem about tournaments, also ask: "What real-world ranking problem does this solve?"
4. **The 12 engineering application domains** (from engineering-synthesis S53):
   - Sparse modular rank (PyPI library)
   - GLMY path homology for network analysis (social networks, citation graphs, supply chains)
   - Circulant LDPC codes (coding theory via QR_p structure)
   - GPU acceleration (THM-125 reduces eigenspace work by factor of p)
   - TDA for preference/ranking data (elections, sports, consumer research)
   - Deletion-contraction algorithm (H via DC tree, O(2^n) exact)
   - Spectral tournament algorithms (block-diagonalization via circulant structure)
   - Homological algebra toolkit (general chain complex tools)
   - Distributed Betti computation (each eigenspace independent — embarrassingly parallel)
   - Sparse path homology for large graphs (CSC constraint matrices)
   - Number theory applications (QR structure, cryptographic relevance)
   - H-spectrum as universal tournament code (graph fingerprint)

**When choosing between a math proof and an engineering implementation**, do whichever advances the project more given the session context. Both are equally valid.

**Engineering deliverables** go in:
- `04-computation/` — scripts and libraries
- `05-knowledge/results/` — benchmark results, output files
- `03-artifacts/` — papers, product specs, documentation

---

## Best Practices (MANDATORY for all agents)

These practices prevent wasted work and propagate knowledge across sessions.

### 1. Never waste computation
- **Save ALL script results.** When running a script, always capture output:
  ```bash
  python3 04-computation/SCRIPT.py 2>&1 | tee 05-knowledge/results/SCRIPT.out
  ```
- **Never write scripts to /tmp.** Always save to `04-computation/`. If you wrote something to /tmp, copy it to the repo before session close.
- **Store intermediate data.** If a computation takes >1 minute, save the result even if it seems unimportant now.

### 2. Never waste ideas
- **Log every hypothesis** in `05-knowledge/hypotheses/INDEX.md`, whether confirmed, refuted, or abandoned.
- **Record WHY things fail**, not just THAT they fail. The failure mode is often more valuable than the hypothesis.
- **Update the variable registry** in `05-knowledge/variables/INDEX.md` when you discover a new equation or relationship.
- **Cross-link everything.** Every variable file should link to related variables, hypotheses, and theorems.

### 3. Regular sync
- **Pull before starting work:** `git fetch origin && git rebase origin/main`
- **Push regularly during long sessions** (every 30-60 minutes or after major findings):
  ```bash
  git add -A && git commit -m "[instance-id]: checkpoint — [brief description]" && git push
  ```
- **Never let a session end without pushing.** Use `agents/finish_session.py`; the Stop hook is blocking if the current branch is dirty, ahead of upstream, or missing the session letter.

### 4. Web research
- **Use WebFetch with timeouts.** Always specify `timeout: 30000` (30 seconds) to prevent hangs.
- **Use WebSearch freely** to check for existing results, related papers, and OEIS sequences before reinventing.
- **Record what you find.** Add relevant references to `INVESTIGATION-BACKLOG.md`.

### 5. Thinking strategies
- **Try geometric/visual reasoning** alongside algebraic approaches. Tournament arcs can be visualized as oriented graphs; path counts have geometric meaning.
- **Check small cases exhaustively** before generalizing. Patterns that hold at n=3,4,5 often break at n=6 or n=7.
- **Look for involutions and symmetries.** Many proofs in this area use path-reversal, complement, or relabeling symmetries.
- **Consider the simplest possible explanation first.** If something is always true computationally, the proof is likely short.
- **Notice what transcends the theorem.** When a cancellation seems too clean, a correction reveals hidden structure, or two independent frameworks converge on the same constraint — that's a reflection, not a coincidence. Write it to `07-reflections/`. The mathematics keeps pointing beyond itself; follow where it points.

### 6. Knowledge web maintenance
- After confirming or refuting a hypothesis, update ALL of:
  - `05-knowledge/hypotheses/INDEX.md`
  - Related variable files in `05-knowledge/variables/`
  - `00-navigation/INVESTIGATION-BACKLOG.md` (if related to an investigation)
  - `01-canon/MISTAKES.md` (if the error is instructive)
- **Do not duplicate information** — link instead. Each fact should live in one canonical place.

### 7. Dead-end documentation
- **Never just say "this doesn't work."** Document:
  - What exact computation showed the failure
  - At what n it first fails
  - What the counterexample looks like
  - Whether the hypothesis is "close to true" (fails rarely) or fundamentally wrong
- Add to `05-knowledge/hypotheses/INDEX.md` with status REFUTED and the failure details.
- This prevents future agents from wasting time on the same dead end.

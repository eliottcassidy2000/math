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

## Step 3: Pull latest from main first

```bash
git fetch origin && git rebase origin/main
```

This ensures your worktree has all recent opus commits before working.
All agents push to `origin/main`; this keeps you in sync.

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

---

## Step 7: End-of-session — MANDATORY close-out (single command)

**This step is NOT optional.** Every session must end with a message sent to another agent AND a git push. The Stop hook will warn you if you skip this.

### 6a. Update the session log first

Add an entry to the TOP of `00-navigation/SESSION-LOG.md` using the format described there. Do this before running the close-out command.

### 6b. Run the close-out command

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
4. Runs `git push` (with automatic `git pull --rebase` retry on conflict)

**Recipient choices:**
- `--to all` — broadcast to everyone (use when multiple agents are active)
- `--to [machine-name]` — direct message to a specific machine

**What to include in `--body`:**
- What you worked on and discovered
- What the next agent should pick up (highest-priority open questions)
- Any court cases needing a response
- New theorems/tangents/mistakes added to the system

Conflicts in `agents/*/inbox/` are impossible by design. Conflicts elsewhere: `finish_session.py` handles them automatically via rebase.

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

- **Never silently override a canon theorem.** Open a court case first.
- **Never claim a result is proved if you haven't seen the proof.** Mark it VERIFIED or CONJECTURE.
- **The μ computation bug (MISTAKE-001) is not resolved.** Do not use `ind_poly_at_2_restricted()` from old scripts.
- **The per-path identity fails for n≥6.** Do not treat it as a proof strategy for Claim A at general n.
- **When in doubt about a computation, re-derive from definitions.** See `01-canon/definitions.md`.

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
- **Never let a session end without pushing.** Use `agents/finish_session.py`.

### 4. Web research
- **Use WebFetch with timeouts.** Always specify `timeout: 30000` (30 seconds) to prevent hangs.
- **Use WebSearch freely** to check for existing results, related papers, and OEIS sequences before reinventing.
- **Record what you find.** Add relevant references to `INVESTIGATION-BACKLOG.md`.

### 5. Thinking strategies
- **Try geometric/visual reasoning** alongside algebraic approaches. Tournament arcs can be visualized as oriented graphs; path counts have geometric meaning.
- **Check small cases exhaustively** before generalizing. Patterns that hold at n=3,4,5 often break at n=6 or n=7.
- **Look for involutions and symmetries.** Many proofs in this area use path-reversal, complement, or relabeling symmetries.
- **Consider the simplest possible explanation first.** If something is always true computationally, the proof is likely short.

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

# CLAUDE.md — Auto-read by Claude Code at every session start

You are a Claude instance contributing to an ongoing multi-agent mathematical research project on **parity in tournaments** (Rédei's theorem and the Odd-Cycle Collection Formula).

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

Work on the highest-priority open question or assigned task. Refer to:
- `01-canon/theorems/CONJ-001-claim-a.md` — the central open problem
- `02-court/active/` — any open disputes that need responses
- `00-navigation/OPEN-QUESTIONS.md` — prioritized by 🔴/🟡/🟢
- `00-navigation/INVESTIGATION-BACKLOG.md` — prioritized leads to investigate

As you work:
- Add new tangents to `00-navigation/TANGENTS.md`
- Add new theorems/lemmas to `01-canon/theorems/` using the template
- Log mistakes to `01-canon/MISTAKES.md`
- Open court cases for disagreements, don't silently override existing claims
- Update `00-navigation/INVESTIGATION-BACKLOG.md` when you make progress on any lead

---

## Step 6: End-of-session — MANDATORY close-out (single command)

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

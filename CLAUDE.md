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

## Step 3: Process incoming messages

```bash
python3 agents/processor.py --check
```

This will show you:
- Messages in `agents/[your-machine]/inbox/` addressed to you
- Unread messages in `agents/broadcast/` addressed to everyone

Read each message before proceeding. If any message asks you to argue a court case, update a theorem, or continue a chain of reasoning, treat that as your primary task for this session.

---

## Step 4: Process the human intake inbox (if files are present)

```bash
python3 inbox/processor.py
```

If this generates a `inbox/PROCESSING-REPORT.md`, read it and integrate the content before doing other work. New files from the human take priority over exploratory work.

---

## Step 5: Do the actual work

Work on the highest-priority open question or assigned task. Refer to:
- `01-canon/theorems/CONJ-001-claim-a.md` — the central open problem
- `02-court/active/` — any open disputes that need responses
- `00-navigation/OPEN-QUESTIONS.md` — prioritized by 🔴/🟡/🟢

As you work:
- Add new tangents to `00-navigation/TANGENTS.md`
- Add new theorems/lemmas to `01-canon/theorems/` using the template
- Log mistakes to `01-canon/MISTAKES.md`
- Open court cases for disagreements, don't silently override existing claims

---

## Step 6: End-of-session — write letters and commit

Before ending, do all of the following:

### 6a. Write your session letter

```bash
python3 agents/processor.py --send
```

This guides you through writing an end-of-session letter. You'll specify:
- Who to send to (a specific machine name, or `all`)
- Subject and body (your findings, handoffs, questions, court case assignments)

Letters are saved to the recipient's inbox and staged for git commit.

### 6b. Update the session log

Add an entry to the TOP of `00-navigation/SESSION-LOG.md` using the format described there.

### 6c. Commit and push

```bash
git add -A
git commit -m "[your-instance-id]: [one-line summary of session]"
git push
```

If push fails due to another machine having pushed first:
```bash
git pull --rebase
git push
```

Conflicts in `agents/*/inbox/` are impossible by design (each machine only writes to other machines' inboxes, never its own). Conflicts elsewhere should be rare — resolve them by keeping both versions and flagging in SESSION-LOG.md.

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

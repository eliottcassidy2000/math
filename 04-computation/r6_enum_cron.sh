#!/bin/bash
# Cron wrapper for the resumable r=6 fine-branch enumeration (death-star-S58).
# Runs one time-budgeted chunk, then commits+pushes the results/progress files.
# Install:  (crontab -l 2>/dev/null; echo "*/33 * * * * /home/claude/math/04-computation/r6_enum_cron.sh >> /home/claude/math/05-knowledge/results/r6_finebranch_enum.cronlog 2>&1") | crontab -
set -u
REPO=/home/claude/math
LOCK=/tmp/r6_enum_cron.lock
# non-blocking lock: if a previous run is still going, exit quietly
exec 9>"$LOCK"
flock -n 9 || { echo "$(date -u +%FT%TZ) skip: previous run still active"; exit 0; }

cd "$REPO" || exit 1
echo "=== $(date -u +%FT%TZ) r6 enum chunk starting ==="
git fetch origin -q 2>&1 | tail -1
GIT_EDITOR=true git rebase origin/main >/dev/null 2>&1 || git rebase --abort 2>/dev/null

# process one ~20-min chunk
python3 04-computation/r6_enum_cron_runner.py 1200

# commit + push just the enum output/progress (small, low-conflict), rebase-retry
git add 05-knowledge/results/r6_finebranch_enum.out 05-knowledge/results/r6_finebranch_enum.progress 2>/dev/null
if ! git diff --cached --quiet; then
  git commit -q -m "death-star-S58 cron: r=6 fine-branch enum chunk (auto)" \
    -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
  for a in 1 2 3 4 5 6; do
    git fetch origin -q 2>&1 | tail -1
    if GIT_EDITOR=true git rebase origin/main >/dev/null 2>&1; then
      if git push origin main 2>&1 | tail -2 | grep -q "main -> main"; then echo "pushed (try $a)"; break; fi
    else
      git checkout --theirs agents/.session-state.json 2>/dev/null
      git add -A; GIT_EDITOR=true git rebase --continue >/dev/null 2>&1
    fi
  done
fi
echo "=== $(date -u +%FT%TZ) chunk done ==="

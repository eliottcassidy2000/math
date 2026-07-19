#!/bin/bash
# Persistent driver for the r=6 fine-branch enumeration (death-star-S58).
# crontab is unavailable to user 'claude' (root-owned spool), so this nohup loop is the
# cron substitute: it repeatedly runs one time-budgeted chunk (via r6_enum_cron.sh, which
# also commits+pushes the results/progress files), sleeping between chunks, until COMPLETE.
# Launch:  nohup /home/claude/math/04-computation/r6_enum_loop.sh >> \
#          /home/claude/math/05-knowledge/results/r6_finebranch_enum.driverlog 2>&1 &
cd /home/claude/math || exit 1
echo "=== $(date -u +%FT%TZ) driver start ==="
while true; do
  /home/claude/math/04-computation/r6_enum_cron.sh
  if grep -q "^COMPLETE" 05-knowledge/results/r6_finebranch_enum.out 2>/dev/null; then
    echo "=== $(date -u +%FT%TZ) enumeration COMPLETE; driver exiting ==="
    break
  fi
  sleep 90
done

#!/bin/bash
# Persistent, git-free driver for the r=6 fine-branch enumeration (death-star-S58).
# crontab is unavailable to user 'claude' (root-owned spool), so this nohup loop is the cron
# substitute. It loops the CHECKPOINTED runner (r6_enum_cron_runner.py), which resumes from
# 05-knowledge/results/r6_finebranch_enum.progress each time -- so it survives interruption/reboot
# (just relaunch this). It does NO git (avoids racing the fleet); results accumulate locally in
# 05-knowledge/results/r6_finebranch_enum.out -- commit them from any session when convenient.
# Launch:  nohup /home/claude/math/04-computation/r6_enum_loop.sh \
#          >> /home/claude/math/05-knowledge/results/r6_finebranch_enum.driverlog 2>&1 &
cd /home/claude/math || exit 1
echo "=== $(date -u +%FT%TZ) git-free driver start ==="
while true; do
  python3 04-computation/r6_enum_cron_runner.py 1200
  if grep -q "^COMPLETE" 05-knowledge/results/r6_finebranch_enum.out 2>/dev/null; then
    echo "=== $(date -u +%FT%TZ) enumeration COMPLETE; driver exiting ==="
    break
  fi
  sleep 30
done

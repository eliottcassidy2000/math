#!/bin/bash
# After d=7 k<=12: run all checks on two lanes (<= 2 cores), then launch k=13 only if every check passed.
cd "$(dirname "$0")/.."
L=logs/post12.out; : > $L
( runs/verify_d7.sh big1 > logs/verify_d7_big1.out 2>&1
  runs/verify_d7.sh small > logs/verify_d7_small.out 2>&1
  nice -n 10 python3 py/validate.py 6 runs/d6_e0 $(seq 1 15) > logs/validate_d6_e0_rerun.out 2>&1
  nice -n 10 python3 py/validate.py 6 runs/d6_e1 14 15 > logs/validate_d6_e1_rerun.out 2>&1
  nice -n 10 python3 py/validate.py 7 runs/d7_e0 $(seq 1 12) > logs/validate_d7_e0.out 2>&1
  nice -n 10 python3 py/leafcount_dp.py 7 13 > logs/dp_d7_k13.out 2>&1
  nice -n 10 python3 py/leafcount_dp.py 7 14 > logs/dp_d7_k14.out 2>&1 ) &
PA=$!
( runs/verify_d7.sh big2 > logs/verify_d7_big2.out 2>&1
  runs/checks.sh ) &
PB=$!
wait $PA $PB
ok=1
for f in logs/verify_d7_big1.out logs/verify_d7_big2.out logs/verify_d7_small.out; do
  grep -q FATAL $f && { echo "FAIL fatal in $f" >> $L; ok=0; }
  grep -q "VERIFY-DONE .* rc=0" $f || { echo "FAIL no clean finish in $f" >> $L; ok=0; }
done
grep -q "XCHECK ALL AGREE" logs/xcheck.out || { echo "FAIL xcheck" >> $L; ok=0; }
[ $(grep -c "noncanonical=2 wrongsize=0 duplicates=1" logs/repcheck_negative.out) -eq 2 ] || { echo "FAIL repcheck negative control" >> $L; ok=0; }
[ $(grep -c "REPCHECK" logs/repcheck.out) -eq 44 ] || { echo "FAIL repcheck count" >> $L; ok=0; }
grep "REPCHECK" logs/repcheck.out | awk '{for(i=1;i<=NF;i++){split($i,a,"="); if((a[1]=="noncanonical"||a[1]=="wrongsize"||a[1]=="duplicates") && a[2]!="0") bad=1}} END{exit bad}' || { echo "FAIL repcheck (see logs/repcheck.out)" >> $L; ok=0; }
grep "NFTEST" logs/nf_test.out | grep -qv "domain_failures=0 same_orbit=True" && { echo "FAIL nf_test" >> $L; ok=0; }
[ $(grep -c "NFTEST" logs/nf_test.out) -eq 5 ] || { echo "FAIL nf_test count" >> $L; ok=0; }
for f in logs/validate_d6_e0_rerun.out logs/validate_d6_e1_rerun.out logs/validate_d7_e0.out; do
  grep -q "ok=False" $f && { echo "FAIL $f" >> $L; ok=0; }
  grep -q "Traceback" $f && { echo "FAIL traceback in $f" >> $L; ok=0; }
done
[ $(grep -c "ok=True" logs/validate_d7_e0.out) -eq 12 ] || { echo "FAIL validate d7 count" >> $L; ok=0; }
if [ $ok = 1 ]; then
  echo "ALL CHECKS PASSED; launching k=13 $(date)" >> $L
  runs/run_d7.sh 13 13
  echo "K13 RUN FINISHED $(date)" >> $L
else
  echo "CHECKS FAILED; k=13 not launched $(date)" >> $L
fi

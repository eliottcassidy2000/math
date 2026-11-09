"""AUDIT G1: my DP + decoding against klein's census witness lines (W/S lines: 'W h code').
Regenerate the inputs first (from this directory):
  D=../../path_spectrum_two_angles_metagraph_20261106/data
  for n in 9 10; do grep "^W \\|^S " $D/census_$n.txt | awk -v n=$n '{print n, $3, $2}' > wit_census_$n.txt; done
  gzip -dc $D/census_11.txt.gz | grep "^W " | awk 'NR%40==1{print 11, $3, $2}' > wit_census_11.txt
Build: cc -O2 -o aud_core aud_core.c (also aud_dec.c, aud_tau.c, aud_bandlemmas.c the same way).
Run:   python3 chk_census_wit.py wit_census_9.txt wit_census_10.txt wit_census_11.txt
Result (2026-11-09): 19,353 witnesses checked, 0 H mismatches.
"""
import subprocess, sys
bad = tot = 0
for fn in sys.argv[1:]:
    lines = [l.split() for l in open(fn)]
    inp = ''.join('%s %s\n' % (a, b) for a, b, c in lines)
    r = subprocess.run(['./aud_core', 'props'], input=inp, capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(r) == len(lines)
    for (n, code, h), row in zip(lines, r):
        tot += 1
        if int(row.split()[1]) != int(h):
            bad += 1
print('census witnesses checked', tot, 'H mismatches', bad)

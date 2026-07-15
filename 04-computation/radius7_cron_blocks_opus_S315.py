# opus-2026-07-15-S315 -- the radius-7 single-cluster cron blocks.
# Runs sweep_single over prefixes NOT covered by the S314 flagship job,
# appending one line per prefix to the progress file. Order: the two
# remaining exceptional prefixes' completion is handled by the S314 job;
# here: all other 787 prefixes, lexicographic, with a per-prefix budget.
import sys, time, itertools
sys.path.insert(0, '04-computation')
from beat_referee_and_sweeps_opus_S314 import sweep_single

FLAGSHIP = [(2,4,6,8,10), (8,9,10,11,12), (2,5,7,10,12), (4,5,6,11,12), (1,2,3,4,5)]
OUT = r'C:/Users/Eliott/AppData/Local/Temp/claude/C--Users-Eliott-Documents-GitHub-ephrepos-math/fdaae644-876b-48c2-8c49-d1dc2f149e53/scratchpad/radius7_cron_progress.txt'

t0 = time.time()
done = 0
with open(OUT, 'a', encoding='utf-8', newline='\n') as f:
    f.write(f"# cron block start {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.flush()
    for Pset in itertools.combinations(range(1, 13), 5):
        if Pset in FLAGSHIP or Pset == (1,2,3,4,6): continue   # (1,2,3,4,6) partially done, requeue at end
        P = list(Pset)
        R = [r for r in range(1, 13) if r not in Pset]
        np_, npk, nfb, nun, complete = sweep_single(P, R, time_budget=300)
        line = (f"P={P}: patterns={np_} packets={npk} fallbacks={nfb} "
                f"unresolved={nun} complete={complete}")
        f.write(line + "\n"); f.flush()
        print(line, flush=True)
        done += 1

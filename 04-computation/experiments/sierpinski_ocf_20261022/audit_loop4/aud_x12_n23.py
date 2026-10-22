# AUDIT: universal parts at n = 23 (HYP-9371 Section 3): X_1 by the THM-4757 closed form and X_2 by the HYP-9367
# engine (loop_engine_fast.XB, a code path independent of loop_engine4), at z = 1, for P_23 and D23_0.
import sys, json
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE)
import loop_engine_fast as LEF
from second_order import X1_of
from drt_switch import code_to_S
cls = json.load(open(f'{HERE}/drt23_classes.json'))
for k in (12, 0):
    S = code_to_S(cls[k]['code']); Sh = S / 23
    print(f"{cls[k]['label']}: X_1(1) closed form {X1_of(Sh, 1.0):+.10f}; X_2(1) HYP-9367 engine {LEF.XB(Sh, 2, 1.0):+.10f}", flush=True)

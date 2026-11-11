#!/usr/bin/env python3
"""Compare provable error budgets, not exact-cache values."""
import json
from pathlib import Path
import mpmath as mp
from predict import logweight,tail_log,choose_cutoff,stirling
rows=[]
for n in [300,1000,10**4,10**6,10**12,10**100]:
    with mp.workdps(2*len(str(n))+90):
        ln10=mp.log(10)
        lead=logweight(n,(3,))
        rows.append({'n':str(n),
          'identity_log10_model_relative_bound':mp.nstr(tail_log(n,3)/ln10,50),
          'three_cycle_only_log10_model_relative_bound':mp.nstr(tail_log(n,5)/ln10,50),
          'three_cycle_relative_to_correction_error_log10_bound':mp.nstr((tail_log(n,5)-lead)/ln10,50),
          'support12_log10_model_relative_bound':mp.nstr(tail_log(n,13)/ln10,50),
          'stirling4_log10count_absolute_error_bound':mp.nstr(stirling(n)[1]/ln10,50),
          'poisson_two3_enhancement':mp.nstr(64*mp.fprod(n-i for i in range(6))/mp.fprod(n-i for i in range(3))**2,50),
          'adaptive':[{'relative_digits':d,'cutoff_and_success':choose_cutoff(n,d)} for d in [50,200,500,1000,2000,10000]]})
Path(__file__).with_name('method_comparisons.json').write_text(json.dumps(rows,indent=2)+'\n')

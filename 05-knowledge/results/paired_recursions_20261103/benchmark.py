"""Single-process timings: exact four normalized bits, not the full count."""
from time import perf_counter
from pathlib import Path
import json,platform
from verify import old,even_four_bits,even_four_bits_crt,need

out=[]
for n in (200,400,800):
    start=perf_counter();generic,stats=old.normalized_residue(n,4);t1=perf_counter()-start
    start=perf_counter();linear=even_four_bits(n);t2=perf_counter()-start
    start=perf_counter();crt=even_four_bits_crt(n);t3=perf_counter()-start
    need(generic==linear==crt,('benchmark correctness',n))
    out.append({'n':n,'normalized_mod16':crt,'partitions_visited':stats['visited'],
                'partition_seconds':t1,'linear_seconds':t2,'CRT_seconds':t3})
large=[{'n':n,'normalized_mod16':even_four_bits_crt(n)} for n in
       (20402,5510,1498,1500,10**6,10**12,2*101**6)]
report={'python':platform.python_version(),'platform':platform.platform(),
        'scope':'one run; factorizations included; four normalized binary digits only',
        'compared':out,'large_exact_formula_evaluations_not_full_counts':large}
(Path(__file__).parent/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

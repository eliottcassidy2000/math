from coprime_kernel import kernel
from time import perf_counter
import json
for N in (300,600,900):
    S=tuple(3*p for p in (1,3,5,7,11,13,17,19,23,29,31,37,41,43))
    t=perf_counter();a,info=kernel(N,S,3);kt=perf_counter()-t
    t=perf_counter();b,_=kernel(N,S,3,gauge=False);bt=perf_counter()-t
    assert a==b
    print(json.dumps(dict(N=N,gauged_seconds=kt,response_seconds=bt,speedup=bt/kt,exact_equal=True,**info)),flush=True)

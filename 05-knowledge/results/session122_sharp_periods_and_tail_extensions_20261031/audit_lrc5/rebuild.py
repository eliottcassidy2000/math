"""Reconstruct all audit records before a producer receipt is available.

This is independent precomputation; later exact comparison is required.
"""
import sys,json
from itertools import combinations
from pathlib import Path
from check import reconstruct,encode
records=[]
for i,C in enumerate(combinations(range(1,14),8)):
    records.append(reconstruct(C))
    if (i+1)%25==0:print('independently rebuilt',i+1,file=sys.stderr,flush=True)
Path(sys.argv[1]).write_bytes(encode(records))
print(json.dumps({'status':'RECONSTRUCTED_NOT_YET_BOUND','core_count':len(records)}))

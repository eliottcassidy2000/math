"""Matched prefix workloads; measurements are not universal speed claims."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import platform
from statistics import median
from time import process_time,perf_counter

from scheduled_counts import count_prefix

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
path=ROOT/'05-knowledge/results/session105_metagraph_recursion_20261029/enumeration/baseline87.py'
spec=importlib.util.spec_from_file_location('baseline87',path)
baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)
rows=[]
for n in (100,150,200):
    records={name:dict(cpu=[],wall=[]) for name in ('baseline87','scheduled')}
    reference=None
    for repeat in range(3):
        methods=[('baseline87',baseline.degree_bucket_dp),('scheduled',count_prefix)]
        if repeat%2:methods.reverse()
        for name,method in methods:
            wall,cpu=perf_counter(),process_time()
            values,stats=method(n)
            records[name]['cpu'].append(process_time()-cpu)
            records[name]['wall'].append(perf_counter()-wall)
            if reference is not None and values!=reference:
                raise ArithmeticError('benchmark outputs differ')
            reference=values
            records[name]['transitions']=stats['transitions']
    for record in records.values():
        record['median_cpu']=median(record['cpu'])
        record['median_wall']=median(record['wall'])
    rows.append(dict(bound=n,methods=records,
                     prefix_sha256=sha256(','.join(hex(v) for v in reference).encode()).hexdigest()))
    print('completed matched prefix',n,flush=True)
result=dict(python=platform.python_version(),platform=platform.platform(),repeat=3,
            baseline_sha256=sha256(path.read_bytes()).hexdigest(),rows=rows)
(HERE/'benchmarks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

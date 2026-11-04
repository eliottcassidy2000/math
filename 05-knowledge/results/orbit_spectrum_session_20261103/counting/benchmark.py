#!/usr/bin/env python3
"""Same requested leading bits, two workloads: full-prefix baseline vs head query.

CPU times are measurements, not correctness data or complexity claims.
The baseline computes a full prefix, so this is not an equal-output-work
comparison and no ratio is advertised. Its returned head must agree exactly.
"""
import importlib.util,json,platform,statistics,time,hashlib
from pathlib import Path
from heads import head,need
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BASE=ROOT/'05-knowledge/results/session87_four_vertices_enumeration_20261028/fast_counts.py'
need(hashlib.sha256(BASE.read_bytes()).hexdigest()=='3844a4a590daaf01410cf99d3584441bfa2f8e5e95371ccaafe534fd152b0ece', 'inherited baseline hash')
spec=importlib.util.spec_from_file_location('s87',BASE)
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)

def run():
    table=[]
    for n in (32,64,100,160,200):
        fullsamples=[];headsamples=[];target=None
        for repeat in range(5):
            order=('full','head') if repeat%2 else ('head','full')
            for method in order:
                start=time.process_time()
                if method=='full':
                    values,_=b.degree_bucket_dp(n)
                    exact=values[-1]
                    this=(exact.bit_length(),exact>>max(0,exact.bit_length()-128))
                    fullsamples.append(time.process_time()-start)
                else:
                    for _ in range(100):
                        rec=head(n,128,max_core_mass=40)
                    this=(rec['bit_length'],int(rec['head_hex'],16))
                    headsamples.append((time.process_time()-start)/100)
                if target is None: target=this
                need(target==this,'benchmark output mismatch')
        table.append({'n':n,'bits':128,
            'full_prefix_cpu_seconds':fullsamples,'head_query_cpu_seconds':headsamples,
            'median_full_prefix_cpu_seconds':statistics.median(fullsamples),
            'median_head_query_cpu_seconds':statistics.median(headsamples),
            'head_hex':hex(target[1]),'bit_length':target[0]})
    large=[]
    for n in (1000,10000):
        samples=[]
        for _ in range(5):
            start=time.process_time();rec=head(n,256,max_core_mass=0)
            samples.append(time.process_time()-start)
        large.append({'n':n,'bits':256,'cpu_seconds':samples,
            'median_cpu_seconds':statistics.median(samples),
            'bit_length':rec['bit_length'],'head_hex':rec['head_hex'],
            'cutoff':rec['cutoff'],'core_types':rec['core_types']})
    result={'status':'PASS','platform':platform.platform(),
        'python':platform.python_version(),'repetitions':5,
        'head_microbenchmark_calls_per_repetition':100,
        'comparison_scope':'Same requested head; baseline computes complete integer prefix, head API returns only a certified prefix. Not a full-count speedup.',
        'table':table,'large_heads':large,
        'source_hashes':{'heads.py':hashlib.sha256((HERE/'heads.py').read_bytes()).hexdigest(),str(BASE.relative_to(ROOT)):hashlib.sha256(BASE.read_bytes()).hexdigest()}}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':run()

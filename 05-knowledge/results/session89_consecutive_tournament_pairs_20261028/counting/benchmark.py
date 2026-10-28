#!/usr/bin/env python3
"""Finite same-output prefix benchmarks: shared graph versus repeated calls."""
import argparse
import gc
from hashlib import sha256
import json
from pathlib import Path
import platform
import statistics
import sys
import time
import tracemalloc

from baseline_sources import load,SOURCES,COMMIT
from channels import scalar_prefix,fused_prefix


def check(condition,message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    return sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def measure(fn):
    gc.collect()
    w=time.perf_counter();c=time.process_time()
    result=fn()
    return result,{'cpu_seconds':time.process_time()-c,'wall_seconds':time.perf_counter()-w}


def memory(fn):
    gc.collect();tracemalloc.start()
    result=fn()
    current,peak=tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result,{'python_peak_bytes':peak,'python_current_bytes':current}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root',type=Path)
    parser.add_argument('--receipt',type=Path)
    parser.add_argument('--repeats',type=int,default=7)
    parser.add_argument('--bounds',type=int,nargs='+',default=[80,120,160])
    parser.add_argument('--no-memory',action='store_true')
    args=parser.parse_args()
    if args.repeats<1:
        raise ValueError('positive repeats required')
    sys.set_int_max_str_digits(0)
    s87=load('session87',args.baseline_root)
    old=load('session22',args.baseline_root)
    rows=[]
    for n in args.bounds:
        expected,stats=fused_prefix(n)
        def repeated_strongest():
            a,astats=s87.degree_bucket_dp(n)
            b,bstats=scalar_prefix(n,2,1)
            check(astats['rows']==bstats['rows']==stats['rows'],'strong graph mismatch')
            return [a,b]
        def repeated_generalized():
            return [scalar_prefix(n,1,0)[0],scalar_prefix(n,2,1)[0]]
        def repeated_historical():
            return [list(map(int,old.channel_dp(n,r,s)[0])) for r,s in ((1,0),(2,1))]
        funcs={'repeated_strongest_s87':repeated_strongest,
               'repeated_generalized_s87':repeated_generalized,
               'repeated_historical_s22':repeated_historical,
               'fused_two_channels':lambda:fused_prefix(n)[0]}
        samples={name:[] for name in funcs}
        names=list(funcs)
        for rep in range(args.repeats):
            order=names[rep%len(names):]+names[:rep%len(names)]
            if rep%2:
                order=order[::-1]
            for name in order:
                output,timing=measure(funcs[name])
                check(output==expected,f'output mismatch {name} n={n}')
                samples[name].append(timing)
        medians={name:{clock:statistics.median(t[clock] for t in ts)
                       for clock in ('cpu_seconds','wall_seconds')}
                 for name,ts in samples.items()}
        mem={}
        if not args.no_memory:
            for name in ('repeated_strongest_s87','fused_two_channels'):
                output,allocation=memory(funcs[name])
                check(output==expected,'memory probe output mismatch')
                mem[name]=allocation
        rows.append({'bound':n,'workload':'entire prefixes F_n(1,0) and F_n(2,1), n=0..N',
                     'output_sha256':digest(expected),'samples':samples,'medians':medians,
                     'cpu_speedup_vs_strongest':medians['repeated_strongest_s87']['cpu_seconds']/medians['fused_two_channels']['cpu_seconds'],
                     'mathematical_graph_per_channel':stats,
                     'physical_graph_traversals':{'repeated':2,'fused':1},
                     'degree_metadata_calculations':{'repeated':2*stats['precomputed_transfers'],'fused':stats['precomputed_transfers']},
                     'memory':mem})
        print(json.dumps({'bound':n,'medians':medians,'memory':mem}),flush=True)
    # Generic batching is useful, but not presumed faster. Keep its loss visible.
    n=80
    batch=((1,0),(1,1),(2,0),(2,1))
    expected,stats=fused_prefix(n,batch)
    funcs={'repeat_generalized_s87':lambda:[scalar_prefix(n,*c)[0] for c in batch],
           'fused_generic_four':lambda:fused_prefix(n,batch)[0]}
    extra={name:[] for name in funcs}
    names=list(funcs)
    for rep in range(args.repeats):
        for name in (names if rep%2==0 else names[::-1]):
            output,timing=measure(funcs[name]);check(output==expected,'generic batch mismatch')
            extra[name].append(timing)
    generic={'bound':n,'channels':batch,'samples':extra,
             'median_cpu_seconds':{name:statistics.median(t['cpu_seconds'] for t in ts) for name,ts in extra.items()},
             'output_sha256':digest(expected)}
    here=Path(__file__).resolve().parent
    receipt={'status':'PASS','python':sys.version,'platform':platform.platform(),
             'baseline_commit':COMMIT,'baseline_sources':SOURCES,
             'sources_sha256':{name:sha256((here/name).read_bytes()).hexdigest()
                               for name in ('channels.py','benchmark.py','baseline_sources.py')},
             'repeats':args.repeats,'rows':rows,'generic_batch':generic,
             'scope':'finite moderate-prefix same-output implementation comparison; not new records, smaller states, or class generation',
             'memory_scope':'untimed tracemalloc peaks, not process RSS; repeated calls retain first output while computing second'}
    if args.receipt:
        args.receipt.write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()

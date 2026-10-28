#!/usr/bin/env python3
"""Same-input, alternating-order finite benchmarks; CPU and wall kept separate."""
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import argparse
import ast
import gc
import json
import platform
import resource
import statistics
import subprocess
import sys
import time
import tracemalloc

from baseline_sources import configure, source_path, SOURCES

from fast_counts import (degree_bucket_dp, packed_integer_dp, count_prefix,
                         count_direct, palette_prefixes)

HERE=Path(__file__).resolve().parent
sys.set_int_max_str_digits(0)


def load(filename,name):
    spec=spec_from_file_location(name,source_path(filename))
    mod=module_from_spec(spec);spec.loader.exec_module(mod)
    return mod


def digest(value):
    return sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def measure(fn):
    gc.collect()
    wall=time.perf_counter();cpu=time.process_time()
    value=fn()
    return value,{'cpu_seconds':time.process_time()-cpu,
                  'wall_seconds':time.perf_counter()-wall}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repeats',type=int,default=5)
    ap.add_argument('--baseline-root',type=Path)
    ap.add_argument('--receipt',type=Path)
    ap.add_argument('--no-memory',action='store_true')
    args=ap.parse_args()
    configure(args.baseline_root)
    old=load('local_scale.py','old')
    channels=load('channels.py','channels')
    # Preserve the original baseline verbatim. This sensitivity control only
    # computes transition telemetry once per source state instead of per edge.
    source=source_path('local_scale.py').read_text()
    modified=source.replace('            for k,nd,base,factor,mask,offset in tables[degree]:',
        "            transitions+=len(tables[degree])\n            for k,nd,base,factor,mask,offset in tables[degree]:")
    modified=modified.replace('                transitions+=1\n','')
    if modified==source:
        raise RuntimeError('failed to construct telemetry sensitivity control')
    light={};exec(compile(modified,'<telemetry-light local_scale>','exec'),light)
    # Import only definitions/imports from the March script: its top-level
    # demo contains unrelated graph benchmarks and historical prose errors.
    turbo_source=source_path('a000568_turbo.py').read_text()
    tree=ast.parse(turbo_source)
    tree.body=[x for x in tree.body if isinstance(x,(ast.FunctionDef,ast.Import,ast.ImportFrom))]
    turbo={};exec(compile(tree,'<March turbo definitions>','exec'),turbo)
    rows=[]
    for N in (100,150,200):
        expected,expected_stats=old.local_dp(N)
        expected=list(map(int,expected))
        funcs={
            'session22_original':lambda:old.local_dp(N),
            'session22_telemetry_light':lambda:light['local_dp'](N),
            'degree_buckets_no_tail_fold':lambda:degree_bucket_dp(N,terminal_fold=False),
            'degree_buckets_tail_fold':lambda:degree_bucket_dp(N),
        }
        samples={name:[] for name in funcs}
        names=list(funcs)
        for rep in range(args.repeats):
            order=names[rep%len(names):]+names[:rep%len(names)]
            if rep%2:order=order[::-1]
            for name in order:
                (values,stats),timing=measure(funcs[name])
                if values!=expected or stats['rows']!=expected_stats['rows']:
                    raise RuntimeError('benchmark answer or transition mismatch')
                samples[name].append(timing)
        medians={name:{clock:statistics.median(s[clock] for s in times)
                       for clock in ('cpu_seconds','wall_seconds')}
                 for name,times in samples.items()}
        rows.append({'bound':N,'workload':'entire prefix 0..N',
            'transitions':expected_stats['transitions'],'peak_states':expected_stats['peak_states'],
            'precomputed_transfers':expected_stats['precomputed_transfers'],
            'output_sha256':digest(expected),'samples':samples,'medians':medians,
            'cpu_speedup_over_original':medians['session22_original']['cpu_seconds']/medians['degree_buckets_tail_fold']['cpu_seconds'],
            'cpu_speedup_over_telemetry_light':medians['session22_telemetry_light']['cpu_seconds']/medians['degree_buckets_tail_fold']['cpu_seconds']})
        print('completed main bound',N,flush=True)
    palettes=[]
    for N,count in ((80,4),(80,16),(100,4),(100,16)):
        cs=[1<<k for k in range(count)]
        funcs={'shared_integer_K3':lambda:palette_prefixes(N,cs)[0],
               'repeated_degree_bucket':lambda:{c:degree_bucket_dp(N,c)[0] for c in cs},
               'repeated_session22_channels':lambda:{c:list(map(int,channels.channel_dp(N,a=1,b=c.bit_length()-1)[0])) for c in cs}}
        samples={k:[] for k in funcs};expected=None
        for rep in range(args.repeats):
            items=list(funcs.items())
            if rep%2:items.reverse()
            for name,fn in items:
                value,timing=measure(fn)
                if expected is None:expected=value
                if value!=expected:raise RuntimeError('palette mismatch')
                samples[name].append(timing)
        palettes.append({'bound':N,'colors':cs,'workload':'all prefix values for each palette; shared build INCLUDED',
                         'output_sha256':digest(expected),'samples':samples,
                         'median_cpu_seconds':{name:statistics.median(s['cpu_seconds'] for s in ss) for name,ss in samples.items()}})
    historical=[]
    for N in (50,80,100):
        expected=degree_bucket_dp(N)[0][N]
        funcs={'March_turbo_single':lambda:turbo['T_tournament'](N),
               'incremental_integer_single':lambda:count_direct(N)[0],
               'degree_bucket_FULL_PREFIX':lambda:degree_bucket_dp(N)[0][N]}
        samples={k:[] for k in funcs}
        for rep in range(3):
            items=list(funcs.items())
            if rep%2:items.reverse()
            for name,fn in items:
                value,timing=measure(fn)
                if value!=expected:raise RuntimeError('historical baseline mismatch')
                samples[name].append(timing)
        historical.append({'n':N,'target':'a(n); DP ALSO computes 0..n',
                           'samples':samples,'median_cpu_seconds':{name:statistics.median(s['cpu_seconds'] for s in ss) for name,ss in samples.items()}})
    memory=[]
    if not args.no_memory:
        for N in (100,150):
            for name,fn in (('session22_original',old.local_dp),('degree_buckets_tail_fold',degree_bucket_dp)):
                gc.collect();tracemalloc.start()
                values,stats=fn(N)
                current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
                memory.append({'bound':N,'method':name,'tracemalloc_current_bytes':current,
                               'tracemalloc_peak_bytes':peak,'note':'Python tracked allocations only; not process RSS; not a timed run'})
    cb=[]
    binary=Path('/tmp/session87-a000568-gmp-v4')
    command=['cc','-O3','-pthread',str(source_path('a000568_gmp_enum_v4.c')),'-lgmp','-o',str(binary)]
    compile_result=subprocess.run(command,text=True,capture_output=True)
    if compile_result.returncode:
        cb={'status':'BLOCKED','stderr':compile_result.stderr,'command':command}
    else:
        for N in (80,100,120):
            expected=degree_bucket_dp(N)[0][N];samples=[]
            for _ in range(3):
                before=resource.getrusage(resource.RUSAGE_CHILDREN)
                wall=time.perf_counter()
                p=subprocess.run([str(binary),str(N),'1'],text=True,capture_output=True,check=True)
                elapsed=time.perf_counter()-wall
                after=resource.getrusage(resource.RUSAGE_CHILDREN)
                value=int(p.stdout.splitlines()[0].split(' = ')[1])
                if value!=expected or p.stderr:raise RuntimeError('C/GMP reference mismatch or diagnostic')
                samples.append({'cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                                'wall_seconds':elapsed,'raw_summary':p.stdout.splitlines()[-1]})
            cb.append({'n':N,'threads':1,'samples':samples,'target':'a(n) only'})
    result={'status':'PASS','environment':{'python':sys.version,'platform':platform.platform()},
            'methodology':'Fresh construction each call; rotating/alternating order; exact equality outside timer; CPU process_time and wall perf_counter; shared host; no universal speed claim.',
            'main':rows,'palettes':palettes,'historical_python':historical,'memory':memory,
            'historical_C_GMP_v4':cb,'telemetry_light_source_sha256':sha256(modified.encode()).hexdigest(),
            'source_sha256':dict({p.name:sha256(p.read_bytes()).hexdigest() for p in [HERE/'fast_counts.py',HERE/'benchmark.py',HERE/'baseline_sources.py']}, **{SOURCES[name][0]:sha256(source_path(name).read_bytes()).hexdigest() for name in ('local_scale.py','channels.py','a000568_turbo.py','a000568_gmp_enum_v4.c')})}
    output=json.dumps(result,indent=2)+'\n'
    if args.receipt:args.receipt.write_text(output)
    print(output)


if __name__=='__main__':main()

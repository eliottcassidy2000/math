"""Compare the recovered Session87 frontier, not only an older packed baseline."""
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
from time import perf_counter
import sys,json,statistics
from bucket_prefix import prefix
HERE=Path(__file__).resolve().parent;sys.set_int_max_str_digits(0)
def main():
    spec=spec_from_file_location('session87',HERE.parent/'session87_four_vertices_enumeration_20261028/fast_counts.py');old=module_from_spec(spec);spec.loader.exec_module(old)
    for n in range(31):
        a=old.degree_bucket_dp(n)[0]
        for tail in ('fold','convolution'):
            if prefix(n,tail)[0]!=a:raise ValueError('small prefix mismatch')
    reports=[]
    for n in (200,256):
        answer=list(map(int,json.loads((HERE/f'a000568_{n}.json').read_text())));samples=[]
        for repeat in range(3):
            for mode in (('session87','flint_fold','flint_convolution') if repeat%2==0 else ('flint_convolution','flint_fold','session87')):
                start=perf_counter()
                if mode=='session87':a,s=old.degree_bucket_dp(n)
                else:a,s=prefix(n,mode.removeprefix('flint_'))
                seconds=perf_counter()-start
                if a!=answer:raise ValueError('complete prefix mismatch')
                samples.append(dict(mode=mode,repeat=repeat,seconds=seconds,stats=s));print(n,mode,seconds,flush=True)
        reports.append(dict(N=n,median_seconds={mode:statistics.median(x['seconds'] for x in samples if x['mode']==mode) for mode in ('session87','flint_fold','flint_convolution')},samples=samples))
    (HERE/'frontier_comparison.json').write_text(json.dumps(dict(status='FINITE_EXACT_SAME_FULL_PREFIX_WORKLOAD',small_prefix_bound=30,benchmarks=reports),indent=2)+'\n')
if __name__=='__main__':main()

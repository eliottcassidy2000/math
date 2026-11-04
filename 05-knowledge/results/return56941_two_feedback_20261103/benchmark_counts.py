"""Same-host, full-prefix comparison; all outputs and exact controls retained."""
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
from time import perf_counter
from hashlib import sha256
import json,sys,statistics
from count_prefix import prefix,need
HERE=Path(__file__).resolve().parent
sys.set_int_max_str_digits(0)
def inherited():
    spec=spec_from_file_location('old_count',HERE.parent/'tournament_enumeration_20261021a/check.py');m=module_from_spec(spec);spec.loader.exec_module(m);return m
def main():
    old=inherited();checks=0
    for n in range(31):
        a,_=prefix(n,audit=True);b,_=prefix(n,'loop',backend='python',audit=True);c,_=old.packed_dp(n,arithmetic='fraction');need(a==b==c,'small rational full prefixes');checks+=1
    reports=[]
    for n in (200,256):
        samples=[];answer=None
        for repeat in range(3):
            for mode in (('old','new') if repeat%2==0 else ('new','old')):
                start=perf_counter()
                if mode=='old':a,stats=old.packed_dp(n);need(all(x.denominator==1 for x in a),'old integral');a=list(map(int,a))
                else:a,stats=prefix(n,audit=True)
                seconds=perf_counter()-start
                need(answer is None or a==answer,'complete prefix comparison');answer=a
                samples.append(dict(mode=mode,repeat=repeat,seconds=seconds,stats=stats));print(n,mode,seconds,flush=True)
        text=json.dumps([str(x) for x in answer],separators=(',',':'))+'\n';(HERE/f'a000568_{n}.json').write_text(text,encoding='utf8',newline='\n')
        reports.append(dict(N=n,coefficients=n+1,decimal_digits_last=len(str(answer[-1])),sha256=sha256(text.encode()).hexdigest(),median_seconds={mode:statistics.median(s['seconds'] for s in samples if s['mode']==mode) for mode in ('old','new')},samples=samples))
    (HERE/'count_comparison.json').write_text(json.dumps(dict(status='FINITE_EXACT_FULL_PREFIX_AND_LOCAL_TIMING',small_rational_prefixes=checks,benchmarks=reports),indent=2)+'\n')
if __name__=='__main__':main()

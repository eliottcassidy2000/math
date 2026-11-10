#!/usr/bin/env python3
import importlib.util,json,time,platform,sys,argparse
from pathlib import Path
from statistics import median
from hashlib import sha256
from engines import integer_response,core_series,direct_series,require
BASE=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('inherited',BASE/'tournament_future_boundary_20261020/verify.py')
inherited=importlib.util.module_from_spec(spec);spec.loader.exec_module(inherited)
parser=argparse.ArgumentParser();parser.add_argument('--orders',type=int,nargs='+',default=[40,60,80]);parser.add_argument('--repeats',type=int,default=3);args=parser.parse_args()
report={'python':sys.version,'platform':platform.platform(),'orders':[]}
known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168,48542114686912,28401423719122304,31021002160355166848]
require(integer_response(15)[0]==known,'known prefix')
for N in args.orders:
    row={'N':N};answers={}
    for name,fn in [('direct_fraction',lambda:direct_series(N)),('inherited_response_fraction',lambda:inherited.profile_dp(N)),('integer_response',lambda:integer_response(N)),('integer_core',lambda:core_series(N))]:
        durations=[]
        for _ in range(args.repeats):
            start=time.perf_counter();res=fn();durations.append(time.perf_counter()-start)
        if isinstance(res,tuple): answer,stats=res;row[name+'_stats']={k:v for k,v in stats.items() if k!='rows'}
        else:answer=res
        answers[name]=[int(x) for x in answer]
        row[name+'_seconds']=median(durations)
    expected=answers['direct_fraction']
    require(all(a==expected for a in answers.values()),'four engine equality')
    for modulus in [2,16,81,1000,1000000007]:
        got,_=integer_response(N,modulus)
        require(got==[a%modulus for a in expected], 'modular engine')
    for bits in [1,4,16,64]:
        got,_=integer_response(N,binary_precision=bits)
        require(got==[a%(1<<bits) for a in expected], 'binary precision engine')
    row['integer_response_speedup_over_inherited']=row['inherited_response_fraction_seconds']/row['integer_response_seconds']
    row['integer_core_speedup_over_direct']=row['direct_fraction_seconds']/row['integer_core_seconds']
    row['sequence_sha256']=sha256(','.join(map(str,expected)).encode()).hexdigest()
    row['last_term']=str(expected[-1]); report['orders'].append(row)
    print(json.dumps(row),flush=True)
Path(__file__).with_name('benchmark.json').write_text(json.dumps(report,indent=2)+'\n')

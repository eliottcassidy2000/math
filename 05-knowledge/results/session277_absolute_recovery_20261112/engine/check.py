from pathlib import Path
import json,sys,time
from fractions import Fraction
from dyadic_completion import residue,complete,recover_from_interval,integer_interval
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[2]
values=list(map(int,json.loads((ROOT/'session270_chanpak_enumeration_20261111/enumeration/prefix.json').read_text())['values']))
controls=0
for n in range(41):
    for bits in (0,1,2,3,5,8,13,21,55,144):
        for mode in (False,True):
            r,_=residue(n,bits,single_inverse=mode)
            assert r==values[n]%(1<<bits),(n,bits,mode)
            controls+=1
    for support in (0,3,6,9):
        result,stats=complete(n,support,backend="partition")
        assert result==values[n],(n,support)
        controls+=1
for lo,hi,r,b in [(0,8,1,3),(2,3,0,2),(5,4,0,2)]:
    try: recover_from_interval(lo,hi,r,b)
    except ValueError: controls+=1
    else: raise AssertionError('hostile input accepted')
rows=[]
for n in (20,30,40,50,64,80,100):
    for s in (0,3,6,9):
        start=time.perf_counter(); x,stats=complete(n,s,backend="partition")
        assert x==values[n]
        total=time.perf_counter()-start
        row=dict(n=n,support=s,count_bits=x.bit_length(),residue_bits=stats['bits'],
                 real_interval_width=str(int(stats['hi'])-int(stats['lo'])),
                 reconstruction_seconds=total,retained_modular_sectors=stats['leaves'],
                 pruned_branches=stats.get('pruned',0),nodes=stats.get('nodes',0))
        rows.append(row)
    print(n,flush=True)
# Compare inversion strategies at same precision, not against the optimized exact engine.
strategies=[]
for n in (40,64,80,100):
    lo,hi=integer_interval(n,3); bits=(hi-lo).bit_length()
    row=dict(n=n,bits=bits)
    for mode in (False,True):
        start=time.perf_counter(); r,stats=residue(n,bits,single_inverse=mode)
        assert r==values[n]%(1<<bits)
        row['single_inverse' if mode else 'termwise_inverse']=time.perf_counter()-start
    strategies.append(row)
out=dict(exact_controls=controls,rows=rows,inversion_strategy=strategies,
         benchmark_scope='Single runs, same process. Cache is independent Session270 exact prefix. No speed claim versus optimized GMP response engine.')
(Path(__file__).parent/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

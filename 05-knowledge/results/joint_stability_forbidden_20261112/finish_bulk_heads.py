"""Close the two remaining binary-head families with the independent search."""
from pathlib import Path
import hashlib,json
import numpy as np
from census_nine import read,save,need,tables
HERE=Path(__file__).resolve().parent

def run():
    d=read('bulk_head_interfaces.json.gz');new=read('head_target_ten.json.gz')['records'];need(len(d['unresolved'])==2 and len(new)==2,'exact two-interface extension')
    heads={4749:[512,129,11,147,103,79,271,373,575,254],9509:[384,97,43,531,519,345,157,62,222,487]}
    for H,rows in heads.items():need(int(tables(np.array(rows,dtype=np.int64))[0])==H,'literal finite remainder')
    extensions=[]
    for o in d['unresolved']:
        G=1;option=next(q for q in o['options'] if q['factor']==1);pair=option['query'];qi=d['queries_by_index'].index(pair)
        r=next(r for r in new if r['defect']==pair[0] and r['U']//(r['U']&-r['U'])==pair[1]);r={k:v for k,v in r.items() if k not in ['kind','catalogue_position']};r['query_index']=qi
        need(not tables(np.array(r['rows'],dtype=np.int64))[3],'the new source has no unit selection at all')
        j=r['dyadic_shift'];off=option['exponent_offset'];start=max(1,j-off+1);finite=[]
        for k in range(1,start):
            H=o['offset']+o['coefficient']*2**k;need(H in heads,'every finite remainder explicitly built')
            finite.append({'k':k,'H':H,'status':'explicit ten-vertex witness','rows':heads[H]})
        answer={k:v for k,v in o.items() if k!='reason'}
        answer.update(factor=1,query_index=qi,dyadic_shift=j,tail_from_parameter=start,output_order_constant=len(r['rows'])+1+off-j,finite_heads=finite)
        d['resolved'].append(answer);d['interfaces'].append(r);extensions.append({'law':[o['A'],o['B']],'s':o['s'],'interface':r,'tail_from_parameter':start,'finite_heads':finite})
    d['scope']='Every one of540 distinct normalized binary coefficient laws from stable birth clocks through10 has all heads and tails realized at every positive parameter. Baseline finite interface compiler extended by two independently found ten-vertex sources with no unit selection. Two finite remainder counts are verified by explicit adjacency rows. This is a finite library of infinite families, not all odd coefficients.'
    d['baseline_input_sha256']=hashlib.sha256((HERE/'bulk_head_interfaces.json.gz').read_bytes()).hexdigest();d['extra_input_sha256']=hashlib.sha256((HERE/'head_target_ten.json.gz').read_bytes()).hexdigest()
    d['baseline_search_retained_interfaces']=d['retained_interfaces'];d['retained_interfaces']=len(d['interfaces']);d['resolved_head_obligations']=len(d['resolved']);d['unresolved_head_obligations']=0;d['unresolved']=[];d['unresolved_laws']=[]
    d['closed_laws']=sorted({(r['odd_coefficient_A'],r['odd_coefficient_B']) for r in d['families']});d['fully_closed_laws']=len(d['closed_laws']);d['extensions']=extensions
    need(d['fully_closed_laws']==540 and d['resolved_head_obligations']==786,'complete scoped closure')
    save('completed_bulk_heads.json.gz',d);print(json.dumps({'fully_closed_laws':540,'resolved_head_obligations':786,'extensions':extensions}),flush=True)
if __name__=='__main__':run()

"""Independent endpoint replay and literal construction of the bulk repairs."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from census_nine import read,save,need,tables,ib
from construct import tail
from complete_rays import join
HERE=Path(__file__).resolve().parent

def run(completed=False):
    filename='completed_bulk_heads.json.gz' if completed else 'bulk_head_interfaces.json.gz'
    data=read(filename);bykey={(r['query_index'],r['dyadic_shift']):r for r in data['interfaces']};seen=set();literal=component=0;controls=[]
    for ix,o in enumerate(data['resolved']):
        key=(o['query_index'],o['dyadic_shift']);r=bykey[key];rows=r['rows'];O=[v for v in range(len(rows)) if r['selected']>>v&1]
        if key not in seen:
            h,u,w,*_=ib.tc.path_tables(rows);need(h==r['H'] and sum(u[v] for v in O)==r['U'],'independent H/U interface')
            need(h-sum(w[v][z] for i,v in enumerate(O) for z in O[i+1:])==r['defect'],'independent signed endpoint defect');seen.add(key)
        option=next(q for q in o['options'] if q['factor']==o['factor'] and list(q['query'])==data['queries_by_index'][o['query_index']])
        factor=data['factor_rows'][str(o['factor'])]
        for k in [o['tail_from_parameter'],o['tail_from_parameter']+1]:
            t=k+option['exponent_offset']-r['dyadic_shift'];need(t>=1,'positive physical tail length');part=tail(rows,O,t)
            target=o['offset']+o['coefficient']*2**k;expected=o['factor']*(r['defect']+r['U']*2**t);need(target==expected,'exact affine carry normalization')
            whole=join(factor,part) if factor else part
            if len(whole)<=18:
                need(int(tables(np.array(whole,dtype=np.int64))[0])==target,'literal complete repaired head');literal+=1
            elif len(part)<=18:
                need(o['factor']*int(tables(np.array(part,dtype=np.int64))[0])==target,'literal affine component and exact product');component+=1
            controls.append({'A':o['A'],'B':o['B'],'s':o['s'],'k':k,'tail_length':t,'H':target,'order':len(whole)})
        if ix%100==0:print('bulk repaired heads audited',ix,flush=True)
    out={'scope':'Independent full Python endpoint replay of every interface selected by the bulk repair compiler. Every constructed affine identity checked at its start and next parameter; literal entire graph counts through order18 and separate component counts where possible. All-length correctness follows from the exact tail and ordered-product laws.','input_sha256':hashlib.sha256((HERE/filename).read_bytes()).hexdigest(),'independent_interfaces':len(seen),'literal_whole_graphs':literal,'literal_components':component,'algebraic_controls':len(controls),'controls':controls}
    save('audit_'+('completed_' if completed else '')+'bulk_heads.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='controls'}),flush=True)
if __name__=='__main__':run('--completed' in sys.argv)

"""Complete small rooted strong current/next-H fibres, not an all-order source bank."""
from pathlib import Path
import sys,json,gzip
from math import gcd
from collections import defaultdict
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ternary_core_grammar_20261113'))
from triangle_probe import bank,need
from atomic_ears import profile

def save(name,obj):
    suffix='' if __debug__ else '.optimized'
    (HERE/(name+suffix+'.json.gz')).write_bytes(gzip.compress(json.dumps(obj,separators=(',',':')).encode(),mtime=0))

def run():
    records=[];pure=defaultdict(list);fibres=defaultdict(list)
    for row in bank(7):
        T=tuple(row['rows']);z=profile(T,len(T)-1);h,s,q,F=z
        rec=dict(n=len(T),rows=T,root=len(T)-1,profile=z,code=row['ordinary_code'],next_H=h+F-s,next_content=gcd(h-s,q,F))
        i=len(records);records.append(rec);pure[len(T),h,s,F].append(i);fibres[h,h+F-s].append(i)
    pure_groups=[ids for ids in pure.values() if len({records[i]['profile'][2] for i in ids})>1]
    repair_groups=[ids for ids in fibres.values() if 1 in {records[i]['next_content'] for i in ids} and any(records[i]['next_content']>1 for i in ids)]
    pure_repairs=[]
    for ids in pure_groups:
        for a in ids:
            for b in ids:
                x,y=records[a],records[b]
                if x['next_content']>1 and y['next_content']==1:
                    pure_repairs.append([a,b])
    obj=dict(scope='Every rooted strong class through7, no arbitrary-order completeness claim',records=records,pure_deletion_groups=pure_groups,primitive_repair_fibres=repair_groups,pure_deletion_repairs=pure_repairs)
    save('fibre_probe',obj)
    print('ROOTED',len(records),'PURE_GROUPS',len(pure_groups),'REPAIR_FIBRES',len(repair_groups),'PURE_REPAIRS',len(pure_repairs))
    shown=set()
    for a,b in sorted(pure_repairs,key=lambda ab:(records[ab[0]]['n'],records[ab[0]]['next_H'],records[ab[0]]['profile'])):
        x,y=records[a],records[b];key=x['profile'],y['profile']
        if key in shown:continue
        shown.add(key)
        if len(shown)<=30:print('PURE',a,b,x['profile'],'->',y['profile'],'target',x['next_H'],'D',x['next_content'],'->1')
    print('DISTINCT_PROFILE_REPAIRS',len(shown))

if __name__=='__main__':run()

"""Complete root-surplus census through eight ordinary tournament vertices."""
from pathlib import Path
from collections import Counter
import gzip,json,sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_context_recursion_20261113'))
from context_head import ordinary_bank
from feedback_probe import all_profiles

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def run():
    bank=ordinary_bank();byh={};rows=[];levels=Counter();nonone=[]
    for (n,code),T in sorted(bank.items()):
        if n<2:continue
        levels[n]+=1
        for r,z in enumerate(all_profiles(np.array(T,dtype=np.int64))):
            h,s,e,q,F,G,P=map(int,z);delta=F-h
            need(delta>=0,'path-cut injection')
            indegree=n-1-int(T[r]).bit_count()
            item=dict(n=n,code=code,rows=T,root=r,h=h,s=s,q=q,F=F,delta=delta,indegree=indegree)
            record=byh.setdefault(h,dict(h=h,min_positive=None,surpluses=set(),witnesses={}))
            record['surpluses'].add(delta)
            if delta>0 and (record['min_positive'] is None or delta<record['min_positive']):record['min_positive']=delta
            if delta in (1,2):
                record['witnesses'].setdefault(delta,item)
                rows.append(item)
                if indegree!=1:nonone.append(item)
    for rec in byh.values():rec['surpluses']=sorted(rec['surpluses'])
    need(dict(levels)==dict(zip(range(2,9),(1,2,4,12,56,456,6880))),'complete ordinary bank sizes')
    head=[byh[h] for h in sorted(byh) if h<=73]
    result=dict(status='FINITE-EXACT',scope='Complete ordinary types2..8, every labelled root occurrence; all-order upgrade only through separately cited strong-order bound and padding analysis',levels=dict(levels),head=head,small_surplus_occurrences=rows,non_indegree_one=nonone)
    suffix='' if __debug__ else '.optimized'
    (HERE/('surplus'+suffix+'.json.gz')).write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode(),mtime=0))
    print('LEVELS',dict(levels),'SMALL_SURPLUS',len(rows),'NON_INDEGREE_ONE',len(nonone),flush=True)
    print('NONONE_EXAMPLES',nonone[:4],flush=True)
    print('HEAD',[(r['h'],r['min_positive'],list(r['witnesses'])) for r in head],flush=True)

if __name__=='__main__':run()

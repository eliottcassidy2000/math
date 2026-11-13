"""Retain actual vertex correspondences while changing the distinguished root."""
import json,gzip
import numpy as np
from fractions import Fraction
from typed_glue import HERE,need,save,diff,endpoint_type,glue,compile_glue
from feedback_probe import all_profiles


def rooted_pairs(rows):
    X,Y=[all_profiles(np.array(T,dtype=np.int64)) for T in rows]
    out=[]
    for r in range(len(X)):
        x,y=tuple(map(int,X[r])),tuple(map(int,Y[r]));d=diff(x,y)
        for mode in endpoint_type(d):
            out.append(dict(root=r,mode=mode,slope=str(Fraction(d[1 if mode=='E' else 2],d[3])),profiles=[x,y],difference=d,current_H=x[0],next_H=x[0]+x[4 if mode=='E' else 5]-x[1 if mode=='E' else 2]))
    return out


def run():
    data=json.loads(gzip.decompress((HERE/'typed_glue.json.gz').read_bytes()))
    contexts=json.loads(gzip.decompress((HERE/'context_head.json.gz').read_bytes()))['first_rules']
    out=[];words=[]
    for p in data['pairs']:
        states=rooted_pairs(p['rows']);out.append(dict(label=p['label'],old_root=p['root'],rows=p['rows'],states=states))
        print('REROOT',p['label'],[(s['root'],s['mode'],s['slope'],s['next_H']) for s in states],flush=True)
        if p['label']!='2/3>1/3':continue
        for s in states:
            if s['mode']!='E' or s['root']==p['root']:continue
            for c in contexts:
                if c['input_slope']!=s['slope'] or c['output_slope'] in ('0','1'):continue
                rows=[glue(c['rows'],c['root'],T,s['root']) for T in p['rows']]
                after=rooted_pairs(rows)
                matching=[v for v in after if v['root']==c['root'] and v['mode']=='E']
                need(len(matching)==1 and matching[0]['slope']==c['output_slope'],'physical reroot/context composition')
                words.append(dict(source=p['label'],reroot=s,context=c,rows=rows,states=after))
                print('WORD',s['slope'],'>',c['output_slope'],'ORDER',len(rows[0]),'H',after[0]['current_H'],'STATES',[(v['root'],v['mode'],v['slope']) for v in after],flush=True)
    save('reroot_probe',dict(status='FINITE-EXACT',scope='Every aligned root of the eleven inherited literal pairs; four nonzero consumers of the newly exposed 1/2 root on the H605 pair, all resulting roots counted directly.',pairs=out,words=words))


if __name__=='__main__':run()

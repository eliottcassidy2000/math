"""Lift numerical neutral directions to actual pairs and actual context words."""
import json,gzip
from fractions import Fraction
from collections import defaultdict
from context_transport import glue
from repair23 import profile
from fibre_probe import HERE,need,save

def run():
    contexts=json.loads(gzip.decompress((HERE/'direction_probe.json.gz').read_bytes()))['examples']
    records=json.loads(gzip.decompress((HERE/'fibre_probe.json.gz').read_bytes()))['records']
    groups=defaultdict(list)
    for i,row in enumerate(records):groups[row['profile'][0],row['next_H']].append(i)
    pairs={}
    for ids in groups.values():
        for a in ids:
            for b in ids:
                x,y=records[a],records[b];v=y['profile'][1]-x['profile'][1];w=y['profile'][2]-x['profile'][2]
                if not w:continue
                key=str(Fraction(v,w))
                if key not in contexts:continue
                rank=(max(x['n'],y['n']),x['n']!=y['n'],x['profile'][0],a,b)
                if key not in pairs or rank<pairs[key][0]:pairs[key]=(rank,a,b)
    outputs=[]
    for key,(_,a,b) in sorted(pairs.items()):
        A=contexts[key];left,right=records[a],records[b]
        T=glue(A['rows'],A['root'],left['rows'],left['root']);U=glue(A['rows'],A['root'],right['rows'],right['root'])
        x,y=profile(T,A['root']),profile(U,A['root']);w=right['profile'][2]-left['profile'][2]
        expected=(0,A['profile'][1]*w,A['profile'][2]*w,A['profile'][1]*w)
        need(tuple(v-u for u,v in zip(x,y))==expected,'physical response direction')
        need(x[0]==y[0] and x[3]-x[1]==y[3]-y[1],'two exact observations')
        outkey=str(Fraction(A['profile'][1],A['profile'][2]))
        obj=dict(input_slope=key,output_slope=outkey,input_indices=[a,b],input_profiles=[left['profile'],right['profile']],context=A,output_profiles=[x,y],output_rows=[T,U],root=A['root'])
        if outkey in contexts and key not in ('0','1'):
            C=contexts[outkey];V=glue(C['rows'],C['root'],T,A['root']);W=glue(C['rows'],C['root'],U,A['root']);vz,wz=profile(V,C['root']),profile(W,C['root'])
            need(vz[0]==wz[0] and vz[3]-vz[1]==wz[3]-wz[1],'literal second neutral transport')
            obj['second']=dict(rows=[V,W],profiles=[vz,wz],root=C['root'],context_slope=outkey,output_slope=str(Fraction(C['profile'][1],C['profile'][2])))
        outputs.append(obj)
        print('LIFT',key,'->',outkey,'H',x[0],'NEXT',x[0]+x[3]-x[1],'INPUT',left['profile'],right['profile'],'OUTPUT',x,y,'SECOND',obj.get('second',{}).get('profiles'))
    save('transport_words',dict(status='PASS',scope='Matches every observed context slope to current/next-count fibres in the complete rooted strong bank through7; selected actual pairs and one/two gluing words literally recounted',missing_slopes=sorted(set(contexts)-set(pairs)),outputs=outputs))
    print('MISSING',sorted(set(contexts)-set(pairs)))

if __name__=='__main__':run()

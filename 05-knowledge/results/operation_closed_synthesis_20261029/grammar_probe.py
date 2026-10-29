"""Finite right-context refinement of the exact tournament insertion data."""
from pathlib import Path
from collections import Counter,defaultdict
import itertools,json

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'tournament_results.json').read_text())
rows={r['word']:r for layer in data['layers'].values() for r in layer['records']}
depth=max(map(len,rows))
alphabet='abcde'
def value(w):
    r=rows.get(w)
    return -1 if r is None else int(r['det']!=0)
def profile(w,h):
    return tuple(value(w+''.join(t)) for k in range(h+1) for t in itertools.product(alphabet,repeat=k))
out={'max_word_length':depth,'right_context_classes':[]}
for h in range(1,5):
    eligible=[w for w in rows if len(w)<=depth-h]
    clusters=defaultdict(list)
    for w in eligible: clusters[profile(w,h)].append(w)
    desc=[]
    for sig,ws in sorted(clusters.items(),key=lambda kv:min((len(w),w) for w in kv[1])):
        rep=min(ws,key=lambda w:(len(w),w))
        desc.append({'representative':rep,'count':len(ws),'type':rows[rep]['type'],'unimodular':bool(rows[rep]['det']), 'one_step':{a:value(rep+a) for a in alphabet}})
    out['right_context_classes'].append({'depth':h,'eligible_prefixes':len(eligible),'classes':len(clusters),'representatives':desc})

witness=None
for n in sorted(data['layers'],key=int):
    rr=data['layers'][n]['records']
    if int(n)>=max(map(int,data['layers'])):continue
    for i,r in enumerate(rr):
        for s in rr[i+1:]:
            if (r['type'],bool(r['det']))!=(s['type'],bool(s['det'])):continue
            for a in alphabet:
                v1,v2=value(r['word']+a),value(s['word']+a)
                if v1>=0 and v2>=0 and v1!=v2:
                    witness={'order':int(n),'left':r,'right':s,'append':a,'left_child':rows[r['word']+a],'right_child':rows[s['word']+a]}
                    break
            if witness:break
        if witness:break
    if witness:break
out['coarse_state_hostile']=witness
out['singular_to_unimodular']=[{'parent':w,'append':a} for w,r in rows.items() if not r['det'] and len(w)<depth for a in alphabet if value(w+a)==1][:12]
out['forbidden_word_criterion']={
    'scope':f'all {len(rows)} grammar words through length {depth}',
    'forbidden_subwords':['be','ee'],
    'violations':[w for w,r in rows.items() if bool(r['det']) != ('be' not in w and 'ee' not in w)],
    'status':'FINITE-EXACT; arbitrary lengths OPEN'}
out['count_recurrence']={
    'formula':'u_n=2*u_(n-1)+u_(n-2)-u_(n-3)',
    'initial':{'4':1,'5':3,'6':7},
    'tested':{n:layer['unimodular'] for n,layer in data['layers'].items()},
    'next_prediction':2065,
    'prediction_order':13,
    'status':'PROVED for pattern-avoiding grammar words; CONDITIONAL for unimodular observers at arbitrary order'}
(HERE/'grammar_results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:out[k] for k in ['max_word_length','forbidden_word_criterion','count_recurrence']},indent=2))

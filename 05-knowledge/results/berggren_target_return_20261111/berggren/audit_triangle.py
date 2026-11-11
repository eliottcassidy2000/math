"""Independent local boundary counts and availability walks for the C2 proof."""
from itertools import product, permutations
from pathlib import Path
import hashlib,json
from check import need

HERE=Path(__file__).resolve().parent


def local(q,fixed=None,endpoint='start'):
    total=0
    for owners in product(range(q),repeat=3):
        bins=[[v for v in range(3) if owners[v]==i] for i in range(q)]
        for orders in product(*(list(permutations(b)) for b in bins)):
            if any(any((a+1)%3!=b for a,b in zip(p,p[1:])) for p in orders):continue
            if fixed is not None:
                p=orders[0]
                if not p or p[0 if endpoint=='start' else -1]!=fixed:continue
            total+=1
    return total


def walks(start,S,E):
    rows=[]
    for word in sorted(set(permutations('S'*S+'E'*E))):
        q=start;counts=[];zero=[];valid=True
        for i,c in enumerate(word):
            if c=='S':q+=1;counts.append(q)
            else:counts.append(q);q-=1
            if q<0:valid=False;break
            if q==0 and i<len(word)-1:zero.append((i,i+1))
        if valid:rows.append({'word':''.join(word),'event_available':counts,'zero_intervals':zero,'end':q})
    return rows


def main():
    counts={str(q):{'unmarked':local(q),'fixed_start':local(q,0),'fixed_end':local(q,0,'end')} for q in(1,2,3)}
    need([counts[str(q)]['fixed_start'] for q in(1,2,3)]==[1,3,7],'exact marked local counts')
    need([counts[str(q)]['unmarked'] for q in(1,2,3)]==[3,12,33],'rotation divisible full local counts')
    need(all(v['fixed_start']==v['fixed_end'] for v in counts.values()),'endpoint reversal')
    four=walks(1,2,2)
    need({r['word'] for r in four}=={'SSEE','SESE','SEES','ESSE','ESES'},'complete four-event walks')
    for r in four:
        if r['word']=='ESES':need(r['zero_intervals']==[(0,1),(2,3)],'two disjoint zero intervals')
        else:need(r['event_available'].count(2)>=2,'at least two local factors of three')
    three_start=walks(0,2,1);three_end=walks(1,1,2)
    for r in three_start+three_end:
        need(2 in r['event_available'] or r['zero_intervals'],'three-event factor or forbidden zero interval')
    source=HERE.parent/'arithmetic/TRIANGLE_C2_PROOF.md'
    receipt={'status':'PASS','scope':'Independent exact local-cover enumeration and complete event walks, supplementing analytic factorization audit; not a full all-graph polynomial rerun.',
             'local_counts':counts,'four_event_walks':four,'three_event_start':three_start,'three_event_end':three_end,
             'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
    (HERE/'triangle_audit.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print('PASS: local counts1,3,7; unmarked counts3,12,33; all four/three-event availability walks')


if __name__=='__main__':main()

"""Independent abstract half-contact budget audit for the all-k proof.

The relaxed model includes non-realizable contact patterns; passing it is a
hostile control, not a substitute for the coefficient factorization proof.
"""
from itertools import product
from pathlib import Path
from collections import Counter
import random,json,hashlib
from check import need
from audit_triangle import local

HERE=Path(__file__).resolve().parent


def valuation3(n):
    out=0
    while n and n%3==0:n//=3;out+=1
    return out


def make_blocks(start,word,cuts,kinds,gaps):
    chunks=[];last=0
    for end in list(cuts)+[len(word)]:chunks.append(word[last:end]);last=end
    q=start;blocks=[]
    for i,(chunk,kind) in enumerate(zip(chunks,kinds)):
        if kind=='S' and (len(chunk)>2 or (len(chunk)==2 and chunk!='SE')):return None
        before=q
        for c in chunk:
            q+=1 if c=='S' else -1
            if q<0:return None
        blocks.append({'type':kind,'m':len(chunk),'before':before,'after':q})
        if i<len(gaps):
            for extra in gaps[i]:
                if q==0:return None
                blocks.append({'type':extra,'m':0,'before':q,'after':q})
    if q not in(0,1):return None
    return blocks


def audit(blocks,start):
    end=blocks[-1]['after'];m=sum(b['m'] for b in blocks)
    endpoint_slack=(1-start)+(1-end)
    need((m+endpoint_slack)%2==0,'integral relaxed birth budget')
    k=(m+endpoint_slack)//2
    contacts=[i for i,b in enumerate(blocks) if b['m']]
    rank={i:j for j,i in enumerate(contacts)}
    regions=[]
    for j in range(len(contacts)+1):
        left=contacts[j-1] if j else None
        right=contacts[j] if j<len(contacts) else None
        q=blocks[left]['after'] if left is not None else start
        interior=list(range(0 if left is None else left+1,len(blocks) if right is None else right))
        regions.append({'left':left,'right':right,'q':q,'interior':interior,'requests':[]})
    loss=0
    for i,b in enumerate(blocks):
        if b['type']!='T':continue
        if b['m']>=2:loss+=1
        elif b['m']==1:
            q=max(b['before'],b['after'])
            need(q*q-q+1>0,'positive one-contact local count')
            if q%3!=2:
                loss+=1
                j=rank[i] if b['before']%3==0 else rank[i]+1
                need(regions[j]['q']%3==0,'request points to divisible region')
                regions[j]['requests'].append(i)
    singleton_charges=Counter();terminal=0;unpaid=[];free_paid=0
    for reg in regions:
        requests=reg['requests'];need(len(requests)<=2,'at most two requests per region')
        if not requests:continue
        if any(blocks[i]['type']=='S' for i in reg['interior']):
            need(reg['q']>0 and reg['q']%3==0,'positive singleton payment')
            free_paid+=len(requests);continue
        endpoints=[i for i in(reg['left'],reg['right']) if i is not None]
        singles=[i for i in endpoints if blocks[i]['type']=='S']
        if singles:
            need(len(requests)==1 and len(singles)==1,'one requester opposite a contact singleton')
            singleton_charges[singles[0]]+=1;continue
        if len(endpoints)<2:
            need(reg['q']==0,'terminal request uses a new endpoint')
            terminal+=len(requests);continue
        unpaid.extend(requests)
    need(all(charges<=blocks[i]['m'] for i,charges in singleton_charges.items()),'singleton contact budget')
    need(terminal<=endpoint_slack,'terminal slack budget')
    runs=[];current=[]
    for i,b in enumerate(blocks):
        if b['type']=='T':current.append(i)
        elif current:runs.append(current);current=[]
    if current:runs.append(current)
    nu=sum(len(run)//2 for run in runs)
    for run in runs:
        requests=sum(i in unpaid for i in run)
        need(requests<=2*(len(run)//2),'odd-run saturation cannot leave one unmatched half-unit')
    free_values=sum(valuation3(b['before']) for b in blocks if b['type']=='S' and b['m']==0)
    need(loss-free_values<=k+nu,'global half-contact lower bound')
    return {'k':k,'nu':nu,'requests':sum(len(r['requests']) for r in regions),
            'free_singleton_paid_halves':free_paid,'contact_singleton_paid_halves':sum(singleton_charges.values()),
            'terminal_paid_halves':terminal,'unpaid_halves':len(unpaid),'loss_minus_free':loss-free_values}


def main():
    local_controls=[]
    for q in range(1,13):
        fixed=local(q,0);unmarked=local(q)
        need(fixed==q*q-q+1 and unmarked==q*(q*q+2),'independent arbitrary-q local cover controls')
        need(valuation3(fixed)==int(q%3==2),'exact one-contact valuation')
        local_controls.append({'q':q,'fixed_endpoint':fixed,'unmarked':unmarked})
    checked=0;payments=Counter()
    for m in range(1,6):
        for chars in product('SE',repeat=m):
            word=''.join(chars)
            for mask in range(1<<(m-1)):
                cuts=[i+1 for i in range(m-1) if mask>>i&1];b=len(cuts)+1
                for kinds in product('TS',repeat=b):
                    for gaps in product(('', 'S'),repeat=b-1):
                        for start in(0,1):
                            blocks=make_blocks(start,word,cuts,kinds,gaps)
                            if blocks is None:continue
                            result=audit(blocks,start);checked+=1
                            for key in('free_singleton_paid_halves','contact_singleton_paid_halves','terminal_paid_halves','unpaid_halves'):
                                payments[key]+=result[key]
    rng=random.Random(1817425);random_checks=0
    for _ in range(20000):
        m=rng.randrange(1,41);start=rng.randrange(2);q=start;word=[]
        for i in range(m):
            c='S' if q==0 or rng.randrange(2) else 'E';word.append(c);q+=1 if c=='S' else -1
        if q not in(0,1):continue
        cuts=[i for i in range(1,m) if rng.randrange(3)==0];b=len(cuts)+1
        kinds=[rng.choice('TTS') for _ in range(b)]
        gaps=[rng.choice(('','T','TT','S','TS','ST')) for _ in range(b-1)]
        blocks=make_blocks(start,''.join(word),cuts,kinds,gaps)
        if blocks is not None:audit(blocks,start);random_checks+=1
    high=make_blocks(1,'SSESESEE',list(range(1,8)),['T']*8,['S']*7)
    high_result=audit(high,1)
    need(high_result['free_singleton_paid_halves']==6,'high-q hostile is paid by three positive gaps')
    source=HERE.parent/'arithmetic/TRIANGLE_ALL_K_PROOF.md'
    data={'status':'PASS','scope':'Exhaustive relaxed patterns through5 contacts, arbitrary contact-block grouping and singleton boundary blocks, binary free-singleton gaps, both terminal states; seeded larger patterns with free triangles and singleton gaps. Not a proof by finite checking.',
          'local_controls':local_controls,'exhaustive_patterns':checked,'payment_totals':dict(payments),'seeded_larger_patterns':random_checks,
          'high_availability_hostile':high_result,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
    (HERE/'all_k_audit.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','exhaustive_patterns':checked,'seeded_larger_patterns':random_checks,'high_q':high_result}))


if __name__=='__main__':main()

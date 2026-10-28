"""A linear-size, operation-closed blue interface for singleton contexts."""
from pathlib import Path
import json
import reconstruct as r

HERE=Path(__file__).resolve().parent
MIXED=2;EMPTY=3

def trit(values):
    s=set(values)
    return EMPTY if not s else next(iter(s)) if len(s)==1 else MIXED

def append_one(x):return 1 if x in (1,EMPTY) else MIXED

def raw(t):
    n=len(t)
    rows=tuple(trit(r.core.edge(t,i,j) for j in range(i+2,n)) for i in range(n))
    cols=tuple(trit(r.core.edge(t,i,j) for i in range(j-1)) for j in range(n))
    prefixes=tuple(r.blue(tuple(row&((1<<k)-1) for row in t[:k])) for k in range(n+1))
    suffixes=tuple(r.blue(tuple(row>>(n-k) for row in t[n-k:])) for k in range(n+1))
    return rows,cols,prefixes,suffixes

def update_raw(profile,left):
    rows,cols,P,S=profile;n=len(rows);N=n+1
    if left:
        rr=((1 if n>=2 else EMPTY),)+rows
        cc=(EMPTY,)+(tuple(append_one(x) if j else EMPTY for j,x in enumerate(cols)))
        pp=tuple(True if k<3 else cols[k-2] in (1,EMPTY) and P[k-2] for k in range(N+1))
        ss=S+(pp[-1],)
    else:
        rr=tuple(append_one(x) if i<n-1 else EMPTY for i,x in enumerate(rows))+(EMPTY,)
        cc=cols+((1 if n>=2 else EMPTY),)
        ss=tuple(True if k<3 else rows[n-k+1] in (1,EMPTY) and S[k-2] for k in range(N+1))
        pp=P+(ss[-1],)
    return rr,cc,pp,ss

def compress(profile):
    rows,cols,P,S=profile;n=len(rows);rr=[];cc=[]
    for x in rows[:max(0,n-2)]:
        if x==MIXED:break
        rr.append(x)
    for x in reversed(cols[2:]):
        if x==MIXED:break
        cc.append(x)
    return n,tuple(rr),tuple(cc),tuple(S[n-d] for d in range(1,len(rr)+1)),tuple(P[n-d] for d in range(1,len(cc)+1)),P[n]

def inflate(profile):
    """Unused entries are formal placeholders, not an asserted tournament."""
    n,rr,cc,plus,minus,b=profile
    rows=[MIXED]*max(0,n-2)+[EMPTY]*min(2,n)
    cols=[EMPTY]*min(2,n)+[MIXED]*max(0,n-2)
    P=[True]*min(3,n+1)+[False]*max(0,n-2);S=P.copy()
    P[n]=S[n]=b
    for i,x in enumerate(rr):rows[i]=x
    for d,x in enumerate(cc,1):cols[n-d]=x
    for d,x in enumerate(plus,1):S[n-d]=x
    for d,x in enumerate(minus,1):P[n-d]=x
    return tuple(rows),tuple(cols),tuple(P),tuple(S)

def signature(t):return compress(raw(t))

def join(profile,left):return compress(update_raw(inflate(profile),left))

def complement(profile):
    n,rr,cc,plus,minus,b=profile
    return n,tuple(1-x for x in rr),tuple(1-x for x in cc),plus,minus,b

def main():
    result=[]
    for n in range(1,8):
        N=1<<((n-1)*(n-2)//2);profiles=set()
        for x in range(N):
            t=r.core.presentation(n,x);full=raw(t);s=compress(full);profiles.add(s)
            r.need(signature(r.core.complement(t))==complement(s),'compressed tau update')
            for left in (False,True):
                literal=raw(r.adjoin(t,left,False))
                r.need(update_raw(full,left)==literal,'raw boundary update')
                r.need(join(s,left)==compress(literal),'compressed boundary update')
        result.append(dict(n=n,presentations=N,boundary_profiles=len(profiles),updates=3*N))
    a,b=(r.core.presentation(7,x) for x in (1072,4176))
    r.need(signature(a)==signature(b)==(7,(),(),(),(),False),'hostile witness has same closed blue interface')
    bank=json.loads((HERE/'rooted-seven.json').read_text());cells=bank['cells']-len(bank['nontrivial'])
    cells+=sum(len({signature(r.core.presentation(7,x)) for x in group}) for group in bank['nontrivial'])
    r.need(cells==32272,'complete all-size sufficient carrier at seven')
    # Inherited all-presentation census through six, plus the five n=6
    # horizon-eight separators, verifies sharpness without assuming it.
    banks,small=r.old.closure(7)
    r.need([x['bounded_attachment_cells'] for x in small]==[1,1,2,8,64,1019,24754],'inherited horizon seven')
    remaining=small[5]['order_eight_residual_probes']
    r.need(len(remaining)==5 and all(x['separating_context'] for x in remaining),'all five remaining pairs separated')
    receipt=dict(status='FINITE-EXACT',profile_updates=result,witness_profile=signature(a),seven_closed_carrier_cells=cells,
                 inherited_sharpness=small,scope='Raw and compressed blue signatures; all presentations through seven, all three one-step operations.')
    (HERE/'boundary-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result));return receipt

if __name__=='__main__':main()

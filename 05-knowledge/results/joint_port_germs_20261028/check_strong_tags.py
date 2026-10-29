"""C3-only selector completeness via intrinsically recoverable strong modules."""
from itertools import combinations,permutations
from pathlib import Path
import json
import port_core as c
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,m):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(m)
def atom(m):
    if m%2!=1 or m<3:raise ValueError('odd strong tag order')
    t=(0,)
    for _ in range((m-1)//2):t=c.blow(t,0)
    return t
def strong(t):
    n=len(t)
    for start in range(n):
        seen=1<<start
        while True:
            more=seen
            for i in range(n):
                if seen>>i&1:more|=t[i]
            if more==seen:break
            seen=more
        if seen!=(1<<n)-1:return False
    return True
def substitute(t,i,j):
    n=len(t);p=n+1+(n%2);q=2*n+1;parts={i:atom(p),j:atom(q)}
    labels=[(a,x) for a in range(n) for x in range(len(parts.get(a,(0,))))]
    out=tuple(sum((c.edge(parts[a],x,y) if a==b else c.edge(t,a,b))<<k
                  for k,(b,y) in enumerate(labels) if (a,x)!=(b,y)) for a,x in labels)
    return out,[a for a,x in labels],p,q
def modules_of_size(t,k):
    n=len(t);full=(1<<n)-1;result=[]
    for chosen in combinations(range(n),k):
        mask=sum(1<<i for i in chosen);outside=full^mask;ok=True
        while outside:
            bit=outside&-outside;outside^=bit;i=bit.bit_length()-1
            if (t[i]&mask) not in (0,mask):ok=False;break
        if ok:result.append(chosen)
    return result
def contract(t,blocks,part):
    part=set(part);remaining=[i for i in range(len(t)) if i not in part]
    representatives=remaining+[min(part)]
    newblocks=[blocks[i] for i in remaining]+[tuple(v for i in sorted(part) for v in blocks[i])]
    return c.induced(t,representatives),newblocks,len(remaining)
def decode(t,n):
    p=n+1+(n%2);q=2*n+1
    blocks=[(i,) for i in range(len(t))]
    big=modules_of_size(t,q);need(len(big)==1,'unique larger strong tag')
    a,blocks,qmark=contract(t,blocks,big[0])
    qlabel=blocks[qmark][0]
    small=modules_of_size(a,p);need(len(small)==1,'unique smaller strong tag after contraction')
    b,blocks,pmark=contract(a,blocks,small[0])
    qmark=next(i for i,block in enumerate(blocks) if qlabel in block)
    return b,(pmark,qmark),blocks
def main():
    rows=[]
    for n in range(2,6):
        size=1<<((n-1)*(n-2)//2)
        masks=range(size) if n<=4 else (0,1,7,size-1)
        pairs=[(i,j) for i in range(n) for j in range(n) if i!=j] if n<=4 else [(0,4),(1,3)]
        tested=0
        for mask in masks:
            t=c.presentation(n,mask)
            for i,j in pairs:
                encoded,labels,p,q=substitute(t,i,j)
                need(strong(atom(p)) and strong(atom(q)),'actual C3-built tags strong')
                shift=(mask+i+j)%len(encoded);order=tuple(range(shift,len(encoded)))+tuple(range(shift))
                actual,marks,blocks=decode(c.relabel(encoded,order),n)
                kinds=[{labels[order[v]] for v in block} for block in blocks]
                need(all(len(k)==1 for k in kinds),'each contracted block has one old owner')
                ids=[next(iter(k)) for k in kinds];need(sorted(ids)==list(range(n)),'all old vertices exactly once')
                need([ids[k] for k in marks]==[i,j],'ordered ports intrinsically decoded')
                need(actual==c.relabel(t,ids),'all old arcs recovered')
                need(len(encoded)==4*n+n%2,'explicit C3 context horizon')
                if n<=3:
                    rooted=sum(g[i]==i and g[j]==j for g in c.oracle.isomorphisms(t,t,True))
                    need(len(c.oracle.isomorphisms(encoded,encoded,True))==9*rooted,'two C3-chain automorphism kernels')
                tested+=1
        rows.append({'order':n,'presentations':len(masks),'port_pairs':len(pairs),'tests':tested,
                     'final_order':4*n+n%2,'scope':'complete' if n<=4 else 'four chosen presentations and two ordered pairs'})
        print(json.dumps(rows[-1]),flush=True)
    out={'status':'FINITE-EXACT','rows':rows,'checks':CHECKS}
    (HERE/'strong_tags.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','checks':CHECKS}))
if __name__=='__main__':main()

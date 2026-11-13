"""Target inversion with an arbitrary marked core; readiness defect is retained."""
from pathlib import Path
import importlib.util,json,math
from itertools import product
import inverse as inv

HERE=Path(__file__).resolve().parent

def labelled(n):
    for mask in range(1<<(n*(n-1)//2)):
        T=[0]*n;i=0
        for a in range(n):
            for b in range(a+1,n):
                u,v=(a,b) if mask>>i&1 else (b,a);T[u]|=1<<v;i+=1
        yield tuple(T)

def bank():
    seen=set();out=[]
    for n in range(3,6):
        for T in labelled(n):
            if not inv.f.a.strong(T):continue
            for root in range(n):
                p,e,s=inv.f.a.bin_values(T,root,3)
                state=(p[1],e[1],s[2],p[2],e[2],s[3]);b=s[1]
                if (state,b) in seen:continue
                seen.add((state,b));out.append(dict(rows=T,root=root,state=state,defect=b))
            inv.f.a.tables.cache_clear()
    return out

def target(source,N):
    h,e,s,p,f,t=source['state'];b=source['defect'];G=s+p
    ns,_=inv.inverse(N,G,h,e,G-e-f-b)
    out=[]
    for n in ns:
        w=tuple((n>>i)&1 for i in range(n.bit_length()))
        h2,e2,s2,p2,F,t2=inv.f.advance(source['state'],w)
        inv.need(h2+F-b==N,'general six-state count')
        # Newborn y has indegree one, but its readiness defect is b.
        newborn=(N,h2-b,h2,N+b)
        # Exchange root and predecessor names: full H is unchanged,
        # and the selected physical root now carries the predecessor state.
        prepared=(N,h,h,2*F-h)
        delta=prepared[3]-N
        d=prepared[3]-h
        exits=[]
        for k in range(1,17):
            tt=(1<<k)*d+h-delta;bb=3*tt+2*delta
            if math.gcd(bb,tt)==1 and tt%3==0 and bb%3:
                exits.append(dict(k=k,t=tt,H=bb))
        out.append(dict(target=N,source=source,tail=n,word=''.join(map(str,w)),newborn_state=newborn,prepared_state=prepared,return_defect=delta,controller_exits=exits))
    return out

def physical(row):
    source=row['source'];w=tuple(map(int,row['word']))
    V,x=inv.f.a.funnel(tuple(source['rows']),(source['root'],),w)
    T,y=inv.f.a.birth(V,x)
    # Rename x<->y, keeping the selected physical name y.
    order=list(range(len(T)));order[x],order[y]=order[y],order[x]
    U=inv.f.a.induced(T,order)
    return T,y,U,y,x

def controls():
    checks=0
    for T in labelled(4):
        if not inv.f.a.strong(T):continue
        for root in range(4):
            p,e,s=inv.f.a.bin_values(T,root,3)
            state=(p[1],e[1],s[2],p[2],e[2],s[3]);b=s[1]
            for w in ((),(1,),(0,1),(1,1)):
                V,x=inv.f.a.funnel(T,(root,),w);U,y=inv.f.a.birth(V,x)
                new=inv.f.advance(state,w)
                inv.need(inv.direct(U)==new[0]+new[4]-b,'arbitrary-core full count')
                inv.need(inv.direct(V,x)==b,'inherited start multiplicity')
                inv.need(inv.direct(V)-inv.direct(U,y)==b,'literal growing defect')
                checks+=1
            inv.f.a.tables.cache_clear()
    return checks

def main():
    checks=controls();sources=bank();rows=[]
    for k in range(1,31):
        N=405*(1<<k)+71;hits=[]
        for source in sources:
            found=target(source,N)
            hits.extend(found)
        hits.sort(key=lambda r:(not bool(r['controller_exits']),r['source']['defect'],len(r['source']['rows'])+len(r['word']),r['tail']))
        selected=hits[0] if hits else None
        if selected:
            T,y,U,r,x=physical(selected)
            selected.update(rows=T,root=y,prepared_rows=U,prepared_root=r,predecessor=x)
            inv.need(inv.f.a.strong(T) and inv.f.a.strong(U),'strong actual sources')
            if len(T)<=18:
                inv.need(inv.direct(T)==N and inv.direct(U)==N,'literal neutral count')
                inv.need(inv.direct(T,y)==selected['newborn_state'][1],'literal newborn start')
                inv.need(inv.direct(U,r)==selected['prepared_state'][1],'literal prepared start')
                selected['literal_checked']=True
            else:selected['literal_checked']=False
        rows.append(dict(k=k,target=N,number_of_hits=len(hits),selected=selected))
        print('k',k,'hits',len(hits),'defect',None if selected is None else selected['source']['defect'],'controller',bool(selected and selected['controller_exits']),flush=True)
    receipt=dict(status='PASS',literal_controls=checks,source_bank=len(sources),targets=rows)
    suffix='' if __debug__ else '.optimized'
    (HERE/('growing_receipt'+suffix+'.json')).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

"""Exact wide-return controls, one-apex1817 inverse, and ready bridge search."""
from pathlib import Path
from itertools import product, combinations
import json, sys
from check import need, state, graph, H, birth, induced, root_profile, strong, reverse_complement_value
from audit_reset import coefficient_polynomial

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'reset'))
import funnel_six_state as funnel


def converse(T):
    return tuple(sum((T[u]>>v&1)<<u for u in range(len(T))) for v in range(len(T)))


def mirror_birth(T,r):
    return converse(birth(converse(T),r))


def rewritten(r,t,z):
    w='O'*r+'I'*t+z;T=graph(w);x=len(T);y=x+1
    G=birth(birth(T,0),x);R=list(G)
    for v in range(r-1):
        R[y]^=1<<v;R[v]|=1<<y
    order=[0,x,*range(1,r-1),y,r-1,*range(r,len(w)),len(w)]
    return w,G,tuple(R),order


def embedding(A,B,forced=None):
    """Exact directed induced embedding, MRV with forward-checking."""
    na,nb=len(A),len(B)
    degrees=[row.bit_count() for row in A]
    domains=[{v for v in range(nb) if B[v].bit_count()>=degrees[u]
              and nb-1-B[v].bit_count()>=na-1-degrees[u]} for u in range(na)]
    forced=forced or {}
    for u,v in forced.items():domains[u]&={v}
    def rec(assigned,dom):
        if len(assigned)==na:return assigned
        u=min((u for u in range(na) if u not in assigned),key=lambda u:len(dom[u]))
        for v in sorted(dom[u]):
            nxt=[set(d) for d in dom];ok=True
            for other in range(na):
                if other in assigned or other==u:continue
                orientation=A[u]>>other&1
                nxt[other]={w for w in nxt[other] if w!=v and (B[v]>>w&1)==orientation}
                if not nxt[other]:ok=False;break
            if ok:
                found=rec({**assigned,u:v},nxt)
                if found is not None:return found
        return None
    if any(not d for d in domains):return None
    result=rec({},domains)
    return None if result is None else [result[i] for i in range(na)]


def main():
    literals=[]
    for r in range(4,10):
        c=3*(1<<(r-3))-1
        for z in ('OO'+'I'*(r-3),'O'*(r-3)+'II'):
            need(state(z)[0]==c,'explicit eligible core at every width')
            for t in range(1,6):
                w='O'*r+'I'*t+z;v='OI'+'O'*(r-2)+'IO'+'I'*t+z
                h=( (1<<r)-1)*(reverse_complement_value(z)+1)*(1<<t)+1-5*(1<<(r-3))
                need(state(w)[0]==h and state(v)[0]==3*h+2,'all-wide scalar identity')
                if r<=6 and t==1:
                    _,G,R,order=rewritten(r,t,z)
                    need(induced(R,order)==graph(v),'literal wide returned order')
                    need(H(G)==H(R)==3*h+2,'independent wide path count')
                    if z.endswith('II'):
                        root=order[-2]
                        cp=root_profile(converse(R),root)
                        need(cp['F']-cp['H']==1,'converse ready root')
                        next_graph=mirror_birth(mirror_birth(R,root),len(R))
                        need(H(next_graph)==9*h+8,'literal mirrored target')
                    literals.append({'r':r,'z':z,'sourceH':h,'targetH':3*h+2,'order':len(R)})
    # Prefix completeness is a theorem; exhaustive small controls include rejected bits/slots.
    prefix_checks=0
    for n in range(2,12):
        for rest in product('OI',repeat=n-2):
            p='OO'+''.join(rest);P=state(p)[1];base=p[0]+'I'+p[1:]
            for position in range(2,len(base)+1):
                for bit in 'OI':
                    q=base[:position]+bit+base[position:]
                    if state(q)[1]==3*P:
                        trimmed=p.rstrip('I');s=len(p)-len(trimmed);r=len(trimmed)
                        need(set(trimmed)=={'O'} and bit=='I','forced run prefix')
                        need(q=='OI'+'O'*(r-2)+'IO'+'I'*s,'forced insertion cut')
                    prefix_checks+=1
    h1817=[];word_count=0
    for m in range(2,13):
        for middle in product('OI',repeat=m-2):
            w='O'+''.join(middle)+'I';word_count+=1
            if state(w)[0]==1817:h1817.append(w)
    need(not h1817,'complete1817 one-apex exclusion')
    banks={1:[funnel.SINGLE]};hits=[]
    for n in range(3,14):
        states=[]
        for m in [1]+list(range(3,n)):
            t=n-m-1
            words=[()] if not t else [tuple((j>>i)&1 for i in range(t-1))+(1,) for j in range(1<<(t-1))]
            for word in words:states.extend(funnel.advance(s,word) for s in banks[m])
        need(len(states)==funnel.a.class_counts(n)[n],'complete inherited ready grammar')
        hits.extend({'n':n,'rank':rank,'state':s} for rank,s in enumerate(states) if s[0]+s[4]-1==1817)
        banks[n]=states
    need([(v['n'],v['rank']) for v in hits]==[(10,344),(10,410)],'complete1817 ready-parent image')
    inflation_rows=[];targets=[]
    for rec in hits:
        P=funnel.a.unrank(rec['n'],rec['rank']);R=birth(P,0);targets.append(R)
        need(H(R)==1817 and root_profile(R,len(R)-1)['F']==1818,'literal ready1817 witness')
        rec['parent_rows']=P;rec['ready_rows']=R
        for module in range(1,len(P)):
            co=coefficient_polynomial(R,module)
            need(not(co[:2]==[1817,945] and not any(co[2:])),'single-module all-t return rejected')
            inflation_rows.append({'rank':rec['rank'],'module':module,'coefficients':co})
    slot_probes=[]
    for t in range(1,16):
        w='OOI'+'O'*t+'IOIIOI';base='OI'+w[1:];target=3*state(w)[0]+2;found=[]
        for pos in range(len(base)+1):
            for bit in 'OI':
                v=base[:pos]+bit+base[pos:]
                if state(v)[0]==target:found.append((pos,bit,v))
        need(not found,'finite mirrored all-slot hostile control')
        slot_probes.append({'t':t,'target':target,'equal_cut_replacements':found})
    # Search an actual higher-core bridge after the mirrored605 ->1817 exit.
    w,G,R,order=rewritten(4,1,'OII');old_apex=len(w)
    source=mirror_birth(mirror_birth(R,order[-2]),len(R))
    need(H(source)==1817 and len(source)==13,'actual mirrored source')
    bridge=None;exhausted=[]
    gauged_targets=[(i,gauge,T if gauge=='forward' else converse(T))
                    for gauge in ('forward','converse') for i,T in enumerate(targets)]
    for k in range(11,6,-1):
        trials=0
        for target_index,gauge,T in gauged_targets:
            for keep in combinations(range(len(T)),k):
                A=induced(T,keep);trials+=1
                mapping=embedding(A,source)
                if mapping is not None:
                    need(induced(T,keep)==induced(source,mapping),'literal retained ready1817 bridge')
                    bridge={'retained':k,'target_index':target_index,'gauge':gauge,'target_kept':keep,'source_kept':mapping,
                            'source_rows':source,'target_rows':T,'old_apex':old_apex,
                            'old_apex_retained':old_apex in mapping}
                    break
            if bridge:break
        if bridge:break
        exhausted.append({'retained':k,'target_subsets_exhausted':trials})
    need(bridge is not None,'some retained ready bridge')
    ready_graph=tuple(bridge['target_rows'])
    need(bridge['gauge']=='forward','selected concrete continuation gauge')
    following=birth(birth(ready_graph,len(ready_graph)-1),len(ready_graph))
    need(H(following)==5453,'literal continued ternary target')
    bridge['following_rows']=following
    bridge['following_H']=5453
    result={'status':'PASS','wide_literal_controls':literals,'prefix_insertion_controls':prefix_checks,
            'H1817_one_apex':{'max_length':12,'strong_words':word_count,'hits':h1817},
            'H1817_ready':hits,'single_module_controls':inflation_rows,
            'mirrored_cut_controls':slot_probes,'ready1817_bridge':bridge,'larger_retention_exhausted':exhausted}
    (HERE/'wide_results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','ready1817_bridge':bridge,'larger_retention_exhausted':exhausted},indent=2))


if __name__=='__main__':main()

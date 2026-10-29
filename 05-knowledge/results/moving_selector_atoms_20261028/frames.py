"""Exact run-order fibers, bases, and the phase/parity boundary."""
from pathlib import Path
from itertools import permutations,combinations,product
from math import factorial
from collections import Counter
import sys,json
import probes as p
r=p.r;core=p.core;need=p.need
HERE=Path(__file__).resolve().parent


def parity(v):return sum(v[i]>v[j] for i in range(len(v)) for j in range(i+1,len(v)))%2


def components(n,F):
    nxt=dict(F);incoming={b for a,b in F};runs=[]
    for a in range(n):
        if a in incoming:continue
        word=[a]
        while word[-1] in nxt:word.append(nxt[word[-1]])
        runs.append(tuple(word))
    need(sum(map(len,runs))==n,'directed linear forest')
    return tuple(runs)


def forest_fibers():
    atlas=json.loads((HERE.parent/'strong_atom_continuation_vitali_20261026/strong-atlas.json').read_text())
    reports=[]
    for entry in atlas:
        if entry['n']>5:continue
        A=tuple(entry['source']);n=len(A);G=r.oracle.isomorphisms(A,A,True);forests=set()
        for word in permutations(range(n)):
            for cuts in range(1<<max(0,n-1)):
                F=tuple(sorted((word[i],word[i+1]) for i in range(n-1) if not cuts>>i&1))
                if all(core.edge(A,a,b) for a,b in F):forests.add(F)
        counts=Counter();orders=0
        for F in forests:
            runs=components(n,F);k=len(runs);H=[g for g in G if tuple(sorted((g[a],g[b]) for a,b in F))==F]
            B=list(core.converse(A))
            for a,b in F:B[a]^=1<<b;B[b]^=1<<a
            need(all(all(core.edge(B,g[a],g[b])==core.edge(B,a,b) for a in range(n) for b in range(n)) for g in H),'joint relation group')
            index={word:i for i,word in enumerate(runs)}
            action=[tuple(index[tuple(g[v] for v in word)] for word in runs) for g in H]
            need(len(set(action))==len(H) and all(parity(g)==0 for g in action),'faithful even component action')
            seen=set()
            for order in permutations(range(k)):
                seen.add(min(tuple(g[i] for i in order) for g in action))
            need(len(seen)*len(H)==factorial(k),'exact component-order orbit count')
            if k>=2:need(Counter(parity(v) for v in seen)=={0:len(seen)//2,1:len(seen)//2},'two equal relative parity sheets')
            counts[k]+=1;orders+=len(seen)
        reports.append(dict(n=n,source=A,forests=len(forests),by_arity=dict(sorted(counts.items())),fiber_orbits_summed=orders))
    return reports


def minimum_base(G,n):
    for k in range(n+1):
        for b in combinations(range(n),k):
            if sum(all(g[v]==v for v in b) for g in G)==1:return b
    raise RuntimeError('nonfaithful vertex action')


def controls():
    A,_=core.substitute(p.cycle3(),(p.cycle3(),(0,),(0,)))
    W,_=core.substitute(p.cycle3(),(p.cycle3(),)*3)
    P=tuple(sum(1<<((v+d)%7) for d in (1,2,4)) for v in range(7))
    objects=[('five',A),('Paley7',P),('wreath9',W)]+[(name,p.controls.build_hostile(v)) for name,v in [('81',False),('27',True)]]
    out=[]
    for name,T in objects:
        G=r.oracle.isomorphisms(T,T,True);base=minimum_base(G,len(T))
        need(all(parity(g)==0 for g in G),'all tournament symmetries even')
        out.append(dict(object=name,order=len(T),group=len(G),base=base,all_singleton_order_orbits=factorial(len(T))//len(G)))
    need([len(v['base']) for v in out]==[1,2,3,3,3],'actual free-frame sizes')
    expected=p.controls.strong_controls()
    need([v['binary_attachment_orbits'] for v in expected[-2:]]==[176,512],'81/27 attachments')
    need([v['transporter_0_to_4'] for v in expected[-2:]]==[27,0],'81/27 transporters')
    for n in (6,8):
        shifts=[tuple((i+k)%n for i in range(n)) for k in range(n)]
        allowed=[k for k,g in enumerate(shifts) if all(c%2 for c in p.controls.old.permutation_cycles(g))]
        need(allowed==([0,2,4] if n==6 else [0]),'normal-basis orientation guard')
    incoming=HERE.parent/'session90_even_graphs_odd_descent_20261029/carriers';sys.path.insert(0,str(incoming))
    import check_transport as ct
    s=(0,2,1,6,8,7,3,5,4);t=(0,2,1,7,6,8,4,3,5);u=(3,5,4,0,2,1,6,8,7)
    h=ct.holonomy(u,t,s)
    need(h!=tuple(range(9)) and parity(h)==0 and ct.compose(h,s)==ct.compose(s,h),'nonzero phase invisible to parity')
    return dict(atoms=out,hostile=expected[-2:],holonomy=h,holonomy_parity=parity(h),F64_shifts=[0,2,4],F256_shifts=[0])


def parity_fiber():
    A,_=core.substitute(p.cycle3(),(p.cycle3(),(0,),(0,)));G=r.oracle.isomorphisms(A,A,True)
    source,_=core.substitute(p.cycle3(),(A,p.transitive(4),p.transitive(4)))
    strong_source,_=core.substitute(p.cycle3(),(A,(2,4,9,3),(2,4,9,3)))
    out=[]
    for middle in permutations(range(3)):
        word=(3,)+middle+(4,);path=[]
        for i,v in enumerate(word):
            path.append(v)
            if i<4:path.extend((5+i,9+i))
        target=core.complement(source,path)
        out.append(dict(middle=middle,parity=parity(middle),H=r.oracle.hamilton(target),
                        all_strong_H=r.oracle.hamilton(core.complement(strong_source,path)),
                        orbit=min(tuple(g[v] for v in word) for g in G)))
    need(len({tuple(v['orbit']) for v in out})==2,'two exact joint classes')
    need({v['H'] for v in out if v['parity']==0}=={2483993} and {v['H'] for v in out if v['parity']==1}=={2491877},'parity complete on this fiber')
    need({v['all_strong_H'] for v in out if v['parity']==0}=={2899675} and {v['all_strong_H'] for v in out if v['parity']==1}=={2906091},'all-strong parity fiber')
    return out


def paley27_phase():
    # F_3[a]/(a^3+2a+1); no root in F_3, hence irreducible.
    def digits(x):return (x%3,(x//3)%3,x//9)
    def add(x,y,sign=1):return sum(((a+sign*b)%3)*3**i for i,(a,b) in enumerate(zip(digits(x),digits(y))))
    def mul(x,y):
        z=[0]*5
        for i,a in enumerate(digits(x)):
            for j,b in enumerate(digits(y)):z[i+j]=(z[i+j]+a*b)%3
        for i in (4,3):
            z[i-2]=(z[i-2]+z[i])%3;z[i-3]=(z[i-3]+2*z[i])%3
        return sum(z[i]*3**i for i in range(3))
    def power(x,k):
        out=1
        for _ in range(k):out=mul(out,x)
        return out
    need(all((x**3+2*x+1)%3 for x in range(3)),'irreducible cubic')
    need(all(power(x,26)==1 for x in range(1,27)),'nonzero field units')
    T=tuple(sum(1<<v for v in range(27) if u!=v and power(add(v,u,-1),13)==1) for u in range(27))
    f=tuple(power(x,3) for x in range(27));f2=tuple(f[f[x]] for x in range(27));G=[tuple(range(27)),f,f2]
    need(all(f[f[f[x]]]==x for x in range(27)) and f!=G[0],'Frobenius order three')
    need(all(r.is_map(T,T,g) for g in G) and core.is_strong(T),'actual strong Paley27 action')
    U=(0,)+tuple(range(2,27))+(1,);V=tuple(5 if v==3 else 3 if v==5 else v for v in U)
    need(all(any(g[a]==b and g[0]==0 and g[1]==1 for g in G) for a,b in zip(U,V)),'Frobenius agrees on every single-port query')
    transposition=tuple(5 if v==3 else 3 if v==5 else v for v in range(27))
    need(not r.is_map(T,T,transposition),'joint ordering not isomorphic')
    need(sum(g[3]==3 for g in G)==1,'one non-prime-field port fixes the supplied Frobenius phase')
    return dict(field_polynomial='x^3+2x+1 over F3',order=27,rows=T,Frobenius=f,
                fixed_ports=[x for x in range(27) if f[x]==x],alpha_orbit=[3,f[3],f2[3]],
                swapped_ports=[3,5],single_port_profiles_agree=True,full_joint_orders_differ=True,
                scope='Displayed C3 subgroup, not a full automorphism-group or modular-primality classification.')


def main():
    result=dict(status='FINITE-EXACT',forests=forest_fibers(),frames=controls(),parity_fiber=parity_fiber(),paley27=paley27_phase())
    (HERE/'frame-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(forests=sum(v['forests'] for v in result['forests']),atoms=result['frames']['atoms'],holonomy=result['frames']['holonomy'])));return result


if __name__=='__main__':main()

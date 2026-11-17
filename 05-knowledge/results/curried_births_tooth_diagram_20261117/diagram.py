from itertools import permutations, product
from functools import lru_cache
from collections import Counter, defaultdict

@lru_cache(None)
def canon(T):
    n=len(T)
    cells=[tuple(i for i in range(n) if T[i].bit_count()==d) for d in range(n)]
    return min(tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        for choice in product(*(permutations(c) for c in cells if c))
        for p in [sum(choice,())])

def delete(T,i):
    p=[v for v in range(len(T)) if v!=i]
    return canon(tuple(sum((T[a]>>b&1)<<k for k,b in enumerate(p)) for a in p))

def opp(T):
    n=len(T)
    return canon(tuple(v^((1<<n)-1)^(1<<i) for i,v in enumerate(T)))

def later(L):
    n=len(next(iter(L))); F=(1<<n)-1
    return {canon(tuple(v|((mask>>i&1)<<n) for i,v in enumerate(T))+(F^mask,))
            for T in L for mask in range(1<<n)}

def ham(T):
    return sum(all(T[a]>>b&1 for a,b in zip(w,w[1:])) for w in permutations(range(len(T))))

def diagram(limit=7):
    L={()}
    for n in range(1,limit+1):
        L=later(L)
        S={T for T in L if T==opp(T)}
        bad=[(T,sorted(set(delete(T,i) for i in range(n)))) for T in sorted(S)
             if all((P:=delete(T,i))!=opp(P) for i in range(n))]
        print(n,len(L),len(S),len(bad),bad[:1],flush=True)
if __name__=='__main__':
    diagram()

@lru_cache(None)
def aut(T):
    n=len(T)
    cells=[tuple(i for i in range(n) if T[i].bit_count()==d) for d in range(n)]
    return [p for C in product(*(permutations(cell) for cell in cells if cell))
        for p in [sum(C,())] if all(T[a]>>b&1==T[p[a]]>>p[b]&1 for a in range(n) for b in range(n))]

def reversers(T):
    return [p for p in permutations(range(len(T)))
            if all(T[a]>>b&1 != T[p[a]]>>p[b]&1 for a in range(len(T)) for b in range(len(T)) if a!=b)]

def svg():
    import html
    layers=[{()}]
    for n in range(6):
        layers.append(later(layers[-1]))
    view=[None]*len(layers)
    P={}
    for n,L in enumerate(layers[1:],1):
        fixed=sorted(T for T in L if T==opp(T))
        chiral=sorted(T for T in L if T<opp(T))
        yy=85+200*(n-1)
        for j,T in enumerate(fixed):
            P[T]=(590,yy+11*(j-(len(fixed)-1)/2))
        for j,T in enumerate(chiral):
            xx=42+22*(j%22)
            Y=yy+14*(j//22)
            P[T]=(590-xx,Y); P[opp(T)]=(590+xx,Y)
        E=[]
        contact_edges=[]
        for T in L:
            rev=reversers(T) if T in fixed else []
            unseen=set(range(n))
            while unseen:
                v=min(unseen)
                O={p[v] for p in aut(T)}
                unseen-=O
                flag=any(p[v]==v for p in rev)
                E.append((T,delete(T,v),v,flag))
        flags=sum(x[-1] for x in E)
        if flags!=(len(fixed) if n%2 else 0):
            raise ValueError(n)
        if n>1:
            # independently count the same incidence set from parent neighborhood orbits
            for U in layers[n-1]:
                masks=set(range(1<<(n-1)))
                while masks:
                    S=min(masks)
                    masks-={sum((S>>i&1)<<p[i] for i in range(n-1)) for p in aut(U)}
                    contact_edges.append(U)
            if Counter(U for T,U,v,f in E)!=Counter(contact_edges):
                raise ValueError(('incidence mismatch',n))
        view[n]=(L,E,flags)
        print('tooth layer',n,'classes',len(L),'pointed incidences',len(E),'fixed incidences',flags)
    bad=(2,4,1,23,39,15)
    ps=sorted({delete(bad,v) for v in range(6)})
    if bad!=opp(bad) or len(ps)!=2 or opp(ps[0])!=ps[1]:raise ValueError('obstruction')
    if len(aut(bad))!=9:raise ValueError('automorphism check')
    a=['<svg xmlns="http://www.w3.org/2000/svg" width="1180" height="1210" viewBox="0 0 1180 1210">',
       '<rect width="1180" height="1210" fill="#f8fafc"/>',
       '<text x="35" y="32" font-size="23" font-family="sans-serif" fill="#0f172a">Tournament tooth diagram: converse is reflection</text>',
       '<line x1="590" y1="45" x2="590" y2="1185" stroke="#9aa7b8" stroke-dasharray="3 5"/>']
    edge_layer={}
    for n in range(2,7):
        for T,U,v,flag in view[n][1]:
            x,y=P[T];u,z=P[U]
            hi=T==bad
            c='#d85b16' if hi else '#176ab1' if flag else '#8695a8'
            edge_layer.setdefault((c,1 if hi else .74 if flag else .24,3 if hi else 1.3 if flag else .85),[]).append(f'M{u},{z}C{u},{(y+z)/2} {x},{(y+z)/2} {x},{y}')
    for (c,o,w),ds in sorted(edge_layer.items(), key=lambda pair:pair[0][1]):
        a.append(f'<path d="{" ".join(ds)}" fill="none" stroke="{c}" stroke-opacity="{o}" stroke-width="{w}"/>')
    for n in range(1,7):
        L,E,f=view[n]; a.append(f'<text x="22" y="{151+200*(n-1)}" font-family="sans-serif" font-size="14" fill="#14213a">n={n} · {len(L)} classes</text>')
        for T in sorted(L):
            x,y=P[T]
            col='#d85b16' if T==bad or T in ps else '#176ab1' if T==opp(T) else '#344866'
            a.append(f'<line x1="{x-7}" y1="{y}" x2="{x+7}" y2="{y}" stroke="{col}" stroke-width="4" stroke-linecap="round"><title>{html.escape(str(T))}; H={ham(T)}; {"self-converse" if T==opp(T) else "chiral"}</title></line>')
    a+=['<text x="26" y="1182" font-family="sans-serif" font-size="14" fill="#14213a">Blue: reflection-fixed tooth or incidence. Orange: a fixed tooth whose two parent classes are swapped.</text>','</svg>']
    from pathlib import Path
    Path(__file__).with_name('tooth_diagram.svg').write_text('\n'.join(a)+'\n')

if __name__=='__main__':
    svg()

def path_profile(T,r):
    n=len(T); full=(1<<n)-1
    H={0:1}
    for S in range(1,full+1):
        U=tuple(i for i in range(n) if S>>i&1)
        H[S]=sum(all(T[x]>>y&1 for x,y in zip(w,w[1:])) for w in permutations(U))
    E=O=0
    for S in range(1,full+1):
        if not S>>r&1:continue
        U=tuple(i for i in range(n) if S>>i&1 and i!=r)
        d=sum(all(T[x]>>y&1 for x,y in zip((r,)+w,w)) for w in permutations(U))*H[full^S]
        if S.bit_count()%2:O+=d
        else:E+=d
    return H[full],H[full^(1<<r)],E,O

def reads(T,r):
    def extend(A,tag):
        n=len(A)
        return tuple(row|((1<<n) if tag=='S' or i!=r else 0) for i,row in enumerate(A))+(0 if tag=='S' else 1<<r,)
    return tuple(ham(finish(T,w,extend)) for w in ('','N','SN','NN'))

def finish(T,w,extend):
    for tag in w:T=extend(T,tag)
    return T

if __name__=='__main__':
    A=(2,4,9,19,7);B=(2,8,3,21,7)
    lhs=(path_profile(A,0),reads(A,0));rhs=(path_profile(B,1),reads(B,1))
    if lhs!=((9,5,4,7),(9,17,31,29)) or rhs!=((9,1,4,5),(9,17,27,33)):raise ValueError((lhs,rhs))
    print('Birth H9 profile/readouts',lhs,rhs)

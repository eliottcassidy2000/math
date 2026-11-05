import importlib.util,itertools,json,math,pathlib,random
ROOT=pathlib.Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
pa=load('prime_audit',ROOT.parent/'session197-prime-atoms-audit'/'audit.py')
fa=load('factor_audit',ROOT.parent/'session197-prime-factor-audit'/'audit.py')
ca=load('census_audit',ROOT.parent/'session197-lowcycle-audit'/'check.py')

def compress(rows,r):
    vs=[v for v in range(len(rows)) if v!=r];S=pa.induced(rows,vs);cs=pa.components(S)
    O={v for v in range(len(S)) if rows[r]&(1<<vs[v])};assert O
    first=next(i for i,C in enumerate(cs) if O.intersection(C))
    keep=[v for i,C in enumerate(cs) if i>=first or len(C)>1 for v in C]
    Q=pa.induced(rows,[vs[v] for v in keep]+[r]);suffix=sum(len(C) for C in cs[first+1:])
    return Q,len(Q)-1,suffix

def check(rows,r,z=None):
    h,e,R,d=z or pa.data(rows)[0][r];U=R-d;f=R+e-h-d
    assert U>0 and h-e>=d and U>=d+f
    Q,q,s=compress(rows,r)
    assert pa.data(Q)[0][q]==(h,e,R,d)
    assert U>=1<<s and len(Q)<=d+U.bit_length()+1
    return Q,q

def main():
    raw=labelled=zero=0;literal_injections=0
    for n in range(1,7):
        for line in (ROOT.parent/'session197-lowcycle-audit'/f'classes{n}.txt').read_text().splitlines():
            key,H,C,st=map(int,line.split());rows=tuple((key>>(8*i))&255 for i in range(n))
            weight=math.factorial(n)//ca.autos(rows)
            for r,z in enumerate(pa.data(rows)[0]):
                if z[2]==z[3]:zero+=weight;continue
                raw+=1;labelled+=weight;check(rows,r,z)
                if n<=5:
                    vs=[v for v in range(n) if v!=r];out=set();deletion=0
                    for p in itertools.permutations(vs):
                        if not all(rows[u]&(1<<v) for u,v in zip(p,p[1:])):continue
                        i=next(i for i,v in enumerate(p) if rows[r]&(1<<v))
                        q=p[:i]+(r,)+p[i:]
                        assert q[-1]!=r and all(rows[u]&(1<<v) for u,v in zip(q,q[1:]))
                        out.add(q);deletion+=1;literal_injections+=1
                    assert len(out)==deletion==z[3]
    # Complete labelled coverage is reconstructed by exact orbit mass, with
    # every root inspected in every independent isomorphism representative.
    source=json.loads((ROOT/'source'/'receipt-sharpening.json').read_text())
    assert labelled==source['insertion_inequality_root_checks_through_order6']
    rng=random.Random(1970291);larger=0
    for n in (7,8,9,10):
        for _ in range(20):
            rows=[0]*n
            for u,v in itertools.combinations(range(n),2):
                a,b=(u,v) if rng.randrange(2) else (v,u);rows[a]|=1<<b
            for r,z in enumerate(pa.data(tuple(rows))[0]):
                if z[2]>z[3]:check(tuple(rows),r,z);larger+=1
    blocks=[(4,1,2),(0,),(0,),(0,),(4,1,2),(0,),(0,),(0,)]
    n=sum(map(len,blocks));S=[];offset=0
    for C in blocks:
        after=sum(1<<v for v in range(offset+len(C),n))
        S.extend((row<<offset)|after for row in C);offset+=len(C)
    adversarial=0
    for O in ((6,),(6,7),(0,7),(3,10)):
        mask=sum(1<<v for v in O)
        T=tuple(row|((1<<n) if v not in O else 0) for v,row in enumerate(S))+(mask,)
        for _ in range(3):
            order=list(range(n+1));rng.shuffle(order)
            check(pa.induced(T,order),order.index(n));adversarial+=1
    sources={};partners={}
    for n in range(1,5):
        for rows in pa.source.labelled(n):
            rows=tuple(rows);pp=pa.data(rows)[0];qq=pa.data(fa.converse(rows))[0]
            for r,z in enumerate(pp):
                if z[2]>z[3]:sources.setdefault(z,(rows,r))
                g,s,L,t=qq[r];partners.setdefault((g,t,L),(rows,r))
    pair_count=actual=0
    for z,(X,r) in sources.items():
        h,e,R,d=z;U=R-d;f=R+e-h-d
        for (g,t,L),(Y,s) in partners.items():
            beta=d*g-f*t
            if beta not in (-1,1):continue
            assert 2*d<=U+beta
            NX,nr=check(X,r,z);assert len(NX)<=(U+beta)//2+U.bit_length()+1
            if L==1:NY=(0,);ns=0
            else:
                assert L-t>0;W,ns=check(fa.converse(Y),s);NY=fa.converse(W)
            gg,ss,LL,tt=pa.data(fa.converse(NY))[0][ns]
            assert (gg,LL,tt)==(g,L,t) and len(NY)<=(1 if L==1 else L+1)
            for m in (1,2,3):
                T,clone_root=fa.star(NX,nr,NY,ns,m);assert pa.data(T)[1][-1]==U*L*(1<<m)+beta;actual+=1
            pair_count+=1
    tables=json.loads((ROOT.parent/'session197-prime-addendum-audit'/'AUDIT.json').read_text())['tables']
    bases={(1,1):[1],(1,-1):[3,5],(3,1):[6],(3,-1):[9,15],(5,1):[11,12,16],(5,-1):[15,39],(7,1):[],(7,-1):[]}
    bases.update({tuple(map(int,k.split(':'))):v['bases'] for k,v in tables.items()})
    allowed=[];unknown=[]
    for p in (3,5,7,11,13,17,19,23,29):
        for beta in (-1,1):
            ds=list(range(1,(p+beta)//2+1,2));hits=[d for d in ds if p in bases.get((d,beta),[])]
            missing=[d for d in ds if (d,beta) not in bases]
            assert hits or not missing
            if hits:allowed.append((p,beta))
            if missing:unknown.append({'p':p,'beta':beta,'d':missing,'already_witnessed':hits})
    assert allowed==[(3,-1),(5,-1),(11,1),(29,1)]
    X=(2,44,49,53,3,17);assert pa.data(X)[0][5]==(23,4,40,11)
    values=[]
    for m in range(1,6):
        T=fa.tail(X,5,m);h=pa.data(T)[1][-1];assert h==29*(1<<m)+1;values.append(h)
    result={'status':'PASS','independent_representative_positive_roots_through6':raw,
            'represented_labelled_positive_roots_through6':labelled,'excluded_zero_U_roots':zero,
            'literal_deletion_path_insertions':literal_injections,'larger_random_compressions':larger,
            'adversarial_multicomponent_shuffled_label_compressions':adversarial,
            'unit_profile_pairs':pair_count,'actual_star_checks':actual,
            'prime_amplitudes_through29':allowed,'unneeded_unknown_rows':unknown,'p29_actual_K_counts':values}
    (ROOT/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()

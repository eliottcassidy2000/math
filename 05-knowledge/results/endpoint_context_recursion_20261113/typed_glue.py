"""Seven-coordinate compiler and exhaustive context search for inherited pairs."""
from pathlib import Path
import gzip,json,sys
from fractions import Fraction
from itertools import product
from collections import Counter
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'surplus_selector_20261113'))
from context_transport import glue,decode
from repair23 import profile,hp2,literal_profile
from inflation_probe import relabel
from context_head import ordinary_bank,need,save


def reverse(T):
    full=(1<<len(T))-1
    return tuple(full^(1<<v)^T[v] for v in range(len(T)))


def seven(T,r):
    h,s,q,F=profile(T,r)
    _,e,_,G=profile(reverse(T),r)
    P=hp2(relabel(T,[v for v in range(len(T)) if v!=r]))[1]
    return h,s,e,q,F,G,P


def compile_glue(a,b):
    """Both rooted operands have at least two vertices. Singleton is identity."""
    h,s,e,q,F,G,P=a
    H,S,E,Q,BF,BG,BP=b
    delta=G-h-q+e
    return (q*H+(e-q)*S-delta*Q+(G-q)*BF,
            s*Q,q*E,q*Q,P*BF+(F-P)*Q,
            q*BG+(G-q)*BP,P*BP)


def diff(a,b):return tuple(y-x for x,y in zip(a,b))


def endpoint_type(d):
    h,s,e,q,F,G,P=d
    return [mode for mode,ok in [('E',F==s),('S',G==e)] if h==0 and q!=0 and ok]


def run():
    checks=literal=0
    # Independent physical construction and permutation controls, all labelled roots.
    for n in range(2,5):
        for code in range(1<<(n*(n-1)//2)):
            A=decode(n,code)
            for r in range(n):
                a=seven(A,r)
                for m in (2,3):
                    for bcode in range(1<<(m*(m-1)//2)):
                        B=decode(m,bcode)
                        for t in range(m):
                            b=seven(B,t);T=glue(A,r,B,t);z=seven(T,r)
                            need(z==compile_glue(a,b),'seven-coordinate physical law')
                            if n<=3 and m==2:
                                h,s,e,q,F,G,P=z
                                need(literal_profile(T,r)==(h,s,q,F),'permutation control')
                                literal+=1
                            checks+=1
    # Explicit physical associativity; the vertex labels also agree in this order.
    assoc=0
    small=[(decode(3,c),r) for c in range(8) for r in range(3)]
    for (A,r),(B,s),(C,t) in product(small[::4],repeat=3):
        AB=glue(A,r,B,s);BC=glue(B,s,C,t)
        X=glue(AB,r,C,t);Y=glue(A,r,BC,s)
        # Labels differ on the right if B's root is not last: compare exact profiles.
        need(seven(X,r)==seven(Y,r)==compile_glue(compile_glue(seven(A,r),seven(B,s)),seven(C,t)),'associative response')
        assoc+=1
    catalog={}
    for (n,code),A in sorted(ordinary_bank().items()):
        if n==1:continue
        for r in range(n):
            z=seven(A,r)
            catalog.setdefault(z,dict(n=n,rows=A,root=r,profile=z))
    print('CATALOG',len(catalog),'PHYSICAL',checks,'PERMUTATIONS',literal,'ASSOCIATIVE',assoc,flush=True)
    raw=json.loads(gzip.decompress((HERE.parent/'surplus_selector_20261113/transport_words.json.gz').read_bytes()))
    pairs=[]
    for row in raw['outputs']:
        pairs.append(dict(label=row['input_slope']+'>'+row['output_slope'],rows=row['output_rows'],root=row['root']))
        if 'second' in row:
            p=row['second'];pairs.append(dict(label=row['input_slope']+'>'+row['output_slope']+'>'+p['output_slope'],rows=p['rows'],root=p['root']))
    outputs=[]
    for pair in pairs:
        x,y=[seven(T,pair['root']) for T in pair['rows']]
        d=diff(x,y);need('E' in endpoint_type(d),'inherited actual E pair')
        hits={};counts=Counter()
        for z,A in catalog.items():
            for side in ('left','right'):
                a,b=(compile_glue(z,x),compile_glue(z,y)) if side=='left' else (compile_glue(x,z),compile_glue(y,z))
                D=diff(a,b)
                for mode in endpoint_type(D):
                    rho=Fraction(D[1 if mode=='E' else 2],D[3])
                    counts[side,mode]+=1
                    key=side,mode,str(rho)
                    rec=dict(side=side,mode=mode,slope=str(rho),context=A,output=[a,b],difference=D,total_orders=[len(T)+A['n']-1 for T in pair['rows']])
                    rank=(A['n'],a[0],z)
                    if key not in hits or rank<hits[key][0]:hits[key]=(rank,rec)
        rec=dict(**pair,profiles=[x,y],difference=d,hits=[v[1] for k,v in sorted(hits.items())],counts=[dict(side=s,mode=m,count=c) for (s,m),c in sorted(counts.items())])
        outputs.append(rec)
        print('PAIR',pair['label'],'ORDER',len(pair['rows'][0]),'DIFFERENCE',d,'TYPED_RULES',[(v[1]['side'],v[1]['mode'],v[1]['slope'],v[1]['context']['n']) for v in hits.values()],flush=True)
    save('typed_glue',dict(status='FINITE-EXACT',scope='All ordinary contexts through8 and all roots; grouped only by full seven-coordinate equality. Ten inherited literal pairs. No exhaustive search over arbitrary input pairs.',physical_checks=checks,permutation_controls=literal,associativity_controls=assoc,context_response_types=len(catalog),pairs=outputs))
    save('context_catalog',dict(status='FINITE-EXACT',records=list(catalog.values())))


if __name__=='__main__':run()

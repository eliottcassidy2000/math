"""Independent actual-piece tests and boundary-cover slack-six table."""
import importlib.util,itertools,json,math,pathlib,random
ROOT=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('prime_audit',ROOT.parent/'session197-prime-atoms-audit'/'audit.py')
pa=importlib.util.module_from_spec(spec);spec.loader.exec_module(pa)

def converse(rows):return tuple(sum(1<<v for v in range(len(rows)) if rows[v]&(1<<u)) for u in range(len(rows)))
def tail(rows,r,j):
    for _ in range(j):
        n=len(rows);rows=tuple(row|((1<<n) if i!=r else 0) for i,row in enumerate(rows))+(1<<r,)
    return rows
def normalize(rows,r):
    C,c,j,u=pa.normalize(rows,r);return tail(C,c,j),c,j
def star(X,r,Y,s,m):
    nx=len(X);B=[v for v in range(len(Y)) if v!=s]
    mapping={s:r,**{v:nx+i for i,v in enumerate(B)}}
    nz=nx+len(B);n=nz+m;rows=list(X)+[0]*(n-nx)
    for u in range(len(Y)):
        for v in range(len(Y)):
            if Y[u]&(1<<v):rows[mapping[u]]|=1<<mapping[v]
    for a in range(nx):
        if a!=r:
            for v in range(nx,n):rows[a]|=1<<v
    for z in range(nz,n):
        rows[z]|=1<<r
        for b in range(nx,nz):rows[z]|=1<<b
        for zz in range(z+1,n):rows[z]|=1<<zz
    return tuple(rows),nz
def k(f):return (1+math.isqrt(1+8*f))//2
def ell(f):return max(1,f).bit_length()-1
def valuation(x):return (x&-x).bit_length()-1

def check_pair(X,r,Y,s,physical=True):
    h,e,R,d=pa.data(X)[0][r];U=R-d;f=R+e-h-d
    g,start,L,t=pa.data(converse(Y))[0][s];beta=d*g-f*t
    assert U>0 and beta in (-1,1)
    assert d+2*f<=2*U-1 and f<=U-1 and 3*d<=2*U+2*beta-1
    if beta==-1:assert f>d
    if beta==1 and f<d:assert (g,t,d-f)==(1,1,1) and L&(L-1)==0
    NX,nr,j=normalize(X,r)
    assert j<=valuation(U) and pa.data(NX)[0][nr]==(h,e,R,d)
    BX=1+(2*U+2*beta-1)//3+k(U-1)+ell(U-1)+valuation(U)
    assert len(NX)<=BX
    if L==1:NY=(0,);ns=0
    else:
        W,ns,jy=normalize(converse(Y),s);NY=converse(W)
        assert L-t>0 and t<=L-1 and jy<=valuation(L-t)<=ell(L-1)
        F=L+start-g-t;assert 0<=F<=L-1
    gg,ss,LL,tt=pa.data(converse(NY))[0][ns]
    assert (gg,LL,tt)==(g,L,t)
    BY=1 if L==1 else L+k(L-1)+2*ell(L-1)
    assert len(NY)<=BY
    if physical:
        for m in (1,2,3):
            T,z=star(NX,nr,NY,ns,m)
            assert len(T)==len(NX)+len(NY)+m-1
            assert pa.data(T)[1][-1]==U*L*(1<<m)+beta
    return U,L,beta,len(NX),len(NY)

def main():
    sources={};partners={};sm={};pm={}
    for n in range(1,5):
        for rows in pa.source.labelled(n):
            rows=tuple(rows);pp=pa.data(rows)[0];qq=pa.data(converse(rows))[0]
            for r,z in enumerate(pp):
                if z[2]>z[3]:sources.setdefault(z,(rows,r));sm[z]=sm.get(z,0)+1
                g,s,L,t=qq[r];key=(g,t,L);partners.setdefault(key,(rows,r));pm[key]=pm.get(key,0)+1
    pairs=[];labelled=0
    for z,(X,r) in sources.items():
        h,e,R,d=z;f=R+e-h-d
        for key,(Y,s) in partners.items():
            g,t,L=key
            if d*g-f*t not in (-1,1):continue
            pairs.append(check_pair(X,r,Y,s));labelled+=sm[z]*pm[key]
    table=pa.table((4,1,2),3,-3)
    # This independent enumeration omits the one zero-outdegree source,
    # which has U=f=0 and cannot enter the positive slack-six table.
    assert table['bases']==[9,15,27,36,60] and table['examined']==753
    for U,w in table['witnesses'].items():
        T=w['rows'];r=w['root'];z=w['profile'];assert pa.literal(T,r)==z
    # Independently inspect every strong H15 class in the complete order<=8
    # census, including the non-nearly-transitive order-six examples.
    H15=[];classes=ROOT.parent/'session197-lowcycle-audit'
    for n in range(1,9):
        for line in (classes/f'classes{n}.txt').read_text().splitlines():
            key,H,C,st=map(int,line.split())
            if H!=15 or not st:continue
            rows=tuple((key>>(8*i))&255 for i in range(n));roots=[]
            for r in range(n):
                # A positive binary root has nonempty O-to-I cross graph with
                # all cross edges sharing a vertex (matching number one).
                O=[v for v in range(n) if rows[r]&(1<<v)]
                I=[v for v in range(n) if rows[v]&(1<<r)]
                edges=[{u,v} for u in O for v in I if rows[u]&(1<<v)]
                binary=bool(edges) and bool(set.intersection(*edges))
                if binary:
                    # Duplicate r with a fresh larger-labelled transitive clone.
                    T=list(rows)+[rows[r]]
                    for v in range(n):
                        if rows[v]&(1<<r):T[v]|=1<<n
                    T[r]|=1<<n
                    h2=pa.data(tuple(T))[1][-1];beta=30-h2
                    assert beta not in (-1,1);roots.append({'root':r,'beta':beta})
            H15.append({'order':n,'rows':rows,'binary_nonunit_roots':roots})
    assert len(H15)==4 and {x['order'] for x in H15}=={5,6}
    report={'status':'PASS','independent_unit_profile_pairs':len(pairs),
            'represented_labelled_unit_pairs':labelled,'actual_star_checks':3*len(pairs),
            'normalized_piece_data':pairs,'independent_d3_f6_table':table,
            'complete_small_strong_H15_controls':H15}
    (ROOT/'AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()

"""Replay literal returns, all-size clone certificates, and transported phases."""
from search import *
from return_program import field,golden
from switch_search import hpaths
from singer import mul,power
from math import comb
import hashlib


def clone(T,q,m):
    addresses=[(v,j) for v in range(len(T)) for j in range(m if v==q else 1)]
    rows=tuple(sum(1<<b for b,(w,k) in enumerate(addresses)
               if (v==w and j<k) or (v!=w and T[v]>>w&1)) for v,j in addresses)
    return rows,addresses


def lift_path(p,addresses):
    return tuple(i for v in p for i,a in enumerate(addresses) if a[0]==v)


def legal(T,p):
    need(sorted(p)==list(range(len(T))) and all(T[a]>>b&1 for a,b in zip(p,p[1:])),'legal named Hamilton path')


def raw(T,root):
    d=profile(T,root);h,x,y=d['observer']
    return [h,y-x,x+y]


def onto(m,r):return sum((-1)**(r-i)*comb(r,i)*i**m for i in range(r+1))
def evaluate(c,m):return sum(x*onto(m,j+1) for j,x in enumerate(c))
def coefficients(a):
    c1=a[0];need((a[1]-c1)%2==0,'two-run integrality');c2=(a[1]-c1)//2
    need((a[2]-c1-6*c2)%6==0,'three-run integrality');c3=(a[2]-c1-6*c2)//6
    need(min(c1,c2,c3)>=0,'actual nonnegative run coefficients')
    return [c1,c2,c3]


def observer(coeff,m):
    h,e,s=[evaluate(c,m) for c in coeff];need((s-e)%2==0,'golden integer guard')
    return h,(s-e)//2,(s+e)//2


def field_frames(v):
    a,b,g=field(v)
    sa,sb=(a+4*b)%7,-b%7
    need(((sb+2*sa)%7,2*sa%7)==(2*a%7,(b+2*a)%7),'normal basis Frobenius swap')
    return dict(phase=(a,b,g),normal_basis=((b+2*a)%7,2*a%7),
                Frobenius=((a+4*b)%7,-b%7,g),
                Q4={q:field(golden.step(v,golden.BANK[q])) for q in ('q0','q1','q2','q4','q5')})


def ancestry(h):
    kind='straight' if h in (1457,1459) else 'revisit'
    # The established ancestral family is at n=3; Singer input was at n=4.
    z4=(4,6,1) if kind=='straight' else (6,4,1)
    z3=mul(power((4,4,1),170),z4)
    need(mul((4,4,1),z3)==z4,'ancestral F343 inverse transport')
    triple=(485,873,1649) if kind=='straight' else (487,889,1631)
    r0=(triple[1]-3*triple[0])%7;r1=(triple[2]-3*triple[1])%7
    return dict(family=kind,n=3,F49_spectral=((r1+4*r0)%7,2*r0%7),F343=z3,
                scope='frozen ancestral phase; clone growth is a different operation')


def run():
    rows=json.loads((HERE/'branches.json').read_text());reset=json.loads((HERE/'reset.json').read_text());literal=[]
    for record in rows:
        T=tuple(record['source']);h=record['H'];need(tables(T)[0][-1]==h,'source H')
        for count,w in record['witnesses'].items():
            p=tuple(w['path']);legal(T,p);U=native.native(T,p)
            W=birth(U,w['mask']);need(W==tuple(w['target']) and strong(W),'literal strong target')
            need(tables(W)[0][-1]==int(count)==sum(1 for _ in hpaths(W)),'independent unpruned path enumeration')
            need(profile(W)['observer']==tuple(w['profile']['observer']),'literal joint root replay')
            literal.append(dict(source=h,target=int(count),moves=1,ancestry=ancestry(h),
                                frames=field_frames(tuple(w['profile']['observer']))));tables.cache_clear()
    T=tuple(reset['source']);U=T
    for name in ('first_path','second_path'):
        p=tuple(reset[name]);legal(U,p);U=native.native(U,p)
    W=birth(U,reset['mask']);need(W==tuple(reset['target']) and strong(W),'two-complement strong reset')
    need(tables(W)[0][-1]==4397==sum(1 for _ in hpaths(W)),'reset independent enumeration')
    literal.append(dict(source=1465,target=4397,moves=2,ancestry=ancestry(1465),
                        frames=field_frames(profile(W)['observer'])));tables.cache_clear()
    need(len(literal)==8,'complete four-source ternary branching')
    census=json.loads((HERE/'census.json').read_text())
    need([(x['source_H'],x['paths'],x['all_hits'],x['strong_hits']) for x in census]==[
        (1457,1457,[10,10],[4,3]),(1459,1459,[12,14],[6,8]),
        (1463,1463,[10,24],[5,12]),(1465,1465,[4,10],[0,10])],'independent complete native/birth census')
    for x in census:
        need(x['paths']==sum(1 for _ in hpaths(bank()[x['source_H']])),'independent source path universe')
        for w in x['weak_targets']:
            W=tuple(w);need(not strong(W) and tables(W)[0][-1] in (3*x['source_H']+2,3*x['source_H']+4),'weak targets stay separate');tables.cache_clear()
    lifts=json.loads((HERE/'lift_certificates.json').read_text());reports=[]
    for d in lifts:
        h=d['source_H'];T=bank()[h];q=d['port'];W=tuple(d['target']);n=len(T)
        p=tuple(d['first_path']);s=tuple(d['second_path']);legal(T,p);U=native.native(T,p);legal(U,s)
        need(birth(native.native(U,s),d['mask'])==W and strong(W),'literal lift target and strength')
        def neighbours(path):
            i=path.index(q);return (path[i-1] if i else None,path[i+1] if i+1<len(path) else None)
        need(neighbours(p)==neighbours(s),'common ordered entry and exit frame')
        degree=T[q].bit_count();need(degree in (1,n-2),'at most two module runs')
        need(bool(d['mask']>>q&1)==(degree==n-2),'birth preserves exceptional in/out neighbour')
        sampled=[[],[]];frames=[]
        for m in range(1,5):
            C,ad=clone(T,q,m);D,wd=clone(W,q,m);pp=lift_path(p,ad);ss=lift_path(s,ad)
            legal(C,pp);V=native.native(C,pp);legal(V,ss);V=native.native(V,ss)
            mask=sum(1<<i for i,a in enumerate(ad) if d['mask']>>a[0]&1)
            need(birth(V,mask)==D,'literal common-frame cloning commutes with both complements and birth')
            inverse=D[:-1];inverse=tuple(row&((1<<len(C))-1) for row in inverse)
            legal(inverse,ss);inverse=native.native(inverse,ss);legal(inverse,pp)
            need(native.native(inverse,pp)==C,'realization-preserving inverse with the named newborn')
            for ix,(X,xd) in enumerate(((C,ad),(D,wd))):
                r=xd.index((2,0));sampled[ix].append(raw(X,r));tables.cache_clear()
            need(sampled[1][-1][0]==3*sampled[0][-1][0]+d['offset'],'literal offset and clone control')
        coeff=[]
        for side in sampled:
            cs=[coefficients([side[m][j] for m in range(3)]) for j in range(3)]
            need(all(evaluate(cs[j],4)==side[3][j] for j in range(3)),'independent fourth clone control')
            need(cs[0][2]==0,'two-run scalar bound')
            coeff.append(cs)
        need(coeff[1][0][0]==3*coeff[0][0][0]+d['offset'] and coeff[1][0][1]==3*coeff[0][0][1],'all-size offset coefficient identity')
        periods=[]
        for cs in coeff:
            phases=[field(observer(cs,m)) for m in range(1,13)]
            periods.append(next(k for k in (1,2,3,6) if all(phases[i]==phases[i+k] for i in range(12-k))))
            # Exact integer observer periods, checked at the first three 7-adic precisions.
            for exponent in range(1,4):
                modulus=7**exponent;period=periods[-1]*7**(exponent-1)
                def reduced(m):return tuple(x%modulus for x in observer(cs,m))
                need(all(reduced(m)==reduced(m+period) for m in (1,2,3)),'precision period upper bound')
                for prime in (2,3,7):
                    if period%prime==0:need(any(reduced(m)!=reduced(m+period//prime) for m in (1,2,3)),'precision period is exact')
        for m in range(1,7):frames.append(dict(m=m,source=field_frames(observer(coeff[0],m)),target=field_frames(observer(coeff[1],m))))
        report=dict(source=h,offset=d['offset'],port=q,source_count=[2*coeff[0][0][1],coeff[0][0][0]-2*coeff[0][0][1]],
                    run_coefficients=coeff,phase_periods=periods,frames=frames,ancestry=ancestry(h))
        reports.append(report);print('proved clone lift',h,d['offset'],q,report['source_count'],periods,flush=True)
    # Full rooted affine return, stronger than matching only the Hamilton scalar.
    for report in reports:
        if report['source']==1457:
            a,b=report['run_coefficients'];c=report['offset']
            for m in range(1,5):
                h,x,y=observer(a,m);need(observer(b,m)==(3*h+c,3*x,3*y+c),'uniform full golden observer return')
        if report['source']==1459 and report['offset']==2:
            a,b=report['run_coefficients']
            for m in range(1,5):
                h,x,y=observer(a,m);need(observer(b,m)==(3*h+2,3*x+2**m,3*y+2**m),'second exact golden return with recoverable dyadic coordinate')
    result=dict(literal=literal,lifts=reports,census=[{k:x[k] for k in ('source_H','paths','all_hits','strong_hits')} for x in census])
    (HERE/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print('audit',hashlib.sha256((HERE/'audit.json').read_bytes()).hexdigest(),flush=True)
    return result


if __name__=='__main__':run()

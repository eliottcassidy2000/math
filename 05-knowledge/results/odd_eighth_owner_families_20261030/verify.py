"""Exact family classification, boundary windows, and signed-cycle checks."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations,product
from math import gcd,comb,factorial,lcm
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
from scout import maximum,owner

HERE=Path(__file__).resolve().parent
GATES=0


def need(ok,message):
    global GATES
    GATES+=1
    if not ok: raise RuntimeError(message)


def safe(V,r):
    current=[(F(0),F(1))]
    for v in sorted(set(V)):
        walls=[((k+r)/v,(k+1-r)/v) for k in range(v)]
        out=[]; i=j=0
        while i<len(current) and j<len(walls):
            lo=max(current[i][0],walls[j][0]); hi=min(current[i][1],walls[j][1])
            if lo<=hi: out.append((lo,hi))
            if current[i][1]<walls[j][1]: i+=1
            else: j+=1
        current=out
        if not current: break
    return current


def measure(intervals): return sum((b-a for a,b in intervals),F())


def distance(x): return min(x%1,1-x%1)


def value(V,t): return min(distance(v*t) for v in V)


def owner_windows(q,E,r):
    c={u:max([u]+[v for v in E if v%q==u]) for u in range(1,q) if gcd(u,q)==1}
    eps=F(1,q)-r
    return [(F(a,q)-eps/c[pow(a,-1,q)],F(a,q)+eps/c[(-pow(a,-1,q))%q])
            for a in range(1,q) if gcd(a,q)==1]


def verify_max(V,r,witnesses=None):
    points=safe(V,r)
    need(bool(points), 'actual closed safe set attains claimed maximum')
    need(all(a==b for a,b in points), 'maximum consists of isolated points')
    gap=F(1,4*max(V)**2+1)
    need(not safe(V,r+gap), 'independent interval exclusion above exact rational maximum')
    if witnesses is not None:
        expected=sorted(set(witnesses)|{1-t for t in witnesses})
        need([a for a,b in points]==expected,'complete actual maximizing times')


def family_checks():
    generic=0; windows=0; weighted=[]
    for q in range(3,10):
        upper_values={}
        for mask in range(1<<(q-1)):
            E=[q+1+i for i in range(q-1) if mask>>i&1]
            core=list(range(1,q))+E
            r0=F(1,q+2) if q+1 in E else F(1,q+1)
            for lam in [F(1,3),F(2,3)]:
                r=r0+lam*(F(1,q)-r0)
                need(safe(core,r)==owner_windows(q,E,r),'uniform two-band owner windows')
                windows+=1
            s=owner(q,E)
            for m in sorted({1,s//2+1,q+1}):
                V=core+[q*m]; actual,times=maximum(V)
                if m==1:
                    missing=next(v for v in range(q+1,2*q+1) if v not in E)
                    need(actual==F(1,missing),'first tail is governed by the first missing speed')
                    expected=[F(a,missing) for a in range(1,missing) if gcd(a,missing)==1]
                    need(sorted(set(times)|{1-t for t in times})==expected,'all primitive first-gap maximizers')
                    eps=F(1,100*(1+max(V))**4)
                    c={u:max(v for v in V if v%missing==u)
                       for u in range(1,missing) if gcd(u,missing)==1}
                    need(measure(safe(V,actual-eps))/eps==2*sum((F(1,v) for v in c.values()),F()),
                         'first-tail terminal mass retains every actual directional owner')
                if 2*m>s:
                    need(actual==F(m,q*m+s),'uniform strict-branch formula')
                    if m>1:
                        a=pow(s,-1,q); h=(a*s-1)//q
                        t=F(a*m+h,q*m+s)
                        need(sorted(set(times)|{1-z for z in times})==sorted({t,1-t}),
                             'two directional-owner maxima on strict non-AP branch')
                else:
                    need(actual<=F(1,q+2),'complete exclusion above common top chamber')
                generic+=1
                if m==q+1: upper_values[mask]=actual
        for model in ['fair','graded','empty','full']:
            probs=[F(1,2) if model=='fair' else F(i+1,q+1) if model=='graded'
                   else F(model=='full') for i in range(q-1)]
            exact=F()
            for mask,r in upper_values.items():
                weight=F(1)
                for i,p in enumerate(probs): weight*=p if mask>>i&1 else 1-p
                exact+=weight*r
            m=q+1; recurrence=F(m,q*m+q+1)
            for u in reversed([u for u in range(1,q) if gcd(u,q)==1]):
                p=probs[u-1]
                recurrence=(1-p)*F(m,q*m+u)+p*recurrence
            need(exact==recurrence,'weighted least-owner recurrence agrees with every native subset')
            weighted.append({'q':q,'m':m,'model':model,'mean':str(exact)})
    # Full finite head; the theorem covers all m>=5.
    rows=[]; profiles={}
    for mask in range(128):
        E=[9+i for i in range(7) if mask>>i&1]; s=owner(8,E)
        head=[]
        for m in range(1,11):
            V=list(range(1,8))+E+[8*m]
            r,times=maximum(V)
            verify_max(V,r,times)
            if m<=4: head.append(str(r))
            if m>=5:
                need(r==F(m,8*m+s),'infinite-family formula matches finite overlap')
                eps=F(1,100*(1+max(V))**4)
                slope=measure(safe(V,r-eps))/eps
                need(slope==2*(F(1,s)+F(1,8*m)), 'actual terminal mass slope')
            rows.append({'mask':mask,'m':m,'s':s,'clearance':str(r),
                         'maximizers':[str(t) for t in sorted(set(times)|{1-t for t in times})]})
        profiles[mask]=tuple(head)+(s,)
        for m in [101,10**6]:
            a=pow(s,-1,8); h=(a*s-1)//8
            t=F(a*m+h,8*m+s); r=F(m,8*m+s)
            need(value(list(range(1,8))+E+[8*m],t)==r,'large exact native witness')
    groups=defaultdict(list)
    for e,profile in profiles.items(): groups[profile].append(e)
    need(len(groups)==14,'exactly fourteen entire clearance functions')
    classes=[]
    for profile,masks in sorted(groups.items(),key=lambda item:item[1][0]):
        pattern=''.join(str(masks[0]>>i&1) if len({e>>i&1 for e in masks})==1 else '*'
                        for i in range(7))
        cube=[e for e in range(128) if all(c=='*' or int(c)==e>>i&1 for i,c in enumerate(pattern))]
        need(cube==masks,'classification row is exactly its Boolean cube')
        classes.append({'pattern_9_through_15':pattern,'count':len(masks),'head':profile[:4],
                        's':profile[4],'masks':masks})
    fullfuture={e:tuple(profiles[e|a] for a in range(128)) for e in range(128)}
    upperfuture={e:tuple(owner(8,[9+i for i in range(7) if (e|a)>>i&1])
                         for a in range(128)) for e in range(128)}
    need(len(set(fullfuture.values()))==128,'all future additions separate every full source')
    need(len(set(upperfuture.values()))==16,'upper-branch future carrier has exactly four bits')
    # A separating continuation for every bit, without relying on a signature count.
    f={(r['mask'],r['m']):F(r['clearance']) for r in rows}
    for e in range(128):
        for bit in range(7):
            context=127^(1<<bit)
            out=e|context
            need(f[out,1]==(F(1,16) if e>>bit&1 else F(1,9+bit)),
                 'source bit reconstructed by first-missing continuation')
    return {'general_core_checks':generic,'literal_window_checks':windows,'weighted_means':weighted,
            'classes':classes,'finite_head_and_overlap':rows,
            'present_classes':14,'eventual_values':5,'upper_future_states':16,'full_future_states':128}


def packet_checks():
    cases=[]
    for A in [[4],[4,28],[4,32],[4,28,32]]:
        V=list(range(1,8))+[8*a for a in A]
        r,times=maximum(V); verify_max(V,r,times)
        cases.append({'packet':A,'maximum':str(r),'maximizers':[str(t) for t in times]})
    need([r['maximum'] for r in cases]==['4/33','4/33','31/257','3/25'],
         'equal first-arrival profiles before attachment diverge afterwards')
    # Check the exact first-arrival feasibility law above the AP top chamber.
    for A in [[1,2],[3,7],[4,32],[4,28,32],[2,3,11]]:
        for r in [F(9,80),F(3,25),F(31,257),F(1,8)]:
            arrivals=safe(A,r); tau=arrivals[0][0] if arrivals else None
            native=safe(list(range(1,8))+[8*a for a in A],r)
            predicted=tau is not None and tau<=1-8*r
            need(bool(native)==predicted,'joint zero-phase attachment keeps first safe arrival')
    horizon=[]
    for V,expected,refuted in [([1,6,7],F(4,13),F(3,8)),([1,2,3,12],F(3,13),F(1,4))]:
        r,times=maximum(V); verify_max(V,r,times)
        need(r==expected and r<refuted,'incoming whole-horizon and joint-fold hostiles reproduced')
        horizon.append({'speeds':V,'maximum':str(r),'refuted_lower_bound':str(refuted)})
    return {'profile_hostile':cases,'incoming_horizon_controls':horizon}


def matmul(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def harmonic_checks():
    cases=0; period_cases=0; samples=[]
    for q in range(3,9):
        n=q+1
        for mask in range(1,1<<(q-1),2):
            E=[q+1+i for i in range(q-1) if mask>>i&1]
            C=[u+q if u+q in E else u for u in range(1,q) if gcd(u,q)==1]
            s=min(C); L=n*lcm(*C)
            def totals(m):
                count=0; weighted=F()
                for e in C:
                    for j in range(1,n):
                        if gcd(j,n)!=1: continue
                        mult=max(0,(m//e-j)//n+1)
                        ell=pow(j,-1,n)
                        count+=2*mult
                        weighted+=2*mult*(F(1,ell)+F(1,n-ell))/q
                return count,weighted # second coordinate is m times slope
            for m in sorted({s-1,s,s+1,2*s+1}):
                V=list(range(1,q))+E+[q*m*j for j in range(1,q+1)]
                actual=safe(V,F(1,n))
                need(bool(actual)==(m>=s),'generalized harmonic threshold is iff')
                if m>=s:
                    points=[]
                    for e in C:
                        a=pow(e,-1,q)
                        for j in range(1,m//e+1):
                            if gcd(j,n)==1:
                                t=F(a,q)-F(j,q*n*m); points.extend([t,1-t])
                    need(len(points)==len(set(points)),'all generalized harmonic addresses are distinct')
                    need(actual==[(t,t) for t in sorted(points)],'complete native harmonic maximizing set')
                    count,mslope=totals(m)
                    need(count==len(points),'floor-period harmonic count')
                    eps=F(1,100*(1+max(V))**4)
                    need(measure(safe(V,F(1,n)-eps))/eps==mslope/m,'physical harmonic terminal slope')
                    later,laterweight=totals(m+L)
                    HE=sum((F(1,e) for e in C),F())
                    HU=sum((F(1,a) for a in range(1,n) if gcd(a,n)==1),F())
                    phi=sum(gcd(a,n)==1 for a in range(1,n))
                    need(later-count==2*L*phi*HE/n,'all-owner count affine-period recurrence')
                    need(laterweight-mslope==4*L*HE*HU/(q*n),'all-owner physical mass affine-period recurrence')
                    period_cases+=1
                    if q==8 and mask in [1,127]:
                        samples.append({'q':q,'mask':mask,'m':m,'count':count,'slope':str(mslope/m),'period':L})
                cases+=1
    return {'native_cases':cases,'exact_period_cases':period_cases,'samples':samples}


def inverse(A):
    n=len(A); M=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if M[i][j])
        M[j],M[pivot]=M[pivot],M[j]; a=M[j][j]
        M[j]=[v/a for v in M[j]]
        for i in range(n):
            if i!=j:
                a=M[i][j]; M[i]=[x-a*y for x,y in zip(M[i],M[j])]
    return [row[n:] for row in M]


def signed_cycle(S):
    n=len(S)
    plus=[[int(i==j)+S[i][j] for j in range(n)] for i in range(n)]
    minus=[[int(i==j)-S[i][j] for j in range(n)] for i in range(n)]
    U=matmul(plus,inverse(minus))
    mono=all(sum(x!=0 for x in row)==1 and all(x in [-1,0,1] for x in row) for row in U)
    return mono,U


def locally_ordered(S):
    n=len(S); D=[1]+S[0][1:]
    R=[[D[i]*S[i][j]*D[j] for j in range(n)] for i in range(n)]
    degrees=[sum(R[i][j]>0 for j in range(1,n)) for i in range(1,n)]
    return sorted(degrees)==list(range(n-1))


def cycle_insertions(S,U):
    n=len(S)
    nxt=[next(i for i in range(n) if U[i][j]) for j in range(n)]
    signs=[int(U[nxt[j]][j]) for j in range(n)]
    generated=set()
    for j in range(n):
        for a in [-1,1]:
            dest=nxt+[nxt[j]]; weights=signs+[-a*signs[j]]
            dest[j]=n; weights[j]=a
            R=[[0]*(n+1) for _ in range(n+1)]
            for start in range(n+1):
                at=start; prod=1
                for k in range(1,n+1):
                    prod*=weights[at]; at=dest[at]
                    R[at][start]=(-1)**(k-1)*prod
            need([row[:n] for row in R[:n]]==S,'signed insertion preserves the actual old tournament')
            need(all(R[i][j]==-R[j][i] for i in range(n+1) for j in range(n+1)),
                 'signed insertion produces a tournament skew matrix')
            need(locally_ordered(R),'every signed insertion is a native local-order child')
            generated.add(tuple(row[-1] for row in R[:-1]))
    native=set()
    for column in product([-1,1],repeat=n):
        R=[row+[x] for row,x in zip(S,column)]+[[-x for x in column]+[0]]
        if locally_ordered(R): native.add(tuple(column))
    need(len(generated)==2*n and generated==native,'all legal native attachments, with no multiplicity loss')


def tournament_checks():
    census=[]
    for n in range(1,7):
        edges=list(combinations(range(n),2)); good=0; spectra=Counter()
        for mask in range(1<<len(edges)):
            S=[[0]*n for _ in range(n)]
            for bit,(i,j) in enumerate(edges):
                S[i][j]=1 if mask>>bit&1 else -1; S[j][i]=-S[i][j]
            defects=0
            for i,j,k,l in combinations(range(n),4):
                pf=S[i][j]*S[k][l]-S[i][k]*S[j][l]+S[i][l]*S[j][k]
                defects+=int(abs(pf)==3)
            ordered=locally_ordered(S); good+=ordered
            need(ordered==(defects==0),'inherited local-order equivalence from actual four-sets')
            S2=matmul(S,S); energy=sum(x*x for row in S2 for x in row)
            need(energy==n*n*(n-1)*(n-1)//2-4*comb(n,4)-32*defects,
                 'exact quartic energy pays thirty-two per vortex')
            if n<=5:
                mono,U=signed_cycle(S)
                need(mono==ordered,'actual rational Cayley matrix is monomial exactly on local orders')
                if ordered: cycle_insertions(S,U)
            spectra[defects]+=1
        need(good==2**(n-1)*factorial(n-1),'inherited local-order count reproduced')
        census.append({'n':n,'labelled_tournaments':1<<len(edges),'local_orders':good,
                       'vortex_histogram':dict(sorted(spectra.items()))})
    for n in range(1,65):
        S=[[0 if i==j else (1 if i<j else -1) for j in range(n)] for i in range(n)]
        W=[[0]*n for _ in range(n)]; W[0][-1]=1
        for j in range(n-1): W[j+1][j]=-1
        minus=[[int(i==j)-S[i][j] for j in range(n)] for i in range(n)]
        plus=[[int(i==j)+S[i][j] for j in range(n)] for i in range(n)]
        need(matmul(minus,W)==plus,'literal signed-cycle Cayley transform at every tested order')
    src=HERE.parent/'regular_parent_symmetry_20261029'/'verify.py'
    spec=importlib.util.spec_from_file_location('native_hostile',src)
    old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    hostile=old.hostile_checks()
    need([r['group_order'] for r in hostile['tournaments']]==[81,27],'native 81/27 control retained')
    return {'census':census,'signed_cycle_identity_orders':[1,64],'hostile_81_27':hostile}


def eighth_checks():
    U=[1,3,5,7]
    chars=[[1]*4,[(-1)**((a-1)//2) for a in U],
           [(-1)**((a*a-1)//8) for a in U]]
    chars.append([a*b for a,b in zip(chars[1],chars[2])])
    need(matmul(chars,list(zip(*chars)))==[[4*int(i==j) for j in range(4)] for i in range(4)],
         'four independent real characters form a Walsh basis')
    for v,a in product(U,repeat=2):
        c=lambda x:(-1)**((x*x-1)//8)
        need(distance(F(v*a,8))==F(1,4)-F(c(v)*c(a),8),'absolute grid distance is one quadratic character')
    # Exact cos/sin numerators: zeta_8^a=(chi_2(a)+i chi_-2(a))/sqrt(2).
    gauss=[(sum(c*x for c,x in zip(chi,chars[2])),sum(c*y for c,y in zip(chi,chars[3]))) for chi in chars]
    need(gauss==[(0,0),(0,0),(4,0),(0,4)],'primitive eighth-root Gauss sums, exact sqrt-two coordinates')
    half=[a%4 for a in U]
    need(Counter(half)==Counter({1:2,3:2}),'phase squaring is two-to-one on primitive eighth roots')
    need({a*a%8 for a in U}=={1},'unit squaring is a different, four-to-one operation')
    coord={a:(chars[2][i],chars[3][i]) for i,a in enumerate(U)}
    for i,a in enumerate(U):
        x,y=coord[7*a%8],coord[5*a%8]
        need((F(x[0]+y[0],2),F(x[1]+y[1],2))==(0,-chars[3][i]),
             'incoming word-address transform sends chi_2 to negative i times chi_-2')
    return {'units':U,'characters':chars,'gauss_numerators_over_sqrt2':gauss}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,default=HERE/'certificate.json')
    args=parser.parse_args()
    data={'family':family_checks(),'packets':packet_checks(),
          'harmonic':harmonic_checks(),'tournaments':tournament_checks(),
          'eighths':eighth_checks(),'gates':GATES}
    raw=json.dumps(data,indent=2)+'\n'; args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps({'gates':GATES,'sha256':hashlib.sha256(raw.encode()).hexdigest()}))


if __name__=='__main__':main()

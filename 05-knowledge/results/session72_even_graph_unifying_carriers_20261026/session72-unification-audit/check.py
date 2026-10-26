#!/usr/bin/env python3
"""Independent Session72 audit. Stdlib-only, no producer imports, explicit gates."""
from itertools import product, combinations, permutations
from fractions import Fraction
from math import gcd,lcm,prod
import argparse,json,hashlib

def require(ok,why):
    if not ok: raise RuntimeError(why)

def rank(matrix,p=None):
    a=[[Fraction(v) if p is None else v%p for v in row] for row in matrix]
    if not a:return 0
    r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        z=1/a[r][j] if p is None else pow(a[r][j],-1,p)
        a[r]=[v*z if p is None else v*z%p for v in a[r]]
        for i in range(r+1,len(a)):
            z=a[i][j]
            a[i]=[v-z*w if p is None else (v-z*w)%p for v,w in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

def determinant(matrix):
    a=[list(row) for row in matrix];n=len(a);prior=1;sign=1
    require(all(len(row)==n for row in a),'square determinant input')
    if n==0:return 1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        value=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*value-a[i][k]*a[k][j]
                require(numerator%prior==0,'exact Bareiss division')
                a[i][j]=numerator//prior
            a[i][k]=0
        prior=value
    return sign*a[-1][-1]

def components(n,edges):
    parts=[{i} for i in range(n)]
    for a,b in edges:
        x=next(s for s in parts if a in s);y=next(s for s in parts if b in s)
        if x is not y:x.update(y);parts.remove(y)
    return parts

def incidence(E,q):
    return [[int(e[side]==v) for e in E] for side in range(2) for v in range(q)]

def response():
    supports=0;fields=[None,2,3,5,7,11]
    for q in range(1,5):
        choices=list(combinations(range(q),2))+[(i,i) for i in range(q)]
        for bits in product(range(2),repeat=len(choices)):
            E=sorted({e for bit,(i,j) in zip(bits,choices) if bit for e in ((i,j),(j,i))})
            if {a for a,b in E}!=set(range(q)):continue
            supports+=1
            B=incidence(E,q);c=len(components(2*q,[(a,q+b) for a,b in E]));r=2*q-c
            require(all(rank(B,p)==r for p in fields),'field-independent incidence rank')
            require((len(E)==r)==all(len([e for e in E if e[0] in C])==len(C)-1 for C in components(2*q,[(a,q+b) for a,b in E])),'forest iff independent')
            # Every edge column has two nonzero coordinates. No other edge has support inside it.
            for e in E:require([f for f in E if f[0]==e[0] and f[1]==e[1]]==[e],'each support ray extreme')
    # Truly full-state reset with duplicate pair fibres, rational nonuniform probabilities.
    X=[(0,0),(0,1),(1,0),(1,1),(0,0),(1,1)]
    tau=[0,2,1,3,4,5];q=[x[0] for x in X];N=len(X)
    base=[[int(a==v) for a,b in X] for v in range(2)]
    T=[[int(j==tau[i]) for j in range(N)] for i in range(N)]
    R=[[Fraction((j+1)*(q[i]+1)+(j%2),sum((k+1)*(q[i]+1)+(k%2) for k in range(N))) for j in range(N)] for i in range(N)]
    countR=[[((j+1)*(q[i]+1)+(j%2)) for j in range(N)] for i in range(N)]
    def expand(rows,M):return [[sum(M[i][j]*v[j] for j in range(N)) for i in range(N)] for v in rows]
    rows=base[:]
    for t in range(5):
        more=rows+expand(rows,T)+expand(rows,R)+expand(rows,countR)
        basis=[]
        for v in more:
            if rank(basis+[v])>len(basis):basis.append(v)
        rows=basis
    require(len(rows)==3,'arbitrary full-state q-dependent reset terminal rank')
    # Same time sample products reveal fourth pair coordinate despite terminal rank three.
    full=[[int(a==i and b==j) for a,b in X] for i,j in product(range(2),repeat=2)]
    require(rank(full)==4,'trajectory observations see joint pair')
    # Drop reset dependence: one terminal operation detects the AND bit.
    bad=[[int(j==(3 if X[i]==(1,1) else 0)) for j in range(N)] for i in range(N)]
    require(rank(rows+expand(base,bad))==4,'nonreset hostile')
    return {'symmetric_full_shore_supports':supports,'fields':['Q',2,3,5,7,11], 'duplicate_state_reset_terminal_rank':3,'trajectory_rank':4,'nonreset_rank':4}

def tournament():
    census=[];all_q={}
    for n in range(1,7):
        nonpath=[(i,j) for i in range(n) for j in range(i+2,n)]
        ps=list(permutations(range(n)));edges=list(combinations(range(n),2))
        labels=[]
        for mask in range(1<<len(nonpath)):
            A=[[0]*n for _ in range(n)]
            for i in range(n-1):A[i][i+1]=1
            for k,(i,j) in enumerate(nonpath):A[i][j]=mask>>k&1;A[j][i]=1-A[i][j]
            code=min(sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(edges)) for p in ps)
            labels.append(code)
        Q=sorted(set(labels));qi={v:i for i,v in enumerate(Q)};E=sorted({(qi[a],qi[labels[i^(len(labels)-1)]]) for i,a in enumerate(labels)})
        c=len(components(2*len(Q),[(a,len(Q)+b) for a,b in E]));r=2*len(Q)-c
        require(rank(incidence(E,len(Q)),2)==r,'native exact mod2 rank independent canonicalization')
        census.append({'n':n,'states':len(labels),'classes':len(Q),'E':len(E),'components':c,'rank':r,'kernel':len(E)-r})
        all_q[n]=labels
    require([v['kernel'] for v in census]==[0,0,0,0,27,621],'first native cycle threshold')
    q=all_q[5];plus=(15,21);minus=(47,5)
    left=sorted(q[x] for x in plus);right=sorted(q[x] for x in minus)
    require(left==right,'native source marginals')
    require(sorted(q[x^63] for x in plus)==sorted(q[x^63] for x in minus),'native target marginals')
    joint=lambda xs:sum((q[x],q[x^63])==(4,10) for x in xs)
    require((joint(plus),joint(minus))==(1,0),'native trajectory distinction')
    return {'census':census,'native_masks':[plus,minus],'joint_event':[4,10],'joint_probabilities':['1/2','0']}

def monoids():
    maps=list(product(range(3),repeat=3));ident=(0,1,2);q=(0,0,1)
    tested=noninvertible=0;monoid_sizes=set()
    for f,g in product(maps,repeat=2):
        M={ident};todo=[ident]
        while todo:
            h=todo.pop()
            for op in (f,g):
                z=tuple(op[h[i]] for i in range(3))
                if z not in M:M.add(z);todo.append(z)
        ordered=sorted(M);monoid_sizes.add(len(M));tested+=1
        noninvertible+=int(len(set(f))<3 or len(set(g))<3)
        signatures={tuple(q[h[x]] for h in ordered) for x in range(3)}
        D=[[int(q[h[x]]==a) for x in range(3)] for h in ordered for a in range(2)]
        # Pullback closure with counted q-dependent reset, over several fields.
        K=[[1+(q[i]+1)*(j+1) for j in range(3)] for i in range(3)]
        for p in (2,3,5):
            base=[[int(q[x]==a) for x in range(3)] for a in range(2)]
            basis=[]
            for v in base:
                if rank(basis+[v],p)>len(basis):basis.append(v)
            while True:
                old=len(basis)
                candidates=[[v[op[x]] for x in range(3)] for v in basis for op in (f,g)]
                candidates += [[sum(K[i][j]*v[j] for j in range(3)) for i in range(3)] for v in basis]
                for v in candidates:
                    if rank(basis+[v],p)>len(basis):basis.append(v)
                if len(basis)==old:break
            require(len(basis)==rank(D,p),'noninvertible monoid marginal-response theorem')
        # The generated point-indicator products separate all realizable signatures.
        indicators=[[int(tuple(q[h[x]] for h in ordered)==s) for x in range(3)] for s in signatures]
        require(rank(indicators,2)==len(signatures),'monoid product algebra rank')
    return {'ordered_generator_pairs':tested,'pairs_with_noninvertible_generator':noninvertible,'monoid_sizes':sorted(monoid_sizes),'fields':[2,3,5]}

def guard_return():
    summary=[];returns=0;walsh=None
    for n in range(1,6):
        nonpath=[(i,j) for i in range(n) for j in range(i+2,n)]
        states=[]
        for mask in range(1<<len(nonpath)):
            A=[[0]*n for _ in range(n)]
            for i in range(n-1):A[i][i+1]=1
            for k,(i,j) in enumerate(nonpath):A[i][j]=mask>>k&1;A[j][i]=1-A[i][j]
            states.append(A)
        guards=[]
        for p in permutations(range(n)):
            inv=tuple(p.index(i) for i in range(n));row=[]
            for A in states:
                legal=all(A[p[i]][p[i+1]] for i in range(n-1));row.append(int(legal))
                if legal:
                    B=[[A[p[i]][p[j]] for j in range(n)] for i in range(n)]
                    require(all(B[inv[i]][inv[i+1]] for i in range(n-1)),'native inverse guard legal after test')
                    C=[[B[inv[i]][inv[j]] for j in range(n)] for i in range(n)]
                    require(C==A,'native test-return recovers complete presentation');returns+=1
            guards.append(row)
        signatures={tuple(row[x] for row in guards) for x in range(len(states))}
        require(len(signatures)==len(states),'native guards separate states')
        raw={str(p):rank(guards,p) for p in (None,2,3,5,7)}
        bitrows=[sum(v<<x for x,v in enumerate(row)) for row in guards]
        pivots={};witness=[]
        for i,a in enumerate(bitrows):
            for j,b in enumerate(bitrows):
                z=a&b
                while z:
                    p=z.bit_length()-1
                    if p in pivots:z^=pivots[p]
                    else:pivots[p]=z;witness.append((i,j));break
        require(len(pivots)==len(states),'two guards suffice through order five')
        selected=[[a*b for a,b in zip(guards[i],guards[j])] for i,j in witness]
        det=determinant(selected)
        require(abs(det)==1,'actual two-guard rows form unimodular basis: every field')
        summary.append({'n':n,'states':len(states),'raw_ranks':raw,'two_guard_rank_every_field':len(pivots),'selected_product_determinant':det,'selected_product_rows':witness if n==5 else []})
        if n==5:
            visible=set()
            for guard in guards:
                transformed=guard[:];step=1
                while step<len(transformed):
                    for start in range(0,len(transformed),2*step):
                        for j in range(step):
                            a,b=transformed[start+j],transformed[start+j+step]
                            transformed[start+j]=a+b;transformed[start+j+step]=a-b
                    step*=2
                visible.update(i for i,v in enumerate(transformed) if v)
            forests=set()
            for mask in range(64):
                E=[e for k,e in enumerate(nonpath) if mask>>k&1]
                deg=[sum(v in e for e in E) for v in range(n)]
                if max(deg)<=2 and len(E)==n-len(components(n,E)):forests.add(mask)
            require(visible==forests and len(visible)==46,'n5 Walsh support exactly linear forests')
            walsh={'visible_characters':46,'absent_characters':18,'Q_rank':raw['None'],'additional_visible_dependencies':46-raw['None'],'free_edge_order':nonpath}
    require(summary[-1]['raw_ranks']=={'None':44,'2':43,'3':44,'5':44,'7':44},'raw n5 characteristic distinction')
    return {'census':summary,'legal_inverse_returns':returns,'law_protocol':'valid X only; failure absorbing; terminal alive','n5_walsh':walsh}

def domain_hostile():
    # f=2x on(0,1),g=3x on(0,2); both formal centres0, but no actual fixed point.
    f_then_g=(Fraction(0),min(Fraction(1),Fraction(2,2)))
    g_then_f=(Fraction(0),min(Fraction(2),Fraction(1,3)))
    require(f_then_g==(0,1) and g_then_f==(0,Fraction(1,3)),'partial composite domains')
    require(f_then_g!=g_then_f,'formal affine commutation does not imply partial commutation')
    require(not(Fraction(0)<0<Fraction(1)),'formal common centre outside domain')
    return {'formal_slopes':[2,3],'formal_common_centre':0,'f_domain':['0','1'],'g_domain':['0','2'],'g_after_f_domain':[str(x) for x in f_then_g],'f_after_g_domain':[str(x) for x in g_then_f],'actual_common_fixed_point':False}

def depth_six():
    n=6;states=[]
    for mask in range(1024):
        A=[[0]*n for _ in range(n)];k=0
        for i in range(n):
            for j in range(i+1,n):
                forward=1 if j==i+1 else mask>>k&1
                if j!=i+1:k+=1
                A[i][j]=forward;A[j][i]=1-forward
        states.append(A)
    guards=[sum(int(all(A[p[i]][p[i+1]] for i in range(5)))<<x for x,A in enumerate(states)) for p in permutations(range(n))]
    basis={};selected=[]
    def add(v,word):
        while v:
            p=v.bit_length()-1
            if p in basis:v^=basis[p]
            else:basis[p]=v;selected.append(word);return
    add((1<<1024)-1,[])
    for i,g in enumerate(guards):add(g,[i])
    raw=len(basis);examined=0
    for i,g in enumerate(guards):
        for j in range(i,len(guards)):
            add(g&guards[j],[i,j]);examined+=1
            if len(basis)==1024:break
        if len(basis)==1024:break
    require(raw==272 and len(basis)==1024 and examined==144190,'native n6 two-guard depth certificate')
    digest=hashlib.sha256(json.dumps(selected,separators=(',',':')).encode()).hexdigest()
    require(digest=='2bb97cd8d267b0e1e11a2274b16ea999e2e7996d0667e28d36c3dcf2817e819e','independent selection agrees with producer certificate')
    return {'n':n,'states':1024,'raw_F2_rank':raw,'two_guard_F2_and_Q_rank':len(basis),'unordered_products_examined':examined,'selection_sha256':digest,'every_field_claim':False}

def cyclic():
    profiles=edgegraphs=delta_states=assignments=0
    for k in range(1,5):
        pairs=list(combinations(range(k),2))
        for ds in product(range(1,5),repeat=k):
            profiles+=1;L=lcm(*ds)
            states=list(product(*(range(d) for d in ds)))
            # Quotient by the explicitly generated diagonal translation group.
            signature=lambda x:tuple((x[j]-x[i])%gcd(ds[i],ds[j]) for i,j in pairs)
            orbits={};seen=set()
            for x in states:
                if x in seen:continue
                orb={tuple((a+t)%d for a,d in zip(x,ds)) for t in range(L)}
                require(not seen&orb and len(orb)==L,'translation orbit length')
                seen|=orb;s=signature(x)
                require(s not in orbits,'pair residues classify orbit');orbits[s]=orb
            require(len(orbits)==prod(ds)//L,'orbit count')
            for values in product(*(range(gcd(ds[i],ds[j])) for i,j in pairs)):
                assignments+=1;r=dict(zip(pairs,values))
                triangular=all((r[i,j]+r[j,h]-r[i,h])%gcd(ds[i],ds[j],ds[h])==0 for i,j,h in combinations(range(k),3))
                require(triangular==(values in orbits),'triangle sufficiency/exact image')
            # Sparse-edge injectivity: count its zero-signature kernel, compare full diagonal.
            for flags in product(range(2),repeat=len(pairs)):
                edgegraphs+=1;E=[e for e,flag in zip(pairs,flags) if flag]
                ker=sum(all((x[j]-x[i])%gcd(ds[i],ds[j])==0 for i,j in E) for x in states)
                delta_states+=len(states)
                layer_condition=True
                for p in (2,3):
                    power=p
                    while power<=max(ds):
                        I={i for i,d in enumerate(ds) if d%power==0}
                        if I:
                            parts=[c&I for c in components(k,[(i,j) for i,j in E if i in I and j in I]) if c&I]
                            if len(parts)!=1:layer_condition=False
                        power*=p
                require((ker==L)==layer_condition,'prime-power layer graph exact criterion')
    # Explicit noncyclic odd-order group obstruction with every pair orbit full.
    orbit={(u,v,(u+v)%3) for u,v in product(range(3),repeat=2)}
    require((0,0,1) not in orbit,'noncyclic triple separation')
    require(all(len({(x[i],x[j]) for x in orbit})==9 for i,j in combinations(range(3),2)),'noncyclic pairs identical')
    return {'profiles':profiles,'candidate_edge_assignments':assignments,'sparse_edge_graphs':edgegraphs,'sparse_delta_states':delta_states,'noncyclic_F3_squared_orbit':len(orbit)}

POLY={1:3,2:7,3:11,4:19,6:67}
def mul(x,y,n):
    r=0
    while y:
        if y&1:r^=x
        y>>=1;x<<=1
        if x>>n:x^=POLY[n]
    return r
def power(x,k,n):
    r=1
    while k:
        if k&1:r=mul(r,x,n)
        x=mul(x,x,n);k>>=1
    return r

def tensor():
    rows=[]
    for ds in [(1,2,3),(2,2),(2,3),(2,2,3),(2,3,3),(4,2)]:
        L=lcm(*ds);D=prod(ds)
        require(all(power(x,2**L-1,L)==1 for x in range(1,2**L)),'ambient binary field invertibility')
        gen=[]
        for d in ds:
            a=next(a for a in range(2**L) if power(a,2**d,L)==a and all(power(a,2**e,L)!=a for e in range(1,d) if d%e==0))
            gen.append(a)
        reps=[];seen=set()
        for s in product(*(range(d) for d in ds)):
            if s in seen:continue
            reps.append(s);seen|={tuple((a+t)%d for a,d in zip(s,ds)) for t in range(L)}
        basis=list(product(*(range(d) for d in ds)))
        def eval(t,s):
            z=1
            for a,e,r in zip(gen,t,s):z=mul(z,power(a,e*2**r,L),L)
            return z
        columns=[[eval(t,s) for s in reps] for t in basis]
        M=[[(col[j]>>bit)&1 for col in columns] for j in range(len(reps)) for bit in range(L)]
        require(len(M)==D and rank(M,2)==D,'all-arity finite-field tensor evaluation bijective')
        # Coordinates of products in each embedded factor power basis.
        coord=[]
        for d,a in zip(ds,gen):
            table={}
            for coefficients in product(range(2),repeat=d):
                z=0
                for e,b in enumerate(coefficients):
                    if b:z^=power(a,e,L)
                table[z]=coefficients
            require(len(table)==2**d,'factor power basis independent')
            coord.append(table)
        for ti,t in enumerate(basis):
            for ui,u in enumerate(basis):
                coeff=[coord[i][power(a,t[i]+u[i],L)] for i,a in enumerate(gen)]
                terms=[v for v in basis if all(coeff[i][v[i]] for i in range(len(ds)))]
                for j,s in enumerate(reps):
                    z=0
                    for v in terms:z^=eval(v,s)
                    require(z==mul(columns[ti][j],columns[ui][j],L),'tensor multiplication preserved')
        rows.append({'degrees':ds,'L':L,'factors':len(reps),'dimension':D,'rank':D})
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--inject',choices=['field','trajectory','sparse']);a=p.parse_args()
    out={'response':response(),'monoids':monoids(),'guard_return':guard_return(),'depth_six':depth_six(),'domain_hostile':domain_hostile(),'tournament':tournament(),'cyclic':cyclic(),'tensor':tensor()}
    if a.inject=='field':require(rank([[1,1],[1,1]],2)==2,'deliberate field-rank corruption')
    if a.inject=='trajectory':require(out['response']['trajectory_rank']==out['response']['duplicate_state_reset_terminal_rank'],'deliberate trajectory conflation')
    if a.inject=='sparse':require(len({(u,v) for u,v,w in product(range(3),repeat=3) if w==(u+v)%3})<9,'deliberate pairwise gluing corruption')
    print(json.dumps({'status':'PASS','checks':out},sort_keys=True,indent=2))
if __name__=='__main__':main()

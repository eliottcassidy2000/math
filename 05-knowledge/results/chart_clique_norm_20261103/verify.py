"""Exact arithmetic, clique-carrier and geometry controls. Checks survive -O."""
from collections import Counter,defaultdict
from itertools import combinations,product
from math import comb,isqrt,gcd
from pathlib import Path
import json,random
from invert import recover

HERE=Path(__file__).resolve().parent

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

def F(q,t):return q**3+1-t**3+3*q*t

def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def sigma_character(q):
    f=factors(q)
    if any(p%3==2 and e%2 for p,e in f.items()):return 0
    ans=1
    for p,e in f.items():
        if p%3==1:ans*=e+1
    return ans

def norm_controls():
    saved={}
    for q in list(range(1,501))+[1729]:
        bd=isqrt(4*q);groups=defaultdict(list)
        for t in range(-bd,bd+1):groups[F(q,t)].append(t)
        collisions={n:ts for n,ts in groups.items() if len(ts)>1}
        sizes=Counter(map(len,collisions.values()));sq=int(isqrt(q)**2==q)
        require(sizes[2]==2*sq and sizes[3]==sigma_character(q)-sq,'complete integer fibre law')
        require(all(len(ts)<=3 for ts in groups.values()),'degree three')
        # An independent norm enumeration checks the classical factor formula.
        norm=sum(a*a-a*b+b*b==q for a in range(-bd,bd+1) for b in range(-bd,bd+1))
        require(norm==6*sigma_character(q),'Eisenstein representation count')
        if q in [1,3,5,7,25,49,55,1729]:saved[q]={'fibres':collisions,'norm_representations':norm}
    precision=[]
    for q in [2,5,8,11,17,23]:
        for m in range(1,7):
            modulus=3**(m+1);values={F(q,t)%modulus for t in range(3**m)}
            residues={F(q,r)%9 for r in range(3)}
            expected={a for a in range(modulus) if a%9 in residues}
            require(values==expected,'sharp ternary precision bijection')
            require(all(recover(q,F(q,t)%modulus,m)==t for t in range(3**m)),'constructive inverse')
            precision.append([q,m,len(values)])
    def vp3(n):
        n=abs(n);require(n!=0,'finite valuation');v=0
        while n%3==0:n//=3;v+=1
        return v
    for q in [2,5,11,17]:
        for t,u in combinations(range(-30,31),2):
            want=0 if (t-u)%3 else vp3(t-u)+1
            require(vp3(F(q,t)-F(q,u))==want,'exact valuation law')
    require(F(1,-2)==F(1,1) and F(3,0)==F(3,3),'nonunit hostiles')
    table=[]
    for a in range(1,10):
        digit=lambda z:1+(z-1)%9
        s=100*a+10*digit(a+6)+digit(a+3);t=a-5
        table.append({'index':a,'solfeggio_integer':s,'trace':t,'count':F(5,t),
                      'mod9':F(5,t)%9,'mod27':F(5,t)%27})
    require(len({r['mod9'] for r in table})==3 and len({r['mod27'] for r in table})==9,'nine states precision')
    square=[[108,148,122],[140,126,112],[130,104,144]]
    lines=square+list(map(list,zip(*square)))+[[square[i][i] for i in range(3)],[square[i][2-i] for i in range(3)]]
    require(all(sum(row)==378 for row in lines),'inherited magic arrangement')
    return {'q_through':500,'selected':saved,'precision':precision,'digit_table':table,
            'solfeggio_sum':sum(r['solfeggio_integer'] for r in table),'magic_square':square}

def keys(n,targets):
    r=len(targets);ans=[('one',)]
    ans += [('edge',c,(a,b)) for a,b in combinations(range(n),2) for c in range(r-1)]
    for c,t in enumerate(targets):
        for k in range(3,min(n,t-1)+1):ans += [('clique',c,s) for s in combinations(range(n),k)]
        if n>=t:ans.append(('total',c,t))
    return ans

def clique(graph,c,s):return all(graph[tuple(sorted(e))]==c for e in combinations(s,2))

def total(graph,n,c,t):return sum(clique(graph,c,s) for s in combinations(range(n),t))

def coordinates(graph,n,targets):
    out=[]
    for key in keys(n,targets):
        if key[0]=='one':out.append(1)
        elif key[0]=='total':out.append(total(graph,n,key[1],key[2]))
        else:out.append(int(clique(graph,key[1],key[2])))
    return out

def dimension(n,targets):
    return (1+(len(targets)-1)*comb(n,2)+
            sum(sum(comb(n,k) for k in range(3,min(n,t-1)+1))+int(n>=t) for t in targets))

def update_matrix(n,N,targets,fixed):
    old=keys(n,targets);idx={k:i for i,k in enumerate(old)};r=len(targets)
    def flag(c,s):
        row=[0]*len(old);a=tuple(v for v in s if v<n)
        if any(fixed[e]!=c for e in combinations(s,2) if e[1]>=n):return row
        if len(a)<2:row[0]=1
        elif len(a)==2:
            if c<r-1:row[idx['edge',c,a]]=1
            else:
                row[0]=1
                for d in range(r-1):row[idx['edge',d,a]]=-1
        else:row[idx['clique',c,a]]=1
        return row
    M=[]
    for key in keys(N,targets):
        if key[0]=='one':M.append([1]+[0]*(len(old)-1))
        elif key[0]=='total':
            c,t=key[1:];row=[0]*len(old)
            if n>=t:row[idx[key]]=1
            for s in combinations(range(N),t):
                if max(s)<n:continue
                term=flag(c,s);row=[a+b for a,b in zip(row,term)]
            M.append(row)
        else:M.append(flag(key[1],key[2]))
    return M

def rank_mod(rows,p):
    basis={}
    for row in rows:
        v=[a%p for a in row]
        for j in sorted(basis):
            if v[j]:
                a=v[j];v=[(x-a*y)%p for x,y in zip(v,basis[j])]
        nz=next((j for j,x in enumerate(v) if x),None)
        if nz is not None:
            inv=pow(v[nz],-1,p);basis[nz]=[x*inv%p for x in v]
    return len(basis)

def colour_graph(n,word):return dict(zip(combinations(range(n),2),word))

def coefficient_minor(n,targets):
    """Independent exact polynomial-coefficient audit of the integral pivot proof."""
    import sympy as sp
    labels=keys(n,targets);last=len(targets)-1;pivots=[]
    for key in labels:
        if key[0]=='one':pivots.append(());continue
        c=key[1];S=tuple(range(key[2])) if key[0]=='total' else key[2]
        es=list(combinations(S,2))
        pivots.append(tuple((e,(1 if j==1 else 0) if c==last else c) for j,e in enumerate(es)))
    def coeff(key,mono):
        if key[0]=='one':return int(not mono)
        c=key[1];edges={e for e,d in mono};vertices={v for e in edges for v in e}
        if key[0]=='total':
            t=key[2]
            if c==last:
                return (-1)**len(mono)*comb(n-len(vertices),t-len(vertices)) if len(vertices)<=t else 0
            return int(len(vertices)==t and edges==set(combinations(sorted(vertices),2)) and all(d==c for e,d in mono))
        S=key[2];full=set(combinations(S,2))
        if c==last:return (-1)**len(mono) if edges<=full else 0
        return int(edges==full and all(d==c for e,d in mono))
    matrix=sp.Matrix([[coeff(key,mono) for mono in pivots] for key in labels])
    det=int(matrix.det(method='domain-ge'))
    require(abs(det)==1,'explicit integral coefficient minor')
    return {'n':n,'targets':targets,'dimension':len(labels),'determinant':det}

def clique_controls():
    target=(3,3,5);ranks={};universes={}
    for n in [3,4]:
        graphs=[colour_graph(n,w) for w in product(range(3),repeat=comb(n,2))]
        universes[n]=graphs;rows=[coordinates(g,n,target) for g in graphs]
        ranks[n]={p:rank_mod(rows,p) for p in [2,3,101]}
        require(set(ranks[n].values())=={dimension(n,target)},'full-universe rank')
        # Actual isolating clique probes, with all original edges kept.
        for key in keys(n,target):
            if key[0] not in ['edge','clique']:continue
            _,c,S=key;t=target[c];N=n+t-len(S)
            for g in graphs:
                ext=dict(g)
                for a,b in combinations(range(N),2):
                    if b<n:continue
                    ext[a,b]=c if a>=n or a in S else (c+1)%3
                require(total(ext,N,c,t)-total(g,n,c,t)==int(clique(g,c,S)),'actual clique probe')
    # Exhaust all 27 attachments on all 27 old n3 colourings.
    checked=0
    for cut in product(range(3),repeat=3):
        fixed={(v,3):c for v,c in enumerate(cut)};M=update_matrix(3,4,target,fixed)
        for g in universes[3]:
            q=coordinates(g,3,target);actual=coordinates(g|fixed,4,target)
            require([sum(a*b for a,b in zip(row,q)) for row in M]==actual,'closed integral update')
            checked+=1
    rng=random.Random(126)
    for _ in range(100):
        g=colour_graph(5,[rng.randrange(3) for _ in range(10)])
        fixed={e:rng.randrange(3) for e in combinations(range(8),2) if e[1]>=5}
        M=update_matrix(5,8,target,fixed);q=coordinates(g,5,target)
        require([sum(a*b for a,b in zip(row,q)) for row in M]==coordinates(g|fixed,8,target),'three-vertex update')
    sample=[coordinates(colour_graph(5,[rng.randrange(3) for _ in range(10)]),5,target) for _ in range(512)]
    for c in range(3):sample.append(coordinates(colour_graph(5,[c]*10),5,target))
    for k in [3,4]:
        for S in combinations(range(5),k):
            g={e:2 if all(v in S for v in e) else 0 for e in combinations(range(5),2)}
            sample.append(coordinates(g,5,target))
    ranks[5]={p:rank_mod(sample,p) for p in [2,3,101]}
    require(set(ranks[5].values())=={39},'order-five sampled lower bound')
    minors=[coefficient_minor(n,t) for n,t in [(3,(3,3,5)),(4,(3,3,5)),(5,(3,3,5)),(6,(3,3,6)),(5,(3,4,5)),(4,(3,3,3,3))]]
    return {'ranks':ranks,'integral_coefficient_minors':minors,
            'n5_scope':'512 seeded plus18 planted sources; all-size theorem supplies exact rank',
            'exhaustive_n3_updates':checked,'sampled_three_birth_updates':100,
            'dimensions':{str(t):{n:dimension(n,t) for n in [0,1,2,3,4,5,6,57,91]} for t in [(3,3,5),(3,3,6)]}}

FRACTIONS=[(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
def step(n):
    for j,(a,b) in enumerate(FRACTIONS):
        if n%b==0:return n//b*a,chr(65+j)
    raise RuntimeError('total program')

def primegame_controls():
    n=1729;seen={};orbit=[];word=[]
    while n not in seen:
        seen[n]=len(orbit);orbit.append(n);n,j=step(n);word.append(j)
    for a in [31,37,31*37,31**3*41]:
        for x in orbit+[1,11,13,17,33,55,65,85]:
            y,j=step(x);z,k=step(a*x)
            require(z==a*y and k==j,'inert spectator equivariance')
    require(step(11)==(13,'K') and step(13)==(11,'J'),'reciprocal pair')
    require(step(33)==(29,'E'),'pair guard hostile')
    blocks={'X0':'JK','X1':'AIMN','X2':'CEFH','X3':'BDGL'};stats={}
    primes=[2,3,5,7,11,13,17,19,23,29]
    for name,letters in blocks.items():
        columns=[];occurrences=0;support=set()
        for letter in letters:
            a,b=FRACTIONS[ord(letter)-65];fa,fb=factors(a),factors(b)
            support|=set(fa)|set(fb);occurrences+=sum(fa.values())+sum(fb.values())
            columns.append([fa.get(p,0)-fb.get(p,0) for p in primes])
        stats[name]={'support':sorted(support),'width':len(support),'factor_occurrences':occurrences,
                     'valuation_rank_mod101':rank_mod(columns,101)}
    require([stats[x]['width'] for x in ['X1','X2','X3']]==[5,6,7],'defined quartet support complexity')
    reciprocal=[(chr(65+i),chr(65+j)) for i,j in combinations(range(14),2)
                if FRACTIONS[i][0]*FRACTIONS[j][0]==FRACTIONS[i][1]*FRACTIONS[j][1]]
    require(reciprocal==[('J','K')],'sole reciprocal pair')
    return {'orbit1729':orbit,'word':''.join(word),'cycle_start':seen[n],
            'factorization':factors(1729),'unit_group_exponent':36,
            'quartets':blocks,'block_statistics':stats,'unique_reciprocal_pair':reciprocal}

def geometry_controls():
    import sympy as s
    r=s.sqrt(13);c=(r-1)/4;sn=s.sqrt((1+r)/8);h=(r-1)*sn/3
    simp=lambda x:s.simplify(s.expand(x))
    u,v=-s.Integer(1),-h;vertices=[(u,v+h)]
    for _ in range(5):
        u,v=simp(c*u+sn*v),simp(-sn*u+c*v)
        vertices.append((u,simp(v+h)))
    require(simp(vertices[-1][0]-1)==0 and simp(vertices[-1][1])==0,'closed unit-chain endpoints')
    for (x,y),(X,Y) in zip(vertices,vertices[1:]):require(simp((X-x)**2+(Y-y)**2-1)==0,'unit fence')
    area=-sum(vertices[i][0]*vertices[(i+1)%6][1]-vertices[(i+1)%6][0]*vertices[i][1] for i in range(6))/2
    area=simp(area);product_area=simp(area**2)
    require(simp(product_area-(89+26*r)/16)==0,'Friedman product exact value')
    best=3*(2+s.sqrt(3));require(float(product_area-best)>0,'strict pseudoforest product barrier')
    # Three inherited boundary preimages, checked against the cancelled map by limits.
    x,y,z=s.symbols('x y z');U=1+3*x*y;G=1-4*x*y-x*x*z
    Fs=[s.cancel((2*U+U**2-3*U**4*G**2)/x**2),s.cancel((1+U-2*U**3*G**2)/x),x*G]
    points=[(0,1,-s.Rational(29,2)),(1,-s.Rational(1,3),s.Rational(7,3)),(-1,1,5)]
    require(all([s.expand(f).subs(dict(zip((x,y,z),p))) for f in Fs]==[0,1,0] for p in points),'quartic boundary replay')
    w,a,q,t=s.symbols('w a q t',nonzero=True)
    P=a-2*a**3;Q=a*a-3*a**4
    require(s.expand(w**4-w*w+2*P*w-Q-(w-a)**2*(w*w+2*a*w+3*a*a-1))==0,'quartic discriminant factorization')
    require(s.simplify(P.subs(a,t/s.sqrt(6*q))-(3*q*t-t**3)/(3*q*s.sqrt(6*q)))==0,'Frobenius cubic to quartic map')
    require(s.simplify(Q.subs(a,t/s.sqrt(6*q))-(2*q*t*t-t**4)/(12*q*q))==0,'retained quartic Q coordinate')
    for tt in range(-4,5):require(tt*tt not in [5,15],'nine trace points avoid cusp and node')
    return {'upper_vertices':[[str(a),str(b)] for a,b in vertices],
            'upper_vertices_numeric':[[float(a),float(b)] for a,b in vertices],
            'face_area':str(area),'product':str(product_area),
            'pseudoforest_product_max_n12':str(best),'product_gain':float(product_area-best),
            'theta_graph':{'vertices':11,'edges':12,'cycle_rank':2},
            'fence_n12':{'pseudoforest_area_max':3,'square_grid_area':4},
            'quartic_trace_transport':{'q':5,'traces':list(range(-4,5)),
                                       'finite_fibre_size_for_each':2,'cusp_test':'t^2=q','node_test':'t^2=3q'},
            'quartic_boundary_points':[[str(a) for a in p] for p in points]}

def main():
    result={'norm':norm_controls(),'clique':clique_controls(),'primegame':primegame_controls(),'geometry':geometry_controls()}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'integer_fibres_at1729':result['norm']['selected'][1729],
                     'clique_dimensions':result['clique']['dimensions'],
                     'primegame1729':result['primegame']['orbit1729'],
                     'product_gain':result['geometry']['product_gain']},indent=2))

if __name__=='__main__':main()

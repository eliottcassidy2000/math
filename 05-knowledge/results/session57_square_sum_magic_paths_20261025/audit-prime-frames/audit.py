import itertools as it, collections as co, json, math

def need(v,msg):
    if not v: raise RuntimeError(msg)

def add(a,b,p):return tuple((x+y)%p for x,y in zip(a,b))
def neg(a,p):return tuple(-x%p for x in a)
def sub(a,b,p):return add(a,neg(b,p),p)
def scale(c,a,p):return tuple(c*x%p for x in a)
def dirs(S,p):
    out=set()
    for a,b in it.combinations(S,2):
        x,y=sub(a,b,p)
        out.add(p if not x else y*pow(x,-1,p)%p)
    return out

def census(p):
    pts=list(it.product(range(p),repeat=2)); hist=co.Counter(); tests=0
    for S in it.combinations(pts,p):
        ds=dirs(S,p);hist[len(ds)]+=1
        need(len(ds)==1 or len(ds)>=(p+3)//2,'direction bound')
        for r in range(p+1):
            vals=[x if r==p else (y-r*x)%p for x,y in S]
            # Exact cyclotomic zero test: coefficient vector constant.
            counts=co.Counter(vals)
            zero=all(counts[i]==counts[0] for i in range(p))
            need(zero==(r not in ds),'Fourier/direction equivalence')
            tests+=1
    need(sum(hist.values())==math.comb(p*p,p),'census total')
    need(hist[1]==p*(p+1),'line total')
    return {'histogram':dict(sorted(hist.items())),'exact_fourier_tests':tests}

def frames(p):
    pts=list(it.product(range(p),repeat=2)); allframes=set(); duplicates=0
    for r in range(p+1):
        v=(0,1) if r==p else (1,r); w=(1,0) if r==p else (0,1)
        L=frozenset(scale(i,v,p) for i in range(p))
        for j in range(p):
            B=frozenset(add(scale(j,w,p),b,p) for b in L)
            for f in it.product(range(p),repeat=p):
                D=frozenset(add(scale(i,w,p),scale(f[i],v,p),p) for i in range(p))
                sums={add(b,d,p) for b in B for d in D}
                need(len(sums)==p*p,'sum bijective')
                need(len({sub(b,d,p) for b in B for d in D})==p*p,'difference bijective')
                key=frozenset((B,D))
                need(len(key)==2,'equal factors impossible')
                if key in allframes:duplicates+=1
                allframes.add(key)
    cp=(p+1)*p**(p-1)-p*(p+1)//2
    need(len(allframes)==cp*p*p,'frame count')
    need(duplicates==p*(p+1)*p*p//2,'double line correction')
    # Canonical orient and form every gauge orbit; test disjointness in ambient n=3
    # for all plane frames, but only once per orbit, embedding into z=0.
    pending=set(allframes); orbits=0; badhist=co.Counter()
    for frame in allframes:
        if frame not in pending:continue
        B,D=tuple(frame); db,dd=dirs(B,p),dirs(D,p)
        if len(db)!=1 or (len(dd)==1 and next(iter(db))>next(iter(dd))): B,D=D,B
        orbit=set()
        for t in pts:
            Bt=frozenset(add(b,t,p) for b in B);Dt=frozenset(sub(d,t,p) for d in D)
            orbit.add(frozenset((Bt,Dt)))
        need(len(orbit)==p*p,'free canonical gauge')
        pending.difference_update(orbit);orbits+=1
        B3=[b+(0,) for b in B];D3=[d+(0,) for d in D]
        bad=0
        for t in it.product(range(p),repeat=3):
            bt={add(b,t,p) for b in B3};dt={sub(d,t,p) for d in D3}
            bad+=bool(bt&dt)
        need(bad==p*p,'bad gauge count')
        badhist[bad]+=1
    need(not pending and orbits==cp,'orbit total')
    return {'frames_in_plane':len(allframes),'both_line_duplicates':duplicates,'gauge_orbits':orbits,'bad_gauge_histogram_n3':dict(badhist),'predicted_disjoint_n3':cp*(p**3-p*p)}

out={'direction_census':{p:census(p) for p in (3,5)},'frame_and_gauge_census':{p:frames(p) for p in (3,5)}}
need(out['direction_census'][5]['histogram']=={1:30,4:1500,5:15600,6:36000},'p5 expected census')
print(json.dumps(out,sort_keys=True,indent=2))

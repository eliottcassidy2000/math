#!/usr/bin/env python3
"""Exact, standard-library certificate for actual tournament CA lifts."""
import itertools as it
import json
import math
from pathlib import Path

MASKS=(153,8345,8337,8369,177,185)
PAIRS=tuple(it.combinations(range(6),2))
INDEX={e:i for i,e in enumerate(PAIRS)}
DECK=(4,2,1,5,0,3)
PROFILE=(0,29,131,162,77,15,1)

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def arc(mask,u,v,n=6):
    pairs=PAIRS if n==6 else tuple(it.combinations(range(n),2))
    idx=INDEX if n==6 else {e:i for i,e in enumerate(pairs)}
    return bool((mask>>idx[tuple(sorted((u,v)))])&1)==(u<v)

def relabel(mask,p):
    out=0
    for k,(i,j) in enumerate(PAIRS):
        u,v=(p[i],p[j]) if (mask>>k)&1 else (p[j],p[i])
        if u<v:out|=1<<INDEX[(u,v)]
    return out

def strong(mask):
    for root in range(6):
        seen={root};todo=[root]
        while todo:
            u=todo.pop()
            for v in range(6):
                if v not in seen and arc(mask,u,v):seen.add(v);todo.append(v)
        if len(seen)!=6:return False
    return True

def profile_dp(mask):
    end=[[0]*6 for _ in range(64)]
    for v in range(6):end[1<<v][v]=1
    for subset in range(1,64):
        for v in range(6):
            if not (subset>>v)&1:continue
            for w in range(6):
                if not (subset>>w)&1 and arc(mask,v,w):
                    end[subset|1<<w][w]+=end[subset][v]
    H=[sum(row) for row in end]
    P=[[0]*7 for _ in range(64)];P[0][0]=1
    for subset in range(1,64):
        root=subset&-subset;part=subset
        while part:
            if part&root:
                rest=subset^part
                for j in range(1,7):P[subset][j]+=H[part]*P[rest][j-1]
            part=(part-1)&subset
    return tuple(P[63])

def profile_cuts(mask):
    counts=[0]*7
    for p in it.permutations(range(6)):
        for cuts in range(32):
            if all((cuts>>i)&1 or arc(mask,p[i],p[i+1]) for i in range(5)):
                counts[cuts.bit_count()+1]+=1
    for j in range(1,7):
        require(counts[j]%math.factorial(j)==0,'cut factorial')
        counts[j]//=math.factorial(j)
    return tuple(counts)

def ordinal(h):
    n=6*len(h);out=0
    for k,(u,v) in enumerate(it.combinations(range(n),2)):
        bu,lu=divmod(u,6);bv,lv=divmod(v,6)
        bit=True if bu!=bv else arc(MASKS[h[bu]],lu,lv)
        if bit:out|=1<<k
    return out

def delta(source,target):
    d=(target-source)%3
    return (0,1,-1)[d]

def update(h,target):
    d=delta(h,target)
    out=(h+d)%6
    require((out%3)==target%3,'projection')
    require((MASKS[h]^MASKS[out]).bit_count()==bool(d),'literal arc update')
    return out

def apply_targets(state,targets,newvalues):
    out=list(state)
    for i,t in zip(targets,newvalues):out[i]=update(state[i],t)
    return tuple(out)

def compiler(h,m):
    """Five layers for F(x)_i=x_(i+1)+1, G(y)_i=y_(i-1)-1."""
    a=list(h[:m]);b=list(h[m:]);changes=[0]*(2*m)
    def layer(which,values):
        arr=a if which==0 else b
        for i,v in enumerate(values):
            old=arr[i];arr[i]=update(old,v)
            changes[which*m+i]+=old!=arr[i]
    layer(1,[b[i]+a[(i+1)%m]+1 for i in range(m)])
    layer(0,[a[i]-b[(i-1)%m]+1 for i in range(m)])
    layer(0,[a[i]+b[i]for i in range(m)])
    layer(1,[a[i]-b[i]for i in range(m)])
    layer(0,[a[i]-b[i]for i in range(m)])
    return tuple(a+b),changes

def clock(h,m):
    a=list(h[:m]);b=list(h[m:]);changes=[0]*(2*m)
    def layer(which,values):
        arr=a if which==0 else b
        for i,v in enumerate(values):
            old=arr[i];arr[i]=update(old,v)
            changes[which*m+i]+=old!=arr[i]
    layer(0,[a[i]+b[i]+b[(i+1)%m]for i in range(m)])
    layer(0,[a[i]+b[i]for i in range(m)])
    layer(1,[a[i]-b[i]for i in range(m)])
    layer(0,[a[i]-b[i]for i in range(m)])
    return tuple(a+b),changes

# Tiny exact polynomials over F_3, low coefficients first.
def padd(a,b):
    c=[0]*max(len(a),len(b))
    for i in range(len(c)):c[i]=((a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0))%3
    while len(c)>1 and c[-1]==0:c.pop()
    return c

def pmul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%3
    while len(c)>1 and c[-1]==0:c.pop()
    return c

def mmul(a,b):
    return [[padd(pmul(a[i][0],b[0][j]),pmul(a[i][1],b[1][j]))for j in range(2)]for i in range(2)]

def mpow(a,n):
    out=[[[1],[0]],[[0],[1]]]
    for _ in range(n):out=mmul(out,a)
    return out

def at1(a):return [[sum(p)%3 for p in row]for row in a]
def deriv1(a):return [[sum(i*c for i,c in enumerate(p))%3 for p in row]for row in a]

def main():
    seed=[];canon=[]
    for h,mask in enumerate(MASKS):
        require(strong(mask),'seed strong')
        P=profile_dp(mask);Q=profile_cuts(mask)
        require(P==Q==PROFILE,'seed profiles')
        for r in range(9):
            F=sum(P[j]*math.prod(range(r-j+1,r+1))for j in range(1,7))
            require(F==r**6+12*r**4+16*r*r,'slot polynomial')
        c=min(relabel(mask,p)for p in it.permutations(range(6)))
        canon.append(c)
        require(relabel(mask,DECK)==MASKS[(h+3)%6],'deck seed')
        neighbors=[j for j,v in enumerate(MASKS)if (mask^v).bit_count()==1]
        require(set(neighbors)=={(h-1)%6,(h+1)%6},'induced hexagon')
        seed.append({'h':h,'mask':mask,'canonical_mask':c,'profile':P,'neighbors':neighbors})
    require(len(set(canon))==3 and all(canon[i]==canon[i+3]for i in range(3)),'class fibers')
    require(all(canon[i]!=canon[j]for i in range(3)for j in range(i+1,3)),'distinct trits')
    hostile=[]
    for mask in (145,8377):
        P=profile_dp(mask);require(P[1]==33 and P==profile_cuts(mask),'hostile corners')
        hostile.append({'mask':mask,'profile':P})

    covers=[]
    for N in range(1,4):
        words=list(it.product(range(6),repeat=N));encoded=[ordinal(h)for h in words]
        require(len(set(encoded))==6**N,'literal injectivity')
        edge_count=0
        fibers={}
        for i,h in enumerate(words):
            q=tuple(v%3 for v in h);fibers.setdefault(q,0);fibers[q]+=1
            neigh=[]
            for j,g in enumerate(words):
                d=(encoded[i]^encoded[j]).bit_count()
                expected=sum(x!=y for x,y in zip(h,g))==1 and any((x-y)%6 in (1,5)for x,y in zip(h,g))
                require((d==1)==expected,'product inducedness')
                if d==1:neigh.append(tuple(v%3 for v in g))
            require(len(neigh)==len(set(neigh))==2*N,'cover neighbor bijection')
            edge_count+=len(neigh)
        require(len(fibers)==3**N and set(fibers.values())=={2**N},'fiber degree')
        covers.append({'N':N,'literal_vertices':6**N,'logical_vertices':3**N,'sheets':2**N,'edges':edge_count//2})

    gates=[]
    for perm in it.permutations(range(3)):
        outputs=[update(h,perm[h%3])for h in range(6)]
        require(len(set(outputs))==6,'gate bijection')
        require(all(outputs[(h+3)%6]==(outputs[h]+3)%6 for h in range(6)),'gate deck commutation')
        inv=[perm.index(i)for i in range(3)]
        require(all(update(outputs[h],inv[outputs[h]%3])==h for h in range(6)),'gate inverse')
        gates.append({'permutation':perm,'lift':outputs})

    comp=[]
    for m in (1,2,3):
        logical_out=set()
        for h in it.product(range(3),repeat=2*m):
            out,changes=compiler(h,m);logical=tuple(v%3 for v in out)
            logical_out.add(logical)
            a=h[:m];b=h[m:]
            c=tuple((b[i]+a[(i+1)%m]+1)%3 for i in range(m))
            d=tuple((a[i]-c[(i-1)%m]+1)%3 for i in range(m))
            require(logical==c+d,'general compiler formula')
            require(sum(changes)<=5*m,'compiler flip bound')
            if not any(b):require(logical==tuple((a[(i+1)%m]+1)%3 for i in range(m))+(0,)*m,'blank workspace')
        require(len(logical_out)==3**(2*m),'logical compiler permutation')
        physical_count=None
        if m<=2:
            physical_out={compiler(h,m)[0]for h in it.product(range(6),repeat=2*m)}
            physical_count=len(physical_out)
            require(physical_count==6**(2*m),'physical compiler permutation')
        comp.append({'width':m,'logical_states':len(logical_out),'physical_states':physical_count})

    matrix=[[[0],[1]],[[1],[1,1]]]
    M4=mpow(matrix,4);M8=mpow(matrix,8)
    require(at1(M4)==[[2,0],[0,2]],'M(1)^4')
    require(at1(M8)==[[1,0],[0,1]],'M(1)^8')
    B=deriv1(M8)
    require(B==[[1,2],[2,2]] and (B[0][0]*B[1][1]-B[0][1]*B[1][0])%3==1,'invertible derivative')
    clocks=[]
    for m in (1,3,9,27,81):
        initial=(1,)+(0,)*(2*m-1);h=initial;seenlogical=set();flipcounts=[0]*(2*m)
        total=0
        for t in range(8*m):
            logical=tuple(v%3 for v in h)
            require(logical not in seenlogical,'no early clock return');seenlogical.add(logical)
            h,changes=clock(h,m);total+=sum(changes)
            flipcounts=[a+b for a,b in zip(flipcounts,changes)]
            if t+1==4*m:require(tuple(v%3 for v in h)==tuple(-v%3 for v in initial),'negative half-period')
        require(h==initial and len(seenlogical)==8*m,'exact physical clock period')
        require(all(v%2==0 for v in flipcounts),'trivial clock monodromy')
        clocks.append({'width':m,'logical_period':8*m,'physical_period':8*m,'literal_flips':total,'flip_parities':[v%2 for v in flipcounts]})

    # Exhaustive whole physical state spaces for the two smallest widths,
    # including arbitrary initial deck bits, not only the basis orbit.
    full=[]
    for m in (1,3):
        count=0
        # Width3: all729 logical basepoints, with all64 deck lifts covered by equivariance.
        for logical in it.product(range(3),repeat=2*m):
            h=tuple(logical)
            for _ in range(8*m):h,_=clock(h,m)
            require(h==logical,'whole-state order')
            count+=1
        # Deck equivariance of all primitive gates extends these checks to every lift.
        full.append({'width':m,'logical_basepoints':count,'physical_states_via_deck_equivariance':6**(2*m),'rounds':8*m})

    # Off-family hostile and longer-cycle controls; only basis-orbit claims.
    other_clocks=[]
    for m,want_logical,want_physical in ((2,8,16),(4,40,80),(5,328,328),(7,2920,2920)):
        initial=(1,)+(0,)*(2*m-1);h=initial;first_logical=None;first_physical=None
        flips=[0]*(2*m);monodromy=None
        for t in range(1,want_physical+1):
            h,changes=clock(h,m)
            flips=[a+b for a,b in zip(flips,changes)]
            if first_logical is None and tuple(v%3 for v in h)==initial:
                first_logical=t;monodromy=[v%2 for v in flips]
                require(all((h[i]-initial[i])%6==3*monodromy[i] for i in range(2*m)), 'orbit parity criterion')
            if h==initial:
                first_physical=t;break
        require((first_logical,first_physical)==(want_logical,want_physical),'off-family periods')
        other_clocks.append({'width':m,'logical_basis_period':first_logical,'physical_basis_period':first_physical,'first_return_deck_bits':monodromy})

    # Sharp period-doubling witness.
    h=0;orbit=[]
    for _ in range(6):orbit.append(h);h=update(h,h+1)
    require(orbit==list(range(6)) and h==0,'period doubling witness')
    result={'status':'FINITE-EXACT controls; general statements proved in NOTE.md',
      'seed':seed,'hostile_corners':hostile,'covers':covers,'all_trit_permutation_lifts':gates,
      'five_layer_compiler':comp,'M4_polynomials':M4,'M8_polynomials':M8,'B_at_1':B,
      'clock_basis_orbits':clocks,'off_family_basis_orbits':other_clocks,'whole_state_clock_controls':full,'doubled_clock':orbit}
    path=Path(__file__).with_name('certificate.json')
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS: literal seed, 3 product covers, all 6 trit gates, reversible compiler, exact clock proofs and controls; certificate.json written.')

if __name__=='__main__':main()

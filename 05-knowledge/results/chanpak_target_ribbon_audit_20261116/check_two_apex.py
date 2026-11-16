#!/usr/bin/env python3
"""Standalone exact controls. No dependencies and no assertions."""
import hashlib
import itertools
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def words(m):
    for bits in itertools.product('01', repeat=m-2):
        yield '1'+''.join(bits)+'0'


def stats(w):
    h,a,p=1,0,0
    for b in map(int,w):
        if b:
            a,p=2*a+1,1
        else:
            h,a,p=h+a+p,2*a,0
    b=sum((c=='0')*(1<<i) for i,c in enumerate(w))
    return h,a,b


def word_rows(w):
    m=len(w)
    r=[sum(1<<j for j in range(i+1,m)) for i in range(m)]+[0]
    for i,b in enumerate(w):
        if b=='1':r[m]|=1<<i
        else:r[i]|=1<<m
    return tuple(r)


def ordinal(*blocks):
    n=sum(map(len,blocks));out=[];off=0
    for block in blocks:
        future=((1<<n)-1)^((1<<(off+len(block)))-1)
        out.extend((r<<off)|future for r in block)
        off+=len(block)
    return tuple(out)


def bridge(u,v,k=0):
    chain=tuple(sum(1<<j for j in range(i+1,k)) for i in range(k))
    r=list(ordinal(word_rows(u),chain,word_rows(v)))
    x=len(u);y=len(r)-1
    require(r[x]>>y&1,'original bridge arc')
    r[x]^=1<<y;r[y]|=1<<x
    return tuple(r),x,y


def hp(r):
    n=len(r);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,weight in enumerate(dp[mask]):
            if not weight:continue
            avail=r[i]&~mask
            while avail:
                bit=avail&-avail;avail-=bit
                dp[mask|bit][bit.bit_length()-1]+=weight
    return sum(dp[-1])


def cycles(r):
    """All directed odd simple cycles, one canonical cyclic rotation."""
    out=[];n=len(r)
    for start in range(n):
        def visit(last,mask,length):
            if length>=3 and length%2 and (r[last]>>start&1):
                out.append(mask)
            avail=r[last]&~mask&~((1<<(start+1))-1)
            while avail:
                bit=avail&-avail;avail-=bit
                visit(bit.bit_length()-1,mask|bit,length+1)
        visit(start,1<<start,1)
    return out


def strong(r):
    n=len(r);reach=[r[i]|1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return all(x==(1<<n)-1 for x in reach)


def is_transitive_on(r,mask):
    verts=[i for i in range(len(r)) if mask>>i&1]
    return sorted((r[i]&mask).bit_count() for i in verts)==list(range(len(verts)))


def polynomial_from_cycles(cc):
    # Called only with a known two-vertex feedback set, hence degree <= 2.
    d=sum(not (a&b) for i,a in enumerate(cc) for b in cc[i+1:])
    return [1,len(cc),d]


def two_apex_rows(w):
    n=len(w)+2;x=n-2;y=n-1
    r=[sum(1<<j for j in range(i+1,n-2)) for i in range(n-2)]+[1<<y,0]
    for i,t in enumerate(w):
        for j,k in ((x,1),(y,2)):
            if t&k:r[j]|=1<<i
            else:r[i]|=1<<j
    return tuple(r),x,y


def finite_compiler(cap):
    L=cap.bit_length()-1
    bank={m:[(w,*stats(w)) for w in words(m)] for m in range(2,L+1)}
    witness={};presentations=0
    for m in bank:
        for n in bank:
            if m+n>L+2:continue
            for u,h,A,_ in bank[m]:
                for v,j,_,B in bank[n]:
                    F=h*j;G=(A+1)*(B+1);N=F+G-1
                    if N<=cap:
                        witness.setdefault(N,[u,v,0]);presentations+=1
                    k=1
                    while F+(G<<k)<=cap:
                        witness.setdefault(F+(G<<k),[u,v,k]);presentations+=1;k+=1
    one={1}
    for m in range(2,2+(cap-1).bit_length()):
        for w in words(m):
            h,_,_=stats(w)
            if h<=cap:one.add(h)
    repairs=sorted(set(witness)-one)
    return witness,repairs,len(one),presentations


def main():
    word_checks=0
    for m in range(2,11):
        for w in words(m):
            h,A,B=stats(w)
            require(h<=A+1<=2*h-3 and h<=B+1<=2*h-3,'port bound '+w)
            require((A+1==h)==(w=='1'+'0'*(m-1)),'left equality '+w)
            require((A+1==2*h-3)==(w=='1'*(m-1)+'0'),'right equality '+w)
            require((B+1==h)==(w=='1'*(m-1)+'0'),'dual left equality '+w)
            require((B+1==2*h-3)==(w=='1'+'0'*(m-1)),'dual right equality '+w)
            require(h>=2**(m-2)+1,'one apex length cutoff')
            if m<=7:require(hp(word_rows(w))==h,'literal one apex DP')
            word_checks+=1
    bridge_checks=0
    for m in range(2,6):
        for n in range(2,6):
            if m+n>7:continue
            for u in words(m):
                for v in words(n):
                    h,A,_=stats(u);j,_,B=stats(v)
                    a=(h-1)//2;b=(j-1)//2;P=A+1;Q=B+1
                    vals=[]
                    for k in range(3):
                        r,x,y=bridge(u,v,k);N=h*j+P*Q-1 if k==0 else h*j+(P*Q<<k)
                        cc=cycles(r);poly=polynomial_from_cycles(cc)
                        c=(P*Q-1)//2 if k==0 else (P*Q<<(k-1))
                        require(poly==[1,a+b+c,a*b], 'literal polynomial')
                        xc=[t for t in cc if t>>x&1 and not(t>>y&1)]
                        yc=[t for t in cc if t>>y&1 and not(t>>x&1)]
                        both=[t for t in cc if t>>x&1 and t>>y&1]
                        require(len(xc)==a and len(yc)==b and len(both)==c,'literal clique counts')
                        require(all(not(s&t) for s in xc for t in yc),'all cross old cycles disjoint')
                        require(hp(r)==N==1+2*poly[1]+4*poly[2],'literal bridge Hamilton DP')
                        require(strong(r),'bridge strongness')
                        full=(1<<len(r))-1
                        require(is_transitive_on(r,full^(1<<x)^(1<<y)),'feedback pair')
                        require(all(not is_transitive_on(r,full^(1<<i)) for i in range(len(r))),'not one apex')
                        vals.append(N);bridge_checks+=1
                    require(vals[1]==2*vals[0]-h*j+2 and vals[2]==2*vals[1]-h*j,'padding discontinuity')
    general_checks=0
    for m in range(5):
        for w in itertools.product(range(4),repeat=m):
            r,x,y=two_apex_rows(w);cc=cycles(r);poly=polynomial_from_cycles(cc)
            a=sum(bool(t>>x&1) and not(t>>y&1) for t in cc)
            b=sum(bool(t>>y&1) and not(t>>x&1) for t in cc)
            c=sum(bool(t>>x&1) and bool(t>>y&1) for t in cc)
            require(a+b+c==poly[1] and 0<=poly[2]<=a*b,'generic two apex coefficients')
            require(poly[1]**2-4*poly[2]>=0,'quadratic real roots')
            require(hp(r)==1+2*poly[1]+4*poly[2],'generic OCF direct DP')
            general_checks+=1
    cap=10000
    witness,repairs,n_one,presentations=finite_compiler(cap)
    examples=[]
    for u,v,N in [('10','110',29),('10','100',35),('100','110',49),('10','1010',59),('10','1000',71)]:
        r,x,y=bridge(u,v);poly=polynomial_from_cycles(cycles(r))
        require(hp(r)==N and N in repairs,'promised exact one-apex repair')
        examples.append({'u':u,'v':v,'H':N,'polynomial':poly,'rows':r})
    require(7 not in witness and 13 not in witness and 21 not in witness,'bridge lower exclusions')
    # Check every compiler witness from its literal word data and formula, not saved profile numbers.
    for N,(u,v,k) in witness.items():
        h,A,_=stats(u);j,_,B=stats(v)
        actual=h*j+(A+1)*(B+1)-1 if not k else h*j+((A+1)*(B+1)<<k)
        require(actual==N,'compiler round trip')
    cert={'cap':cap,'witnesses':{str(N):witness[N] for N in sorted(witness)},'one_apex_absent':repairs}
    encoded=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()
    with open('finite_certificate.json','wb') as f:f.write(encoded+b'\n')
    print(json.dumps({'status':'PASS','word_bound_checks':word_checks,'literal_bridge_controls':bridge_checks,
        'generic_two_apex_controls':general_checks,'cap':cap,'bridge_values':len(witness),
        'bridge_presentations':presentations,'one_apex_values':n_one,'exact_new_to_one_apex':len(repairs),
        'repairs_at_most_100':[n for n in repairs if n<=100],
        'finite_certificate_sha256':hashlib.sha256(encoded+b'\n').hexdigest(),
        'examples':examples},sort_keys=True,indent=2))

if __name__=='__main__':main()

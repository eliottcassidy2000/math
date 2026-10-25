#!/usr/bin/env python3
"""Independent exact finite-field Parker census; standard library only."""
import argparse, hashlib, itertools, json, math, pathlib, time

class Field:
    def __init__(self,p,n=1,modulus=None):
        self.p,self.n,self.q=p,n,p**n
        self.modulus=modulus or (self.find_modulus() if n>1 else [0,1])
        self.digits=[self.decode(i) for i in range(self.q)] if n>1 else None
        self.negatives=[self.neg(i) for i in range(self.q)]
        self.squares={self.mul(i,i) for i in range(self.q)}
    def decode(self,a):
        out=[]
        for _ in range(self.n): out.append(a%self.p); a//=self.p
        return out
    def encode(self,a):
        out=0
        for x in reversed(a):out=out*self.p+x%self.p
        return out
    def add(self,a,b):
        if self.n==1:return (a+b)%self.p
        return self.encode([x+y for x,y in zip(self.digits[a],self.digits[b])])
    def neg(self,a):
        if self.n==1:return -a%self.p
        return self.encode([-x for x in self.digits[a]])
    def sub(self,a,b):return self.add(a,self.negatives[b])
    def mul(self,a,b):
        if self.n==1:return a*b%self.p
        z=[0]*(2*self.n-1)
        for i,x in enumerate(self.digits[a]):
            for j,y in enumerate(self.digits[b]):z[i+j]=(z[i+j]+x*y)%self.p
        for k in range(len(z)-1,self.n-1,-1):
            for j in range(self.n):z[k-self.n+j]=(z[k-self.n+j]-z[k]*self.modulus[j])%self.p
        return self.encode(z[:self.n])
    def find_modulus(self):
        # Trial division by EVERY monic polynomial of degree <= floor(n/2).
        # A reducible monic polynomial has such a factor, so this is complete.
        divisors=[]
        for d in range(1,self.n//2+1):
            divisors.extend([list(c)+[1] for c in itertools.product(range(self.p),repeat=d)])
        for c in itertools.product(range(self.p),repeat=self.n):
            if not c[0]:continue
            f=list(c)+[1]
            if all(any(poly_rem(f,g,self.p)) for g in divisors):return f
        raise AssertionError('no irreducible polynomial')

def poly_rem(f,g,p):
    r=f[:]
    for k in range(len(f)-1,len(g)-2,-1):
        for j in range(len(g)-1):r[k-len(g)+1+j]=(r[k-len(g)+1+j]-r[k]*g[j])%p
        r[k]=0
    return r[:len(g)-1]

def magic(F,c,u,v):
    a,s=F.add,F.sub
    return [a(c,u),s(s(c,u),v),a(c,v),a(s(c,u),v),c,s(a(c,u),v),s(c,v),a(a(c,u),v),s(c,u)]

def valid(F,m):
    if len(set(m))!=9 or not set(m)<=F.squares:return False
    lines=[[0,1,2],[3,4,5],[6,7,8],[0,3,6],[1,4,7],[2,5,8],[0,4,8],[2,4,6]]
    sums=[F.add(F.add(m[a],m[b]),m[c]) for a,b,c in lines]
    return len(set(sums))==1

def normalized(F,full=False):
    # ALL labeled center-one squares, and center-zero squares with u=1.
    # Square scaling maps any nonzero center to 1; if center=0, u is a
    # nonzero square and its inverse scales u to 1.
    count0=count1=0; witness=None
    for v in (range(F.q) if F.negatives[1] in F.squares else []):
        m=magic(F,0,1,v)
        if valid(F,m):
            count0+=1;witness=witness or m
            if not full:return {'witness':m,'method':'center_zero'}
    C=[u for u in range(F.q) if F.add(1,u) in F.squares and F.sub(1,u) in F.squares]
    for v in C:
        for u in C:
            m=magic(F,1,u,v)
            if valid(F,m):
                count1+=1;witness=witness or m
                if not full:return {'witness':m,'method':'center_one'}
    return {'witness':witness,'method':'exhaustive_normalized','center_zero_u1_count':count0,'center_one_count':count1}

def prime_powers(limit):
    primes=[p for p in range(3,limit,2) if all(p%d for d in range(3,math.isqrt(p)+1,2))]
    rows=[]
    for p in primes:
        q=p;n=1
        while q<limit:rows.append((q,p,n));q*=p;n+=1
    return sorted(rows)

def census(limit,full_limit):
    rows=[]
    for q,p,n in prime_powers(limit):
        F=Field(p,n)
        r={'q':q,'p':p,'n':n,'modulus_low_first':F.modulus}
        r.update(normalized(F,full=q<=full_limit))
        if r['witness'] is not None:assert valid(F,r['witness'])
        rows.append(r)
    return {'schema':'finite-field-parker-v1','universe':f'all odd prime powers q < {limit}', 'limit_exclusive':limit,'full_counts_through':full_limit,'parker_orders':[r['q'] for r in rows if r['witness'] is None],'rows':rows}

def check_certificate(path):
    data=json.loads(pathlib.Path(path).read_text()); rows=data['rows']
    assert [(r['q'],r['p'],r['n']) for r in rows]==prime_powers(data['limit_exclusive'])
    for r in rows:
        F=Field(r['p'],r['n'])
        assert F.modulus==r['modulus_low_first']
        if r['q']<=data['full_counts_through'] or r['witness'] is None:
            exact=normalized(F,full=True)
            for key in ('center_zero_u1_count','center_one_count'):
                assert exact[key]==r[key],(r['q'],key)
            assert (exact['witness'] is None)==(r['witness'] is None)
        if r['witness'] is not None:assert valid(F,r['witness'])
    assert data['parker_orders']==[r['q'] for r in rows if r['witness'] is None]
    print('PASS: complete universe, irreducible field moduli, direct eight-line witnesses, exhaustive negative searches and reported bounded counts')

def controls():
    # Direct brute triple path independently tests normalization and orbit mass.
    for p in [3,5,7,11,17,29]:
        F=Field(p);r=normalized(F,True)
        brute=sum(valid(F,magic(F,c,u,v)) for c in F.squares for u in range(p) for v in range(p))
        expected=(p-1)//2*(r['center_zero_u1_count']+r['center_one_count'])
        assert brute==expected,(p,brute,expected)
    # Finite field 27 differs from residue ring 27; zero divisors do not occur here.
    F=Field(3,3)
    assert len(F.squares)==14 and all(F.mul(a,b) for a in range(1,27) for b in range(1,27))
    assert normalized(F,True)['witness'] is None
    print('PASS: brute triple controls, normalization mass, F27 field and Parker control')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int,default=14400);ap.add_argument('--full-counts-through',type=int,default=300);ap.add_argument('--output',default='census.json');ap.add_argument('--verify');ap.add_argument('--controls',action='store_true');args=ap.parse_args()
    if args.controls:controls()
    elif args.verify:check_certificate(args.verify)
    else:
        t=time.time();data=census(args.limit,args.full_counts_through);pathlib.Path(args.output).write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'fields':len(data['rows']),'parker':data['parker_orders'],'seconds':time.time()-t,'sha256':hashlib.sha256(pathlib.Path(args.output).read_bytes()).hexdigest()}))

#!/usr/bin/env python3
"""Canonical #P witnesses for tournament determinant gaps; no third-party imports.

Run normally for exhaustive n<=5 checks and adversarial rejection controls.
Verifier inputs are (adjacency matrix, mode, fixed-length bit string).
The enumerator is diagnostic only: verification itself is polynomial-time.
"""
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial
from hashlib import sha256
from pathlib import Path

DELTA, E, MARGIN, GAMMA = range(4)
SKIP, RECURSE, LEAF = range(3)
NAMES = ('delta', 'e', 'H-P', 'Gamma')

def require(condition, detail='check failed'):
    if not condition:raise AssertionError(detail)


def determinant(a):
    a=[list(r) for r in a]; n=len(a)
    if not n:return 1
    sign=1; previous=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if a[i][k]),None)
        if p is None:return 0
        if p!=k:a[k],a[p]=a[p],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                quotient,remainder=divmod(numerator,previous)
                if remainder:raise ArithmeticError('Bareiss division was not exact')
                a[i][j]=quotient
        for i in range(k+1,n):a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]


def valid_tournament(a):
    n=len(a)
    return (all(len(r)==n for r in a)
            and all(a[i][i]==0 for i in range(n))
            and all(a[i][j] in (0,1) and a[j][i] in (0,1)
                    and a[i][j]+a[j][i]==1 for i,j in combinations(range(n),2)))


class GapSystem:
    def __init__(self,a):
        self.a=tuple(tuple(r) for r in a); self.n=len(a)
        self.good=valid_tournament(a)
        n=self.n
        self.V=tuple(range(n))
        self.L=max(1,n.bit_length()) # ceil(log2(n+1)), also one bit at empty
        self.B=(2*factorial(n)).bit_length() # ceil(log2(2 n!+1))
        self.widths=(2,2,self.L)+ (self.L,)*n +(1,self.B)
        self.record_bits=sum(self.widths)
        self.zero=(0,)*(n+5)
        self.uv_cache={}
        if self.good:
            self.epsilon=-1 if self.uv(self.V)[1]>0 else 1
        else:self.epsilon=1

    def uv(self,U):
        if U not in self.uv_cache:
            a=self.a;k=len(U)
            u=determinant([[a[i][j]+int(i==j) for j in U] for i in U])
            v=(-1)**k*determinant([[a[i][j] for j in U] for i in U])
            self.uv_cache[U]=(u,v)
        return self.uv_cache[U]

    def leaf_bound(self,mode,R):
        u,v=self.uv(R)
        if mode==DELTA:return u
        if mode==E:return u-1
        if mode==MARGIN:return u-self.epsilon*v
        if mode==GAMMA:return u-self.epsilon*v-int(bool(R))
        raise ValueError('unknown mode')

    def record(self,mode,tag,P=(),color=0,label=0):
        return (mode,tag,len(P))+tuple(i+1 for i in P)+(0,)*(self.n-len(P))+(color,label)

    def pack_records(self,records):
        if len(records)!=self.n:raise ValueError('record count')
        pieces=[]
        for row in records:
            if len(row)!=len(self.widths):raise ValueError('record shape')
            for value,width in zip(row,self.widths):
                if not isinstance(value,int) or not 0<=value<1<<width:
                    raise ValueError('out-of-range field')
                pieces.append(format(value,'0%db'%width))
        return ''.join(pieces)

    def pack_trace(self,trace):
        if len(trace)>self.n:raise ValueError('trace too long')
        return self.pack_records(tuple(trace)+(self.zero,)*(self.n-len(trace)))

    def unpack(self,witness):
        if (not isinstance(witness,str) or len(witness)!=self.n*self.record_bits
            or any(c not in '01' for c in witness)):return None
        rows=[];p=0
        for _ in range(self.n):
            row=[]
            for width in self.widths:
                row.append(int(witness[p:p+width],2));p+=width
            rows.append(tuple(row))
        return tuple(rows)

    def verify(self,initial_mode,witness):
        """Polynomial-time verifier. All four functions are zero at empty/invalid."""
        if not self.good or not self.n or initial_mode not in range(4):return False
        rows=self.unpack(witness)
        if rows is None:return False
        U=self.V;mode=initial_mode
        for pos,row in enumerate(rows):
            if not U:return False
            rm,tag,length=row[:3]
            slots=row[3:3+self.n];color,label=row[-2:]
            if rm!=mode:return False
            if tag==SKIP:
                if length or any(slots) or color or label:return False
                U=U[1:]
                mode={MARGIN:DELTA,GAMMA:E}.get(mode,mode)
                continue
            if tag not in (RECURSE,LEAF):return False
            if not 2<=length<=len(U)-1:return False
            if any(slots[length:]):return False
            P=tuple(z-1 for z in slots[:length]);x=U[0]
            if len(set(P))!=length or any(i not in U or i==x for i in P):return False
            C=(x,)+P
            if any(not self.a[i][j] for i,j in zip(C,C[1:]+C[:1])):return False
            Cset=set(C);R=tuple(i for i in U if i not in Cset)
            if tag==RECURSE:
                if len(C)%2!=1 or label:return False
                # color may be 0 or 1: the two genuine copies in 2 delta(R).
                U=R;mode=DELTA
                continue
            if color:return False
            if any(other!=self.zero for other in rows[pos+1:]):return False
            bound=self.leaf_bound(mode,R)
            return 0<=label<bound
        return False # a skip/recurse chain without a leaf is never accepting

    def traces(self,initial_mode):
        """Generate each accepting trace once; exponential diagnostic counter."""
        if not self.good or not self.n or initial_mode not in range(4):return ()
        @lru_cache(None)
        def walk(mode,U):
            if not U:return ()
            x=U[0];out=[]
            next_mode={MARGIN:DELTA,GAMMA:E}.get(mode,mode)
            for tail in walk(next_mode,U[1:]):
                out.append((self.record(mode,SKIP),)+tail)
            for length in range(2,len(U)):
                for P in permutations(U[1:],length):
                    C=(x,)+P
                    if any(not self.a[i][j] for i,j in zip(C,C[1:]+C[:1])):continue
                    Cset=set(C);R=tuple(i for i in U if i not in Cset)
                    if len(C)%2:
                        for color in (0,1):
                            for tail in walk(DELTA,R):
                                out.append((self.record(mode,RECURSE,P,color),)+tail)
                    bound=self.leaf_bound(mode,R)
                    if bound<0:raise AssertionError(('negative bound',mode,U,P,bound))
                    for label in range(bound):out.append((self.record(mode,LEAF,P,0,label),))
            return tuple(out)
        return walk(initial_mode,self.V)

    def certificates(self,mode):
        return tuple(self.pack_trace(trace) for trace in self.traces(mode))


def literal_expected(a):
    """Independent literal words, minimum-rooted cycles, direct determinants."""
    n=len(a)
    if not n:return (0,0,0,0)
    V=tuple(range(n));H=0
    for p in permutations(V):
        H+=all(a[i][j] for i,j in zip(p,p[1:]))
    call=cproper=0
    for k in range(3,n+1):
        for S in combinations(V,k):
            for tail in permutations(S[1:]):
                C=(S[0],)+tail
                if all(a[i][j] for i,j in zip(C,C[1:]+C[:1])):
                    call+=1;cproper+=int(k<n)
    u=determinant([[a[i][j]+int(i==j) for j in V] for i in V])
    v=(-1)**n*determinant(a)
    delta=H-u
    return delta,delta-call,delta-abs(v),delta-abs(v)-cproper


def tournament(n,mask):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(mask>>k)&1;a[j][i]=1-a[i][j]
    return a


def check_one(a):
    S=GapSystem(a);expect=literal_expected(a);total=0
    for mode in range(4):
        certs=S.certificates(mode)
        require(len(certs)==len(set(certs)),('duplicate certificate',mode,a))
        require(len(certs)==expect[mode],('wrong count',mode,a,len(certs),expect))
        require(all(S.verify(mode,w) for w in certs),('valid certificate rejected',mode,a))
        total+=len(certs)
    return total


def rejection_controls():
    # C3 => C3: all four functions positive, so all four verifier modes tested.
    a=tournament(6,(1<<15)-1)
    for x,y in [(2,0),(5,3)]:a[x][y]=1;a[y][x]=0
    require(literal_expected(a)==(5,3,4,2))
    require(check_one(a)==14) # includes a nonzero Gamma=2 and colored recursion
    S=GapSystem(a);controls=0
    for mode in range(4):
        trace=S.traces(mode)[0];w=S.pack_trace(trace)
        require(S.verify(mode,w))
        rows=[list(r) for r in S.unpack(w)]
        leaf=next(i for i,r in enumerate(rows) if r[1]==LEAF)
        cyc=next(i for i,r in enumerate(rows) if r[1] in (LEAF,RECURSE))
        def reject(mutator):
            nonlocal controls
            changed=[r[:] for r in rows];mutator(changed)
            require(not S.verify(mode,S.pack_records(changed)),('bad certificate accepted',mode))
            controls+=1
        reject(lambda r:r[0].__setitem__(0,r[0][0]^1)) # wrong active mode
        reject(lambda r:r[cyc].__setitem__(4,r[cyc][3])) # repeated path vertex
        reject(lambda r:r[leaf].__setitem__(-1,(1<<S.B)-1)) # out-of-range label
        reject(lambda r:r[leaf].__setitem__(-2,1)) # inactive color ambiguity
        reject(lambda r:r[cyc].__setitem__(3+r[cyc][2],1)) # nonzero path padding
        reject(lambda r:r[-1].__setitem__(-1,1)) # data after terminal leaf
        reject(lambda r:r[cyc].__setitem__(1,3)) # reserved action tag
        require(not S.verify(mode,w+'0'));controls+=1
        require(not S.verify(mode,w[:-1]));controls+=1
        require(not S.verify(mode,'x'+w[1:]));controls+=1
        # Same certificate cannot be relabelled as another initial function.
        require(not S.verify(mode^1,w));controls+=1
    for bad in [[[0,1],[1,0]],[[1]],[[0],[0,0]]]:
        Sbad=GapSystem(bad)
        require(not Sbad.verify(DELTA,'0'*(Sbad.n*Sbad.record_bits)));controls+=1
    require(not GapSystem([]).verify(DELTA,''));controls+=1
    return controls


def main():
    inputs=accepted=0
    for n in range(6):
        for mask in range(1<<(n*(n-1)//2)):
            accepted+=check_one(tournament(n,mask));inputs+=1
        print('exhaustive n=%d passed; inputs=%d accepting certificates=%d'%(n,inputs,accepted),flush=True)
    # A source before C3=>C3 forces Gamma's skip-to-e transition with e>0.
    a=tournament(7,(1<<21)-1)
    for x,y in [(3,1),(6,4)]:a[x][y]=1;a[y][x]=0
    require(literal_expected(a)==(5,3,5,3))
    require(check_one(a)==16)
    print('source=>C3=>C3 passed; gaps=(5,3,5,3), certificates=16',flush=True)
    print('adversarial rejection controls passed:',rejection_controls(),flush=True)
    print('PASS exact canonical #P certificate counts for delta, e, H-P, Gamma',flush=True)
    print('script_sha256',sha256(Path(__file__).read_bytes()).hexdigest(),flush=True)

if __name__=='__main__':main()

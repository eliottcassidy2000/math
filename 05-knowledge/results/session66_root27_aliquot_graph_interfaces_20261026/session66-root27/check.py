"""Exact source-fixed shape-switch cocycle and hostile controls; standard library only."""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb
from collections import Counter
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def val(n, p=2):
    n=abs(n)
    require(n != 0, 'finite valuation input')
    k=0
    while n % p == 0:
        n//=p; k+=1
    return k


def qval(x):
    return val(x.numerator)-val(x.denominator)


def carry(w):
    c=0
    for i,b in enumerate(w):
        c=3*c+int(b)*(1<<i)
    return c


def direct_carry(w):
    return sum(int(b)*(1<<i)*3**(len(w)-1-i) for i,b in enumerate(w))


def param(w):
    n=len(w); a=3**n; b=1<<n
    return n,a,b,a-b,carry(w)


def fuel(w,z):
    n,a,b,d,c=param(w)
    return val(d*z+c)


def center(w):
    n,a,b,d,c=param(w)
    return -Fraction(c,d)


def transform(w,z):
    n,a,b,d,c=param(w)
    return Fraction(a*z+c,b)


def literal(z,n,q=1):
    w=''; values=[z]; credits=[q]
    for _ in range(n):
        bit=z%2
        require(q>0, 'invalid native credit')
        require(not (q==1 and bit==0), 'attempted blocked native macro')
        w+=str(bit); q+=3*bit-1; z=(3*z+bit)//2
        values.append(z); credits.append(q)
    return w,values,credits


def residue(w):
    mod=1<<len(w)
    return (-carry(w)*pow(3**len(w),-1,mod))%mod


def balanced(w):
    q=0
    for b in w:
        q+=3*int(b)-1
        if q<0:return False
    return q==0


@lru_cache(None)
def trees(e):
    if e==0:return ('0',)
    return tuple('1'+a+b+c for i in range(e) for j in range(e-i)
                 for a in trees(i) for b in trees(j) for c in trees(e-1-i-j))


def primitive(w):
    for n in range(1,len(w)+1):
        if len(w)%n==0 and w[:n]*(len(w)//n)==w:return w[:n]
    raise RuntimeError('primitive root missing')


def main():
    stats=Counter(); blocks=[]
    for e in range(1,5):
        ts=trees(e)
        require(len(ts)==comb(3*e,e)//(2*e+1),'complete tree count')
        blocks.extend(w[:-1] for w in ts)
    require(len(blocks)==71 and all(balanced(w) for w in blocks),'balanced bank')
    # All binary words through length9, including the carry-cylinder injection.
    for n in range(10):
        residues=[]; carries=[]
        for tup in product('01',repeat=n):
            w=''.join(tup);c=carry(w)
            require(c==direct_carry(w),'independent carry formula')
            residues.append(residue(w));carries.append(c)
            stats['binary_words']+=1
        require(len(set(residues))==len(set(carries))==1<<n,'word injection')
    # All pairs in the bounded balanced-tree-prefix bank.
    common_centers=0
    for p in blocks:
        lp,ap,bp,dp,cp=param(p)
        for q in blocks:
            lq,aq,bq,dq,cq=param(q)
            delta=cq*dp-cp*dq
            require(carry(p+q)==aq*cp+bp*cq,'concatenation')
            require(transform(q,transform(p,7))-transform(p,transform(q,7))
                    ==-Fraction(delta,bp*bq),'commutator defect')
            require((delta==0)==(p+q==q+p)==(primitive(p)==primitive(q)),
                    'common center iff common primitive')
            common_centers+=delta==0
            # A positive integer in the p cylinder makes p an exact legal block.
            for offset in (1,2):
                z=residue(p)+offset*bp;y=transform(p,z)
                require(y.denominator==1,'exact source cylinder')
                y=int(y)
                np=dp*z+cp;nq=dq*y+cq
                require(bp*dp*nq==ap*dq*np+bp*delta,'switch cocycle')
                hp=fuel(p,z);hq=fuel(q,y)
                require(hp>=lp,'source guard for block')
                if delta==0:
                    require(hq==hp-lp,'common-center consumption')
                else:
                    t=val(delta)
                    if hp-lp!=t:require(hq==min(hp-lp,t),'unequal-depth no cancellation')
                    else:require(hq>=t+1,'equal-depth cancellation gain')
                require(literal(z,lp)[0]==p,'literal balanced execution')
                stats['pair_source_cocycles']+=1
            stats['block_pairs']+=1
    # Fixed-source pullbacks, with every chronological guard checked.
    w,values,credits=literal(7,110)
    require(min(credits)>=1,'marked-source full-prefix guards')
    for n in range(111):
        prefix=w[:n];z=values[n]
        require(z==transform(prefix,7),'marked exact affine endpoint')
        for p in blocks:
            rho=((1<<n)*center(p)-carry(prefix))/3**n
            require(qval(Fraction(7)-rho)==n+fuel(p,z),'source-fixed precision')
            # Pulling center through its own block changes no source center.
            rho2=((1<<(n+len(p)))*center(p)-carry(prefix+p))/3**(n+len(p))
            require(rho2==rho,'same-shape pulled center invariant')
            stats['marked_pullback_centers']+=1
    p='100';q='101000';lp,ap,bp,dp,cp=param(p);lq,aq,bq,dq,cq=param(q)
    delta=cq*dp-cp*dq
    require((delta,val(delta))==(684,2),'exact shape-pair separation')
    n=86;z=values[n];y=values[n+lp]
    require(w[n:n+lp+lq]==p+q,'marked switch word')
    actual={'macro_index':n,'z':z,'pending_slots':credits[n],
            'P':p,'Q':q,'endpoint_after_P':y,'hP_before':fuel(p,z),
            'hP_after':fuel(p,y),'hQ_after':fuel(q,y),'center_separation':val(delta)}
    require(actual=={'macro_index':86,'z':10952205219944621,'pending_slots':47,
                     'P':'100','Q':'101000','endpoint_after_P':36963692617313097,
                     'hP_before':5,'hP_after':2,'hQ_after':9,'center_separation':2},
            'actual marked source refill hostile')
    # Arbitrarily large refills: finite controls of the symbolic all-m construction.
    refill_rows=[]
    for m in (1,2,3,4,5,8,16,32,64):
        word=p+q*m+'0';z=residue(word)+(1<<len(word));y=int(transform(p,z))
        got,vs,qs=literal(z,len(word)-1)
        require(got==word[:-1] and min(qs)>=1,
                'all proper tree-prefix macros executed legally')
        require(fuel(p,z)==5 and fuel(q,y)==6*m,'unbounded refill formula')
        # Read the blocked leaf, but never apply its arithmetic update.
        require(qs[-1]==1 and vs[-1]%2==0,'last leaf is blocked')
        try:
            literal(vs[-1],1,qs[-1])
        except RuntimeError as error:
            require(str(error)=='attempted blocked native macro',
                    'specific native guard rejection')
        else:
            raise RuntimeError('blocked native macro was executed')
        stats['blocked_leaf_rejections']+=1
        stats['refill_family_executed_native_macros']+=len(got)
        refill_rows.append({'repeats':m,'seed':str(z),'hP':5,'hQ_after_P':6*m,
                            'tree_length':len(word),'executed_native_macros':len(word)-1})
    # Exact clock/source clone hostile: P and P^2 share center but not action.
    require(center(p)==center(p+p) and transform(p,7)!=transform(p+p,7),
            'same periodic center is not same clocked affine action')
    # Reconstruct the actual old-section reset hostile and corrected indexing.
    x=27;a=6;returned_v3=0;old_hostile=None
    require(a==6,'empty returned-state v3 sum')
    for j in range(20):
        k=val(3*x-1);y=3**(k-1)*((3*x-1)>>k);newa=a+2*k-3
        returned_v3+=val(y,3)
        require(newa==6+2*returned_v3-(j+1),'correct returned-state indexing')
        require(a>0 and newa>0,'old-section finite guards')
        if j==19:
            old_hostile={'index':j,'x':x,'a':a,'k':k,'v3_x_minus_1':val(x-1,3),
                         'charge_before':a+val(x-1,3),'charge_after':newa+val(y-1,3)}
        x,a=y,newa
    require(old_hostile=={'index':19,'x':57971215,'a':21,'k':2,'v3_x_minus_1':4,
                         'charge_before':25,'charge_after':22},'old-root charge hostile')
    return {'status':'PASS','scope':'PROVED_ELEMENTARY_COCYCLE_WITH_FINITE_EXACT_CONTROLS',
            'counts':dict(sorted(stats.items())),'common_center_ordered_pairs':common_centers,
            'actual_marked_source_refill':actual,'unbounded_family_controls':refill_rows,
            'old_root27_charge_hostile':old_hostile,
            'open':'All-prefix root27 survival and arbitrary-shape global rank remain unproved'}


if __name__=='__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))

#!/usr/bin/env python3
"""Exact checks of the forward-polynomial integer-span theorem."""
import importlib.util,pathlib,json,math,itertools
parent=pathlib.Path(__file__).resolve().parent.parent
candidates=[parent/name/'check.py' for name in ('session43-cyclotomic','cyclotomic')]
source=next((p for p in candidates if p.is_file()),None)
if source is None:raise RuntimeError('Missing companion cyclotomic/check.py or session43-cyclotomic/check.py beside this packet')
spec=importlib.util.spec_from_file_location('base',source);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
R=b.require

def generators(n):
    out=[b.eulerian(n)]
    for s in range(1,(n-1)//2+1):
        c=[0]*n
        for k,v in enumerate(b.eulerian(n-2*s)):
            for j in range(2*s+1):c[k+j]+=2*v*math.comb(2*s,j)*(-1)**j
        out.append(c)
    return out

def coordinates(f,g):
    rem=[x-y for x,y in zip(f,g[0])];coef=[1];n=len(f)
    for s,c in enumerate(g[1:],1):
        num=sum(math.comb(k,2*s)*rem[k] for k in range(2*s,n));den=2*math.factorial(n-2*s)
        R(num%den==0,'integral span coordinate');a=num//den;R(a>=0,'nonnegative cycle-excess coordinate');coef.append(a)
        rem=[x-a*y for x,y in zip(rem,c)]
    R(not any(rem),'span reconstruction')
    return coef

def block_poly(n,s):
    pairs=list(itertools.combinations(range(n),2));code=(1<<len(pairs))-1
    if s:code^=1<<pairs.index((0,2*s))
    return b.polynomial(n,code)

def run():
    checks=0;families=0;order5=set()
    for n in range(1,6):
        gs=generators(n)
        for code in range(1<<(n*(n-1)//2)):
            f=b.polynomial(n,code);coordinates(f,gs);checks+=1
            if n==5:order5.add(tuple(f))
    for n in range(1,11):
        gs=generators(n)
        for s in range((n-1)//2+1):
            f=block_poly(n,s);expected=list(gs[0])
            for j in range(1,s+1):expected=[x+math.comb(2*s-1,2*j-1)*y for x,y in zip(expected,gs[j])]
            R(f==expected,'unitriangular actual block family');families+=1
    hostile=[7,38,30,38,7]
    R(tuple(hostile) not in order5,'nonrealizable positive affine hostile')
    R(coordinates(hostile,generators(5))==[1,3,0],'hostile coordinates')
    R(all(c>=0 for c in hostile) and sum(hostile)==120 and hostile==hostile[::-1],'hostile elementary gates')
    R(hostile==[3*x-2*y for x,y in zip(block_poly(5,1),b.eulerian(5))],'hostile affine certificate')
    table=[];strict=[]
    for n in range(1,102):
        gs=generators(n);row={'n':n,'rank':len(gs),'primes':{}}
        for p in (2,3,5,7,11):
            vals=[b.valuation(b.reduce_cyclotomic(c,p),p) for c in gs]
            nonzero=[v for v in vals if v is not None];best=min(nonzero) if nonzero else None
            generic=b.bound(n,p);generic=generic if p==2 else 2*((generic+1)//2)
            R(best is None or best>=generic,'generic bound follows')
            if p==2 and n%2:R(best==b.binary_bound(n),'sharp binary formula')
            if p==2 and n%2==0:R(best is None,'even binary zero')
            if best is not None and best>generic:strict.append({'n':n,'p':p,'generic':generic,'exact':best})
            # Every block valuation computed from the unitriangular family.
            blockvals=[]
            for s in range(len(gs)):
                f=list(gs[0])
                for j in range(1,s+1):f=[x+math.comb(2*s-1,2*j-1)*y for x,y in zip(f,gs[j])]
                blockvals.append(b.valuation(b.reduce_cyclotomic(f,p),p))
            R((min(v for v in blockvals if v is not None) if any(v is not None for v in blockvals) else None)==best,'explicit family attains exact ideal bound')
            row['primes'][p]={'exact_universal_pi_valuation':best,'first_extremal_block_s':next((s for s,v in enumerate(blockvals) if v==best),None),'generic_lower_bound':generic}
        table.append(row)
    print(json.dumps({'exhaustive_coordinate_checks':checks,'actual_block_family_checks':families,'all_order_formula_table_through':101,'strict_improvements_over_generic':strict,'hostile_order5_polynomial':hostile,'table':table},indent=2))
if __name__=='__main__':run()

"""Exhaustive Frobenius-form and subfield-line controls in F16,F64,F256."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path


def need(ok,msg):
    if not ok: raise RuntimeError(msg)


def rank(rows):
    piv={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in piv: x ^= piv[p]
            else: piv[p]=x; break
    return len(piv)


def kernel(rows,n):
    a=list(rows); j=0; piv=[]
    for k in range(n):
        p=next((i for i in range(j,len(a)) if a[i]>>k&1),None)
        if p is None: continue
        a[j],a[p]=a[p],a[j]
        for i in range(len(a)):
            if i!=j and a[i]>>k&1: a[i] ^= a[j]
        piv.append(k); j+=1
    result=[]
    for k in range(n):
        if k not in piv:
            x=1<<k
            for i,p in enumerate(piv):
                if a[i]>>k&1: x |= 1<<p
            result.append(x)
    need(all(not ((x&r).bit_count()%2) for x in result for r in rows),'nullspace')
    return result


def span(basis):
    out=[0]
    for b in basis: out += [x^b for x in out]
    return out


class Field:
    def __init__(self,r,poly): self.r=r; self.size=1<<r; self.poly=poly
    def mul(self,a,b):
        z=0
        while b:
            if b&1: z ^= a
            b >>= 1; a <<= 1
            if a>>self.r: a ^= self.poly
        return z
    def power(self,a,n):
        z=1
        while n:
            if n&1: z=self.mul(z,a)
            a=self.mul(a,a); n >>= 1
        return z
    def trace(self,a):
        z=0
        for _ in range(self.r): z ^= a; a=self.mul(a,a)
        need(z in (0,1),'trace field'); return z
    def orbit(self,a):
        o=[a]; x=self.mul(a,a)
        while x!=a:
            o.append(x); x=self.mul(x,x)
            need(len(o)<=self.r,'Frobenius orbit')
        return o


def analyze(r,poly):
    f=Field(r,poly); q=1<<(r//2)
    need(all(f.power(a,f.size-1)==1 for a in range(1,f.size)),'field inverse test')
    pairs=list(combinations(range(r),2)); ambient=len(pairs)
    def coeff(a,b):
        return sum((((a>>i&1)*(b>>j&1)) ^ ((a>>j&1)*(b>>i&1)))<<k for k,(i,j) in enumerate(pairs))
    def bilinear(mask,a,b): return (mask & coeff(a,b)).bit_count()%2
    eq=[coeff(f.mul(1<<i,1<<i),f.mul(1<<j,1<<j))^coeff(1<<i,1<<j) for i,j in pairs]
    ibasis=kernel(eq,ambient); forms=span(ibasis)
    nondeg=[z for z in forms if rank([sum(bilinear(z,1<<i,1<<j)<<j for j in range(r)) for i in range(r)])==r]
    target=sum(f.trace(f.mul(1<<i,f.power(1<<j,q)))<<k for k,(i,j) in enumerate(pairs))
    need(target in nondeg,'trace form')
    sub=[a for a in range(1,f.size) if f.power(a,q)==a]
    lines=sorted({tuple(sorted(f.mul(a,x) for a in sub)) for x in range(1,f.size)})
    need(len(lines)==q+1,'spread size')
    def line_basis(line):
        b=[]
        for a in line:
            if rank(b+[a])>len(b): b.append(a)
        need(len(b)==r//2,'line dimension'); return b
    rows=[]; survivors=[]; line_data=[]
    for line in lines:
        basis=line_basis(line)
        constraints=[sum(bilinear(z,a,b)<<k for k,z in enumerate(ibasis)) for a,b in combinations(basis,2)]
        allowed=[z for z in nondeg if all(bilinear(z,a,b)==0 for a,b in combinations(basis,2))]
        orbit={line}; nxt=tuple(sorted(f.mul(a,a) for a in line))
        while nxt!=line:
            orbit.add(nxt); nxt=tuple(sorted(f.mul(a,a) for a in nxt))
        rows.append(constraints); survivors.append(set(allowed))
        line_data.append(dict(points=line,basis=basis,Frobenius_period=len(orbit),constraint_rank=rank(constraints),
                              nondegenerate_survivors=len(allowed),surviving_masks=allowed))
    def minimal(kind):
        for k in range(len(lines)+1):
            good=[]
            for ids in combinations(range(len(lines)),k):
                if kind=='linear': passed=rank([row for i in ids for row in rows[i]])==len(ibasis)-1
                else:
                    keep=set(nondeg)
                    for i in ids: keep &= survivors[i]
                    passed=keep=={target}
                if passed: good.append(ids)
            if good: return dict(minimum=k,number_of_minimizers=len(good),first=good[0])
        raise RuntimeError('all lines must work')
    allrows=[row for l in rows for row in l]
    need(rank(allrows)==len(ibasis)-1,'full-spread reconstruction')
    normal=next(f.orbit(a) for a in range(1,f.size) if len(f.orbit(a))==r and rank(f.orbit(a))==r)
    element_orbits=[]; seen=set()
    for a in range(1,f.size):
        if a not in seen:
            o=f.orbit(a); element_orbits.append(o); seen.update(o)
    fixed_tournaments=[]
    for power in range(r):
        perm=[f.power(a,1<<power) for a in range(1,f.size)]
        seen=set(); lengths=[]
        for a in range(1,f.size):
            if a in seen: continue
            x=a; cycle=[]
            while x not in seen: seen.add(x); cycle.append(x); x=perm[x-1]
            lengths.append(len(cycle))
        if any(l%2==0 for l in lengths): exponent=None
        else:
            from math import gcd
            exponent=sum((l-1)//2 for l in lengths)+sum(gcd(a,b) for a,b in combinations(lengths,2))
        fixed_tournaments.append(dict(Frobenius_power=power,cycle_lengths=dict(sorted(Counter(lengths).items())),
                                      fixed_tournament_log2=exponent))
    return dict(degree=r,polynomial_bits=poly,normal_basis=normal,
                invariant_form_dimension=len(ibasis),invariant_form_basis=ibasis,nondegenerate_invariant_forms=nondeg,
                trace_form_mask=target,spread_size=len(lines),line_data=line_data,
                minimum_lines_nondegenerate=minimal('nondeg'),minimum_lines_linear=minimal('linear'),
                line_summary=[dict(period=p,rank=k,survivors=s,count=c) for (p,k,s),c in sorted(Counter((d['Frobenius_period'],d['constraint_rank'],d['nondegenerate_survivors']) for d in line_data).items())],
                element_orbit_lengths=dict(sorted(Counter(map(len,element_orbits)).items())),
                tournament_symmetry_controls=fixed_tournaments)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output'); args=parser.parse_args()
    result={str(1<<r):analyze(r,p) for r,p in [(4,19),(6,67),(8,283)]}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')


if __name__=='__main__':main()

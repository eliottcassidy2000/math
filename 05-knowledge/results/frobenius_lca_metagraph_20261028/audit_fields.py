"""Independent circulant coordinates and oriented-pair orbit audit."""
import argparse
from itertools import combinations
import json
from pathlib import Path
from fields import Field, need


def remainder(a,b):
    while a and a.bit_length()>=b.bit_length(): a ^= b<<(a.bit_length()-b.bit_length())
    return a


def polynomial_gcd(a,b):
    while b: a,b=b,remainder(a,b)
    return a


def audit(d):
    r=d['degree']; m=r//2; f=Field(r,d['polynomial_bits']); normal=d['normal_basis']
    # Build the inverse normal-coordinate map by exhaustive vector addition.
    coordinate={}
    for word in range(1<<r):
        v=0
        for i,b in enumerate(normal):
            if word>>i&1: v ^= b
        need(v not in coordinate,'normal basis independent'); coordinate[v]=word
    pairs=list(combinations(range(r),2)); nondeg=[]; coefficients={}
    for a in range(1<<m):
        c=0
        for i in range(1,m):
            if a>>(i-1)&1: c ^= (1<<i)|(1<<(r-i))
        if a>>(m-1)&1: c ^= 1<<m
        if polynomial_gcd(c,(1<<r)|1)!=1: continue
        def B(x,y):
            u=coordinate[x]; v=coordinate[y]; value=0
            for i in range(r):
                if u>>i&1:
                    for j in range(r):
                        if v>>j&1: value ^= (c>>((j-i)%r))&1
            return value
        mask=sum(B(1<<i,1<<j)<<k for k,(i,j) in enumerate(pairs))
        nondeg.append(mask); coefficients[str(mask)]=a
        for line in d['line_data']:
            isotropic=all(B(x,y)==0 for x,y in combinations(line['basis'],2))
            need(isotropic==(mask in line['surviving_masks']),'line constraint cross-check')
    need(sorted(nondeg)==sorted(d['nondegenerate_invariant_forms']),'circulant enumeration')
    orbit_controls=[]
    for item in d['tournament_symmetry_controls']:
        k=item['Frobenius_power']; perm={a:f.power(a,1<<k) for a in range(1,1<<r)}
        seen=set(); count=0; obstruction=False
        for a,b in combinations(range(1,1<<r),2):
            if (a,b) in seen: continue
            count+=1; x,y=a,b
            while True:
                seen.add(tuple(sorted((x,y))))
                x,y=perm[x],perm[y]
                if (x,y)==(b,a): obstruction=True
                if (x,y)==(a,b): break
        exponent=None if obstruction else count
        need(exponent==item['fixed_tournament_log2'],'oriented-pair orbit consistency')
        orbit_controls.append(dict(power=k,unordered_pair_orbits=count,orientation_obstruction=obstruction))
    return dict(degree=r,circulant_nondegenerate_coefficients=coefficients,
                trace_coefficients=coefficients[str(d['trace_form_mask'])],oriented_pair_audit=orbit_controls)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output'); args=parser.parse_args()
    data=json.loads(Path(__file__).with_name('fields-normal.json').read_text())
    result=[audit(data[str(1<<r)]) for r in (4,6,8)]
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')


if __name__=='__main__':main()

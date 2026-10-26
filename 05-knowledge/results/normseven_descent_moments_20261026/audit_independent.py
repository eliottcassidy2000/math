"""No producer imports: finite generated-group audit and integer inverse audit."""
import itertools
import json
from pathlib import Path
from math import gcd


def need(ok,label):
    if not ok: raise RuntimeError(label)


def perm_det(a):
    total=0
    for p in itertools.permutations(range(len(a))):
        inversions=sum(p[i]>p[j] for i in range(len(a)) for j in range(i+1,len(a)))
        term=(-1)**inversions
        for i,j in enumerate(p): term*=a[i][j]
        total+=term
    return total


def main():
    base=Path(__file__).resolve().parent
    producer=json.loads((base/'lifts.normal.json').read_text())
    rows=[]
    for m in range(2,13):
        seen={(1,0,0,1)}; todo=list(seen)
        while todo:
            a,b,c,d=todo.pop()
            for t in [((a+c)%m,(b+d)%m,c,d),(a,b,(c+a)%m,(d+b)%m),(c,d,a,b)]:
                if t not in seen: seen.add(t); todo.append(t)
        invertible={v for v in itertools.product(range(m),repeat=4) if gcd(v[0]*v[3]-v[1]*v[2],m)==1}
        predicted={v for v in invertible if (v[0]*v[3]-v[1]*v[2])%m in {1,m-1}}
        need(seen==predicted,'BFS image equals determinant-sign subgroup')
        row=producer['matrix_census'][m-2]
        need(len(seen)==row['same_rank_lifts'] and len(invertible)==row['invertible_2_by_2'],'independent census agreement')
        rows.append({'m':m,'generated':len(seen),'all_invertible':len(invertible)})
    need(perm_det(producer['six_to_seven_mod7_lift'])==1,'Leibniz 7x7 determinant')
    need(perm_det(producer['two_coordinate_mod7_lift'])==1,'Leibniz 2x2 determinant')
    # Composite unit rings, without assuming cyclic unit groups.
    cases=0
    for N in range(2,43):
        units=[u for u in range(N) if gcd(u,N)==1]
        for b in units:
            original={u for u in units if pow(u,7,N)==b}
            target={(w,z) for w in units for z in units
                    if pow(w,7,N)==b*b%N and pow(z,7,N)==pow(b,7,N)
                    and pow(pow(w,7,N),-1,N)*z*z%N==1}
            image={(u*u%N,pow(u,7,N)) for u in original}
            need(target==image,'marked section over noncyclic unit rings')
            for w,z in target:
                u=pow(w,25,N)*pow(pow(z,7,N),-1,N)%N
                need(u in original,'ring unit inverse')
            cases+=1
    # Independent 35 determinant list and characteristic-three dependencies.
    lines=[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1)]
    normal=[]; mod3=[]
    for ids in itertools.combinations(range(7),3):
        d=perm_det([lines[i] for i in ids])
        if d==0: normal.append(list(ids))
        if d%3==0: mod3.append(list(ids))
    need(normal==[[0,1,2],[0,3,4],[0,5,6]],'normal triple list')
    need(mod3==normal+[[1,3,6],[2,4,5]],'char-three triple list')
    return {'group_census':rows,'unit_ring_cases':cases,'ring_moduli':[2,42],
            'normal_triples':normal,'char3_triples':mod3,'seven_rank_lift_determinant':1}


if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))

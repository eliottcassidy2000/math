"""Exact fruit-curve trees and their positivity/torsion obstructions.

Standard library only. All checks survive -O. Full Mordell--Weil generation
is a cited input in the note, not inferred from this finite computation.
"""
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from math import gcd, isqrt, lcm
from pathlib import Path
import json

CHECKS = 0
G = (F(-4), F(28))
T = (F(56), F(728))
O = None
A = 154476802108746166441951315019919837485664325669565431700026634898253202035277999
B = 368751317941299998271978115652254748254929799689719709962831374716372246340555790
C = 43736126779286972578612526023713901528165375581616136186214379933784234677720360
BERGGREN = (
    ((1,-2,2),(2,-1,2),(2,-2,3)),
    ((1,2,2),(2,1,2),(2,2,3)),
    ((-1,2,2),(-2,1,2),(-2,2,3)),
)


def check(value, label):
    global CHECKS
    if not value:
        raise RuntimeError(label)
    CHECKS += 1


def neg(p):
    return None if p is None else (p[0],-p[1])


def add(p,q,a2=109,a4=224):
    if p is None:
        return q
    if q is None:
        return p
    x,y=p; u,v=q
    if x==u and y==-v:
        return None
    slope=(v-y)/(u-x) if x!=u else (3*x*x+2*a2*x+a4)/(2*y)
    z=slope*slope-a2-x-u
    return z,slope*(x-z)-y


def mul(n,p,a2=109,a4=224):
    if n<0:
        return mul(-n,neg(p),a2,a4)
    result=None
    while n:
        if n%2:
            result=add(result,p,a2,a4)
        p=add(p,p,a2,a4)
        n//=2
    return result


def on_curve(p,a2=109,a4=224,a6=0):
    return p is None or p[1]**2==p[0]**3+a2*p[0]**2+a4*p[0]+a6


def normalize(v):
    v=tuple(F(x) for x in v)
    den=lcm(*(x.denominator for x in v))
    out=tuple(int(x*den) for x in v)
    content=gcd(gcd(out[0],out[1]),out[2])
    check(content!=0,'nonzero projective vector')
    out=tuple(x//content for x in out)
    return out if next(x for x in out if x)>0 else tuple(-x for x in out)


def from_point(p):
    if p is None:
        return (1,-1,0)
    x,y=p
    return normalize((56-x+y,56-x-y,-56-12*x))


def to_point(v):
    a,b,c=v
    x,y,z=-28*(a+b+2*c),364*(a-b),6*(a+b)-c
    return None if z==0 else (F(x,z),F(y,z))


def cubic(v):
    a,b,c=v; s=a+b+c
    return s**3-6*s*(a*b+a*c+b*c)+7*a*b*c


def fruit(v):
    a,b,c=v
    return F(a,b+c)+F(b,a+c)+F(c,a+b)


def positive(v):
    return all(x>0 for x in v) or all(x<0 for x in v)


def point_summary(p):
    v=from_point(p)
    return {'positive':positive(v),'signs':[1 if x>0 else -1 if x<0 else 0 for x in v],
            'coordinate_digits':[len(str(abs(x))) for x in v],
            'primitive_triple_sha256':sha256(json.dumps(v).encode()).hexdigest()}


def matrix_vector(matrix,v):
    return tuple(sum(row[j]*v[j] for j in range(3)) for row in matrix)


def coefficient_parent(n):
    if n==1:
        return None
    p=(n+1)//3
    return p,n-3*p


def inverse_three(n,k):
    if n%3 or k%3:
        return []
    return [(n//3,j) for j in range(6) if 3*j%6==k]


def inverse_tree(n,k):
    seen=set(); queue=deque([(n,k)])
    while queue:
        p=queue.popleft()
        check(p not in seen,'inverse tripling has no repeated vertex')
        seen.add(p)
        queue.extend(inverse_three(*p))
    return seen


def count_mod(p,a2,a4,a6):
    brute=1+sum(y*y%p==(x**3+a2*x*x+a4*x+a6)%p for x in range(p) for y in range(p))
    characters=1
    for x in range(p):
        z=(x**3+a2*x*x+a4*x+a6)%p
        characters+=1 if z==0 else 2 if pow(z,(p-1)//2,p)==1 else 0
    check(brute==characters,'independent finite-field point count')
    return brute


def cube_root(n):
    low,high=0,1<<((n.bit_length()+2)//3)
    while low<high:
        middle=(low+high+1)//2
        if middle**3<=n:
            low=middle
        else:
            high=middle-1
    return low


def main():
    global CHECKS
    CHECKS=0
    check(B%10==C%10==0,'exact decimal divisions')
    repaired=(A,B//10,C//10)
    check(cubic((A,B,C))!=0 and fruit((A,B,C))!=4,'literal input fails')
    check(cubic(repaired)==0 and fruit(repaired)==4 and positive(repaired),'repaired positive solution')
    check(normalize(repaired)==repaired,'repaired primitive content')
    torsion=[mul(k,T) for k in range(6)]
    check(len(set(torsion))==6 and mul(6,T) is None,'exact torsion order six')
    check(mul(2,T)==(F(4),F(52)) and mul(3,T)==(F(0),F(0)),'torsion labels')
    check(G not in torsion and on_curve(G),'infinite-order control after inherited torsion bound')
    check(count_mod(11,109,224,0)==12 and count_mod(17,109,224,0)==18,'torsion reduction bound inputs')
    points={0:None}
    current=None; positive_indices=[]
    for n in range(1,41):
        current=add(current,G)
        points[n]=current
        check(on_curve(current) and current==mul(n,G),'sequential versus binary multiplication')
        v=from_point(current)
        check(to_point(v)==current and cubic(v)==0 and fruit(v)==4,'curve fruit round trip')
        if positive(v):
            positive_indices.append(n)
    check(positive_indices==[9,17],'positive multiple census through forty')
    H=points[9]
    check(from_point(H)==repaired,'large solution exactly nine G')
    children=[points[18],points[27],points[36]]
    check(all(not positive(from_point(p)) for p in children),'all three seeded children fail positivity')
    check(positive(from_point(points[17])) and 17%9!=0,'positive point outside seed cyclic subgroup')

    # The whole level-eight graph is checked in exact integer coefficients and
    # independently in primitive Pythagorean side coordinates, not giant E points.
    queue=deque([(1,(3,4,5),())]); seen={}; triangle_seen=set(); level_counts={}
    depth_cap=8
    while queue:
        n,triangle,word=queue.popleft()
        check(n not in seen and triangle not in triangle_seen,'tree injectivity in both carriers')
        seen[n]=word; triangle_seen.add(triangle)
        a,b,c=triangle
        check(a*a+b*b==c*c and gcd(gcd(a,b),c)==1 and a%2==1 and b%2==0,'primitive triangle preserved')
        depth=len(word); level_counts[depth]=level_counts.get(depth,0)+1
        if word:
            parent,digit=coefficient_parent(n)
            check(seen[parent]==word[:-1] and digit==word[-1],'unique parent with branch label')
        if depth<depth_cap:
            for digit,matrix in zip((-1,0,1),BERGGREN):
                queue.append((3*n+digit,matrix_vector(matrix,triangle),word+(digit,)))
    check(set(seen)==set(range(1,(3**(depth_cap+1)+1)//2)),'all positive coefficient indices to level cap')
    check(all(level_counts[d]==3**d for d in range(depth_cap+1)),'full ternary level sizes')
    for n in range(1,14):
        for digit in (-1,0,1):
            actual=add(mul(3,points[n]),mul(digit,G))
            check(actual==points[3*n+digit],'actual elliptic tree branches through depth three')
    q=mul(2,T)
    for digit in (-1,0,1):
        check(add(mul(3,G),mul(digit,G))==add(mul(3,add(G,q)),mul(digit,G)),
              'rational three-torsion collision off the cyclic tree')

    group_orbit={add(mul(sign,H),r) for sign in (-1,1) for r in torsion}
    check(len(group_orbit)==12,'finite symmetry orbit twelve')
    fruit_orbit={from_point(p) for p in group_orbit}
    permuted={normalize(v) for v in permutations(repaired)}
    check(len(permuted)==6 and permuted=={v for v in fruit_orbit if positive(v)},'exactly six positive permutation images')
    check(all(cubic(v)==0 and fruit(v)==4 for v in fruit_orbit),'all symmetry images valid signed solutions')
    berggren_hostiles=[]
    for name,matrix in zip(('L','M','R'),BERGGREN):
        small=matrix_vector(matrix,(11,4,-1))
        large=matrix_vector(matrix,repaired)
        check(cubic(small)!=0 and cubic(large)!=0,'literal Berggren matrix does not preserve fruit cubic')
        berggren_hostiles.append({'branch':name,'small_image':small,'small_residual':cubic(small),
                                  'large_residual_sign':1 if cubic(large)>0 else -1})
    discriminants=[]
    for i in range(3):
        a,b,c=(11,4,-1)[i],(11,4,-1)[(i+1)%3],(11,4,-1)[(i+2)%3]
        disc=-3*a*a+6*a*(b+c)+5*b*b+30*b*c+5*c*c
        check(disc<0 or isqrt(disc)**2!=disc,'coordinate cubic has no other rational root')
        discriminants.append(disc)

    # An independent forward map on the group coefficients checks inverse guards.
    preimages={}
    for n in range(-81,82):
        if n:
            for k in range(6):
                preimages.setdefault((3*n,3*k%6),[]).append((n,k))
    inverse_cases=0
    for n in range(-81,82):
        if n:
            for k in range(6):
                check(inverse_three(n,k)==preimages.get((n,k),[]),'inverse tripling full finite coefficient universe')
                inverse_cases+=1
    inverse_records=[]
    for r in range(9):
        vertices=inverse_tree(3**r,0)
        check(len(vertices)==1+3*r,'inverse tripling linear size law')
        inverse_records.append({'root_coefficient':3**r,'r':r,'vertices':len(vertices)})
    seven=inverse_tree(9,0)
    check(seven=={(9,0),(3,0),(3,2),(3,4),(1,0),(1,2),(1,4)},'nine G inverse tree seven vertices')
    for n,k in seven-{(9,0)}:
        p=add(points[n],torsion[k])
        target=add(points[3*n],torsion[3*k%6])
        check(mul(3,p)==target,'actual rational inverse tree edge')

    check(2**10==10**3+24,'exact power gap twenty-four')
    mordell=(F(10),F(32))
    double=mul(2,mordell,0,0)
    check(on_curve(mordell,0,0,24) and on_curve(double,0,0,24),'Mordell gap point and double')
    check(double==(F(505,256),F(23053,4096)),'Mordell double exact')
    check(count_mod(19,109,224,0)==18 and count_mod(19,0,0,24)==27,'good reduction distinguishes isogeny classes')
    power_gaps=[]
    for k in range(1,201):
        target=2**k-24
        if target>0:
            n=cube_root(target)
            if n**3==target:
                power_gaps.append([k,n])
    check(power_gaps==[[5,2],[10,10]],'finite power-gap census only')

    source=Path(__file__)
    seed_summaries={str(n):point_summary(points[n]) for n in (18,27,36)}
    outside_summary={'n':17,**point_summary(points[17])}
    result={'status':'PASS','checks_passed':CHECKS,
        'source_sha256_lf':sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        'scope':'Written subgroup/tree/positivity proofs; finite exact controls. Full rational inverse-tree completeness uses cited Mordell--Weil generation; no rank or global-minimality computation.',
        'universe':{'elliptic_multiples':[1,40],'actual_elliptic_tree_depth':3,
                    'integer_coefficient_and_Berggren_tree_depth':8,'tree_vertices':len(seen),
                    'inverse_coefficient_cases':inverse_cases,'inverse_roots_3_power_range':[0,8],
                    'power_gap_exponents':[1,200]},
        'literal_input':{'A':str(A),'B':str(B),'C':str(C),'fruit_sum':str(fruit((A,B,C)))},
        'repaired_input':[str(x) for x in repaired],
        'repaired_digits':[len(str(x)) for x in repaired],
        'positive_nG_through_40':positive_indices,
        'seed_children':seed_summaries,
        'positive_outside_seed':outside_summary,
        'ternary_level_counts':level_counts,
        'symmetry_orbit':{'total':len(group_orbit),'positive':len(permuted)},
        'literal_Berggren_hostiles':berggren_hostiles,
        'coordinate_replacement_discriminants':discriminants,
        'inverse_tripling_counts':inverse_records,'inverse_nine_G_vertices':sorted(seven),
        'mordell_gap_curve':{'point':['10','32'],'double':[str(x) for x in double],
                             'counts_mod19':{'fruit':18,'gap':27},'finite_power_gap_solutions':power_gaps}}
    source.with_suffix('.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':'PASS','checks_passed':CHECKS,'tree_vertices':len(seen),'output':source.with_suffix('.json').name}))


if __name__=='__main__':
    main()

"""Independent root-deletion checks. No producer implementation imported."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import sys


def demand(test, why):
    if not test:
        raise RuntimeError(why)


def rank(matrix):
    a = [[Fraction(v) for v in row] for row in matrix]
    r = 0
    for j in range(len(a[0]) if a else 0):
        k = next((k for k in range(r, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[r], a[k] = a[k], a[r]
        f = a[r][j]
        a[r] = [v / f for v in a[r]]
        for k in range(len(a)):
            if k != r:
                f = a[k][j]
                a[k] = [u - f * v for u, v in zip(a[k], a[r])]
        r += 1
    return r


def det(a, b, c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


forms = ((1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1),(1,1,-1),(1,-1,1))
triples = {t for t in combinations(range(9), 3) if det(*(forms[i] for i in t)) == 0}
demand(len(triples) == 8, "global magic incidence")
hist = Counter()
for seven in combinations(range(9), 7):
    own = [t for t in triples if set(t) <= set(seven)]
    for omitted in seven:
        t = sum(omitted in line for line in own)
        hist[t] += 1
demand(sum(hist.values()) == 252, "deletion universe")

models = []
for t in range(4):
    vertices = tuple(product((0,1), repeat=t))
    count = len(vertices)
    degree = 32 // count
    branch_points = 6-t
    # Hurwitz, computed in rational arithmetic before extracting genus.
    genus = 1 + Fraction(degree * (branch_points-4), 4)
    demand(genus.denominator == 1, "integral genus")
    g = int(genus)
    mult = 16 // count
    matrix = [[2*g-2 if a == b else mult*int(sum(x!=y for x,y in zip(a,b)) == 1) for b in vertices] for a in vertices]
    observed_rank = rank(matrix)
    demand(observed_rank == count-t, "independent rational matrix rank")
    null_vectors = [[1 if v[j] == 0 else -1 for v in vertices] for j in range(t)]
    demand(rank(null_vectors) == t, "hyperplane relation independence")
    demand(all(all(sum(row[k]*u[k] for k in range(count)) == 0 for row in matrix) for u in null_vectors), "hyperplane relations span nullspace")
    eigs = Counter()
    for mask in product((0,1), repeat=t):
        vector = [(-1)**sum(a*b for a,b in zip(mask,v)) for v in vertices]
        value = 32 // count * (1-sum(mask))
        demand(all(sum(matrix[i][j]*vector[j] for j in range(count)) == value*vector[i] for i in range(count)), "Walsh action")
        eigs[value] += 1
    demand(sum(v for v in eigs.values()) == count, "complete Walsh basis")
    demand(sum(x for row in matrix for x in row) == 32, "branch self-intersection")
    demand(count*g + 8*t-count+1 == 17, "arithmetic genus17")
    demand(2*(24+8*t)-count*(2-2*g) == 80, "resolved quotient Euler80")
    models.append({"t":t,"components":count,"genus":g,"matrix_rank":observed_rank,"eigenvalues":dict(sorted(eigs.items()))})

# Recheck the explicit central curves over a DIFFERENT field from the producer.
q = 41
iunit = 9
demand(iunit*iunit % q == q-1, "sqrt minus one")
conic_points = set()
for x,y in product(range(q), repeat=2):
    if (1+x*x-y*y) % q == 0:
        conic_points.add((1,x,y))
for y in range(q):
    if (1-y*y) % q == 0:
        conic_points.add((0,1,y))
demand(len(conic_points) == q+1, "central conic enumeration")
conics=[]
for signs in product((1,-1), repeat=3):
    ps=set()
    for a,b,c in conic_points:
        ps.add((a,iunit*signs[0]*a % q,b,iunit*signs[1]*b % q,c,iunit*signs[2]*c % q))
    conics.append(ps)
intersections=Counter()
nodes=set()
for a,b in combinations(range(8),2):
    common=conics[a]&conics[b]
    want=2 if (a^b).bit_count()==1 else 0
    demand(len(common)==want,"F41 actual conic intersections")
    intersections[len(common)] += 1
    nodes |= common
demand(len(nodes)==24 and all(sum(p in c for c in conics)==2 for p in nodes), "ordinary two-component nodes")

# All 5040 labellings verify the common simplex/cube criterion, independently
# of the root producer: pair differences sum to the central binary vector.
from itertools import permutations
joint=Counter()
for v in permutations(range(1,8)):
    d=(v[1]^v[2],v[3]^v[4],v[5]^v[6])
    demand(d[0]^d[1]^d[2]==v[0],"pair difference sum")
    orbit={0}
    for x in d:
        orbit |= {y^x for y in tuple(orbit)}
    binary_rank=(len(orbit)).bit_length()-1
    bad=sum(x==v[0] for x in d)
    demand((bad,binary_rank) in ((0,3),(1,2),(3,1)), "cube/Fano shared rank criterion")
    joint[(bad,binary_rank)] += 1

root=Path(__file__).resolve().parents[1]
producer=root/'root_deletion'/'deletion_receipt.json'
prior=json.loads(producer.read_text(encoding='utf-8'))
demand(prior['t_histogram']=={str(k):v for k,v in hist.items()},"producer incidence histogram")
data={"status":"FINITE-EXACT independent controls; accompanying audit owns proof judgment",
      "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      "producer_proof_sha256":hashlib.sha256((root/'root_deletion'/'DELETION_AND_CUBE.md').read_bytes()).hexdigest(),
      "producer_script_sha256":hashlib.sha256((root/'root_deletion'/'check_deletion.py').read_bytes()).hexdigest(),
      "cases":252,"incident_triple_histogram":dict(sorted(hist.items())),"cube_models":models,
      "central_F41":{"points_per_conic":len(conic_points),"intersection_size_counts":dict(sorted(intersections.items())),"nodes":len(nodes)},
      "joint_code_cube_labellings":[{"bad_triples":r,"translation_rank":s,"labelings":n} for (r,s),n in sorted(joint.items())]}
encoded=json.dumps(data,indent=2,sort_keys=True)+'\n'
if len(sys.argv)>1:
    Path(sys.argv[1]).write_text(encoded,encoding='utf-8',newline='\n')
else:
    print(encoded,end='')

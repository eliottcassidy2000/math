from itertools import combinations,permutations
p=list(combinations(range(7),2))
def H(code):
    arcs={e if code>>i&1 else e[::-1] for i,e in enumerate(p)}
    return sum(all(e in arcs for e in zip(w,w[1:])) for w in permutations(range(7)))
D=(0,24,2,768,24^2,24^768,2^768,24^2^768)
vals=[H(863667^i) for i in D]
if vals!=[79,133,77,15,131,69,13,65]:
    raise RuntimeError(vals)
if [vals[k]-vals[i]-vals[j]+vals[0] for k,i,j in ((4,1,2),(5,1,3),(6,2,3))]!=[0]*3:
    raise RuntimeError("nonzero pair at base")
if vals[7]-sum(vals[4:7])+sum(vals[1:4])-vals[0]!=-2:
    raise RuntimeError("third-order defect")
print("literal 7-vertex cube: all base pairs vanish, third dividend −2")

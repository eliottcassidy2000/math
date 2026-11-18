from collections import Counter
from itertools import permutations

# A is the run of a variable-size transitive block. B and C are two
# separate (fixed) vertices of the former B block.
order = "0ABC456"
BASE = tuple(order)
fwd = {(x,y) for i,x in enumerate(order) for y in order[i+1:]}
def flip(edges,x,y):
    edges.remove((x,y))
    edges.add((y,x))
for x in "BC": flip(fwd,"0",x)
hi=fwd.copy()
flip(hi,"0","4")
flip(hi,"4","6")

def templates(edges,k):
    return sorted({w for w in permutations(order+"A"*(k-1))
                   if all(e in edges for e in zip(w,w[1:]))})

if __name__=="__main__":
    for k in (1,2):
        a,b=templates(fwd,k),templates(hi,k)
        c=set(a)&set(b)
        print(k,len(a),len(b),len(c),''.join(BASE) in {''.join(w) for w in c})
        print("low:"," ".join(map("".join,a)))
        print("shared:"," ".join(sorted(map("".join,c))))
        print("high:"," ".join(map("".join,b)))

def deletions(edges, word):
    one = set(templates(edges, 1))
    return tuple(word[:i] + word[i+1:] if word[:i] + word[i+1:] in one else None
                 for i, x in enumerate(word) if x == "A")

def face_histogram(edges, base):
    return Counter(sum(v not in (None,base) for v in deletions(edges,w))
                   for w in templates(edges,2))

for E, expected, alt in ((fwd, (2,1), {2,4,6}), (hi, (16,5), {30,32,34})):
    for w in templates(E,2):
        print(''.join(w),'→',tuple(''.join(v) if v else 'invalid' for v in deletions(E,w)))
    hist = face_histogram(E,BASE)
    assert (hist[1],hist[2]) == expected and hist[0]==0
    assert {2*face_histogram(E,b)[1] for b in templates(E,1)} == alt
    print("two-run template faces (one survives, both survive):",expected)
    print("width-two, exactly one nonbase face, all choices:",sorted(alt))

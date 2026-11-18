from itertools import combinations
A=[252,189,504,128,8,24,58,112,251]
B=[244,253,504,129,8,24,56,112,251]
def h(out):
 n=len(out); D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,v in enumerate(row):
   if not v: continue
   todo=out[i]&~S
   while todo:
    j=todo&-todo; D[S|j][j.bit_length()-1]+=v; todo^=j
 return sum(D[-1])
def edits(out,arcs):
 out=out.copy()
 for i,j in arcs:out[i]^=1<<j;out[j]^=1<<i
 return out
es=list(combinations(range(9),2));es.remove((1,6));es.remove((0,3))

assert [h(A),h(B)]==[131,133]
hits={}
for k in (1,2):
    hits[k]={}
    for arcs in combinations(es,k):
        x,y=h(edits(A,arcs)),h(edits(B,arcs))
        if y-x==2 and (x,y)!=(131,133): hits[k][arcs]=(x,y)
assert not hits[1]
assert hits[2]=={((0,5),(2,6)):(325,327),((0,5),(5,8)):(441,443),((1,2),(3,4)):(57,59),((3,8),(6,8)):(391,393)}
for (e,f),xy in hits[2].items():
    l=[h(edits(A,S)) for S in ([e],[f],[e,f])]
    r=[h(edits(B,S)) for S in ([e],[f],[e,f])]
    ja=l[2]-l[0]-l[1]+131
    jb=r[2]-r[0]-r[1]+133
    g=(r[0]-l[0],r[1]-l[1])
    assert g[0]+g[1]-4+jb-ja==0
    print(e,f,g,(ja,jb),xy)
n=0
for E,F in combinations(hits[2],2):
    if set(E)&set(F):continue
    n+=1
    x,y=h(edits(A,E+F)),h(edits(B,E+F))
    assert y-x!=2
    print("noncomposing",E,F,x,y)
assert n==5

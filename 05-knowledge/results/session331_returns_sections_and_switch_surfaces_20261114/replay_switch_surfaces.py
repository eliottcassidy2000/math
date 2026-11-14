from itertools import combinations,product
from math import isqrt

C=(1,63,81,88,56,65,35,86,58,42,79,90,54,67,33,48,96,73,8,92,77,44,37,27,94,50,14,22,99,70,51,93,28,72,97,47,53,91,78,66,34,30,19,17,83,38,11,25,75,46,98,71,29,52,69,31,18,82,87,57,64,36,85,84,16,9,7,2,62,59,41,80,89,55,45,76,68,32,49,15,10,39,61,20,5,95,74,26,23,13,12,24,40,60,4,21,43,6,3)
R=(1,80,89,32,68,76,93,51,49,95,74,70,99,97,72,28,53,91,30,34,47,17,64,36,45,19,81,63,37,84,85,15,66,55,26,38,83,86,58,42,79,65,35,29,92,52,69,75,94,6,43,78,22,27,9,40,41,8,73,96,48,16,20,61,60,21,4,77,44,56,88,12,24,57,87,13,23,98,71,50,14,2,7,18,82,62,59,5,11,25,39,10,90,31,33,67,54,46,3)
# Inherited odd Gerbicz return: singleton tokens or (offset, reverse).
G=(1,(-1,0),(1,0),(-7,1),(6,0),(-6,0),(0,1),11,(-5,1),5,4,12,
   (-12,0),(12,0),(7,1),(-8,0),(2,1),(-3,0),9,7,(4,0),(-4,0),10,6,
   (5,0),(-11,1),2,(-2,0),8,(3,0),(-9,1),(9,1),(-10,0),(10,0),(11,0),(8,1),3)

def req(ok,msg):
    if not ok:raise RuntimeError(msg)
def sq(n):
    return isqrt(n)**2==n
def edge(a,b):
    return frozenset((a,b))
def edges(P):
    return {edge(a,b) for a,b in zip(P,P[1:]+P[:1])}
def create(P,i,j,twist):
    a,b,c,d=P[i],P[i+1],P[j],P[j+1]
    req(len({a,b,c,d})==4,"non-disjoint cuts")
    return {edge(a,b),edge(c,d)},({edge(a,c),edge(b,d)} if twist else {edge(a,d),edge(b,c)})
def test(path,moves,pick):
    E=edges(path)
    req(len(E)==len(path) and all(sq(sum(e)) for e in E),"invalid base spanning cycle")
    for j in pick:
        A,B=moves[j]
        req(A<=E and not B&E,"wrong toggle")
        E=(E-A)|B
    N=len(path)
    v={i:set() for i in range(1,N+1)}
    for x in E:
        a,b=x
        req(a in v and b in v and sq(a+b),"not a square-sum edge")
        v[a].add(b);v[b].add(a)
    req(all(len(s)==2 for s in v.values()),"not a spanning 2-factor")
    unseen=set(v);cs=0
    while unseen:
        cs+=1;todo=[unseen.pop()]
        while todo:
            nxt=v[todo.pop()]&unseen
            unseen-=nxt;todo.extend(nxt)
    return cs

def lift(P,moves):
    O=[];M=[]
    for token in G:
        if isinstance(token,int): O.append(token); continue
        c,rev=token
        f={x:25*x+(-1)**i*c for i,x in enumerate(P)}
        q=[f[x] for x in P]
        O.extend(q[::-1] if rev else q)
        for A,B in moves:
            conv=lambda E:{frozenset(f[x] for x in e) for e in E}
            M.append((conv(A),conv(B)))
    req(O[0]==1 and O[-1]==3 and set(O)==set(range(1,len(O)+1)),"wrong lift coverage")
    return tuple(O),M

def report(P,ij,twist,expected):
    pairs=[create(P,*k,twist) for k in ij]
    req(P[0]==1 and P[-1]==3 and set(P)==set(range(1,100)),"bad base")
    S=set(P[::2])
    K=set().union(*(set().union(*x) for x,_ in pairs))
    req(len(K)==8,"supports overlap")
    for old,new in pairs:
        req(edge(1,3) not in old and all(len(e&S)==1 for e in new),"port or parity lost")
    actual=[test(P,pairs,t) for t in ((),(0,),(1,),(0,1))]
    req(actual==expected,"wrong base circuit counts")
    Q,M=lift(P,pairs)
    pick=[0,1,2,len(M)-1]
    for S in ([],[0],[1],[0,1],pick,list(range(len(M))),list(range(0,len(M),2))):
        q=1+sum(expected[int(i in S)+2*int(i+1 in S)]-1 for i in range(0,len(M),2))
        req(test(Q,M,S)==q,"one-level factor law failed")
    print("99 base",ij,"circuit vector",actual)
    print("2487",len(M),"moves; single, double, mixed and all-instance composite toggles replayed")
    print("label pairs",[(sorted(map(sum,A)),sorted(map(sum,B))) for A,B in pairs])

def pairings(xs):
    if not xs:
        yield ()
        return
    for j in range(1,len(xs)):
        for x in pairings(xs[1:j]+xs[j+1:]):
            yield ((xs[0],xs[j]),)+x

def rank_binary(rows):
    a=list(rows)
    r=0
    while a:
        x=a.pop()
        if x:
            r+=1
            z=1<<(x.bit_length()-1)
            a=[y^x if y&z else y for y in a]
    return r

def abstract_control():
    counts=0
    for k in range(1,5):
        n=4*k;base=tuple(range(n))
        for chords in pairings(tuple(range(2*k))):
            for signs in product((0,1),repeat=k):
                rows=[sum(1<<j for j,(c,d) in enumerate(chords)
                      if (i==j and signs[i]) or (i!=j and (a<c<b<d or c<a<d<b)))
                      for i,(a,b) in enumerate(chords)]
                E=edges(base)
                for (i,j),sgn in zip(chords,signs):
                    U,V=create(base,2*i,2*j,sgn)
                    E=(E-U)|V
                D=list(range(n))
                def find(x):
                    while D[x]!=x:x=D[x]
                    return x
                for a,b in E:D[find(a)]=find(b)
                req(len({find(x) for x in range(n)})==1+k-rank_binary(rows),
                    "abstract interlacement mismatch")
                counts+=1
    print("abstract disk and band circuit formula verified on",counts,"signed pairings through four bands")

if __name__=="__main__":
    report(C,((9,66),(48,97)),1,[1,1,1,2])
    report(R,((0,56),(39,81)),0,[1,2,2,1])
    abstract_control()

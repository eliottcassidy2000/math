"""Wave 24 / THM-5162 (A): exact check of the block shape lemma and of the invertibility characterization.
Block J w on (H,E,lam) (blocks11.py convention). For w = O v Z with Z-positions S (1-based, 1 notin S, n in S):
  P = sum_{s in S} 2^(s-1), b = 1+P, R = sum_{s in S} 2^(n-s), Q = sum_{t in S} (sum_{s in S, s<t} 2^(t-1-s) + [t-1 in S]),
  M_w = [[b-Q, -1, b], [1, 0, 0], [2^n - R, 0, 2^n]],  det = 2^n.
J w is invertible iff w is empty or w = O v Z (w ending in O: E' = 0; w starting with Z: J Z factors through (H, lam))."""
import itertools
def step(s,m):
    H,E,L,h=s
    if m=='O': return (H,0,2*L,h)
    if m=='Z': return (H+L-E,h,2*L-h,h)
    return (H+L-E,H-E,H+L,H)
def blockmat(w):
    cols=[]
    for e in ((1,0,0),(0,1,0),(0,0,1)):
        s=step((e[0],e[1],e[2],0),'J')
        for ch in w: s=step(s,ch)
        cols.append(s[:3])
    return [[cols[c][r] for c in range(3)] for r in range(3)]
def det3(M):
    a,b,c=M[0];d,e,f=M[1];g,h,i=M[2]
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
nshape=0; ninv=0; nsing=0
for n in range(2,13):
    for v in itertools.product('OZ',repeat=n-2):
        w='O'+''.join(v)+'Z'; S=[i+1 for i,ch in enumerate(w) if ch=='Z']
        P=sum(2**(s-1) for s in S); b=1+P; R=sum(2**(n-s) for s in S)
        Q=sum(sum(2**(t-1-s) for s in S if s<t)+(1 if (t-1) in S else 0) for t in S)
        assert blockmat(w)==[[b-Q,-1,b],[1,0,0],[2**n-R,0,2**n]], w
        nshape+=1
    bs=sorted(1+sum(2**i for i,ch in enumerate('O'+''.join(v)+'Z') if ch=='Z') for v in itertools.product('OZ',repeat=n-2))
    assert bs==list(range(2**(n-1)+1,2**n,2)), n
for n in range(0,11):
    for w in itertools.product('OZ',repeat=n):
        w=''.join(w); d=det3(blockmat(w))
        inv=(n==0) or (w[0]=='O' and w[-1]=='Z')
        assert (d!=0)==inv, (w,d)
        ninv+=inv; nsing+=(not inv)
print("shape lemma verified for", nshape, "blocks J O v Z with |w| <= 12; b runs over the odd numbers of (2^(n-1), 2^n) once each")
print("invertibility characterization verified for all", ninv+nsing, "blocks J w with |w| <= 10:", ninv, "invertible,", nsing, "singular")

from replay_ready_overlaps import req,end,converse,strong,profile,pair,triple
from replay_binary_seam import glue

def K(b):
    if not b:return 0
    c=b//2
    return K(c)+(b%2)*(c+c%2)

def right(b):
    bits=list(map(int,bin(b)[2:]));n=len(bits)
    # Mark 0; bits record the in-cut backwards.
    return (sum((1-x)<<(i+1) for i,x in enumerate(bits[::-1])),)+tuple(
        (((1<<(n+1))-1)^((1<<(i+2))-1))|x for i,x in enumerate(bits[::-1]))

for A,v,H in pair+[(T,5,135+2*i) for i,T in enumerate(triple)]:
    h,L=profile(A,v)[1:]
    E=end(A)[-1][v];req(E==(3 if len(A)==5 else 15))
    if len(A)==7:req(end(converse(A))[-1][v]==30)
    for b in (2,4,6,10,12):
        B=right(b);req(strong(B) and end(B)[-1][0]==1)
        T=glue(A,v,B,0);req(strong(T))
        n=b.bit_length(); rev=int(bin(b)[2:][::-1],2)
        req(profile(T,v)==(H+L*b-h*K(b),h,(L<<n)-h*rev))
        req(end(T)[-1][v]==h)
    req(glue(glue(A,v,right(6),0),v,right(12),0)==glue(A,v,right((12<<3)+6),0))
print('all-cut ready corridor, complete probes, literal concatenation and end refresh pass')

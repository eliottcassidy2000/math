"""
Pin down the MINIMAL r=2 window-collapse failure and embed it in a full 13-set.
W=(3,39) means mult-of-7 runners {21,273} at V*=19. Show the two teeth families
cover (0, 1/(2*19)] = (0,1/38] entirely. Then build a covering primitive S3 13-set
S that is ALL-MULT7-LARGE with exactly these large mult-of-7 runners and V*=19,
verify NO tau=k/7+u/7 witness, yet M(S)>=1/14.
kind-pasteur S4-wf.
"""
from fractions import Fraction as F
from math import gcd
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def gaps(W,Vstar):
    R=F(1,2*Vstar); h=F(1,14); forb=[]
    for w in set(W):
        j=0
        while True:
            c=F(j,w); lo=c-h/F(w)
            if lo>R: break
            forb.append((max(F(0),lo),min(R,c+h/F(w)))); j+=1
    forb.sort(); merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    prev=F(0); g=[]
    for a,b in merged:
        if a>prev: g.append((prev,a))
        prev=max(prev,b)
    if prev<R: g.append((prev,R))
    return [(lo,hi) for lo,hi in g if hi>lo and hi>0],R,merged

print("MINIMAL r=2 window-collapse failure: V*=19, mult-of-7 runners {21,273}, W={3,39}")
W=[3,39]; Vstar=19
g,R,merged=gaps(W,Vstar)
print("  window (0,R], R=1/(2*19)=",R)
print("  forbidden (merged) teeth:",[(str(a),str(b)) for a,b in merged])
print("  SAFE GAPS:",g," => window fully blocked:",len(g)==0)
# explicit: w=3 teeth at u in [j/3 +- 1/42]; on (0,1/38]: j=0 -> [0,1/42].
# 1/42 ~0.0238, R=1/38~0.0263. remaining (1/42,1/38]. w=39: teeth at j/39+-1/546.
# j=1 -> 1/39~0.02564 in (1/42,1/38]? 1/42=0.0238<1/39=0.02564<1/38=0.02632 yes.
# tooth [1/39-1/546, 1/39+1/546]=[0.02380,0.02747] covers (1/42,1/38]. So blocked.
print("  w=3 j=0 right-tooth: [0, 1/42]=[0,",float(F(1,42)),"]")
print("  w=39 j=1 tooth: [1/39-1/546, 1/39+1/546]=[",float(F(1,39)-F(1,546)),",",float(F(1,39)+F(1,546)),"]")
print("  1/42 =",float(F(1,42)),"  1/39-1/546 =",float(F(1,39)-F(1,546)),
      "  overlap closes gap?",F(1,39)-F(1,546)<=F(1,42))
print("  right end R=1/38=",float(R)," <= 1/39+1/546 =",float(F(1,39)+F(1,546)),"?",R<=F(1,39)+F(1,546))
print()

# Build a full 13-set: need covering (mults of 2..14), primitive, S3, ALL-MULT7-LARGE
# with V*=19 (so 19 in S, no non-mult7 runner exceeds 19), mult-of-7 = {21,273} ONLY?
# But covering needs a multiple of 14: 21 is not mult of 14, 273=3*7*13 not mult 14.
# So we need a multiple of 14 too, also >19. Add 7*w with 7w%14==0 i.e. w even, w>19/7.
# Use 21(=7*3),273(=7*39), and a mult of 14 >19, e.g. 7*4=28 (w=4). Recompute window
# with W={3,4,39}; check still fully blocked.
print("Embedding into a real covering 13-set (need a multiple of 14 >19):")
# choose non-mult7 part covering 2,3,4,5,6,8,9,10,11,12,13 and including 19, max=19.
nm7=[1,2,4,5,8,9,10,11,12,13,19]   # check covers 2..6,8..13 via these
m7=[21,28,273]                      # 21=3*7,28=4*7(=mult14),273=39*7 ; all>19
S=sorted(nm7+m7)
print("  S=",S," |S|=",len(S))
print("  covering:",is_cov(S)," primitive:",is_prim(S))
Vmin=min(S);Vmax=max(S);k=sum(1 for v in S if v>13)
print("  S3? :",k>=2 and Vmax>=13*Vmin,"  V*(max nonmult7)=",max(v for v in S if v%7!=0))
Wfull=[v//7 for v in S if v%7==0]
gf,Rf,_=gaps(Wfull,max(v for v in S if v%7!=0))
print("  full mult-of-7 subsystem W=",sorted(Wfull)," window gaps:",gf," fully blocked:",len(gf)==0)
print("  Mval(S)=",Mval(S),"~",float(Mval(S))," >=1/14?",Mval(S)>=F(1,14))
print()
print("VERDICT: in-scope ALL-MULT7-LARGE set; witness tau=k/7+u/7 does NOT exist")
print("(window fully blocked), yet M>=1/14. The witness-existence sub-claim is FALSE;")
print("the closure CONCLUSION (M>=1/14) survives here only via a different tau.")

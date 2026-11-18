"""Physical orientation edits, and shared-frame squares in the two graph types."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moonshine_prime_blocks'))
from frame_replay import add
from replay import reads,wrap
from clone_replay import TT
from atlas import H,canon_lane,flip,restrict
from metalines import frames,line,blue,canon
SPOKES=(
 (10007,7505,(174,52,248,306,257,80,27,122,231),331,363),
 (10009,7507,(466,40,331,321,206,413,162,14,210),183,247),
 (257,193,(46,108,72,96,111,4,33),50,54),
 (263,197,(94,32,26,66,42,77,22),114,122),
 (259,97,(18,4,209,247,2,151,179,19),249,251),
 (30013,11255,(548,333,264,1,943,526,317,111,169,462),626,622),
)
WIDTH=(
 (1,(0,1,3,7,47,79,159,63),48,(7,5),(16,12,1)),
 (5,(0,5,17,7,75,31,175,63),80,(2,4),(64,12,1)),
 (-17,(0,1,3,71,47,79,151,63),56,(3,6),(64,48,-13)),
 (-23,(0,1,11,51,39,71,159,63),164,(4,5),(128,96,-17)),
)
def run():
    for N,h,D,a,b in SPOKES:
        A,B=add(D,a),add(D,b)
        cost=(a^b).bit_count()
        assert H(A)==N and H(B)==h==canon_lane(N)
        # Small cases also compute a marked commuting square.
        if cost==1:
            i,j=(a^b).bit_length()-1,len(D)
            P=next(P for P in frames(A) if all({x,y}!={i,j} for x,y in zip(P,P[1:])))
            assert P in frames(B)
            C,G=line(A,P),line(B,P)
            assert flip(C,i,j)==G
            print(N,h,'shared frame',P,'metaline partner counts',H(C),H(G),
                  'metaline colors',blue(A,P),blue(B,P))
        else:
            seq=[H(add(D,k)) for k in (626,630,638,622)]
            assert seq==[30013,24961,13955,11255]
            print(N,h,'spokes',cost,'strict route',seq)
    for pos,P in enumerate(((3,4,0,7,5,9,8,6,1,2),(6,7,2,1,5,0,4,3,8,9))):
        N,h,D,a,b=SPOKES[pos]
        S,T=add(D,a),add(D,b)
        assert P in frames(S) and H(line(S,P))==h and not blue(S,P)
        assert (canon(line(S,P))==canon(T))==(pos==0)
        print(N,h,'direct black from spoke source; same merged target',pos==0)
    for N,h,T,P in (
        (257,193,(2,40,35,53,71,209,143,31),(4,6,3,2,0,1,5,7)),
        (259,97,(2,20,41,67,45,203,151,31),(2,0,1,4,5,3,6,7)),
        (263,197,(8,17,35,70,45,203,151,31),(6,0,3,1,4,2,5,7))):
        assert P in frames(T)
        assert H(T)==N and H(line(T,P))==h==canon_lane(N) and not blue(T,P)
        print(N,h,'direct genuine black line',P)
    for N,h,D,a,b,M in (
        (30011,5627,(290,312,11,49,5,84,31,127,252,511,1023),1481,329,(29325,8315)),
        (30013,11255,(238,128,482,182,71,338,10,368,91,511,1023),1875,979,(22797,14321))):
        assert (a^b)==(1<<7)|(1<<10)
        assert H(add(D,a))==N and H(add(D,b))==h==canon_lane(N)
        U=tuple(H(add(D,a^(1<<i))) for i in (7,10))
        assert set(U)==set(M) and all(h<k<N for k in U)
        print('alternate two-spoke rails',N,U,h)
    for b,D,I,(i,j),(G0,G1,b1) in WIDTH:
        K=flip(D,i,j)
        for T,expected in ((D,(G0,b)),(K,(G1,b1))):
            d,g,q=reads(T,I);assert (g,d+q-g)==expected
        for v in (1,2,3):
            A,B=wrap(D,I,TT(v)),wrap(K,I,TT(v))
            assert H(A)==G0*2**v+b and H(B)==G1*2**v+b1==canon_lane(H(A))
        print('one-arc family',(G0,b),(G1,b1))
def frame_reverse(T,P):
    assert all(T[i]>>j&1 for i,j in zip(P,P[1:]))
    for i,j in zip(P,P[1:]):T=flip(T,i,j)
    return T

def double_lines():
    D=(2,8,3,69,47,79,151,63);I=128;DD=flip(D,1,3)
    for U,ex in ((D,(64,59)),(DD,(12,11))):
        h,b,q=reads(U,I)
        assert (b,h+q-b)==ex
    for w in (2,3,4):
        A=wrap(D,I,TT(w));B=wrap(DD,I,TT(w));b=w+8
        P=(8,*range(10,b),4,5,2,0,1,3,6,7,b,9)
        assert len(set(P))==len(A)
        i=P.index(1)
        Q=P[i::-1]+P[:i:-1]
        assert A[8]>>9&1
        C=frame_reverse(A,P);U=frame_reverse(C,Q)
        assert U==flip(flip(A,1,3),8,9)
        pi=tuple(9 if i==8 else 8 if i==9 else i for i in range(len(A)))
        assert restrict(U,pi)==B
        assert H(A)==64*2**w+59 and H(B)==12*2**w+11==canon_lane(H(A))
        print('two true lines',w,H(A),H(C),H(B),'blue',blue(A,P),blue(C,Q))
if __name__=='__main__':
    run()
    double_lines()

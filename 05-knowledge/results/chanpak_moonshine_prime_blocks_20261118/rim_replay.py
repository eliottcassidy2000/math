"""Universal cyclic lift replay for count-changing one-arc edits."""
from physical_edges import SPOKES,add,frames,frame_reverse,restrict,TT,H
from clone_replay import inflate
from replay import profile
def universal_lift():
    for N,h,D,a,b in SPOKES[:2]:
        X,Y=add(D,a),add(D,b)
        i,j=(a^b).bit_length()-1,len(D)
        P=next((P for P in frames(X) if any({e,f}=={i,j} for e,f in zip(P,P[1:]))),None)
        if P is None:
            X,Y=Y,X
            P=next(P for P in frames(X) if any({e,f}=={i,j} for e,f in zip(P,P[1:])))
        for w in (2,3):
            A=inflate((2,4,1),{0:TT(w),1:X});B=inflate((2,4,1),{0:TT(w),1:Y})
            R=(0,*range(2,w),*(w+j for j in P),len(A)-1,1)
            l=next(l for l,(e,f) in enumerate(zip(R,R[1:])) if {e,f}=={w+i,w+j})
            S=R[l::-1]+R[:l:-1]
            C=frame_reverse(A,R);T=frame_reverse(C,S)
            pi=(1,0)+tuple(range(2,len(A)))
            assert restrict(T,pi)==B
            gx,hy=profile(X),profile(Y)
            x=2*gx[0]+(2**w-1)*(gx[1]-gx[0])
            y=2*hy[0]+(2**w-1)*(hy[1]-hy[0])
            assert H(A)==x and H(B)==y
            print('universal padded route for',N,h,'width',w,x,H(C),y)
if __name__=='__main__':universal_lift()

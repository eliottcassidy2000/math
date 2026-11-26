"""Uniform representative-free lift of the published mod-16 guarded descent."""
import sys

def clocks(max_depth,N=1<<71):
    d=p=1; j=0; maximum=0; argmax=None; active=0
    for k in range(1,max_depth+1):
        d<<=1
        if p*3<d:
            p*=3;j+=1
        # p=largest power of 3 strictly below d; exact first crossing implies p>=d/2.
        assert p<d<=3*p
        if 2*p<d:continue
        active+=1
        den=3*(d-p)
        assert j*p < den*N
        threshold=j*p//den
        if threshold>maximum:maximum=threshold;argmax=(k,j)
    return active,maximum,argmax

def residues(k):
    good=[]
    for r in range(1<<k):
        p,x=1,r
        for i in range(1,k+1):
            if x&1:p*=3;x=(3*x+1)//2
            else:x//=2
            if p<1<<i:break
        else:good.append(r)
    return good

if __name__=='__main__':
    w=residues(4)
    assert w==[7,11,15]
    print("mod 16 unresolved by a coefficient drop:",w)
    for depth in ([int(sys.argv[1])] if len(sys.argv)>1 else [100000]):
        print("depth, possible first-crossing clocks, largest floor threshold, attaining clock",
              depth,*clocks(depth))

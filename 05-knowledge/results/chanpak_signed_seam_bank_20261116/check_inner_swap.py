from check_seam_cube import actual, coordinates, tournament, value
def check():
    A=tournament(7,2096861);B=tournament(6,149)
    def cover(C):
        s={0}
        for _ in C:s|={v for u in s for v in range(len(C)) if C[u][v]}
        return len(s)==len(C)
    assert all(cover(D) and cover(tuple(zip(*D))) for D in (A,B))
    for i,(v,w) in enumerate(((110,61),(86,53))):
        p=coordinates(A,v);q=coordinates(B,w)
        h1,h2=p[0],q[0]
        assert (h1,h2)==(55,31)
        F=(2*p[1]-p[3],2*q[1]+q[2])
        assert F==((1107,880),(1188,820))[i]
        assert F[0]*F[1]==974160
        d=(2*p[2]-p[6],-2*q[3]-q[6])
        assert d==((846,-702),(548,-534))[i]
        assert h1*h2+h2*d[0]+h1*d[1]==-10679+2*i
        for t in (1,2):
            H=974160*2**t-10679+2*i
            assert value(p,q,2,2,2**t,0)[0]==H
            assert actual([A,B],[v,w],(1,t,1))==H
            print(t,v,w,F,d,H)
    assert (110^86).bit_count()+(61^53).bit_count()==4
if __name__=='__main__':check()

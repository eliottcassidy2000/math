"""Literal small regimes of the full offset-59 signed lane."""
from atlas import H,canon_lane
from metalines import line,blue,canon,frames
from physical_edges import double_lines
D=(
 (1,23,(2,64,3,7,47,79,157,63),(0,1,6,7,4,5,3,2)),
 (6,23,(2,16,11,35,77,87,143,63),(0,1,4,6,7,5,2,3)),
 (7,35,(4,17,34,7,109,203,143,31),(0,2,5,3,1,4,6,7)),
)
detours=(
 (2,171,47,(16,5,17,7,74,31,175,63),(7,5,1,2,0,4,6,3),
 (0,65,11,67,15,31,181,63),(3,1,6,7,5,4,2,0)),
 (3,105,25,(0,5,25,67,43,79,151,63),(1,2,4,3,6,7,5,0),
 (0,1,19,23,99,79,143,63),(3,2,4,6,7,5,1,0)),
 (5,65,17,(2,4,9,51,71,87,143,63),(4,6,7,0,1,2,3,5),
 (0,1,11,19,39,79,159,63),(2,3,4,5,6,7,1,0)),
)
def run():
    for m,y,T,P in D:
        assert P in frames(T)
        assert H(T)==2**m+59 and H(line(T,P))==y==canon_lane(H(T))
        assert not blue(T,P)
        print(m,H(T),y)
    for m,middle,y,A,P,B,Q in detours:
        assert P in frames(A) and Q in frames(B)
        assert H(A)==2**m+59 and H(B)==y==canon_lane(H(A))
        U,V=line(A,P),line(B,Q)
        assert H(U)==middle and canon(U)==canon(V)
        assert not blue(A,P) and not blue(B,Q)
        print(m,H(A),middle,y)
    assert canon_lane(75)==7
    double_lines()
if __name__=='__main__':run()

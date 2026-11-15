from reference import profiles, tt, join, close

B=[18,8,3,5,14]; O=6; P=3
assert profiles(B,O)==((5,32,13,6),104)
assert profiles(B,P)==((7,32,13,6),104)
for t in range(1,10):
    x=profiles(close(B,O,tt(t)),0)
    y=profiles(close(B,P,tt(t)),0)
    assert (x[0][2],y[0][2])==(5+(32<<t),7+(32<<t))
    assert x[1]-y[1]==2*(-4*2**t+3*3**t-3)
C=[2,4,1]; A=join(C,B); D=join(A,C)
for mark,r in ((O,5),(P,7)):
    assert profiles(D,mark<<3)==((9*r,32*3*12,117,0),104*12**2)
print("same-source fork; nine literal closings; one cyclic context")
# Nontransitive tail, with ordered nonempty three-cover mass 6.
qx=profiles(close(B,O,C),0)
qy=profiles(close(B,P,C),0)
assert (qx[0][2],qy[0][2])==(5*3+32*12,7*3+32*12)
assert qx[1]-qy[1]==2*(-4*12+3*33-3*3)
print("nontransitive common tail also changes quotient")
for marks,(a,b,c) in ((O,(76,132,44)),(P,(80,129,41))):
    for t in range(1,5):
        assert profiles(close(B,marks,tt(t)),0)[1]==2*(a*2**t+b*3**t-c)

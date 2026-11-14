from fractions import Fraction as F
suffixes=['','O','OO','OZ','OOO','OOZ','OZO','OZZ','OOOZ','OZOZ','OZZZ']
def step(x,w):
    A,E,L=x; H=A+E; h=H
    H,E,L=H+L-E,H-E,H+L
    for c in w:
        if c=='O': H,E,L=H,0,2*L
        else:H,E,L=H+L-E,h,2*L-h
    return (H-E,E,L)
def det(M,la=0):
    a,b,c,d,e,f,g,h,i=M;a-=la;e-=la;i-=la
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
upper=[434,154,154,40,154,11,21,13,3,4,4]
for w,u in zip(suffixes,upper):
    cols=[step(x,w) for x in [(1,0,0),(0,1,0),(0,0,1)]]
    mat=tuple(cols[j][i] for i in range(3) for j in range(3))
    assert all(x>=0 for x in mat) and det(mat,2)!=0
    for (a,e,l),k in zip(cols,[0,1,2]):
        assert a>=int(k==2) and l>=1 and a+e>=int(k!=1) and l>=a+e
    a,e,l=step((1,0,F(7,4)),w); c=min(a,l/F(7,4))
    assert c>abs(det(mat))
    assert F(u,1000)**20*c**37>1
    print('J'+w, 'c =',c,'det =',det(mat))
assert sum(upper)==992
print('All eleven positive cone, growth, and determinant gates verified.')

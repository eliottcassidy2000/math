from verify_physical_refill import require, responses, strong, prime

M=[6,12,80,37,11,23,59]
N=[2,12,145,165,203,151,175,3]
L=[2,20,9,35,45,7]
R=[2,4,9,51,135,215,159,15]
h,p,w=responses(M)
j,q,z=responses(N)
a,b,c=responses(L)
x,y,v=responses(R)
require((h,w,j,q,a,c,x,y)==(101,8823,131,2864,29,1845,59,1712))
require(all(strong(T) for T in (M,N,L,R)))
count=[60*c*y*2**(5-k)*3**k+a*x for k in range(6)]
require(count[0]==60*w*q*4+h*j==6064590511==7*866370073)
require(prime(866370073))
require(all(2*v==3*u-1711 for u,v in zip(count,count[1:])))
require((q*4&-(q*4),y*32&-(y*32))==(64,512))

def graph(k):
    tail=[20,21,22,23,24]
    order=[30]+list(range(8))+tail[:5-k]+[31]+tail[5-k:]+list(range(8,14))+[32,33,34,35,36]
    arcs={(a,b) for i,a in enumerate(order) for b in order[i+1:]}
    for base,rows in ((0,R),(8,L)):
        for i in range(len(rows)):
            for j in range(i+1,len(rows)):
                if rows[j]>>i&1:
                    arcs.remove((base+i,base+j))
                    arcs.add((base+j,base+i))
    for a,b in ((30,32),(30,33),(31,35),(31,36)):
        arcs.remove((a,b));arcs.add((b,a))
    return arcs

for k in range(5):
    old,new=graph(k),graph(k+1)
    require(old-new=={(24-k,31)} and new-old=={(31,24-k)})
print('PASS: four actual modules, 7p rebase, five physical transfers')

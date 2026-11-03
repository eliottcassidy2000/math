def require(x):
 if not x: raise AssertionError("strong-family check failed")
M=75
Dh=9*M; Dg=225*M; hp=450; gp=2250
require(pow(4,hp,Dh)==1 and 6*hp%Dh==0)
require(pow(4,gp,Dg)==1 and pow(9,gp+2,Dg)==pow(9,2,Dg))
require(150*gp%Dg==0 and 60*gp%Dg==0)
H={((10*pow(4,n,Dh)+6*n-1)%Dh)//9 for n in range(1,hp+1)}
G={((25+150*n+(668+60*n)*pow(4,n,Dg)+432*pow(9,n,Dg))%Dg)//225 for n in range(1,gp+2)}
missing=sorted(set(range(M))-(H|G))
require(missing==[9,12,24,27,54,57,69,72])
print('h residues',len(H),'g residues',len(G),'missing union',missing)

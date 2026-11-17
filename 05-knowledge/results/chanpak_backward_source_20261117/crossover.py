from reader.replay import pair,decode,literal,phi
def ck(z):
    if not z:raise ArithmeticError("crossover")
h,p,j=197,3728,25
ck(pow(2,j,h)*p%h==2)
ck(pow(2,196,h)==1 and all(pow(2,196//l,h)!=1 for l in (2,7)))
ck(all(pow(2,i,h)*p%h!=2 for i in range(j)))
N,u=pair(h,p,j,0)
ck((N,u)==(158860293255395167777,42612739620732928))
ck(u.bit_length()+9==65 and phi(u+(1<<j),h,p)==N+2)
ck((j,0) in decode(h,p,N))
print("PASS original source to new reader crossover and its exact eligibility period")

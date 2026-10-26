\\ Independent modular data for the audit (PARI/GP mf package; no eta/theta products used).
\\ f = 8.4.a.a (weight 4, level 8, trivial character): the unique newform of S_4^new(Gamma0(8)).
\\ g = 32.3.c.a (weight 3, level 32, character chi_-4 = Conrey 32.31): one form of the Galois orbit; y^2 = -1.
default(parisize, 200000000);
N = 2200;
mf8 = mfinit([8,4],0); B8 = mfeigenbasis(mf8);
if(#B8 != 1, error("8.4 new space not 1-dim"));
c8 = mfcoefs(B8[1], N);
mf32 = mfinit([32,3,Mod(31,32)],0); B32 = mfeigenbasis(mf32);
if(#B32 != 1, error("32.3 chi_-4 new space not a single orbit"));
c32 = mfcoefs(B32[1], N);
\\ check chi(31,.) is chi_-4 on small n
for(n=1,31, if(n%2==1, if(chareval(znstar(32,1),znconreylog(znstar(32,1),31),n,[-1,2]) != kronecker(-4,n), error("char")) ));
print("dims: S4new(8) = ", mfdim(mf8), ", S3new(32,chi_-4) = ", mfdim(mf32));
forprime(p=3, N, a = c8[p+1]; b = lift(c32[p+1]); re = polcoeff(b,0); im = polcoeff(b,1); print("p=", p, " a8=", a, " bre=", re, " bim=", im));
quit;

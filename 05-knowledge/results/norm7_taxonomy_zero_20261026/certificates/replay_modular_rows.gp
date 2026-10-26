\\ Independent fresh PARI reconstruction of the13-prime modular columns.
\\ Emits q,a84,a86,a46,norm_b32; no stored coefficient table is read.
N=103;
A=mfcoefs(mfeigenbasis(mfinit([8,4,1],0))[1],N);
C=mfcoefs(mfeigenbasis(mfinit([8,6,1],0))[1],N);
D=mfcoefs(mfeigenbasis(mfinit([4,6,1],0))[1],N);
B=mfcoefs(mfeigenbasis(mfinit([32,3,-4],0))[1],N);
P=[7,11,19,23,31,43,47,59,67,71,79,83,103];
for(j=1,#P,p=P[j];b=lift(B[p+1]);b2=if(type(b)=="t_POL",polcoef(b,0)^2+polcoef(b,1)^2,b^2);print(p," ",A[p+1]," ",C[p+1]," ",D[p+1]," ",b2));
quit;

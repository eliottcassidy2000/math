/* [a+bw]/[c+dw], coefficients -2..2, nonzero determinant; exact finite probe. */
f=y^3+2*y^2-4*y-4; nf=nfinit(f); w=Mod(y,f);
seen=Map(); tested=0; found=0; nontrivial=0; by=matrix(4,2);
{
for(a=-2,2,
 for(b=-2,2,
  for(c=-2,2,
   for(d=-2,2,
    if(a*d-b*c!=0,
     r=(a+b*w)/(c+d*w);
     if(!mapisdefined(seen,r),
      mapput(seen,r,1);
      for(ti=1,4,
       t=ti/4; DD=t^2+1-2*t*(r^2-11)/(r^2+11);
       for(z=1,2,
        disc=DD*(4-DD)/if(z==1,11,55); rr=nfroots(nf,x^2-disc); tested++;
        if(#rr,
         if(#rr!=2,error("unexpected zero or inseparable contact"));
         for(ii=1,2,if(rr[ii]^2!=disc,error("unverified square root")));
         found++;by[ti,z]++;
         if(ti!=4 || z!=1,nontrivial++)
        )
       )
      )
     )
    )
   )
  )
 )
);
}
if(tested!=1792 || found!=224 || nontrivial!=0,error("finite probe changed"));
print("{\"status\":\"FINITE-EXACT; bounded Mobius bank only\",\"distinct_parameters\":",#seen,",\"tests\":",tested,",\"hits\":",found,",\"nontrivial_hits\":",nontrivial,",\"hits_by_separation_and_class\":[",by[1,],",",by[2,],",",by[3,],",",by[4,],"]}");

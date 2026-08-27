\\ Independent PARI/GP audit for the W=0 CM ideal conventions.
\\ Run: gp -fq 04-computation/jc23_w0_prime_ideals_independent_audit_thm4241.gp

x='x;
pol=x^4-x^2+1;
K=bnfinit(pol);
z=Mod(x,pol);
w=z^4;
T=z^5;
p2=1-z^3;
p3=z+z^3;
sqrt3=z+z^-1;
trace_q_1_plus_T=trace((3+sqrt3)*(1+T)*(1+z^-5))/2;
if(lift(sqrt3^2)!=3,error("sqrt3 embedding mismatch"));
if(trace_q_1_plus_T!=6,error("3+sqrt3 trace convention mismatch"));

print("W=0 PARI ideal audit");
print("polynomial=",pol);
print("discriminant=",K.disc);
print("class_number=",K.no);
print("omega=",lift(w)," T=",lift(T)," T^2=",lift(T^2));
print("sqrt3=",lift(sqrt3)," scalar_trace_q(1+T)=",trace_q_1_plus_T);
print("pi2=",lift(p2)," norm=",norm(p2)," pi2^2/2=",lift(p2^2/2));
print("pi3=",lift(p3)," norm=",norm(p3)," pi3^2/3=",lift(p3^2/3));
print("prime_decomposition_2=",idealprimedec(K,2));
print("prime_decomposition_3=",idealprimedec(K,3));
print("minkowski_bound_18_over_pi2=",18/Pi^2);
print("result=PASS");
quit;

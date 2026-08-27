\\ Independent PARI/GP theta audit for the normalized full Hom lattice.
\\ qfrep counts nonzero vectors modulo the pair {x,-x}; double its entries.
\\ Run: gp -fq 04-computation/jc23_w0_full_theta_independent_audit_thm4241.gp

Q=[4,-2,0,0,0,0,0,0;-2,4,0,0,0,0,0,0;0,0,6,-3,-3,0,-3,3;0,0,-3,6,3,-3,0,-3;0,0,-3,3,6,-3,3,-3;0,0,0,-3,-3,6,0,3;0,0,-3,0,3,0,4,-2;0,0,3,-3,-3,3,-2,4];

half_counts=qfrep(Q,42);
raw_counts=2*Vec(half_counts);
expected=[0,0,0,60,0,72,0,324,0,864,0,276,0,2592,0,2868,0,1152,0,5832,0,9504,0,2556,0,15552,0,15456,0,6480,0,22356,0,36288,0,7836,0,49248,0,44280,0,16992];

if(raw_counts!=expected,error("theta mismatch"));
if(matdet(Q)!=2916,error("quadratic determinant mismatch"));

print("W=0 independent full-theta audit");
print("quadratic_matrix_determinant=",matdet(Q));
print("qfrep_pair_counts_through_42=",half_counts);
print("raw_counts_through_42=",raw_counts);
print("degree34_raw=",raw_counts[34]," degree42_raw=",raw_counts[42]);
print("result=PASS");
quit;

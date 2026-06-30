print("(B) q = n^2-n+1 = Phi_6(n) (6th cyclotomic poly); n^3 = -1 mod q, n has ORDER 6:")
for n in [4,7,14]:
    q=n*n-n+1; p=[pow(n,i,q) for i in range(7)]
    print(f"   n={n:>2}: q={q}=Phi6({n}); n^0..n^6 mod q = {p};  n^3=-{q-p[3]}, order6={p[6]==1}")
print("   6 powers {1,n,n-1,-1,-n,1-n} = the runners' landing at witness t=n/q (the 6=phi(14) binding).")
print("\n(C) key congruence (n-1)n*n = -n mod Phi_6(n) (the killer binds at distance n/q):")
for n in [4,7,14,20,100]:
    q=n*n-n+1; lhs=((n-1)*n*n)%q
    print(f"   n={n:>3}: (n-1)n^2 mod {q} = {lhs} = -{q-lhs}; equals -n: {lhs==(q-n)%q}")

"""
a(x)=x/2 (DIVISION), b(x)=x+1 (ADDITION). Parity-routed:
  f(x) = b(x) if x even else x        = the ODD one of {x, x+1}
  g(x) = a(x) if x even else a(b(x))  = HALF the EVEN one of {x,x+1} = ceil(x/2)
=> f*g = x(x+1)/2 = T_x (triangular). The '2' in T_x is absorbed by whichever of x,x+1 is even -- parity routing.
Think recursively: g is the 2-adic descent (apply repeatedly -> 1); f is the odd core; b is the observer +1.
"""
a=lambda x:x//2
b=lambda x:x+1
f=lambda x: b(x) if x%2==0 else x
g=lambda x: a(x) if x%2==0 else a(b(x))
print("x : f(x)=odd one  g(x)=ceil(x/2)  f*g   T_x=x(x+1)/2")
for x in range(1,11):
    print(f"{x:>2} :   {f(x):>3}          {g(x):>3}        {f(x)*g(x):>4}   {x*(x+1)//2:>4}   match={f(x)*g(x)==x*(x+1)//2}")
print(f"\nf*g = T_x for x=1..40: {all(f(x)*g(x)==x*(x+1)//2 for x in range(1,41))}")
print()
# RECURSION: g is the descent. Apply g repeatedly = halve to 1; collect the f's (odd cores) = 2-adic structure
def descent(x):
    cores=[]; cur=x
    while cur>=1:
        cores.append(("odd-core f",f(cur)) if False else f(cur))
        if cur==1: break
        cur=g(cur)
    return cores
print("RECURSIVE descent g^k(x) (halve to 1), collecting f = odd one at each level:")
for x in [13,14,20]:
    seq=[]; cur=x
    while True:
        seq.append((cur,f(cur)))
        if cur<=1: break
        cur=g(cur)
    print(f"  x={x}: (value, f=odd) chain = {seq}")
print("  => g = the 2-adic DESCENT (renormalization S->E/2); f = the odd core at each level; depth ~ log2(x).")
print()
# the three dualities as the operation group <a,b> and inverses
print("THE THREE DUALITIES (the generators + inverses = the dyadic affine group):")
print("  EVEN/ODD      = parity = which of {x,x+1} carries the factor 2 (the f/g routing); the descent level.")
print("  ADD/DIVISION  = b (x+1) vs a (x/2): the observer's +1 (baseline) vs the 2-adic halving (descent).")
print("  POSITIVE/NEG  = b (x+1) vs b^{-1} (x-1): the antipode; killer = -1 mod Phi6, complement = reversal.")
print("  <a,b> = the dyadic/Stern-Brocot tree -- the SAME tree as the observer's Farey escape [0;n-1,n].")
print()
print("THE PROJECT, FUNCTIONALLY DECOMPOSED:")
print("  TRIANGLE (staircase delta) = f*g = T_x  -- 'everything is the triangle' = the parity-routed product.")
print("  DESCENT  (renormalization) = recursive g = a on the even part -- the 2-adic flow cusp->doublet.")
print("  OBSERVER (the +1 baseline) = b -- Redei's inshat=1+2..., LRC escape n/(n(n-1)+1), the Farey hair.")
print("  Phi_6(n) = n(n-1)+1 = b(a^{-1}...) : the antipodal +1 (b) on the pronic n(n-1) (=lcm killer).")
print("  => a (divide) builds the descent; b (add 1) builds the observer; parity routes; product = the triangle.")

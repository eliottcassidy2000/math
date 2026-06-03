# Apex-lift certificate sheaf for n=14 (q=7), concrete over F_7.
# SITE: the lane set = residues mod 7 = {0,1,..,6}. Antipodal involution sigma: a -> -a mod 7.
#   nonzero lanes pair up {1,6},{2,5},{3,4}; lane 0 is sigma-FIXED = the apex seam.
# STALK at lane a: a runner forbids an affine line in the plane K^2 (covector (a_i,b_i), target c_i).
# A LOCAL CERTIFICATE over a sub-site U = a point of K^2 avoiding every forbidden line indexed by U.
# SHEAF GLUING:  Cert(U union V) = Cert(U) intersect Cert(V).  Global section = nonempty total intersection.
q=7
K=list(range(q))
P=[(x,y) for x in K for y in K]               # the plane A^2 over F_7
def line(a,b,c): return {(x,y) for (x,y) in P if (a*x+b*y)%q==c}   # forbidden set (a line, or all/none)
def cert(forbidden_sets):                      # certificate locus = avoid every forbidden set
    bad=set().union(*forbidden_sets) if forbidden_sets else set()
    return [p for p in P if p not in bad]

# --- one runner per nonzero lane a, tight tuple: covector (a,1), target c=-a (slope collapses to -1) ---
transverse=[ (a,1,(-a)%q) for a in range(1,q) ]          # 6 transverse lanes
F_trans=[line(a,b,c) for (a,b,c) in transverse]
print("STALK sizes (transverse lanes):", [len(line(a,b,c)) for (a,b,c) in transverse], " each = q =",q)
print("certificate count per single transverse lane:", len(cert([F_trans[0]])), "= q^2 - q =", q*q-q)

# --- GLUING law check: Cert(U cup V) == Cert(U) cap Cert(V) ---
import itertools, random
random.seed(0); ok=True
for _ in range(200):
    idx=list(range(6)); random.shuffle(idx); k=random.randint(1,5)
    U,V=idx[:k],idx[k:]
    cu=set(cert([F_trans[i] for i in U])); cv=set(cert([F_trans[i] for i in V]))
    cuv=set(cert([F_trans[i] for i in U+V]))
    if cuv!=cu&cv: ok=False; break
print("GLUING law  Cert(U cup V)=Cert(U) cap Cert(V):", ok)

# --- global section over the 6 transverse lanes (no apex yet) ---
glob_trans=cert(F_trans)
print("global section over 6 transverse lanes: size", len(glob_trans), "(nonempty:", len(glob_trans)>0, ")")

# --- THE APEX. Add the apex runner. Its speed is a multiple of q => covector (0,0), target 0 => forbids WHOLE plane
apex_unlifted=line(0,0,0)
print("\nAPEX stalk in A^2: forbidden-set size =", len(apex_unlifted), "= whole plane", q*q, "-> local cert set EMPTY:", len(cert([apex_unlifted]))==0)
glob_with_apex=cert(F_trans+[apex_unlifted])
print("GLOBAL section WITH unlifted apex: size", len(glob_with_apex), "-> sheaf has NO global section (gluing obstructed at apex seam)")

# --- THE LIFT. r/p lift: add a 3rd coordinate u with nonzero apex coeff d. Apex (0,0,d), d!=0, forbids a HYPERPLANE.
P3=[(x,y,u) for x in K for y in K for u in K]
def line3(a,b,d,c): return {(x,y,u) for (x,y,u) in P3 if (a*x+b*y+d*u)%q==c}
def cert3(fs):
    bad=set().union(*fs) if fs else set()
    return [p for p in P3 if p not in bad]
F3_trans=[line3(a,b,0,c) for (a,b,c) in transverse]      # transverse lanes lifted (u-coeff 0)
apex_lifted=line3(0,0,1,0)                                # d=1 != 0
print("\nAPEX stalk in A^3 (lifted, d=1): forbidden size =", len(apex_lifted), "= q^2 =", q*q, "(codim 1, NOT whole space", q**3,") -> local cert nonempty:", len(cert3([apex_lifted]))>0)
glob_lift=cert3(F3_trans+[apex_lifted])
print("GLOBAL section in the LIFT (6 transverse + lifted apex): size", len(glob_lift), "-> RESTORED, nonempty:", len(glob_lift)>0)

# --- ANTIPODAL sigma-equivariance: sigma(x,y)=(-x,-y). Cert locus is sigma-stable; apex lane is the sigma-fixed lane.
def sig(p): return ((-p[0])%q,(-p[1])%q)
stable = set(map(sig,glob_trans))==set(glob_trans)
print("\nANTIPODAL sigma-equivariance: transverse global section is sigma-stable:", stable)
print("apex lane a=0 is sigma-fixed (a=-a mod 7):", (0)==((-0)%q), "; nonzero lanes pair:", [(a,(-a)%q) for a in range(1,q)])

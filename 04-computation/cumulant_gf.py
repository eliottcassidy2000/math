#!/usr/bin/env python3
"""
cumulant_gf.py - Explore generating function for tournament cumulants.

For the transitive tournament, the CGF is related to Bernoulli numbers:
  K_trans(t) = sum_{k>=1} kappa_{2k}(trans) * t^{2k} / (2k)!
             = sum_{k>=1} (-1)^{k-1}(n+1)|B_{2k}|/(2k) * t^{2k}/(2k)!

The generating function of Bernoulli numbers is:
  t/(e^t - 1) = sum_{k>=0} B_k * t^k / k!
  => sum B_{2k} * t^{2k} / (2k)! = t/(e^t-1) - 1 + t/2

So K_trans(t) involves a modified version of t/(e^t-1).

For a general tournament T, the CGF is:
  K_T(t) = K_trans(t) + delta_K(t)
where delta_K encodes the cycle corrections.

The KEY question: is delta_K(t) itself a nice function of the tournament invariants?

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction
from math import comb

# The cumulant corrections are:
# delta_kappa_2 = (2/C(n,2)) * t3
# delta_kappa_4 = (2/C(n,4)) * (t5 + 2*a2) - 3*(2*t3/C(n,2))^2
# delta_kappa_6 = (2/C(n,6)) * t7 + cross_terms

# The LEADING new-cycle term in delta_kappa_{2k} is (2/C(n,2k)) * t_{2k+1}.
# Let me define: tau_{2k+1} = 2*t_{2k+1} / C(n, 2k).
# Then:
#   delta_kappa_2 = tau_3
#   delta_kappa_4 = tau_5 + 2*tau_alpha_2 - 3*tau_3^2
#   delta_kappa_6 = tau_7 + cross_terms

# If we write the full CGF correction as:
#   delta_K(t) = sum_{k>=1} delta_kappa_{2k} * t^{2k} / (2k)!
# The leading term is:
#   sum_{k>=1} tau_{2k+1} * t^{2k} / (2k)!
#   = sum_{k>=1} (2*t_{2k+1}/C(n,2k)) * t^{2k} / (2k)!

# Now 1/C(n,2k) = (2k)!(n-2k)!/n!, so:
# tau_{2k+1} * t^{2k} / (2k)! = 2*t_{2k+1} * (n-2k)!/n! * t^{2k} / 1

# Hmm, this doesn't simplify to a standard generating function easily.

# ALTERNATIVE APPROACH: Think of the cumulants as moments of a different object.
# The CGF of fwd(sigma) under T is:
#   K_T(t) = log(sum_{k=0}^{n-1} F_k(T) * e^{kt} / n!)
# where F_k(T) = #{sigma : fwd(sigma)=k}.

# The Worpitzky expansion gives:
#   sum_{m>=0} F(T,m) * x^m = W*(T,x) / (1-x)^n
# where W*(T,x) is related to the Worpitzky polynomial.

# The moment generating function is:
#   M_T(t) = E[e^{t*fwd}] = sum_{k>=0} F_k * e^{kt} / n!
# And K_T(t) = log M_T(t).

# The question is whether there's a SIMPLER representation.

# Actually, let me think about this differently.
# The fwd polynomial F(T,x) evaluated at x=e^t gives:
#   F(T, e^t) = sum_k F_k * e^{kt} = n! * M_T(t)
# So: K_T(t) = log(F(T, e^t) / n!)

# For the transitive tournament, F(T,x) = A_n(x) (Eulerian polynomial).
# K_trans(t) = log(A_n(e^t)/n!)
# This is a well-studied function!

# The Worpitzky identity gives: A_n(x) = sum_j w_j * x^j * (1-x)^{n-1-j}
# Actually the standard Worpitzky: sum_j A(n,j) * C(m+n-1-j, n-1) = m^n
# So: sum_j A(n,j) * C(m+n-1-j, n-1) * x^m = x^n/(1-x)^{n+1} * A_n(1/(1-x))
# Hmm this is getting complicated.

# Let me instead look at the EXPONENTIAL generating function for cumulants
# across different tournament invariants.

# Define the "cumulant Dirichlet series":
#   C(s) = sum_{k>=1} kappa_{2k} * s^{2k}

# For the transitive tournament:
#   C_trans(s) = sum_{k>=1} (-1)^{k-1}(n+1)|B_{2k}|/(2k) * s^{2k}

# This is related to: (n+1) * sum_{k>=1} (-1)^{k-1}|B_{2k}|/(2k) * s^{2k}
# = (n+1) * [log(s/sinh(s)) - s^2/6 + ...]  -- not quite standard.

# Actually: sum_{k>=1} B_{2k} * s^{2k} / (2k) = log(s / (e^s - 1)) + s/2
# (from the expansion of log(1 - sum B_{2k} s^{2k} / (2k)!) * (2k)!/(2k) ...)
# This is getting complicated. Let me just note the empirical structure.

print("CUMULANT GENERATING FUNCTION STRUCTURE")
print("=" * 60)
print()
print("For tournament T on n vertices:")
print()
print("CGF: K_T(t) = log(F(T, e^t) / n!)")
print("     = sum_{k>=1} kappa_{2k}(T) * t^{2k} / (2k)!")
print("     (only even cumulants by symmetry)")
print()
print("Decomposition: K_T(t) = K_trans(t) + delta_K(t)")
print()
print("K_trans(t) = log(A_n(e^t)/n!)")
print("  = (n+1) * sum_{k>=1} (-1)^{k-1} |B_{2k}|/(2k) * t^{2k}/(2k)!")
print()
print("delta_K(t) encodes tournament cycle structure:")
print("  Leading terms: delta_kappa_{2k} contains 2*t_{2k+1}/C(n,2k)")
print("  Cross terms: products of lower cycle corrections")
print()
print("OPEN: Is there a closed form for delta_K(t) in terms of")
print("      the tournament adjacency matrix A?")
print()

# An interesting possibility: delta_K might be expressible in terms of
# the logarithm of the Ihara zeta function of the tournament.
# Ihara: log zeta_T(u) = sum_{k>=1} tr(A^k)/k * u^k = -log det(I - uA)
# (for digraphs without backtracking, which tournaments satisfy)

# The tr(A^k) = sum of ALL directed closed walks of length k, including non-simple.
# For k=3: tr(A^3) = 6*t3 (all 3-walks are simple in a tournament).
# For k=5: tr(A^5) includes non-simple 5-walks.
# tr(A^5) = 10*t5 + correction_from_non_simple_walks

# But the cumulant hierarchy involves SIMPLE cycle counts.
# The Ihara zeta involves ALL walks.
# So the connection is not direct.

# However, for tournaments there's a nice relationship between
# the adjacency matrix eigenvalues and cycle counts.
# A tournament on n vertices has adjacency matrix A with:
# - all eigenvalue norms <= (n-1)/2
# - tr(A) = 0 (diagonal zero)
# - rank = n
# The spectral radius is at most (n-1)/2 (equality for regular tournaments).

# The determinant det(I - uA) = Ihara determinant for tournaments
# (no backtracking correction needed since tournaments have no 2-cycles).

# This means: the POLYNOMIAL det(I - uA) encodes all cycle information.
# And the cumulant hierarchy is a FILTERED version of this information,
# seen through the lens of the forward-edge moment expansion.

print("SPECTRAL CONNECTION:")
print("  det(I - uA) = Ihara zeta polynomial for tournament T")
print("  log det(I - uA) = -sum tr(A^k)/k * u^k")
print("  tr(A^{2k+1}) = (2k+1) * #simple (2k+1)-walks + correction")
print("  kappa_{2k} probes t_{2k+1} via cluster expansion")
print()
print("The cumulant hierarchy is a STATISTICAL SHADOW of the Ihara zeta:")
print("  both encode the same cycle data, but cumulants use binomial")
print("  normalization while Ihara uses exponential.")

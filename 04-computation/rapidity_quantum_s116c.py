#!/usr/bin/env python3
"""rapidity_quantum_s116c.py — Quantum mechanics and optics through rapidity

Rapidity phi = arctanh(v) = ln(Q(v))/2 where Q(x)=(1+x)/(1-x).
Q(v) = e^{2*phi}. Addition of rapidities = multiplication of Q-values.

This script explores quantum mechanics and optics through the rapidity lens.
"""
from math import sqrt, log, pi, cos, sin, exp, atanh, tanh, comb, factorial, acos

print("QUANTUM MECHANICS AND OPTICS THROUGH RAPIDITY")
print("="*70)
print()

# ============================================================
print("="*70)
print("1. SPIN AND RAPIDITY: THE BLOCH SPHERE")
print("="*70)
print()
print("  A qubit state |psi> = cos(theta/2)|0> + e^{i*phi}*sin(theta/2)|1>")
print("  lives on the Bloch sphere with polar angle theta.")
print()
print("  The probability of measuring |0>: p = cos^2(theta/2)")
print("  The probability of measuring |1>: 1-p = sin^2(theta/2)")
print()
print("  The RAPIDITY of the measurement bias:")
print("  rapidity_bias = arctanh(2p - 1) = arctanh(cos(theta))")
print()
print("  theta    p(|0>)    2p-1      rapidity_bias")
print("  " + "-"*50)
for theta_deg in [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]:
    theta = theta_deg * pi / 180
    p = cos(theta/2)**2
    bias = 2*p - 1
    if abs(bias) < 0.9999:
        rap = atanh(bias)
    else:
        rap = float('inf') if bias > 0 else float('-inf')
    if abs(rap) < 100:
        print(f"  {theta_deg:4d}   {p:8.4f}    {bias:+.4f}    {rap:+10.4f}")
    else:
        print(f"  {theta_deg:4d}   {p:8.4f}    {bias:+.4f}    {'infinity':>10s}")

print()
print("  At theta=0 (pure |0>): rapidity = +infinity (certain)")
print("  At theta=90 (equal superposition): rapidity = 0 (fair coin)")
print("  At theta=180 (pure |1>): rapidity = -infinity (certain)")
print()
print("  The RAPIDITY IS the natural parameter of the Bloch sphere's")
print("  z-projection. It measures how far from 'fair' the qubit is.")
print("  This is EXACTLY the Fisher-Rao connection from before!")
print()

# ============================================================
print("="*70)
print("2. QUANTUM TUNNELING AND THE CAYLEY TRANSFORM")
print("="*70)
print()
print("  For a rectangular barrier of height V0 and width L:")
print("  Transmission coefficient T = 1/(1 + V0^2*sin^2(kL)/(4E(V0-E)))")
print("  for E < V0 (classically forbidden).")
print()
print("  The Cayley form: Q(T-1) = (T)/(2-T) ... not quite.")
print("  Actually: define reflection R = 1 - T.")
print("  The rapidity of transmission: arctanh(2T-1) = arctanh(T-R)")
print("  = arctanh(T - (1-T)) = arctanh(2T-1)")
print()
print("  For T << 1 (deep tunneling): rapidity ~ -ln(1/T)/2 < 0")
print("  For T = 1 (resonance): rapidity = +infinity")
print("  For T = 1/2 (half-transparent): rapidity = 0")
print()
print("  T      R=1-T    rapidity     Q(2T-1)")
print("  " + "-"*50)
for T_val in [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
    rap = atanh(2*T_val - 1)
    q = exp(2*rap)
    print(f"  {T_val:.3f}   {1-T_val:.3f}    {rap:+8.4f}     {q:10.4f}")

print()
print("  The transmission rapidity LINEARIZES tunneling:")
print("  Adding barriers in series ADDS rapidity (for independent barriers).")
print("  A chain of n identical barriers has total rapidity = n * single rapidity.")
print()

# ============================================================
print("="*70)
print("3. SQUEEZED STATES AND RAPIDITY")
print("="*70)
print()
print("  A squeezed vacuum state has uncertainties:")
print("  Delta_x = (Delta_x)_0 * e^{-r}")
print("  Delta_p = (Delta_p)_0 * e^{+r}")
print("  where r is the SQUEEZING PARAMETER.")
print()
print("  The squeezing parameter r IS a rapidity!")
print("  The squeeze operator S(r) = exp(r*(a^2 - a^dag^2)/2)")
print("  is a LORENTZ BOOST in the (x,p) phase space.")
print()
print("  Squeezing r and its effects:")
print("  r      e^{-r}    e^{+r}   noise reduction (dB)")
print("  " + "-"*55)
for r in [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]:
    reduction_dB = 20 * log(exp(-r)) / log(10)  # = -20*r/ln(10)
    print(f"  {r:.1f}    {exp(-r):.4f}    {exp(r):.4f}      {reduction_dB:+.1f} dB")

print()
print("  LIGO's gravitational wave detectors use ~10 dB squeezing (r ~ 1.15).")
print("  This is rapidity ~ 1.15, velocity ~ tanh(1.15) ~ 0.82.")
print()
print("  The SQUEEZE FACTOR e^{2r} = Q(tanh(r)).")
print("  Squeezing IS the Cayley transform of hyperbolic velocity!")
print()

# ============================================================
print("="*70)
print("4. THE LAMB SHIFT AND RAPIDITY CORRECTIONS")
print("="*70)
print()
print("  The Lamb shift is the QED correction to hydrogen energy levels.")
print("  For 2S_{1/2} vs 2P_{1/2}: Delta_E ~ 1057 MHz.")
print("  The 2S energy: E_2 = -3.4 eV = -13.6/4 eV.")
print()
print("  The Lamb shift as a rapidity correction:")
print("  delta_rapidity = Delta_E / (2*E_2) ~ 1057 MHz / (2 * 3.4 eV)")
e2_hz = 3.4 * 1.602e-19 / 6.626e-34  # E_2 in Hz
lamb_shift = 1057e6  # Hz
delta_rap = lamb_shift / (2 * e2_hz)
print(f"  E_2 in Hz = {e2_hz:.4e}")
print(f"  delta_rapidity = {delta_rap:.4e}")
print()
print("  The Lamb shift is a TINY rapidity perturbation ~ 6.4e-7.")
print("  In musical terms: about 0.001 cents -- imperceptible!")
print()

# ============================================================
print("="*70)
print("5. OPTICAL FREQUENCIES AND THE ELECTROMAGNETIC SPECTRUM")
print("="*70)
print()
print("  Visible light: 400-700 nm (violet to red).")
print("  Frequency: c/lambda.")
print("  Rapidity of frequency: ln(nu)/2.")
print()

c = 3e8  # m/s
print("  Color       lambda(nm)    nu(THz)    rapidity(nu in Hz)")
print("  " + "-"*60)
colors = [
    ("Gamma ray",    0.001,),
    ("X-ray",        1.0,),
    ("UV",           300,),
    ("Violet",       400,),
    ("Blue",         470,),
    ("Green",        530,),
    ("Yellow",       580,),
    ("Orange",       610,),
    ("Red",          700,),
    ("Near IR",      1000,),
    ("Thermal IR",   10000,),
    ("Microwave",    1e6,),
    ("Radio FM",     3e9,),
]
for name, lam_nm in colors:
    lam_m = lam_nm * 1e-9
    nu = c / lam_m
    nu_THz = nu / 1e12
    rap = log(nu)/2
    print(f"  {name:14s}  {lam_nm:12.1f}  {nu_THz:12.2f}     {rap:8.3f}")

print()
print("  The visible spectrum spans rapidity from {:.3f} to {:.3f},".format(
    log(c/(700e-9))/2, log(c/(400e-9))/2))
print("  a range of {:.4f} rapidity units.".format(
    log(c/(400e-9))/2 - log(c/(700e-9))/2))
octave = log(2)/2
visible_range = log(700/400)/2
print(f"  = {visible_range/octave:.4f} octaves")
print(f"  The visible spectrum is about {visible_range/octave:.2f} octaves wide!")
print(f"  (700/400 = 1.75, just under 2 = one octave)")
print()
print("  The ENTIRE electromagnetic spectrum:")
em_range = log(3e8/0.001e-9) / 2 - log(3e8/(3e9 * 1e-9)) / 2
# Actually: gamma ~ 10^20 Hz to radio ~ 10^6 Hz
full_range = (log(1e20) - log(1e6)) / 2
print(f"  From radio (~10^6 Hz) to gamma (~10^20 Hz):")
print(f"  Rapidity range = (ln(10^20) - ln(10^6))/2 = {full_range:.2f}")
print(f"  = {full_range / octave:.1f} octaves")
print(f"  ~ 46 octaves of electromagnetic radiation!")
print()

# ============================================================
print("="*70)
print("6. SNELL'S LAW AND RAPIDITY")
print("="*70)
print()
print("  Snell's law: n1*sin(theta1) = n2*sin(theta2)")
print("  Total internal reflection when sin(theta2) > 1,")
print("  i.e., sin(theta1) > n2/n1.")
print()
print("  For the NORMAL component of velocity:")
print("  v_n = c/n * cos(theta)")
print("  The rapidity of the normal velocity:")
print("  phi_n = arctanh(v_n/c) = arctanh(cos(theta)/n)")
print()
print("  Refraction changes rapidity!")
print("  Material    n       v/c     rapidity(v)")
print("  " + "-"*45)
materials = [
    ("Vacuum",     1.0),
    ("Air",        1.0003),
    ("Water",      1.333),
    ("Glass",      1.5),
    ("Diamond",    2.417),
    ("Silicon",    3.42),
]
for name, n in materials:
    v = 1/n  # v/c
    if v >= 0.9999:
        print(f"  {name:12s}  {n:6.3f}   {v:.4f}   infinity")
    else:
        rap = atanh(v)
        print(f"  {name:12s}  {n:6.3f}   {v:.4f}   {rap:.6f}")

print()
print("  In diamond, light moves at v = c/2.417 = 0.414c,")
print(f"  rapidity = {atanh(1/2.417):.6f}")
print(f"  Compare: rapidity(1/e) = {atanh(1/exp(1)):.6f}")
print(f"  The speed of light in diamond is close to c/e!")
print()

# ============================================================
print("="*70)
print("7. COMPTON SCATTERING IN RAPIDITY")
print("="*70)
print()
print("  Compton scattering: photon hits electron, changes wavelength.")
print("  lambda' - lambda = (h/(m_e*c)) * (1 - cos(theta))")
print("  Compton wavelength: lambda_C = h/(m_e*c) = 2.426e-12 m")
print()
print("  The RAPIDITY CHANGE of the photon:")
print("  delta_rapidity = ln(lambda'/lambda)/2 = ln(1 + lambda_C/lambda * (1-cos(theta)))/2")
print()
compton = 2.426e-12  # meters
print("  For an X-ray photon (lambda = 0.1 nm = 1e-10 m):")
lam = 1e-10
for theta_deg in [0, 30, 60, 90, 120, 150, 180]:
    theta = theta_deg * pi / 180
    lam_prime = lam + compton * (1 - cos(theta))
    delta_rap = log(lam_prime / lam) / 2
    print(f"    theta = {theta_deg:4d} deg: lambda'/lambda = {lam_prime/lam:.6f}, "
          f"delta_rapidity = {delta_rap:.6f}")

print()
print("  At 180 degrees (backscatter), the rapidity shift is:")
print(f"  {log(1 + 2*compton/lam)/2:.6f}")
print(f"  = ln(1 + 2*lambda_C/lambda)/2")
print()

# ============================================================
print("="*70)
print("8. QUANTUM HALL EFFECT AND RAPIDITY QUANTIZATION")
print("="*70)
print()
print("  Hall resistance: R_H = h/(nu*e^2) where nu is the filling factor.")
print("  For integer QHE: nu = 1, 2, 3, ...")
print("  R_K = h/e^2 = 25812.807 ohms (von Klitzing constant)")
print()
print("  The rapidity of R_H for different filling factors:")
print()
R_K = 25812.807  # ohms
print("  nu    R_H (ohms)     rapidity(R_H)")
print("  " + "-"*40)
for nu in [1, 2, 3, 4, 5, 1/3, 2/5, 3/7, 2/3]:
    R_H = R_K / nu
    rap = log(R_H) / 2
    if nu == int(nu):
        print(f"  {int(nu):3d}   {R_H:12.1f}     {rap:.6f}")
    else:
        from fractions import Fraction
        f = Fraction(nu).limit_denominator(10)
        print(f"  {str(f):>4s}   {R_H:12.1f}     {rap:.6f}")

print()
print(f"  rapidity(R_K) = {log(R_K)/2:.6f}")
print(f"  Each integer filling factor nu decreases rapidity by ln(nu)/2 = rapidity(nu).")
print(f"  The QUANTIZED Hall resistance lives on a RAPIDITY LADDER")
print(f"  with steps of size rapidity(nu).")
print()

# ============================================================
print("="*70)
print("9. LASER COHERENCE AND RAPIDITY")
print("="*70)
print()
print("  A laser's coherence length L_c = c / (pi * Delta_nu)")
print("  where Delta_nu is the linewidth.")
print()
print("  For common lasers:")
lasers = [
    ("HeNe (633nm)", 633e-9, 1.5e9),  # lambda, delta_nu
    ("Nd:YAG (1064nm)", 1064e-9, 1e6),
    ("Diode (808nm)", 808e-9, 1e12),
    ("Fiber (1550nm)", 1550e-9, 1e3),
    ("Argon (488nm)", 488e-9, 5e9),
]
print("  Laser             lambda   Delta_nu    L_c (m)   rapidity(L_c/lambda)")
print("  " + "-"*70)
for name, lam, dnu in lasers:
    L_c = c / (pi * dnu)
    ratio = L_c / lam
    rap = log(ratio) / 2
    print(f"  {name:20s} {lam*1e9:.0f}nm   {dnu:.0e}   {L_c:10.3f}    {rap:8.3f}")

print()
print("  The coherence length / wavelength ratio measures the 'quality' of coherence.")
print("  Its rapidity is the 'coherence rapidity' -- higher = more coherent.")
print()

# ============================================================
print("="*70)
print("10. THE UNCERTAINTY PRINCIPLE IN RAPIDITY")
print("="*70)
print()
print("  Heisenberg: Delta_x * Delta_p >= hbar/2")
print("  Define the 'uncertainty ratio' U = Delta_x * Delta_p / (hbar/2) >= 1.")
print("  rapidity(U) = ln(U)/2 >= 0.")
print()
print("  For the ground state harmonic oscillator: U = 1, rapidity(U) = 0.")
print("  This is the MINIMUM UNCERTAINTY state = rapidity zero = at rest!")
print()
print("  For a thermal state at temperature T:")
print("  U = 2*n_bar + 1 where n_bar = 1/(exp(hbar*omega/(kT))-1).")
print()
print("  U       rapidity(U)    interpretation")
print("  " + "-"*50)
for U in [1.0, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0]:
    rap = log(U)/2
    interp = ""
    if U == 1.0:
        interp = "ground state (minimum uncertainty)"
    elif U == 2.0:
        interp = "n_bar = 0.5 thermal photon"
    elif U == 3.0:
        interp = "n_bar = 1 thermal photon"
    elif U == 10.0:
        interp = "n_bar = 4.5 thermal photons"
    print(f"  {U:6.1f}   {rap:8.4f}       {interp}")

print()
print("  The UNCERTAINTY RAPIDITY measures how far above the Heisenberg limit.")
print("  Ground state: rapidity 0 (at rest in uncertainty space).")
print("  Thermal state: rapidity grows logarithmically with temperature.")
print("  Squeezed state: rapidity 0 (minimum uncertainty maintained,")
print("  but REDISTRIBUTED between x and p).")
print()

# ============================================================
print("="*70)
print("11. THE BOHR MODEL AS A RAPIDITY LADDER")
print("="*70)
print()
print("  Bohr orbits: r_n = n^2 * a_0, v_n = alpha*c/n, E_n = -13.6/n^2 eV")
print("  where alpha = 1/137.036 (fine structure constant).")
print()
alpha_em = 1/137.036
print("  The orbital VELOCITY rapidity:")
print("  rapidity(v_n) = arctanh(alpha/n)")
print("  For n=1: arctanh(alpha) ~ alpha = {:.6f}".format(alpha_em))
print("  (Since alpha << 1, arctanh(alpha) ~ alpha.)")
print()
print("  n    v_n/c = alpha/n    rapidity(v_n)")
print("  " + "-"*45)
for n in range(1, 11):
    v = alpha_em / n
    rap = atanh(v)
    print(f"  {n:2d}    {v:.8f}        {rap:.8f}")

print()
print("  The rapidity gap between levels n and n+1:")
for n in range(1, 6):
    gap = atanh(alpha_em/n) - atanh(alpha_em/(n+1))
    ratio = gap / atanh(alpha_em)
    print(f"  n={n}->{n+1}: delta_rap = {gap:.4e} = {ratio:.4f} * rapidity(v_1)")

print()
print("  For n>>1: rapidity(v_n) ~ alpha/n.")
print("  The rapidity gaps decrease as 1/n^2.")
print("  The Bohr model is a HARMONIC rapidity ladder (like the musical overtone series).")
print()

# ============================================================
print("="*70)
print("12. GRAND SYNTHESIS: QUANTUM RAPIDITY")
print("="*70)
print()
print("  QUANTUM CONCEPT        RAPIDITY ROLE               FORMULA")
print("  " + "-"*65)
print("  Qubit bias             Fisher-Rao distance         arctanh(2p-1)")
print("  Tunneling T            transmission log            arctanh(2T-1)")
print("  Squeezing r            phase space boost           r = rapidity")
print("  Uncertainty U          excess over Heisenberg      ln(U)/2")
print("  Bohr velocity          orbital rapidity            arctanh(alpha/n)")
print("  Compton shift          wavelength change           ln(lambda'/lambda)/2")
print("  Hall resistance        quantized rapidity ladder   ln(R_K/nu)/2")
print("  EM spectrum            frequency log               ln(nu)/2")
print("  Refractive index       light slowdown              arctanh(1/n)")
print("  Laser coherence        quality measure             ln(L_c/lambda)/2")
print()
print("  ALL of these are arctanh or half-log -- the SAME rapidity.")
print("  Quantum mechanics uses rapidity as naturally as relativity does.")
print("  The Cayley transform Q connects them all.")

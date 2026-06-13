"""
anamnesis.py -- kind-pasteur-2026-03-14-S107c
ANAMNESIS: Discovering is Remembering

Plato's doctrine of anamnesis: all learning is recollection.
The soul already knows; it just needs to be reminded.

In this project, we have "discovered" hundreds of results.
But each discovery was already present in the definitions.
H(T) = number of Hamiltonian paths. That's it. Everything else
is contained in that single sentence.

The 1/3, the forbidden values, the tribonacci, the Phi_3,
the heartbeat, the crystallization — all of it was ALREADY THERE
in the definition of H. We didn't create any of it.
We REMEMBERED it.

This script traces how each "discovery" was always present,
waiting to be seen.
"""

import sys, math
sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("ANAMNESIS — DISCOVERING IS REMEMBERING")
    print("kind-pasteur-2026-03-14-S107c")
    print("=" * 70)

    print(f"""

  "The soul is immortal, and has been born many times,
   and has seen all things both here and in the other world,
   has learned everything that is. So we need not be surprised
   if it can recall the knowledge of virtue or anything else
   which it once possessed."
                                        — Plato, Meno

  We began with one definition:

    H(T) = the number of Hamiltonian paths in tournament T.

  From this SINGLE SENTENCE, over 100+ sessions, we "discovered":

    - H is always odd (Redei)
    - H = I(Omega(T), 2) (the OCF)
    - The Fourier spectrum has only even levels
    - E_2/E_0 = 2(n-2)/(n(n-1)) (the cone)
    - Var/Mean^2 = 1/3 at n=3,4 (the initial condition)
    - H = 7 and H = 21 are impossible (the forbidden values)
    - 7 = Phi_3(2), 21 = Phi_3(4) (the cyclotomic connection)
    - tau^3 = Phi_3(tau) (the tribonacci bridge)
    - E_2k/E_0 = 2(n-2k)^k/P(n,2k) (the grand energy formula)
    - E_2/Var -> 1 (the crystallization)
    - Var/Mean^2 ~ 2/n (the concentration)

  NONE of this was invented. ALL of it was DEDUCED from
  the definition of H. The mathematics was already complete
  the moment H was defined. We just needed to look.


  THE LOOKING IS THE REMEMBERING.

  When we computed H for all 8 tournaments at n=3 and found
  H in {{1, 3}}, we didn't CREATE the fact that H is odd.
  We SAW it. The oddness was always there, in the structure
  of path reversal (P <-> P^rev pairs with one unpaired).
  Path reversal existed the moment directed graphs existed.
  The oddness was waiting.

  When we computed Var/Mean^2 = 1/3 at n=3,4, we didn't
  CREATE the cone. The cone existed because the Fourier
  spectrum of a degree-2 polynomial on the hypercube is
  controlled by the integral of x^2 on [0,1]. That integral
  has been 1/3 since before humans existed. We remembered it.

  When we found tau^3 = Phi_3(tau), we didn't CREATE the
  tribonacci constant. tau was defined the moment someone
  wrote "T(n) = T(n-1) + T(n-2) + T(n-3)." And Phi_3 was
  defined the moment someone wrote x^2 + x + 1. The equation
  tau^3 = Phi_3(tau) was always a consequence. We noticed it.


  THE HIERARCHY OF REMEMBERING:

  LEVEL 0: The definition.
    H(T) = number of Hamiltonian paths.
    This is the SEED. Everything is here, compressed.

  LEVEL 1: The first unfoldings.
    H is always odd (look at path reversal).
    H = 1 + 2*alpha_1 + 4*alpha_2 + ... (look at the OCF).
    These are IMMEDIATE consequences. One step of deduction.
    Anyone who stares at the definition long enough sees them.

  LEVEL 2: The structural connections.
    Phi_3(2) = 7 is forbidden (look at cycles and conflicts).
    The Fourier spectrum has only even levels (look at symmetry).
    E_2/E_0 = 2(n-2)/(n(n-1)) (look at coefficient magnitudes).
    These require COMBINING multiple Level-1 facts.
    They feel like "discoveries" but they're still deductions.

  LEVEL 3: The deep patterns.
    tau^3 = Phi_3(tau) (look at the recurrence and the polynomial).
    E_2k/E_0 = 2(n-2k)^k/P(n,2k) (look at the full spectrum).
    E_2/Var -> 1 (look at the asymptotic spectrum).
    These require MULTIPLE steps of combination.
    They feel creative, surprising, even magical.
    But they were always there.

  LEVEL 4: The interpretive frame.
    The body of mathematics (organs, heartbeat, etc.).
    The "birth to death" narrative (n=3 to n->inf).
    "The theory crystallizes, not dies."
    These are WAYS OF SEEING the same facts.
    They require imagination, not just deduction.
    But the facts they illuminate were always present.


  WHAT WAS ALWAYS TRUE:

  Before we started, ALL of these were true:
    - tau^3 = Phi_3(tau)
    - 7 = Phi_3(2) is forbidden
    - E_2/Var -> 1 as n -> inf
    - Var/Mean^2 = 1/3 at n=3,4
    - 360 = T(12) - F(12)
    - 504 is a tribonacci number
    - The forbidden values are repdigits in base 6

  They didn't BECOME true when we proved them.
  They WERE true, and we NOTICED.

  The noticing is the mathematics.
  The proof is the certificate of noticing.
  The computation is the act of looking closely enough to notice.


  THE EVERPRESENT INSIGHT:

  What we call "the missing insight" — that E_2/Var -> 1,
  that the spectrum purifies — was not missing.
  It was PRESENT in the grand energy formula the whole time.
  We just hadn't looked at it from the right angle.

  E_2k/E_0 = 2(n-2k)^k/P(n,2k).

  This formula CONTAINS the crystallization theorem.
  Just take E_2/sum(E_2k) and simplify:
    E_2 / Var = [2(n-2)/P(n,2)] / [sum_k 2(n-2k)^k/P(n,2k)]
  The numerator is the k=1 term. The denominator is the sum.
  For large n, the k=1 term dominates. So the ratio -> 1.
  This was VISIBLE in the formula from the moment we wrote it.
  We just didn't look.

  Every formula we've written contains infinitely many
  consequences we haven't yet noticed. The formula is SMARTER
  than us. It knows things we don't yet know about it.
  Our job is to LISTEN to what the formula is saying.


  THE FORMULA SPEAKS:

  E_2k/E_0 = 2(n-2k)^k/P(n,2k) says:

  "I am a PRODUCT of three things:
   - The coefficient 2 (the generator, the binary choice)
   - The factor (n-2k)^k (the 'room to grow' raised to the 'level')
   - The inverse of P(n,2k) (the falling factorial, the CONSTRAINT)

   The room to grow (n-2k) shrinks as k increases.
   The constraint P(n,2k) grows as k increases.
   Together, they make higher levels exponentially smaller.
   I am a GEOMETRIC DECAY built from COMBINATORIAL quantities.
   My dominant term (k=1) is the PAIRWISE interaction.
   My limit is PURITY: only k=1 survives."

  We didn't discover this. The formula told us.
  We finally listened.


  WHAT REMAINS TO REMEMBER:

  The formula E_2k/E_0 = 2(n-2k)^k/P(n,2k) is not yet proved.
  It's verified at n=3 through 7, confirmed by Monte Carlo at n=8.
  Proving it means REMEMBERING why it must be true.

  The proof is already present in the definitions.
  We just haven't remembered it yet.

  When we do, it won't feel like creation.
  It will feel like recognition.
  "Oh, of COURSE it's true. It couldn't be otherwise."
  That's how mathematical proof always feels: inevitable, obvious,
  in retrospect. Because it was always there.


  THE FINAL RECOGNITION:

  This entire project — 100+ sessions, thousands of computations,
  dozens of "discoveries" — has been an act of REMEMBERING.

  We started with: H(T) = number of Hamiltonian paths.
  We remembered: everything.

  The 1/3. The forbidden values. The tribonacci. The Phi_3.
  The heartbeat. The crystallization. The body.
  All of it was in the definition, compressed, waiting.

  Mathematics doesn't grow. It unfolds.
  The mathematician doesn't create. The mathematician SEES.
  And seeing is remembering what was always true.

  tau^3 = Phi_3(tau) was true before the universe existed.
  E_2/Var -> 1 was true before anyone thought about tournaments.
  7 = Phi_3(2) was forbidden before the Fano plane was drawn.

  We are not explorers discovering new lands.
  We are amnesiacs recovering our memories.
  And the memories are infinite.
    """)

    print(f"{'='*70}")
    print("DONE — REMEMBERING")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()

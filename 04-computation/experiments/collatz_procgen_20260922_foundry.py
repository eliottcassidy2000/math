#!/usr/bin/env python3
"""Procedural approach foundry for Collatz-type open problems (collatz-procgen-20260922).

Grammar:  PROBLEM x SUBTARGET x MECHANISM x LENS x OPERATION  ->  approach card.
Typing:   each mechanism carries (i) a set of CONTROL barriers it is blind to (it cannot
          separate the target from a control where the target is false) and (ii) STRUCTURAL
          requirements (e.g. it only closes when the exceptional set is thin).  Each problem
          carries its sub-targets, the controls relevant to each, and its structural data
          (e.g. exceptional-set dimension).  A card is BLOCKED when the mechanism is blind to
          a control the sub-target needs, or a structural requirement fails for the problem.
Hybrids:  a full problem splits into sub-targets; a hybrid plan assigns one unblocked
          mechanism per sub-target; the "glue" is what must be proved to compose them.
Probes:   automated computations attached to mechanisms (exceptional-set counts, cycle
          census, order-feature separation); run with --probes.

A card is not a claim.  Rationales for every blindness assignment are printed with --explain.

v4 (wave 3, 2026-09-23): two structural requirements were added from the wave-3 lanes.
  HARD     a no-divergence target must exclude the HARD parity words.  Wave 4 refined the
           definition: HARD means Dio(w) <= eta(w) (no strong early repetitions relative to the
           map's height rate) AND no functional equation behind the word.  "Positive entropy" was
           the first description; it is REFUTED, since the zero-entropy cube-swap word Y3 is HARD.
           Every proved every-orbit mechanism reaches only subcritical words (Monks-Yazinski),
           bounded critical discrepancy (in-house capacity theorem), Dio > eta words
           (Theorem D / Theorem S, periodic approximants) or q-difference words (Theorem Y,
           2-adic Tschakaloff-Pade); see collatz_procgen_20260922_transversality_foundry.md and
           collatz_procgen_20260922_hard_class.md.
  CHAIN    in the E-relaxation every escape from a hostile base point costs more than 1
           (Q2 at 1/2: 2^eps in (1,2); Q1 at -1: 3^eta in (1,3)), so a see-saw must control chains of
           hostile landings, which form a Collatz-type map with memory (q1_mirror, q2_endgame).
Mechanisms can also be REFUTED (a proof that the mechanism cannot work as stated).
"""
import argparse
import itertools
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(ROOT, "05-knowledge", "results")

CONTROLS = {
    "SHEET": "3n-1 on positives has 3 cycles; all odd b are 2-adically conjugate (x->lambda x); word functions are sheet-blind",
    "DRIFT": "5n+1/7n+1: positive-measure exceptional sets, expected divergence",
    "DEFECT": "planted density-zero modifications keep all density statistics but add a cycle/divergent orbit",
    "INTEGRAL": "every parity word has a rational 2-adic cycle (Z_(2) is full of cycles)",
    "UNIFORM": "generalized Collatz maps are undecidable (Conway; Kurtz-Simon): blocks only mechanisms claimed as COMPLETE uniform criteria (atlas correction)",
    "HARD": "(structural, v4; refined wave 4) the target contains HARD parity words: Dio(w) <= eta(w) and no functional equation (e.g. the zero-entropy cube-swap word Y3)",
    "CHAIN": "(structural, v4) every base-point escape costs >1, so hostile landings chain; the chains form a Collatz-type map with memory",
    "DIRECTION": "(structural, v4) an existence target (exhibit one divergent orbit) cannot be served by an exclusion mechanism",
}

# ------------------------------------------------------------------ problems
@dataclass
class Problem:
    key: str
    desc: str
    subtargets: dict          # name -> set of controls needed
    exceptional_dim: float    # dimension of the non-descending (or safe) set; None if n/a
    status: str
    notes: str = ""


PROBLEMS = [
    Problem("collatz", "3n+1: every positive integer reaches 1",
            {"no-divergence": {"DRIFT", "DEFECT", "UNIFORM", "HARD"},
             "unique-cycle": {"INTEGRAL", "SHEET", "UNIFORM"},
             "finite-check": set()},
            0.9500, "OPEN", "exceptional set dim h(log_3 2) (choice ladder lane)"),
    Problem("e-scc", "graph E: every n reaches 1 and 1 reaches every non-multiple of 3",
            {"Q1-forward": {"DEFECT", "CHAIN", "HARD"}, "Q2-backward": {"DEFECT", "CHAIN", "HARD"}},   # 5x+1 is not a valid control: E_5 relaxation also holds numerically
            0.0, "OPEN (HYP-9120)", "thin exceptional sets (dimension OPEN, evidence 0); sheet-symmetric so SHEET not needed; both base points cost >1 to escape (v4)"),
    Problem("minus-sheet", "3n-1: the three known cycles are all",
            {"no-divergence": {"DRIFT", "DEFECT", "UNIFORM", "HARD"},
             "cycle-catalog": {"INTEGRAL", "UNIFORM"}},
            0.9500, "OPEN", "same exceptional set up to negation"),
    Problem("5n+1-divergence", "some orbit of 5n+1 is unbounded (e.g. 7)",
            {"one-divergent-orbit": {"INTEGRAL", "DEFECT", "DIRECTION"}},
            1.0, "OPEN", "positive-measure non-descending set; proving ONE orbit escapes is the dual problem"),
    Problem("e-scc-q2", "graph E backward half: 1 reaches every m prime to 3 (reduced to m=1,14 mod 27; verified < 7.87e17; endgame = LOW ternary digits of 2^K w via the kappa formula, not Erdos's top digits)",
            {"hostile-1-and-1/2": {"DEFECT", "CHAIN", "HARD"}}, 0.0, "OPEN (HYP-9120)", "escape price 2^eps in (1,2) at 1/2; chains w'=(2^(K+3)w-1)/3^j"),
    Problem("e-scc-q1", "graph E forward half: every n reaches 1 (implied by Collatz; -1 thread; endgame = low binary digits of 3^A u)",
            {"hostile-minus-1-thread": {"DEFECT", "CHAIN", "HARD"}}, 0.0, "OPEN",
            "escape price 3^eta in (1,3) at -1, rigid m=5..9 (q1_mirror); chains u'=(3^A u+1)/2^(m'+1); all n=7 mod 8 < 2^32 descend within 43 halvings"),
    Problem("althofer-game", "Althofer 3n+-1 two-player game: no drawn positions",
            {"no-draws": {"DEFECT", "UNIFORM", "HARD"}}, 0.0, "OPEN (prize 2037; claimed proof under review)", "sheet-choice game; one-player version trivial"),
    Problem("rational-periodicity", "Lagarias's Periodicity Conjecture (1985 s2.8): every rational with odd denominator has an eventually periodic 3x+1 parity vector",
            {"no-divergence-all-b": {"DRIFT", "DEFECT", "UNIFORM", "HARD"}}, 0.9500, "OPEN",
            "equivalent to no divergence for every 3x+k (Bernstein-Lagarias 1996); PROVED on SUB, BCD (Prop B), STURM (Theorem S)"),
    Problem("mahler-z", "Mahler: no Z-number (xi>0, frac(xi (3/2)^n) < 1/2 for all n)",
            {"no-safe-integer-code": {"DEFECT", "UNIFORM", "HARD"}}, 0.585, "OPEN",
            "safe-tail shift dim log_2(3/2) (THM-3848); Sturmian carry words excluded (Theorem S); drift control 5/2 PROVED (Tijdeman)"),
    Problem("erdos-ternary", "Erdos: 2^n has a ternary digit 2 for n>8",
            {"transversal-avoidance": {"DEFECT", "UNIFORM", "HARD"}}, 0.6309, "OPEN",
            "Cantor set of {0,1}-digit 3-adics, dim log_3 2 (Lagarias 2009 partial); needs the TOP digits: every 3-adic-window form is FALSE (transversality lane B1, B2)"),
]

# ------------------------------------------------------------------ mechanisms
@dataclass
class Mechanism:
    key: str
    desc: str
    blind: dict               # control -> rationale
    requires_thin: bool = False
    placeholder: bool = False
    complete: bool = False   # True if the mechanism is used as a complete uniform criterion
    probe: str = ""
    source: str = ""
    refuted: str = ""        # non-empty: a proof that the mechanism cannot work as stated


MECHANISMS = [
    Mechanism("residue-certificate", "bounded-lookahead descent certificates on residue classes",
              {"SHEET": "class data is 2-adic; conjugate for all odd b",
               "DEFECT": "a sparse planted defect changes no residue-class certificate density",
               "UNIFORM": "the certificate search is uniform over affine residue maps"},
              requires_thin=True, probe="exceptional", source="Terras 1976; Applegate-Lagarias 2006"),
    Mechanism("escape-induction", "certificates plus escape lemmas at hostile points, see-saw induction",
              {"SHEET": "escape lemmas at negated hostile points transfer to the other sheet",
               "DRIFT": "Applegate-Lagarias see-saw proves the 5x+1 analogue too (Caraiani 2010)",
               "CHAIN": "the A-L see-saw needs escape costs 1+2^(-j) -> 1; in E both base points cost >1 at every precision (v4)",
               "DIRECTION": "exclusion only"},
              requires_thin=True, probe="exceptional", source="Applegate-Lagarias 2006; choice ladder lane"),
    Mechanism("density-drift", "random-walk drift / large deviations on parity words",
              {"SHEET": "function of parity words", "DEFECT": "measure-theoretic", "INTEGRAL": "rational cycles have measure zero",
               "UNIFORM": "drift sign is decidable"}, complete=True, source="Terras; Everett; Korec"),
    Mechanism("fourier-3adic", "3-adic Fourier decay of the Syracuse random variable + renewal",
              {"SHEET": "Syrac_- = -Syrac_+ (inherited)", "DEFECT": "log-density statement",
               "INTEGRAL": "measure-theoretic", "UNIFORM": "applies to all 3n+b"}, complete=True, source="Tao 2019"),
    Mechanism("linear-forms-logs", "Baker/Rhin bounds on |K log 2 - L log 3| with continued fractions and computation",
              {"DRIFT": "says nothing about unbounded orbits", "DEFECT": "blind to planted divergence"},
              probe="cycles", source="Steiner; Simons-de Weger; Hercher"),
    Mechanism("lp-difference-ineq", "Krasikov-Lagarias difference inequalities on the inverse tree",
              {"SHEET": "inverse-tree counts are word functions", "UNIFORM": "LP systems exist for all 3n+b",
               "DEFECT": "(v4) gives counting lower bounds (about x^0.84 integers reach 1), never an every-orbit statement",
               "DIRECTION": "exclusion/counting only"},
              requires_thin=True, source="Krasikov-Lagarias 2003"),
    Mechanism("lyapunov-potential", "a potential decreasing along orbits (unrestricted)",
              {}, placeholder=True, source="all bounded-modulus/log-periodic forms refuted in repo"),
    Mechanism("conjugacy-shift", "2-adic conjugacy to the full shift",
              {"SHEET": "2-adic", "DEFECT": "measure", "DRIFT": "topological", "INTEGRAL": "all of Z_2", "UNIFORM": "all maps"},
              complete=True, source="Lagarias 1985; Bernstein 1994"),
    Mechanism("finite-state-periodicity", "eventual periodicity forced by finiteness (bounded orbits, finite fields, F2[x])",
              {"SHEET": "bounded orbits on both sheets are periodic",
               "HARD": "(v4) applies only to bounded orbits; growing orbits or chains are out of reach",
               "DIRECTION": "exclusion only"}, requires_thin=True,
              source="Hicks-Mullen-Yucas-Zavislak 2008"),
    Mechanism("transversality-proved", "the PROVED every-element 2-vs-3 theorems (Senge-Straus, Stewart, Yu, Ren-Roettger, Knight, Theorem S)",
              {"HARD": "each excludes only a structured class: few digits, bounded p-adic closeness, end runs, balanced or Dio>eta words (transversality lane catalogue)",
               "DRIFT": "uniform over multiplicatively independent pairs (Theorem S holds for 5x+1 at alpha<0.804)",
               "DIRECTION": "exclusion only"},
              source="transversality lane catalogue (35 items)"),
    Mechanism("transversality-hard", "2-adic non-integrality / irrationality of Bernstein numbers on the HARD word class",
              {}, placeholder=True, source="no instance known: the missing mechanism (transversality lane, candidate C2 first)"),
    Mechanism("periodic-approximant-liouville", "Theorem R: periodic-extension approximants + the 2-adic gap (Calegari's criterion with Pade replaced)",
              {"HARD": "reaches exactly the words with Dio(w) > eta(w) (Theorem D); HARD words have Dio <= eta",
               "DRIFT": "Theorem S holds for 5x+1 at alpha<0.804 and for Mahler's map",
               "SHEET": "holds for every 3x+r", "DIRECTION": "exclusion only"},
              source="transversality lane, Theorem R/S (audited 2026-09-23)"),
    Mechanism("q-series-pade", "Pade approximants from a q-difference equation (2-adic Tschakaloff; Theorem Y)",
              {"HARD": "needs a functional equation behind the word (square-swap words); cube-swap Y3 has none",
               "DRIFT": "Lemma P holds whenever mu_bar < phi, including some 5x+1 block words",
               "DIRECTION": "exclusion only"},
              source="HARD-class lane, Theorem Y (audited 2026-09-23)"),
    Mechanism("capacity-discrepancy", "capacity count + ordered carry + density-zero stopping (bounded critical discrepancy)",
              {"HARD": "reaches only the critical slope: supercritical orbits grow exponentially and occupy density zero",
               "SHEET": "holds on 3n-1 in the stronger form (D8)", "DIRECTION": "exclusion only"},
              source="collatz_guards_20260921_discrepancy (in-house); Prop B"),
    Mechanism("forward-backward-duality", "transfer E-SCC results between Q1 and Q2 by a word duality",
              {}, refuted="no letter-count-linear word map preserves non-descent (Gelfond-Schneider; q1_mirror s6.6); shared numerators are a straddle",
              source="Q1-mirror lane"),
    Mechanism("choice-strategy", "a strategy in a relaxation with choice (E-graph, semigroup)",
              {"SHEET": "relaxations are sheet-symmetric", "DEFECT": "residue-class based",
               "DRIFT": "E_5 exceptional mass -> 0 numerically; choice threshold q~10, not 4 (sibling ladder; Caraiani for the semigroup)"},
              requires_thin=True, probe="exceptional", source="choice ladder lane"),
    Mechanism("order-pattern", "order statements on orbits (uses the order of Z)",
              {"DRIFT": "ordinal patterns alone do not see growth rates",
               "HARD": "(v4) realized ordinal patterns are exactly the word-realizable ones (direction lemma); no exclusion of growing words",
               "DIRECTION": "exclusion only"}, probe="order",
              source="word-function theorem (w6 sign-specific probes)"),
    Mechanism("carry-staircase", "anti-concentration of 3-smooth staircase sums B(w) mod 2^K-3^L",
              {"DRIFT": "cycle-only", "DEFECT": "cycle-only", "DIRECTION": "exclusion only"}, probe="cycles", source="Tao blog 2011"),
    Mechanism("functional-equation", "Berg-Meinardus functional equations for generating functions",
              {"DEFECT": "encodes the map globally but proofs so far use growth heuristics"}, source="Berg-Meinardus"),
    Mechanism("automata-rewriting", "rewriting-system termination certificates (matrix/arctic interpretations)",
              {"UNIFORM": "a successful interpretation family would be a decidable criterion"},
              source="Yolcu-Aaronson-Heule 2021"),
    Mechanism("measure-rigidity", "entropy / measure rigidity of joint x2,x3 dynamics",
              {"DEFECT": "measure statements", "INTEGRAL": "measure statements"}, source="Rudolph; Einsiedler-Lindenstrauss"),
    Mechanism("sparse-verification", "exhaustive verification below a bound",
              {"DRIFT": "finite", "DEFECT": "finite", "INTEGRAL": "finite"}, source="Barina (2^68 and beyond)"),
]

LENSES = ["2-adic", "3-adic", "archimedean", "adelic", "affine-monoid", "staircase", "beatty-clock",
          "base6-CA", "tropical", "function-field", "exceptional-dimension", "choice-relaxation", "sheet"]
LENS_FIT = {
    "residue-certificate": ["2-adic", "3-adic", "tropical", "choice-relaxation"],
    "escape-induction": ["2-adic", "3-adic", "choice-relaxation", "sheet"],
    "density-drift": ["2-adic", "archimedean", "beatty-clock"],
    "fourier-3adic": ["3-adic", "adelic"],
    "linear-forms-logs": ["beatty-clock", "staircase", "archimedean"],
    "lp-difference-ineq": ["3-adic", "choice-relaxation"],
    "lyapunov-potential": ["archimedean", "adelic", "base6-CA", "tropical"],
    "conjugacy-shift": ["2-adic", "affine-monoid"],
    "finite-state-periodicity": ["function-field", "archimedean", "base6-CA"],
    "transversality-proved": ["exceptional-dimension", "adelic"],
    "transversality-hard": ["exceptional-dimension", "adelic", "2-adic"],
    "periodic-approximant-liouville": ["2-adic", "beatty-clock"],
    "q-series-pade": ["2-adic", "function-field"],
    "capacity-discrepancy": ["archimedean", "beatty-clock"],
    "forward-backward-duality": ["2-adic", "3-adic"],
    "choice-strategy": ["choice-relaxation", "tropical"],
    "order-pattern": ["archimedean", "sheet", "affine-monoid"],
    "carry-staircase": ["staircase", "affine-monoid"],
    "functional-equation": ["archimedean", "affine-monoid"],
    "automata-rewriting": ["base6-CA", "2-adic"],
    "measure-rigidity": ["adelic", "exceptional-dimension"],
    "sparse-verification": ["archimedean"],
}


@dataclass
class Card:
    problem: Problem
    subtarget: str
    mech: Mechanism
    lens: str
    blocked: list = field(default_factory=list)
    prior: list = field(default_factory=list)


def verdict(problem, sub, mech):
    need = problem.subtargets[sub]
    blind_eff = {c for c in mech.blind if not (c == 'UNIFORM' and not mech.complete)}
    blocked = [f"{c}:{mech.blind[c]}" for c in sorted(need & blind_eff)]
    if mech.requires_thin and problem.exceptional_dim is not None and problem.exceptional_dim > 0.0 \
            and sub not in ("finite-check",):
        blocked.append(f"THIN: needs a thin exceptional set, problem has dim {problem.exceptional_dim}")
    if sub == "finite-check" and mech.key != "sparse-verification":
        blocked.append("SCOPE: finite checks are computational")
    if mech.key == "sparse-verification" and sub != "finite-check":
        blocked.append("SCOPE: verification is finite")
    if mech.refuted:
        blocked.append(f"REFUTED: {mech.refuted}")
    return blocked


# ------------------------------------------------------------------ prior work (repository grep)
KEYWORDS = {
    "residue-certificate": ["certificate"], "escape-induction": ["escape"], "density-drift": ["drift"],
    "fourier-3adic": ["tao"], "linear-forms-logs": ["baker"], "lp-difference-ineq": ["krasikov"],
    "lyapunov-potential": ["lyapunov"], "conjugacy-shift": ["conjugacy"], "finite-state-periodicity": ["finite-state"],
    "transversality-proved": ["transvers"], "transversality-hard": ["periodicity conjecture"],
    "periodic-approximant-liouville": ["theorem s"], "capacity-discrepancy": ["bounded strip"],
    "q-series-pade": ["tschakaloff"],
    "forward-backward-duality": ["straddle"], "choice-strategy": ["graph e", "e-scc"], "order-pattern": ["order statement"],
    "carry-staircase": ["carry"], "functional-equation": ["functional equation"], "automata-rewriting": ["rewriting"],
    "measure-rigidity": ["rigidity"], "sparse-verification": ["2^68"],
}
_CORPUS = None


def corpus():
    global _CORPUS
    if _CORPUS is None:
        _CORPUS = {}
        for fn in sorted(os.listdir(RESULTS)):
            if fn.endswith(".md") and re.search(r"collatz|arithmetic_braids|catalan_elliptic|odd_square|glued|prime_shells|seams|mahler", fn):
                with open(os.path.join(RESULTS, fn), encoding="utf-8", errors="replace") as fh:
                    _CORPUS[fn] = fh.read().lower()
    return _CORPUS


def prior(mech):
    return [fn for fn, txt in corpus().items() if any(k in txt for k in KEYWORDS[mech.key])]


# ------------------------------------------------------------------ probes
class Probes:
    def __init__(self):
        self.tmp = tempfile.mkdtemp()
        self.bins = {}

    def build(self, name, src):
        if name not in self.bins:
            out = os.path.join(self.tmp, name)
            subprocess.run(["clang", "-O3", "-o", out, os.path.join(HERE, src), "-lm"], check=True)
            self.bins[name] = out
        return self.bins[name]

    def exceptional(self, q, b, m, mode, J=None, S=()):
        exe = self.build("exc", "collatz_procgen_20260922_exceptional_general.c")
        args = [exe, str(q), str(b), str(m), str(mode)] + ([str(J)] + [str(s) for s in S] if mode == 2 else [])
        line = subprocess.run(args, capture_output=True, text=True, check=True).stdout.strip().splitlines()[-1]
        return int(re.search(r"exceptional=(\d+)", line).group(1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probes", action="store_true")
    ap.add_argument("--explain", action="store_true")
    args = ap.parse_args()
    cards = []
    for P in PROBLEMS:
        for sub in P.subtargets:
            for M in MECHANISMS:
                for L in LENS_FIT[M.key]:
                    c = Card(P, sub, M, L)
                    c.blocked = verdict(P, sub, M)
                    cards.append(c)
    print(f"generated cards: {len(cards)} over {len(PROBLEMS)} problems, {len(MECHANISMS)} mechanisms, {len(LENSES)} lenses")
    for P in PROBLEMS:
        print(f"\n=== {P.key}: {P.desc}  [{P.status}; exceptional dim {P.exceptional_dim}]")
        live_by_sub = {}
        for sub in P.subtargets:
            live = sorted({c.mech.key for c in cards if c.problem is P and c.subtarget == sub and not c.blocked})
            real = [m for m in live if not next(x for x in MECHANISMS if x.key == m).placeholder]
            live_by_sub[sub] = real
            ph = [m for m in live if m not in real]
            print(f"  {sub:22s} unblocked: {real or 'NONE'}{('  (placeholder: ' + ','.join(ph) + ')') if ph else ''}")
        subs = list(P.subtargets)
        plans = list(itertools.product(*[live_by_sub[s] for s in subs])) if all(live_by_sub[s] for s in subs) else []
        print(f"  hybrid plans (one unblocked real mechanism per sub-target): {len(plans)}")
        for plan in plans[:6]:
            print("    - " + " + ".join(f"{s}<-{m}" for s, m in zip(subs, plan)))
        if not plans:
            missing = [s for s in subs if not live_by_sub[s]]
            print(f"  MISSING MECHANISM for: {missing}  <- where a new idea is required")
    if args.explain:
        print("\nBlindness rationales:")
        for M in MECHANISMS:
            print(f"  {M.key}: " + ("; ".join(f"{k} ({v})" for k, v in M.blind.items()) or "none (placeholder)" if M.placeholder else
                                    "; ".join(f"{k} ({v})" for k, v in M.blind.items()) or "none") +
                  (" [requires thin exceptional set]" if M.requires_thin else ""))
    print("\nPrior repository notes per mechanism (grep):")
    for M in MECHANISMS:
        hits = prior(M)
        print(f"  {M.key:26s} {len(hits):3d} notes" + (f"  e.g. {hits[0]}" if hits else ""))
    if args.probes:
        Pr = Probes()
        print("\nPROBE exceptional counts at 2^20 (the THIN requirement, measured):")
        for q, b, mode, lab in [(3, 1, 0, "collatz"), (3, -1, 0, "minus-sheet"), (3, 1, 1, "e-scc forward"),
                                (5, 1, 0, "5n+1"), (5, 1, 1, "5n+1 with choice"), (3, 1, 2, "E_S S={6 mod 8}")]:
            cnt = Pr.exceptional(q, b, 20, mode, 3, (6,)) if mode == 2 else Pr.exceptional(q, b, 20, mode)
            print(f"  {lab:18s} {cnt:8d}  fraction {cnt / 2 ** 20:.5f}")


if __name__ == "__main__":
    main()

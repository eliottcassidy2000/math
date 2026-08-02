#!/usr/bin/env python3
"""Exact depth-six non-resurrection certificate for THM-3169.

For support ``(1,3)``, bank ``I2``, the companion reconstructs every one of
the 1,189 multiplicity-valid physical prefix states of depth at most six.
It first shows that the THM-3158 nine-row separator is crossed by exactly 24
of the 507 new depth-six states.  It then verifies a primitive positive
eleven-row combination which is strictly negative on the complete bank.

The floating separation oracle used to discover the eleven rows is not part
of the proof.  All promoted checks below use integers and ``Fraction`` only.
"""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations_with_replacement
from math import gcd
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"
UPSTREAM_BANK = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
UPSTREAM_LF_SHA256 = (
    "151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7")
UPSTREAM_BANK_LF_SHA256 = (
    "15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f")


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


def load_upstream_prefix(maximum_degree):
    actual = hashlib.sha256(lf_bytes(UPSTREAM)).hexdigest()
    require(actual == UPSTREAM_LF_SHA256,
            ("pole-prefix helper hash drift", actual))
    bank_actual = hashlib.sha256(lf_bytes(UPSTREAM_BANK)).hexdigest()
    require(bank_actual == UPSTREAM_BANK_LF_SHA256,
            ("signed-bank helper hash drift", bank_actual))
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "MAXIMUM_DEGREE"
                        for target in node.targets)):
            node.value = ast.Constant(maximum_degree)
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "UNIVERSE"
                        for target in node.targets)):
            break
    module = ast.fix_missing_locations(
        ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(13)
# The upstream routines are pure exact maps.  A fixed prefix state invokes
# them repeatedly on the same residual-root alphabets, so memoizing these two
# maps changes no arithmetic and keeps the immutable replay inexpensive.
UP["residual_roots"] = lru_cache(maxsize=None)(UP["residual_roots"])
UP["all_monomial_values"] = lru_cache(maxsize=None)(
    UP["all_monomial_values"])
partitions = UP["partitions"]
coarsens = UP["coarsens"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
BANK = UP["BANKS"][1]

POLES, _ = reduced_poles(1, BANK, 1, 3)
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
BY_DEPTH = tuple(
    tuple(
        state
        for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    )
    for depth in range(1, 7)
)
STATES = tuple(state for layer in BY_DEPTH for state in layer)
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,3) pole multiset drift")
require(tuple(map(len, BY_DEPTH)) == (8, 33, 93, 200, 348, 507)
        and len(STATES) == 1189,
        "depth-six physical-state census drift")

VECTORS = {
    state: coefficient_vectors(1, BANK, 1, 3, state)
    for state in STATES
}


def complement_upset(degree, excluded):
    return frozenset(shape for shape in partitions(degree)
                     if shape not in excluded)


def generated_upset(degree, generators):
    return frozenset(
        shape for shape in partitions(degree)
        if any(coarsens(generator, shape) for generator in generators)
    )


def response_row(degree, upset):
    row = tuple(
        sum(VECTORS[state][degree][shape] for shape in upset)
        for state in STATES
    )
    require(all(value.denominator == 1 for value in row),
            "response row lost integrality")
    return tuple(value.numerator for value in row)


def combine(rows, coefficients):
    return tuple(
        sum(coefficient * row[index]
            for coefficient, row in zip(coefficients, rows))
        for index in range(len(STATES))
    )


# ---------------------------------------------------------------------------
# 1. The old depth-five wall is genuinely crossed at depth six.
# ---------------------------------------------------------------------------

OLD_R0 = complement_upset(8, {(1,) * 8})
OLD_R1 = complement_upset(10, {(1,) * 10})
OLD_R2 = complement_upset(11, {(1,) * 11})
OLD_R3 = complement_upset(12, {(1,) * 12})
OLD_R4 = complement_upset(13, {(1,) * 13})
OLD_R5 = generated_upset(11, ((2, 2) + (1,) * 7,))
OLD_R6 = generated_upset(13, (
    (4, 2) + (1,) * 7,
    (3, 3) + (1,) * 7,
    (3, 2, 2) + (1,) * 6,
    (2, 2, 2, 2) + (1,) * 5,
))
OLD_R7 = generated_upset(13, (
    (4, 2) + (1,) * 7,
    (3, 3, 3) + (1,) * 4,
    (3, 2, 2, 2) + (1,) * 4,
    (2, 2, 2, 2, 2) + (1,) * 3,
))
OLD_R8 = generated_upset(13, (
    (4, 2) + (1,) * 7,
    (3, 3) + (1,) * 7,
    (3, 2, 2) + (1,) * 6,
    (2, 2, 2, 2, 2) + (1,) * 3,
))
OLD_ROW_SPECS = (
    (8, OLD_R0), (10, OLD_R1), (11, OLD_R2), (12, OLD_R3),
    (13, OLD_R4), (11, OLD_R5), (13, OLD_R6), (13, OLD_R7),
    (13, OLD_R8),
)
OLD_ROWS = tuple(response_row(degree, upset)
                 for degree, upset in OLD_ROW_SPECS)
OLD_COEFFICIENTS = (
    79966346203432495238210892467836051822404864748631537035246352913301789115999332193913950905464025191074830103726371948258801337873186194518,
    4385616886032821959435837177631104375265963397917293088514990249982327757237457317872391016149793462027134777492436937280742150390259271553,
    323165012303885417971182559526221438425147324616172737787579961666908713240511160571617207173787027117769558974661863206105407364575577287,
    7744440757659416200591308784108711415468181697333643558792830072051161414809397905536549772356510263760667669854566940871988150246747864,
    64198555842605394277362557914437771970209622650006876621643073243158046817590826310245420481922857881786832819412787764831035747707567,
    1238919946181093630207021295990664443370069758642989075525265110060651909301010930905141007031425295949856900214021045838758589520118799,
    57173834483208673247111573734889531022707614332675275902916730146155429982624862452119240877347978257760901011795049436871545379049,
    3220669545142804044646293595442679616378197077725758366899484010296430542402560592047500094580557804747642540554956255672701117951,
    43322631946675983735399838865716397620757342948083795752004906378365718688259231686082777317064579418881565313974450622142635728100,
)
OLD_COMBINED = combine(OLD_ROWS, OLD_COEFFICIENTS)
require(all(value < 0 for value in OLD_COMBINED[:682]),
        "THM-3158 separator stopped killing its original bank")
OLD_DEPTH6 = OLD_COMBINED[682:]
OLD_POSITIVE_STATES = tuple(
    state for state, value in zip(BY_DEPTH[5], OLD_DEPTH6) if value > 0
)
EXPECTED_OLD_POSITIVE_STATES = (
    (1, 1, 1, 1, 2, 2),
    (1, 1, 1, 1, 2, 3),
    (1, 1, 1, 1, 3, 3),
    (1, 1, 1, 2, 2, 2),
    (1, 1, 1, 2, 2, 3),
    (1, 1, 1, 2, 2, 4),
    (1, 1, 1, 2, 3, 3),
    (1, 1, 2, 2, 2, 3),
    (1, 1, 2, 2, 3, 3),
    (1, 4, 5, 5, 6, 7),
    (1, 4, 5, 5, 6, 8),
    (2, 3, 5, 5, 6, 7),
    (2, 3, 5, 5, 6, 8),
    (2, 4, 5, 5, 6, 7),
    (2, 4, 5, 5, 6, 8),
    (3, 3, 5, 5, 6, 7),
    (3, 3, 5, 5, 6, 8),
    (3, 3, 5, 6, 7, 8),
    (3, 4, 4, 5, 6, 7),
    (3, 4, 4, 5, 6, 8),
    (3, 4, 5, 5, 6, 7),
    (3, 4, 5, 5, 6, 8),
    (4, 4, 5, 5, 6, 7),
    (4, 4, 5, 5, 6, 8),
)
require(OLD_POSITIVE_STATES == EXPECTED_OLD_POSITIVE_STATES
        and sum(value < 0 for value in OLD_DEPTH6) == 483
        and not any(value == 0 for value in OLD_DEPTH6),
        "old-wall depth-six sign census drift")


# ---------------------------------------------------------------------------
# 2. Eleven exact facets kill the complete depth-at-most-six bank.
# ---------------------------------------------------------------------------

S0 = frozenset({(13,)})
S1 = complement_upset(8, {(1,) * 8})
S2 = complement_upset(10, {(1,) * 10})
S3 = complement_upset(11, {(1,) * 11})
S4 = complement_upset(12, {(1,) * 12})
S5 = generated_upset(13, (
    (3,) + (1,) * 10,
    (2, 2) + (1,) * 9,
))
S6 = generated_upset(11, (
    (4,) + (1,) * 7,
    (3, 2) + (1,) * 6,
    (2, 2, 2, 2) + (1,) * 3,
))
S7 = generated_upset(11, ((2, 2) + (1,) * 7,))
S8 = generated_upset(12, (
    (3,) + (1,) * 9,
    (2, 2) + (1,) * 8,
))
S9 = generated_upset(13, (
    (4,) + (1,) * 9,
    (3, 2) + (1,) * 8,
    (2, 2, 2, 2, 2) + (1,) * 3,
))
S10 = complement_upset(13, {(1,) * 13})

ROW_SPECS = (
    (13, S0), (8, S1), (10, S2), (11, S3), (12, S4),
    (13, S5), (11, S6), (11, S7), (12, S8), (13, S9),
    (13, S10),
)
require(tuple(len(upset) for _, upset in ROW_SPECS)
        == (1, 21, 41, 55, 76, 99, 51, 53, 75, 95, 100),
        "mutated-wall upset-size census drift")
ROWS = tuple(response_row(degree, upset) for degree, upset in ROW_SPECS)
COEFFICIENTS = (
    436764920456040688124702805158394740410824823162360475870244123960010343376945223591925129315884594186165234915730957865748549194770202817403205624731015963923247146653830852889963929924015390519626718929519399552,
    9551799852454867057515847755276355878567471861993028124266405320805490862429008184204335626760543022099449401625305839756542506132381832835093111668856644690943854835592937780230792620027369940322405467055637355683687721755645,
    181690090865129024070142414040531237570039269619373274125889177961365196824025571003184918811598833112964833122639047586840276287574158242300310380572491280528499893595355716663719070139814157591594512992823869285034951539800,
    8982429103668787766639028057728600231554118212172089403046640439328603704399954637742515027736592766933989277549028987692528837240204299388230216981978371505594126607961393707856202627672206179408517195477493538765015048000,
    79172248975625840555721593329182474339490552128605885135094961368642170982122964749488326020538841671268762353157399451759233890364503032244985634983832840133942866436883304364285557867449975508780819119751393291123232000,
    341519987408195470267273533132603237707964318115050465525824336395982551789901713041685008426272887965508080649523677493305739754358978768621800766271924094111910096412323831077601169773870405239145023021545592884699860,
    114731005179118650713176772902537204580415667472311389212031597481812363757709170392642149190081054416440270842153326431885136975662388879751896315959655072703574261329702039607683508968715245153119324431472466113640000,
    71391041616705886039084637642084093269363743638689454233367169528826501176607680169242703293113681229424147786671866977250603898359432443771039270940570447807991113912883285652070291593893636668965150795924503875902792000,
    2081112935854219705773075390541783716634305704402824165084445274723470514894540052681996006435085206075749691101154818806875735856921190739006289672516946189456899827935813506516400990344807267086414169219880690395397000,
    1046850574416691070619177800001976239240110498075960382009355742605843054593189678673898588642858351388755112350614072139212283860792851195414238751962827539414511051692680293104797499413452308265848802621268224112640,
    1388626266511533582492503837397424937821351557825202440089217505537118806918191540873803356752805333632793398896386130623849997122180832308648736748029161872426636091643364892961513617960409455522758235919339673413959040,
)
require(all(value > 0 for value in COEFFICIENTS)
        and reduce(gcd, COEFFICIENTS) == 1,
        "depth-six Farkas coefficients lost positivity/primitivity")
COMBINED = combine(ROWS, COEFFICIENTS)
COMBINED_RANGE = (min(COMBINED), max(COMBINED))
MIN_STATES = tuple(
    state for state, value in zip(STATES, COMBINED)
    if value == COMBINED_RANGE[0]
)
MAX_STATES = tuple(
    state for state, value in zip(STATES, COMBINED)
    if value == COMBINED_RANGE[1]
)
require(all(value < 0 for value in COMBINED),
        "degree-thirteen depth-six separator lost strict negativity")
require(COMBINED_RANGE == (
    -47332088032127152443471443179776917270582014627631861565338703674497140122346501038549032754323400731113994666515396848189060964522146817633983301402670857070313337034266756845387023940839558373061259579883655633703709902609007717030259031488,
    -333069732042673888589973785127001015895756948560433115539988367179103639879784946812755225140971938122209275099364169743436382909125214539183879574694029514651992331418795357906555966338863788271102786523483721177978066454717946469619104,
), "depth-six separator range drift")
require(MIN_STATES == ((8,),)
        and MAX_STATES == ((1, 1, 1, 1, 2, 2),),
        "depth-six separator extreme-state drift")


print("THM-3169 depth-six degree-thirteen non-resurrection")
print("pole_prefix_dependency_sha256=" + UPSTREAM_LF_SHA256)
print("signed_bank_dependency_sha256=" + UPSTREAM_BANK_LF_SHA256)
print("poles=" + repr(POLES))
print("state_counts_depth1_depth6_total="
      + repr(tuple(map(len, BY_DEPTH)) + (len(STATES),)))
print("old_separator_depth6_signs="
      + repr((sum(value < 0 for value in OLD_DEPTH6),
              sum(value == 0 for value in OLD_DEPTH6),
              sum(value > 0 for value in OLD_DEPTH6))))
print("old_separator_positive_states=" + repr(OLD_POSITIVE_STATES))
print("mutated_upset_sizes="
      + repr(tuple(len(upset) for _, upset in ROW_SPECS)))
print("mutated_Farkas_coefficients=" + repr(COEFFICIENTS))
print("mutated_Farkas_range=" + repr(COMBINED_RANGE))
print("mutated_Farkas_min_states=" + repr(MIN_STATES))
print("mutated_Farkas_max_states=" + repr(MAX_STATES))
print("barcode=(C12_depth5_nonempty,C13_depth5_empty,C13_depth6_empty)")
print("scope=fixed_support_bank_averaged_virtual_prefix_currents")
print("all_exact_checks=PASS")

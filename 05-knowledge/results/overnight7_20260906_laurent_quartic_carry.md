# Endpoint 27: five first channels with the complete lower carry

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
[Independent proof and exact audit: PASS](overnight7_20260906_laurent_independent_audit.md).
For every integer `g>=14` with `gcd(g,27)=1`, and every nonzero complex
coefficient triple, a Laurent polynomial on support

```text
(-27,2g-27,3g-27)
```

has first nonzero constant-term moment exactly `g` or `2g`. Both values
are attained on every support. At every one of the four first-row scalar
cancellation phases, the doubled moment divided by the square of the first
anchor monomial is strictly negative. This is one unbounded endpoint-27
family, not all supports of smaller endpoint at most 27.

## 1. Inheritance and coverage

The method continues the independently audited
[endpoint-21 cubic family](overnight6_20260906_laurent_cubic_carry.md).
The immediate earlier mechanism is the
[endpoint-15 lower-carry theorem](trinomial_width15_empty_core_returns_sep06.md).
The canonical source for first-row roots is
[THM-4436, complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
The source tuple here is `(A,B,h,r,z,x)=(2,3,4,0,1,g-13)`.
It satisfies every positive-count hypothesis.

Before deriving, I read incoming `6f450d8cb`'s
`05-knowledge/results/nc2_hadamard_transport_overnight_hexagon_sep05.md`
from `origin/main` in full. It proves virtual negativity, a midpoint
coefficient inclusion, and a stronger transport theorem for the already
closed endpoint-15 family; its general signed transport remains OPEN.
Consequently it is not a proved dependency for the new family below.
Its specific missing-midpoint object motivates the separate
[three-group transport note](overnight7_20260906_laurent_midpoint_transport.md).

Targeted exact support/parameter searches found no existing endpoint-27
family theorem. This extends recovered repository coverage, not external
priority. All supports have smaller extreme endpoint at least 15, and
their unequal support gaps `2g,g` exclude the AP sector. The first row has
five channels, so THM-4432's two-first-channel theorem does not apply.

The corrected near miss is omission of the lower carry. Here, as in the
endpoint-15 and endpoint-21 families, the uncorrected canonical doubled
polynomial is positive at every first root, while the correctly anchored
doubled moment is negative. The retained board is: full fibers, lower
carry, quotient spectrum, shifted coefficient signs, parameter degree,
and same-constraint signed transport.

## 2. Literal source and complete moment identities

Write

```text
f(u)=alpha*u^-27+beta*u^(2g-27)+gamma*u^(3g-27),
tau=alpha*gamma^2/beta^3,
X=alpha^(g-13)*beta^12*gamma.
```

The charge equation for a multiplicity vector of mass `m` is
`g(2y+3z)=27m`. The gcd condition forces `g|m`. All channels at `g` and
`2g` are, respectively,

```text
(g-13+j,12-3j,1+2j),       j=0,...,4,
(2g-27+j,27-3j,2j),       j=0,...,9.
```

The first coordinates are positive for `g>=14`. Thus the first support
return is `g`, there are no intermediate return masses before `2g`, and
the second anchor is twice the first minus `(1,-3,2)`. The lower carry is
one; the upper carry is zero.

Set

```text
P(t)=(g-9)_4+220(g-9)_3*t+5544(g-9)_2*t^2
     +15840(g-9)*t^3+1320t^4,
L=(g)_9/12!,                      K=(2g)_18,
Qbar(t)=sum_(j=0)^9 (2g-18)_(9-j)*t^j /[(27-3j)!(2j)!].
```

Literal multinomial coefficients give

```text
CT(f^g)=X L P(tau),
CT(f^(2g))=X^2 K tau^-1 Qbar(tau).                 (1)
```

Both scale factors `L,K` are strictly positive. By the explicitly eligible
THM-4436 tuple above, `P` has four distinct strictly negative roots at every
integer parameter under consideration. That canonical theorem supplies
the real spectrum; positive coefficients alone would not establish it.

## 3. Four positive characteristic coefficients

Let `p=P/1320`, and form the quotient algebra `Q(g)[t]/(p)`. Define

```text
R=t^-1 Qbar mod p=sum_(j=0)^3 R_j t^j,
C=the monic companion matrix for multiplication by t,
V=sum_(j=0)^3 R_j C^j,
det(zI-V)=z^4+B1(g)z^3+B2(g)z^2+B3(g)z+B4(g).     (2)
```

The only apparent variable denominator in `R` cancels exactly:

```text
Qbar(0)/p(0)
 =42240(g-13)(2g-19)(2g-21)(2g-23)(2g-25)/27!.
```

After reduction, `deg_g R_j<=8-j`. The four characteristic coefficients
have degrees `8,16,24,32`. The exact output publishes every coefficient of
`B_k(s+14)` in a positive-denominator integer convention. There are
`9+17+25+33=84` coefficients, and every one is strictly positive.
Appendix A repeats these complete certificate arrays.

This is an identity certificate, not a finite parameter census. The
producer constructs the complete symbolic quotient over `Q(g)` and checks
the characteristic coefficients independently through Newton's identities

```text
k B_k=-sum_(r=1)^k B_(k-r) Tr(V^r),       B_0=1.
```

Consequently the polynomial in `(2)` is positive at every real `z>=0`
when `g>=14`. At an admissible integral `g`, evaluation at the four distinct
negative first roots conjugates `V` to their four real response values
`R(lambda_i)`. None can be nonnegative. Hence

```text
R(lambda_i)<0,       i=1,2,3,4.                    (3)
```

The signs in the coefficient arrays hold in the real parameter extension;
the first-return and canonical real-root inputs are invoked only for the
stated integer source domain. There is no hidden real-parameter claim about
Laurent moments.

For an independent identity audit, assigning weight one to both `g,t`
shows that reduction by monic `p` preserves weighted degree. Every term of
`Qbar/t` has weight at most eight after the displayed inverse cancellation.
Newton sums therefore give `deg B_k<=8k`. Thirty-three distinct rational
parameter values suffice to certify all four identities by interpolation
once these degree bounds are retained.

## 4. Exact consumer, attainment and hostiles

Combining `(1)` and `(3)` gives `CT(f^(2g))/X^2<0` at every zero of the
first scalar moment. This quotient is real even when `X^2` is complex.
The first nonzero moment is therefore exactly `g` off the four first roots
and `2g` at each of them. Positive coefficients realize the first case;
`beta=gamma=1, alpha=lambda_i` realizes each cancellation case. The worst
order `2g` is below total width `3g`, so this family does not saturate the
linear-width conjecture.

There is no equality or common-root case in the domain. An anchor change
by `k(1,-3,2)` multiplies the normalized doubled response by
`tau^(-2k)>0` at first roots, so the sign is anchor invariant. Omitting the
lower carry reverses it: `Qbar(lambda)=lambda R(lambda)>0`.

Dropping the gcd condition changes the first-return interpretation.
At `g=15`, the support `(-27,3,18)` has a return at mass five, not fifteen.
Its balanced channel is `(1,3,1)`. The displayed rows at masses 15 and 30
remain correct identities, but cease to be the first two support rows.

THM-4440's ordinary real-rooted-core theorem does not already cover the
new family's entire phase locus. At `g=14`, `P(-59)>0` and `P(-58)<0`,
so there is a first root in `(-59,-58)`. For the real ordinary cubic
`1+beta*s^2+s^3` with `beta^3=1/tau`, its discriminant is
`-4/tau-27<0` there. A nonreal conjugate pair remains outside the real-core
SOS source domain. This is a non-subsumption witness, not an SOS hostile.

## 5. Verification and remaining general mechanism

The companion imports no repository producer. Six named parameters
`14,16,17,19,23,29` are reconstructed by both the original charge equation
and repeated Laurent multiplication retaining the gamma exponent. They
verify both entire coefficient rows and every earlier/intermediate empty
return, specialized carry inversion, four-root reality and complete
negative-definite weighted trace forms. The generic coefficient identities
are symbolic; these finite controls test the source-to-consumer map.

Normal and optimized outputs match with **98 explicit gates**:

```text
python -B 04-computation/overnight7_20260906_laurent_quartic_carry.py
python -B -O 04-computation/overnight7_20260906_laurent_quartic_carry.py
```

```text
source 3ab5a563fbb0a8f1dfceb35f348ff9d8e2175265f41f3e3fd2afdc614257e2e3
output 316e1bff9b8b4a3a445eb6d6c80570049bea0b1f9686851b70b4f89e89f39674
semantic d6f03cfe59f5a7374558b9bede55afc3c337f7c69f7112d99950655b245f8e3b
```

There is a general algebraic family, rather than a need to rediscover the
normalization at each endpoint. For `a=6h+3`, `h>=1`, `g>=3h+2`, the first
and second counts have `h+1` and `2h+2` channels. After removing positive
contents, `deg R_j<=2h-j` and characteristic degree bounds are `2hk`.
These are proved in the transport companion. They do not prove positivity
uniformly in `h`. The new three-group midpoint identity and its explicit
same-root-condition hostile give the next structural object to pursue;
no larger fixed-h census was run.

## Appendix A. Complete positive-shift identity certificate

For each `k`, the displayed descending coefficient list divided by its
denominator is exactly `B_k(s+14)`. The producer uses an integer polynomial
and the least common multiple of its rational coefficient denominators;
all listed coefficients are integers and strictly positive.

```text
B1 denominator = 74860977471626171105280000000
coefficients = [
  19451826205848692368573,
  787455838384465942339138,
  13946297132868416991900238,
  141137507095522233828485428,
  892678153196839055968496413,
  3613402243640892028666948810,
  9141274352770790905525181400,
  13214422476190843062998040000,
  8357123003336257784266080000
]
```

```text
B2 denominator = 4981480842673174326979650300930531179416780800000000000000
coefficients = [
  15947543841494344912477006491749447,
  1214895534270251462895467773479034724,
  43364125314359075001171767629306940036,
  962653992244482089826819681153721861328,
  14875995824127975414090465233998621716410,
  169678385767526930702641093500723168531368,
  1477716286368299762086280206178553483071716,
  10023215245198401855787991093857585891290304,
  53513282061037716135112128223874061167875687,
  225629587459270423219783194867765913051175380,
  748795147850295596337228207992483160547708000,
  1935394089105442687949088227874884047111766800,
  3819282637155943644559841202715787609380846800,
  5562816856930741175382235980876727809598992000,
  5639672221356385939369360127574307517915040000,
  3555698107969668265521614023487939746214400000,
  1050268635132476223834534591824133203712000000
]
```

```text
B3 denominator = 29833482011095508645861359680923986616910739336740661153545238609920000000000000000000
coefficients = [
  142895213635567973672728962843696494038505,
  15351496454492439459911058637608625469299158,
  789662023476822504821184243479456928612041234,
  25882252366891552498389205467946474300287132756,
  606852435722302327713517776520049805488039079287,
  10832349179155285297986426682734169678058013719586,
  152952902889532410442015796926709672615540960589812,
  1752319253169565461415204486973341942754978556727752,
  16576640208430136110710973962764890605310990973022167,
  131079995679500434502305541406906965704746298558526874,
  873821628046796431930648845675192932384611619178727010,
  4938332287775313664902832431318460845308456102736342228,
  23734812048588870444620948192912762050059369807939928409,
  97114081678668406831530732646856711923290046438140703182,
  337920226224064224073909652116254757093653361019596711544,
  996764336225168642210985865195071775888602727209335607504,
  2478448856124771984464618616039528266717044432250479897232,
  5150794428936517936797932009990195173737599684455793921760,
  8837588086321068050551058109105959947380964825528576320000,
  12300944984647840977604451581184648345383655877145211904000,
  13540994185313302806276782979839971152237343900143682560000,
  11345402196887128629651504762613080412506709465201075200000,
  6798215508043276119745699946461782380983275022906368000000,
  2594849712203075520428415491267284507097084609372160000000,
  474072997824590055722621175826266672444748016025600000000
]
```

```text
B4 denominator = 107201453987173704440272878544104053485493191964989010053189624956661145230232396984274124800000000000000000000000000
coefficients = [
  24326701156875917876749764204876467618869586903963,
  3288033209170910304416296354087715126247661034370248,
  215009153192899043333167832951734894473736926882717736,
  9059940928564451123492622710302562719734418985910620608,
  276442176297320527280043734613331126764400519991363490900,
  6507269103696828944492319843248598845808210435416603925120,
  122935143032509208444233935725745953587033267147269371567256,
  1914551778636322018227334418137920956372256943756416253066080,
  25053958605846022443035032422157707623949138261751149046634466,
  279408779284789832337213782930989884735042985384264579801113936,
  2684050410283126518512435649153713636484449543734964173827034008,
  22390466507906819335016193606532762423573640707896669897621833184,
  163213098197578004160416018860981684028242991317908852419467796532,
  1044475403116173888584729755698730264809014123634552988926624285504,
  5888087392104331966259936190684528713567646883990476846507690593448,
  29308307555995938930256613764571042331176303338810594897035046309152,
  128982353947709691772247118905159733660324042154204490739837493109787,
  502085441933110531159478747634005224396510632919436807063413331453352,
  1727938336489587112519774888153579903338002352331213371724330572710272,
  5250246458158332795393054314105529047548473039544795170063012331534816,
  14050908088959558522124914231724508258930584172925354732497201214235232,
  33006597288844988251386544911553976393213127375428826608060987490320000,
  67735523174538017906817107413081159896609582552223286264364634544102400,
  120679481330099030988172018164197562141232581043323962197978448998464000,
  185138049533909768613747979812249468335730035246695340122639209182304000,
  241963217076078459648074044712228610044715065360733223827819337019904000,
  265604989249537488386222249496876970853848081639581839698625018188800000,
  240231192547473177722292878493380948402266994306659858144196301004800000,
  174290877217420211907191029029288330279544596318249651450836007936000000,
  97494619626451907552412858502850330275086810165473440763900559360000000,
  39465476877277215134742329485762906213469679126610005476209459200000000,
  10287167158359218876463719396080861344095563095812196750852096000000000,
  1296309473866388503919864896633677326166561438301533626695680000000000
]
```

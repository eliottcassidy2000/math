import ProcgenSelfieEdim.ListLemmas
import ProcgenSelfieEdim.Hypercube
import ProcgenSelfieEdim.EdimQ6
import ProcgenSelfieEdim.CountingBound
import ProcgenSelfieEdim.Stabilizer
import ProcgenSelfieEdim.KeyBlocks
import ProcgenSelfieEdim.EdimQ7Data
import ProcgenSelfieEdim.EdimQ7P0
import ProcgenSelfieEdim.EdimQ7P1
import ProcgenSelfieEdim.EdimQ7
import ProcgenSelfieEdim.EdimQ8Data
import ProcgenSelfieEdim.EdimQ8P0
import ProcgenSelfieEdim.EdimQ8P1
import ProcgenSelfieEdim.EdimQ8P2
import ProcgenSelfieEdim.EdimQ8P3
import ProcgenSelfieEdim.EdimQ8
import ProcgenSelfieEdim.EdimQ9Data
import ProcgenSelfieEdim.EdimQ9P0
import ProcgenSelfieEdim.EdimQ9P1
import ProcgenSelfieEdim.EdimQ9P2
import ProcgenSelfieEdim.EdimQ9P3
import ProcgenSelfieEdim.EdimQ9P4
import ProcgenSelfieEdim.EdimQ9P5
import ProcgenSelfieEdim.EdimQ9P6
import ProcgenSelfieEdim.EdimQ9P7
import ProcgenSelfieEdim.EdimQ9P8
import ProcgenSelfieEdim.EdimQ9P9
import ProcgenSelfieEdim.EdimQ9P10
import ProcgenSelfieEdim.EdimQ9P11
import ProcgenSelfieEdim.EdimQ9
import ProcgenSelfieEdim.Tournament
import ProcgenSelfieEdim.HamPath
import ProcgenSelfieEdim.ArcParity
import ProcgenSelfieEdim.AltSum
import ProcgenSelfieEdim.TournamentCode
import ProcgenSelfieEdim.SelfieFinite
import ProcgenSelfieEdim.Shaved
import ProcgenSelfieEdim.CycleCount
import ProcgenSelfieEdim.Circulant
import ProcgenSelfieEdim.HPExist
import ProcgenSelfieEdim.SwitchSum
import ProcgenSelfieEdim.Redei
import ProcgenSelfieEdim.ConstantH
import ProcgenSelfieEdim.ParityBreak
import ProcgenSelfieEdim.CollatzDrop
import ProcgenSelfieEdim.PolyNorm
import ProcgenSelfieEdim.DropZ
import ProcgenSelfieEdim.DropPair
import ProcgenSelfieEdim.DropUnits
import ProcgenSelfieEdim.AntiAut
import ProcgenSelfieEdim.CayleyAbelian
import ProcgenSelfieEdim.AntiCirculant

/-!
Explicit library root: every theorem of the package is reachable from here.

THM-4525 (edge multiset dimension): `Hypercube`, `EdimQ6`, `CountingBound`, `Stabilizer`,
`KeyBlocks`, `EdimQ7*`, `EdimQ8*`, `EdimQ9*` (explicit upper bounds for `Q_7`, `Q_8`, `Q_9`), `AltSum` (L3).
THM-4524 (selfie tournaments, arc-HP parity): `Tournament`, `HamPath`, `ArcParity`,
`TournamentCode`, `SelfieFinite`, `Circulant`, `HPExist`, `SwitchSum`, `Redei` (Rédei's theorem and its corollaries), `ConstantH`, `ParityBreak`.
THM-4526 (shaved tournaments): `Shaved`, `CycleCount`. THM-4527 / S15 Theorems 1-2: `CollatzDrop`.
Round 2. THM-4530 (drop multiplicities): `PolyNorm` (a small `Int` polynomial normaliser),
`DropZ` (two copies over Z, Prop. 2.6), `DropPair` (pair-injectivity, Theorem 4.5, both sheets),
`DropUnits` (unit branches of `qx + 1`, Theorem 5.2). THM-4529 (Paley coordinates,
anti-circulants): `AntiAut` (anti-automorphisms, the reflection lemma, `N ≡ 2 mod 4`),
`CayleyAbelian` (Theorem 5.1 for every finite abelian group of odd order), `AntiCirculant`
(Theorem 6.1, antipodal arcs odd).
-/

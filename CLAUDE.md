# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

**Atlas Orbitalis** (Atlas) is a Fortran toolkit for satellite orbit computation and analysis by Shawn Hu (Nanjing University, MIT-licensed): orbit propagation (RKF78 numerical, analytic two-body), observation simulation (angles, ranging, Doppler, GNSS), IOD (General Laplace), POD (batch least-squares, EKF, sequential filter), multi-body dynamics (CRTBP, Halo, Lissajous, Lambert) and Kalman filtering.

**User manual:** `manual/atlas_manual.md` (Chinese, PDF alongside; `data/doc/` holds git-ignored algorithm notes). **The manual lags the code** — trust the source; the known drift is listed under [Known gaps](#known-gaps-traps-and-stale-artifacts).

---

## Build

```bash
make                     # library + every app and test (auto-discovered from app/, test/)
make build/libatlas.a    # library only
make bin/test_CRTBP      # one target (same as bin/test_CRTBP.exe on MinGW)
make clean               # remove build/ and bin/
```

`gfortran -O2 -std=gnu -fall-intrinsics -fPIC -fopenmp`; gfortran ≥ 9 or Intel `ifx` ≥ 2021. Use `-g -Og` for debugging.

- **Run every executable from the repository root.** All runtime data paths are relative (`data/AtlasData/…`); there is no path resolution or `chdir`, so running from `bin/` silently fails to find the ephemeris, gravity or EOP files.
- **SOFA is built from source** (`sofa/`, 249 F77 files → `build/libsofa.a`, linked by every target); first build takes ~50 s, then it is cached. `lib/libsofa.{a,lib}` are dead binaries (wrong object format / Intel symbol mangling) — delete them once the VS build no longer needs them.
- On MinGW the Makefile names targets *with* `.exe`; the extensionless `make bin/xxx` forwards to `bin/xxx.exe`. Do not "fix" this — without the suffix every `make` relinks everything.
- `-fopenmp` is **inert**: all OpenMP calls are commented out, and the only `use OMP_LIB` is commented out too (this machine's MinGW gfortran has `libgomp.a` but no `omp_lib.mod`).
- Visual Studio + Intel Fortran (`Atlas.sln` / `Atlas.vfproj`) is a **stale secondary** build: it omits `module_shadow.f90` and registers only one test.

### ⚠️ Changing a type that other modules `use` requires `make clean`

Library objects carry **no `.mod` dependency edges** (`build/%.o: src/%.f90` and nothing else); when `src/` changes, the Makefile recompiles all app/test sources but *not* the library modules above the change. So a module embedding a changed type keeps a stale `.mod`: the next build either fails on an unrelated target (`Mismatch in components of derived type … expecting 'x', but got 'y'`) or links objects that disagree about the layout. After touching `type_access`, `type_clock` (embedded by value in `type_satellite`/`type_station`), `type_satellite`, `type_ODor`, …: run `make clean` first. Adding `-MMD -MP` plus `-include $(OBJ_FILES:.o=.d)` would fix this properly.

---

## Running tests

No framework: every `test/*.f90` is a standalone program that prints results for human inspection (some also write `data/output/` files). Build with `make`, run from the repo root (`./bin/test_elements.exe`). All 16 tests build and pass as-is.

| Test | Covers |
|------|--------|
| `test_CRTBP`, `test_halo`, `test_LissajousOrbit2` | CRTBP systems/Lagrange points, 3rd-order Halo, Lissajous (writes the synodic ephemeris) |
| `test_fastPrediction` | fast propagator vs `NumOP`: speed, accuracy, STM |
| `test_c2t06a`, `test_sp3`, `test_OS` | IAU2006 CIO transform (SOFA), SP3 parsing, OS detection |
| `test_moon_coordinate_transformation`(2), `test_lunar_gravity` | LCF↔LCI round trip, Moon-centred frames, lunar gravity field |
| `test_DXDR_nsg2` | NSG2 Jacobians: mutual inverses + finite-difference verification |
| `test_elements`, `test_thruster` | `type_orbit`, thruster ledger |
| `test_jpl_fallback` | JPL ephemeris failure → analytical Sun/Moon fallback, warning-once, km/km/s contract |
| `BDS_orbitPrediction`, `orbit_prediction_ECOM` | BeiDou MEO 1-day prediction vs STK; ECOM SRP model |

`test_thruster` needs neither a scenario nor the integrator. Six tests and the shadow plot script were deleted in Oct 2026 (`test_JPLEPH`, `test_LissajousOrbit`, `test_accessIO`, `test_clockIO`, `test_shadow`, `test_verify_DXDR_nsg2`, `visualize_shadow.py`) — recover from git history if needed; `module_shadow` currently has **no** test exercising it.

---

## Architecture

### Module dependency hierarchy (7 levels)

Defined by `LEVEL0`…`LEVEL6` in the [Makefile](Makefile); a module may only `use` modules from lower levels. `module_Atlas` (Level 6) is the public API: `use module_ODor, module_IODor, module_sp3` only — everything below arrives transitively.

| Level | Role | Modules |
|-------|------|---------|
| **0** | Infrastructure | `module_const`, `module_toolKits`, `module_ErrorHandler`, `module_progressBar`, `module_rotator`, `module_randum`, `module_String`, `module_FileDef`, `module_elements`, `module_MathsMatrix`, `module_sp3`, `module_EarthRadiationPressure`, `module_shadow` |
| **1** | Dynamics & time | `module_CRTBP`, `module_lambert`, `module_time`, `module_EOPParam` |
| **2** | Ephemerides & special orbits | `module_matrix`, `module_planets`, `module_halo`, `module_LissajousOrbit`, `module_timeSeries`, `module_clock` |
| **3** | Coordinates & attitude | `module_space`, `module_SLE`, `module_attitude`, `module_cartesian` |
| **4** | Forces & integrator | `module_Gravity`, `module_atmosphere`, `module_perturbationPN/SRP/ThirdBody`, `module_thruster`, `module_orbitPredictionNumerical` |
| **5** | Domain objects | `module_station`, `module_satellite`, `module_access`, `module_residual`, `module_satelliteScenario`, `module_IODor`, `module_ODor` |
| **6** | Public API | `module_Atlas` |

`module_EarthRadiationPressure` and `module_shadow` sit at Level 0 because both are self-contained (no Atlas `use` at all). `module_STM` (analytic STM, NSG2 only) has no dependencies, is **not in the Makefile** and is never compiled — add it to a `LEVEL` before using it.

Legacy Fortran 77, compiled separately: `PLEPH.for` (JPL ephemeris reader), `nrlmsise00_sub.for` + `nrlmsise00_driver.for` (atmosphere).

#### Module names that differ from their file names

`use` names the *module*, not the file — these are the traps:

| File | Module |
|------|--------|
| `module_lambert.f90` | `module_LambertSolution` |
| `module_EOPParam.f90` | `EOPParam` |
| `module_EarthRadiationPressure.f90` | `EarthRadiationPressure` |
| `module_ErrorHandler.f90` | `ErrorHandler` |
| `module_FileDef.f90` | `FileDef` |
| `module_MathsMatrix.f90` | `MathsMatrix` |
| `module_String.f90` | `String` |
| `module_progressBar.f90` | `CMD_ProgressBar` |

### Unit system

Internal computation uses **normalized units** from Earth parameters: length = 6378.1363 km, time = 806.8109913067326 s, mass = 5.9742272407074334×10²⁴ kg, acceleration = 9.79828762253515 m/s². `module_const` defines the SI magnitudes (`LENGTH_UNIT`, `TIME_UNIT`, `MASS_UNIT`, `ACC_UNIT`, `FORCE_UNIT`) and the conversion factors `u_km`, `u_meter`, `u_sec`, `u_day`, `u_kg`, `u_m_s2`, `u_Newton`, …

> `u_X` converts **from unit X into normalized units** (`u_km = 1/6378.1363`), so `r_km * u_km` is normalized. Manual appendix A.3 has this inverted.

File I/O converts at the boundary (users work in km, km/s, degrees). Exceptions to keep straight: `type_cartesian` stores **SI**; `type_access%dataTable` is **normalized** (`readData`/`saveData` convert); perturbation routines mix units internally (e.g. `solar_pressure_perturbation` returns `F = F*acc_unit` while the force accumulator ends with `force = force/acc_unit`) — **check every F's unit at each call site**. CRTBP uses its own unit system (length = primary-secondary distance, time = period/2π).

---

## Key types and conventions

### Time: `type_t` ([module_time.f90](src/module_time.f90), Level 1)

MJD day + seconds of day; fields `ttype`, `day`, `sec`, `val` (`TIME_TYPE_UTC`=1, `TT`=2, `UT1`=3, `BDT`=4, `GPST`=5).

```fortran
call t%init(2025d0, 1d0, 1d0, 12d0, 0d0, 0d0)   ! calendar (UTC by default)
call t%init(51544.5d0)                           ! MJD input
call t%init(51544d0, 43200d0)                    ! MJD day + seconds
t%toTT(); t%toUTC(); t%toUT1(); t%toGPST(week, sow); t%toBDT(week, sow)
t%toDate()                                       ! "YYYY-MM-DD HH:MM:SS.sss" (UTC)
dt = timeMinus(t1, t2)                           ! days
t2 = timeAdd(t1, dt_days)
```

**Two independent leap-second tables** — update both: `module_time%UTC2TT_daySec` is hard-coded and ends at MJD 57754 (2017-01-01, ΔLS = 37; later dates raise an error), while `module_EOPParam` reads `data/AtlasData/leap.txt`. UT1 conversions print an "unfinished part" warning.

### Orbital elements ([module_elements.f90](src/module_elements.f90), Level 0)

`STATE_TYPE_KEPLERIAN`=1 `(a,e,i,Ω,ω,M)`, `STATE_TYPE_CARTESIAN`=2 `(x,y,z,vx,vy,vz)`, `STATE_TYPE_NONSINGULAR_1`=3 NSG1 `(a,i,Ω,ξ,η,λ)` (integrator default), `STATE_TYPE_NONSINGULAR_2`=4 NSG2 `(a,h,k,p,q,λ)` (needed for the analytic STM). All interconvert through Keplerian (`kpl2cts`, `kpl2nsg1`, `kpl2nsg2`, `nsg12nsg2`, …); conversions are normalized (μ = 1), `*_SI2normal` / `*_normal2SI` handle the boundary. `OP2Body` propagates two-body in NSG1.

Prefer the OOP wrapper: `type(type_orbit)` holds `elements(6)` + `eleType` and stores normalized units — `call o%init(x0, STATE_TYPE_KEPLERIAN)` (input assumed SI unless `isNormalized=.true.`), `o%toCartesian()/toNSG1()/toNSG2()/toKeplerian()`, `o%toSI()/toNormalized()`, `o%get_rEf(E,f)`, `o%get_rut(u,utilde)`, `o%propagateTwoBody(dt)`, `o%show()`.

### Cartesian: `type_cartesian` ([module_cartesian.f90](src/module_cartesian.f90), Level 3)

Position, velocity, frame, epoch — in **SI**. `init(position, velocity, COORD_ICRS|COORD_ITRS|COORD_HORIZON, t)`, `ICRS2ITRS()/ITRS2ICRS()` (velocity included), `ITRS2Horizon(lon,lat)`, `getReferenceName()`. Frame mismatches are checked; velocity transforms cost ~2× a position transform and are skipped when velocity is zero; horizon-frame velocity is **not** implemented (zeroed).

### Satellite: `type_satellite` ([module_satellite.f90](src/module_satellite.f90), Level 5)

```fortran
call sat%init(t0, x0, number=1, name="MySat", mass=500d0, eleType=1, propagator="numerical", degree=8, order=8)
sat%opMod%key_moonGravity = 1     ! switches default to 0 — see "Force model system"
call sat%generateTrajectory(epochs)            ! or generateTrajectoryTwobody(epochs)
call sat%savetrajectory(STATE_TYPE_CARTESIAN, "out.txt")   ! readtrajectory(...) to load back
pv = sat%getPV(t);  state = sat%getStates(t)   ! interpolated J2000 PV / elements
call sat%setTrajectory(eleType, nData, mat);  call sat%freeTrajectory();  call sat%show()
```

`propagator`: `'numerical'` (RKF78 + perturbations), `'two-body'` (analytic), `'SGP4'` (placeholder). `eleType` defaults to NSG1 (3). `mat` for `setTrajectory` is `(nData, 8)` = MJD-day, MJD-sec, 6 elements.

### Station: `type_station` ([module_station.f90](src/module_station.f90), Level 5)

```fortran
call sta%init(number=1, name="Beijing", longitude=116.4d0, latitude=39.9d0, altitude=50d0)
!                          degrees, degrees, metres → rad / normalized internally
call sta%generateTrajectory(ts)     ! over a type_timeSeries
pv = sta%getPV(t)                   ! normalized J2000 PV — calls stop if t is outside the span
```

Geodetic↔ECF helpers in the same module: `Geodetic`, `geo2pos`, `pos2geo` (WGS-84). `sta%clock` is a `type_clock` that nothing fills automatically (assign it yourself); `calObsData1` applies it as the receiver clock, and an unset clock contributes exactly zero.

### Visibility & observations: `type_access` ([module_access.f90](src/module_access.f90), Level 5)

Holds the link identity (`number`, `fromType/fromNum`, `toType/toNum`, `obsType`), the generated data (`timeTable(:)` of `type_t`, `dataTable(:,:)` in **normalized** units, `nData`, `npe`), per-column `bias(10)`/`errStd(10)`, the constraints (`minEle`/`maxEle`, `solarDistance`, `lunarDistance`, `shadow`, `EarthBlockage`, `MoonBlockage`, `pointTo`, `t_coverage`, `beamHalfAngle`) and the informational device names `fromName` (transmitter) / `toName` (receiver), which `type_satelliteScenario%observationAccess` fills from the scenario.

```fortran
call ac%init(number=1, name="ranging1", fromType=OBJECT_TYPE_SATELLITE, fromNum=1, &
             toType=OBJECT_TYPE_STATION, toNum=1, obsType=OBS_TYPE_RANGING, minEle=10d0)
call ac%saveData(ts, fn)    ! ts is unused (vestigial)      call ac%readData(fn)
call ac%show()
```

| `OBS_TYPE_*` | Value | Observable | npe |
|---|---|---|---|
| `RA_DEC` | 11 | RA, Dec | 2 |
| `AZIMUTH_ELEVATION` | 12 | Az, El | 2 |
| `RANGING` | 13 | Range | 1 |
| `RANGERATE` | 14 | Range-rate | 1 |
| `POSITION` | 15 | x, y, z | 3 |
| `POSITION_VELOCITY` | 16 | Position + velocity | 6 |

`setNpe()` derives `npe` from `obsType`. Generation functions (called by the scenario): `calObsData1(ts, sat, station, access, force)`, `calObsData2(ts, sat1, sat2, …)`, `calObsData3(ts, sat, …)`. Correction flags — mind the defaults, they are **not** all off: `key_light_time_delay` = `.true.`, `key_atmospherical_refraction` = `.false.`, `key_clock_error` = `.true.`. `AzimuthElevation2RaDec(station, access)` converts Az/El to RA/Dec.

For RANGING, `key_clock_error` makes the observation follow the receiver-timed model `rho_obs = rho_true + (clk_R - clk_T)·c` and `t_obs = t_true + clk_R`, taking `clk_R`/`clk_T` from `sta%clock`/`sat%clock` (satellite–satellite links use both satellites' clocks) — so an *unset* clock contributes exactly zero. Both sat–station (`calObsData1`) and sat–sat (`calObsData2`) implement it; a station link whose receiver clock is zero keeps the ideal epochs. **No measurement noise is ever added** (`bias`/`errStd` are informational — `calW` in `module_ODor` returns the identity).

**Observation file format** (`saveData`/`readData`): a self-describing `#` header followed by `MJD_day  MJD_sec  obs1 … obsN` lines (SI in the file).

```
# generator   : Atlas (Atlas Orbitalis)     # transmitter : satellite #1 (LEO-600km)
# access      : 1 (ranging1)                # receiver    : station #1 (station1)
# data type   : RANGING (obsType = 13)      # data format : MJD_day  MJD_sec  Range/km
```

`readData` skips `#`/`*`/`!`/`！` and blank lines and **adopts the header's `obsType`** (fixing `npe`; a conflict with the object is reported and the header wins), plus the device names. The header's type name is checked against the integer: a stale name/int pair (pre-redesign `GNSS` = 14, `AZIMUTH_ELEVATION` = 16, `RADAR` = 12) is reported and ignored rather than silently reinterpreting the old columns. Reading is fault tolerant: `readDataMatrix` reports and **skips** malformed data lines (never zero-fills) and returns the count in its optional `nRead`, which `readAccessData` and `%readtrajectory` use to shrink their tables. Headerless legacy files still read (lines skipped one by one, `obsType` kept from the object).

### Scenario: `type_satelliteScenario` ([module_satelliteScenario.f90](src/module_satelliteScenario.f90), Level 5)

```fortran
call ss%init(t0, te, step, eleType, introduction)   ! builds its own type_timeSeries — NOT init(ts)
call ss%addSatellite(sat);  call ss%addStation(sta);  call ss%addAccess(acc)
call ss%observation(accIdx);  call ss%reportAccess();  call ss%freeSS()
```

There is **no `generateObservations()`** (manual only) — the entry point is `observation(accIdx)`.

### Orbit determination (Level 5)

`type_ODor` ([module_ODor.f90](src/module_ODor.f90)) drives `OD_METHOD_BATCH` = 1 (batch least-squares), `OD_METHOD_EKF` = 2, `OD_METHOD_SQUENTIAL` = 3: `solve(maxIter, convergence)`, `invoke(...)`, `saveEph(fn)`, `SS2ODor(scenario)`, `ODor2SS(scenario)`. `type_varGroup%label` selects the parameter set: 1 states (6), 2 + `amrDrag` (7), 3 + `amrSRP` (7), 4 + both (8), 5 + the **first `m = nVar − 6` ECOM parameters** (so `nVar` alone decides how many are estimated: 11 → 6 states + ECOM(1:5); 15 = 6 + 9 is just the maximum, `idx(7:nVar)` is what selects them) — for labels 2–3 the estimated quantity is the **area-to-mass ratio**, not the `Cd` the manual names. Label 5 requires `SRP_Model = 2`: under cannonball the ECOM parameters do not enter the force model, so `numericalSTM_ECOM`'s finite differences give identically zero columns. Labels **11/12/13** are **clock groups** — the first `nVar = 1/2/3` coefficients `(bias0, drift, drift2)` of one object's clock polynomial, referenced to `type_ODor%t`, with `type_varGroup%objType` (`OBJECT_TYPE_SATELLITE`/`_STATION`) + `origin` naming the object. They carry no orbit state, so `calB` short-circuits them (`calB_clock`) *before* the satellite-participation test and fills direct columns `±SPEED_OF_LIGHT*u_meter·dt^k` (receiver +, transmitter −) instead of going through the STM; `buildSTMseries` deliberately gives them no series. Because range data only sees `clk_R − clk_T`, one reference clock must stay fixed by the application (the library adds no constraint); `navSat_OD` fixes the station clock at zero and estimates the three satellite clocks.

`type_ODor%showProgress` (default `.false.`) turns on an iteration progress display: the trajectory phase reuses `sat%opMod%showProgressBar` (temporarily enabled, restored at the end) and the observation sweep reports processed epochs through the global `Progress` bar, followed by a one-line per-iteration summary (epochs + residual RMS + elapsed). `CMD_ProgressBar%Put` gained an optional trailing-text argument for that. The per-link `RMS` lines still go to `data/output/log_batchOD.txt` unconditionally but are suppressed on screen **only** while the display is on, so other OD apps keep their old output. Both progress bars refresh every ~1 % (a per-epoch refresh is a 2.4 MB wall of text per navSat run — 6.6 k updates per iteration — for a 40-cell bar that shows nothing new; `generateTrajectoryNumberical` throttles with `nStep = nEpochs/100` and always redraws the last epoch). Cost split measured for a navSat iteration (~1.9 s): trajectory generation ~1.7 s — `generateTrajectoryNumberical` calls `NumOP` once per *output epoch* (2881 per satellite), so per-call overhead dominates, not the integration (a single 1-day MEO call is only 25 ms) — observation sweep ~0.15 s, STM series ~0.05 s. `numericalSTM_ECOM` is expensive (15 finite-difference propagations per call, and the batch solve used to call it once per *observation* — hours per iteration for a 1-day 30 s arc), so `type_ODor` carries `STMseries(:)`: one `type_STM` per group holding `dX(t_k)/dp` on a grid of `stmStepSec` (default 1800 s), rebuilt at the top of every iteration and read back by 4-point Lagrange interpolation in `calB` (measured interpolation error 7e-4 relative; the underlying finite difference itself is good to 0.3%). `calB` also now handles a group whose satellite is the *receiver* of a link (sat-sat ranging) by swapping the two ends in `calDYDR` — it used to `stop` with "unfinished part".

Residual RMS is printed in **display units**, not normalized ones: `resUnitCode(obsType, jc)` maps each observation *component* to 1 = arcsec, 2 = m, 3 = m/s, 4 = raw (normalized), and `resRMSstr(nAcc, str)` pools `res%acc(1:nAcc)` into one RMS per unit group, e.g. `RMS = 6.99E+00 m, 6.40E-03 m/s` for POSITION_VELOCITY (position and velocity are never added into a single number — that sum has no physical meaning) or `RMS = 6.96E-02 arcsec` for Az-El. The convergence test still uses the normalized `res%RMS` from `calRMS`, untouched. Note the log's per-access line is *cumulative* over the accesses processed so far (pre-existing behaviour); the last one is the overall RMS shown in the iteration summary.

Batch OD is already **per-access type-agnostic**: the `orbitDetermination_batch` access loop reads each access's own `npe`, reallocates `B`/`Y`/`W` when it changes, calls `calB`/`calDYDR` for that type, and accumulates all rows into the same `SLE`; `res%acc(iac)` and `saveData`/`calRes` are per access. Mixed types therefore need no structural change (only the label-3 SRP column was missing for Cartesian states). Caveats: `calW` still returns identity (`errStd` is informational), so mixed observations are equally weighted in normalized units, and `res%RMS` is a type-agnostic scalar stopping metric — `resRMSstr` is the physical, per-unit display. `saveResiduals` writes `_accN.txt` when there are several accesses.

`type_ODor%logFile` (default `data/output/log_batchOD.txt`) redirects the batch iteration log, and `nIter` reports the iterations performed; `computeResiduals(ss)` recomputes final O−C residuals at the estimated state (the solve loop leaves them one update stale) and `saveResiduals(fn)` writes them per access. `stm_series_init` now builds the numerical **amrSRP column** for a Cartesian label-3 group (`dX/dx7`, via an `amrSRP` finite difference); it used to be zero, which made that parameter rank-deficient. `module_residual%calRes` wraps the periodic first residual component for RA/Dec and Az/El.

`type_IODor` ([module_IODor.f90](src/module_IODor.f90)) — General Laplace: `init(maxDataLength)`, `solve(initialState)` (J2000 Cartesian, normalized), `show()`, `saveIODor(fn)`.

`type_SLE` ([module_SLE.f90](src/module_SLE.f90), Level 3) — "System of Linear Equations", a generic solver (**not** Satellite Laser Ranging): `initialize(n)`, `addData(m, n, B, Y, W)`, `solve()`, `propagate(Phi, Q)`, `dataFusion(...)`, `freeSLE()`; `solve` zeroes `BTB`/`BTY` after extracting `x` and `P`. `module_STM` — analytic state-transition matrices, valid **only for NSG2**.

### Gravity field: `GravityField` ([module_Gravity.f90](src/module_Gravity.f90), Level 4)

```fortran
call gf%init("EarthGravity", "Earth", 360, 360)
gf%GM = 3.986004415d14;  gf%referenceRadius = 6378136.3d0        ! metres!
call gf%read_gravity_coefficients_4col("data/AtlasData/EGM2008.txt")
call gf%acceleration(position, degree, order, acceleration)      ! METRES in, m/s² out
```

**Unit trap:** `acceleration` divides by `referenceRadius`, which is in **metres** — unlike most of Atlas, which speaks km.

### Other types

| Type | Module | Purpose |
|------|--------|---------|
| `type_CRTBP` | `module_CRTBP` | CRTBP system + Lagrange points L1–L5 |
| `type_halo`, `type_LissajousOrbit` | same-named modules | 3rd-order Halo (Richardson 1980), Lissajous design; `init`, `compute_state` |
| `type_TimeSeries` | `module_timeSeries` | `n`, `interval` (s), `t0`, `tEnd`, `duration` (days), `time(:)`; `init(t0, step, n)` |
| `type_residual` / `type_acc` | `module_residual` | O−C residuals, outlier detection, weighted RMS |
| `type_sp3` | `module_sp3` | SP3 precise-ephemeris parsing |
| *(no derived type)* | `module_lambert` | Lambert solver, module-level routines `(r₁, r₂, Δt)` |
| `type_thruster`, `type_clock` | `module_thruster`, `module_clock` | see below |

### Thruster: lazy, epoch-driven ledger ([module_thruster.f90](src/module_thruster.f90))

Two ledgers, accrued **only when the caller flushes**: `current_*` (burn in progress; zeroed by `deactivate`/`set_fault`) and `total_*` (lifetime; only `init()` resets). Every state-changing method takes a UTC `type_t`.

```fortran
call thr%init(..., initial_fuel_mass=100d0)  ! kg; omit or 0 = fuel budget NOT modelled
call thr%activate(t);  call thr%deactivate(t);  call thr%toggle(t);  call thr%update_state(t)
q = thr%get_total_impulse(t);  mf = thr%get_total_fuel_used(t)   ! intent(inout): flush to t, then return
call thr%get_thrust_ECI(rs, vs)   ! argument-free, intent(in) — on the RKF78 hot path, must stay so
```

- A burn is accounted only when flushed; the flush applies the **current** `current_thrust` to the whole interval (settle it first, or pass `t` to `set_current_thrust`). `module_satellite` flushes once per propagated epoch.
- Fuel is drawn down and clamped at 0; a burn that would empty the tank stops at the flame-out epoch. `initial_fuel_mass = 0` = unlimited (the previous behaviour). Thrust acceleration is `F / (mass - total_fuel_used)`, both taken as a snapshot at the start of the propagation.
- `coordinate_frame` 2 (SBCS) is unimplemented (warns, zero thrust). Epochs must be **UTC** — `timeMinus` ignores `ttype`, and a BDT/GPST epoch would inject tens of thousands of spurious days.

### Satellite clock: `type_clock` ([module_clock.f90](src/module_clock.f90))

Held as `sat%clock` / `sta%clock`; all quantities in **seconds**, bias = 卫星钟面时间 − 参考时间.

```fortran
call clk%init(t0, bias0=1d-4, drift=1d-11, drift2=1d-18, noise=1d-10)   ! resets model + table, keeps name
!   bias(t) = bias0 + drift*dt + drift2*dt^2      (drift2 is the dt² coefficient, NOT halved)
clk%calcMethod = CLOCK_METHOD_INTERPOLATION   ! 2; default CLOCK_METHOD_FORMULA = 1
call clk%generateTable(epochs);  call clk%setTable(n, times, bias)   ! both replace the table wholesale
bias = clk%get_clock_state(t);   rate = clk%get_clock_drift(t)       ! the only query entry points
call clk%saveData("clock.txt", name="LEO-600km");  call clk%readData("clock.txt")
```

- Formula mode uses the analytic model (drift from its derivative); interpolation mode uses 4-point Lagrange on `table_clock`. Interpolation without a table warns and falls back to formula. Getter diagnostics go to **stderr** (`write(0,*)`), never `print*` — see the recursive-I/O trap.
- `saveData`/`readData` mirror the `type_access` convention: a `#` header (generator, object, data type, calc method, the four model parameters, `t0`, layout) then `MJD_day  MJD_sec  Bias/s` lines. Reading restores everything the header declares (file wins) and is fault tolerant; a headerless foreign file still yields its table (and switches `calcMethod` to interpolation, or the table would go unused). `name` is free-form metadata — pass `name=` (or set `clk%name`).
- **A query outside `table_clock` `stop`s** (`clock_locate`, no tolerance, same as `type_satellite%getStates`); formula mode has no such limit and evaluates any epoch. The trap to watch: a RANGING observation is tagged with the *receiver's clock-face time* `t + clk_R`, so a table-backed clock generated over the arc alone cannot be queried at its own observation epochs.

---

## Force model system

### Switches: named `integer` 0/1 attributes on `sat%opMod` (no switch array)

`key_moonGravity`, `key_moonTide`, `key_sunGravity`, `key_sunTide`, `key_srp`, `key_drag`, `key_postNewton` — **all default 0**, so a satellite with no switch assignment gets no luni-solar, SRP, drag or post-Newtonian term. The Earth gravity field is **not** one of these switches: `opMod%degree` and `opMod%order` control it directly (defaults 2 and 2, i.e. the full 2×2 field with J2; `order` = 0 gives zonal terms only; `degree` ≤ 1 disables the whole non-spherical term). `sat%init` takes no switches (the old `K(8)` argument is gone); `opMod%resetSwitches()` zeroes exactly these seven (not `degree`/`order`, `key_ERP`, `orderMoon`, the planet flags or `usePV`). `NumOP`'s `K_(7)` is internal: `module_satellite` assembles it via the private `opMod%forceSwitchArray()` (order: moonGravity, moonTide, sunGravity, sunTide, srp, drag, postNewton).

Common sets: Earth gravity only = `degree`/`order` alone (all switches 0); +luni-solar = `+key_moonGravity, key_sunGravity`; LEO full = `+key_srp, key_drag`; GNSS high precision = all except `key_postNewton`.

### `type_opMod` defaults

`degree`=2 (Earth gravity degree n), `order`=2 (gravity order m; the field is summed over n = 2..degree, m = 0..min(n, order)), `orderMoon`=0 (0 = point mass), `Cd`=2.2, `amrDrag`=0.02 m²/kg, `F107`=130 sfu, `Ap`=15, `reflectivity`=0, `amrSRP`=0.02 m²/kg, `SRP_Model`=1 (1 = cannonball, 2 = ECOM — `ECOM(1:9)` is read as `[D0, Y0, B0, Bc, Bs, Dc, Ds, Yc, Ys]` and gives `F = Σ_{X=D,Y,B} (X0 + Xc·cos u + Xs·sin u)·e_X`, in normalized accel units; **not** the literature ECOM1/ECOM2 parameter set — there are no X0/Z0 terms and the periodic terms are at 1×, not 2×/4×, orbital frequency), `key_ERP`=0 with `reflectivityERP`=0, `amrERP`=0.02, `usePV`=0 (0 = integrate NSG1, 1 = Cartesian), `step`=0 (fixed step in time units; 0 = adaptive), `showProgressBar`=.false., `Mercury`…`Pluto`=0.

- **ERP**: `module_EarthRadiationPressure` (Knocke/Ries/Tapley 1988), gated by `key_ERP` and independent of the eight switches; grids are computed analytically (the `GL01.DAT`/`GLM1.DAT` files are legacy), resolution 1°–10° by altitude.
- **Shadow**: two unconnected implementations — `shadowModel` in `module_perturbationSRP` (1 = cylindrical, 2 = conical, what SRP uses) and the standalone `module_shadow` (`inShadow(...)` → 0 lit / 1 umbra / 2 penumbra), which nothing in the repo exercises since `test_shadow` was deleted. Body radii differ per module: `module_shadow%SUN_RADIUS` = 696340 km, `module_planets%radius_sun` = 696000, `radius_earth` = 6371, `module_const%LENGTH_UNIT` = 6378.1363.
- **Drag**: `module_atmosphere`, `atmosphereModel` 1 = tabular exponential (default, grids in `data/AtlasData/nrlmsise00_cheng/`) or 2 = NRLMSISE-00 via the legacy `.for` (interface incomplete). `-½·Cd·(A/m)·ρ·|v_rel|·v_rel`, atmosphere co-rotating with the Earth.

### Integrator

**RKF78** — 13-stage embedded 7/8 method, error tolerance 1.5×10⁻¹³ (normalized), adaptive step (double below tolerance/100, halve and redo above); `step ≠ 0` selects fixed step. The step loop **aborts on non-finite values** (checks `ieee_is_finite` on `deviation` and the state, writes to stderr, `stop`s; second guard after 60 halvings) — keep that diagnostic on unit 0. `usePV` selects NSG1 or Cartesian integration; `NumOP` also accepts `eleType` 1/2/4, and Keplerian/Cartesian perigee below 100 km aborts propagation with a warning. All contributions are summed in SI and converted once at the end: `F_total = F_gravity + F_drag + F_SRP + F_ERP + F_thirdbody + F_PN + F_thrust`.

---

## Coordinate systems

`COORD_ICRS`=1 (J2000, inertial), `COORD_ITRS`=2 (ECF), `COORD_HORIZON`=3 (topocentric).

- **ECI↔ECF**: built-in IAU1980 (`HG_get`: IAU1976 precession + 20 nutation terms + Earth rotation + polar motion, needs `syncEOP(t)` for polar motion and DUT1) or IAU2006A via SOFA (`HG_get_sofa_06a`, `iau_C2T06A`). Velocity: `v_ECF = HG·v_ECI + ḢG·r_ECI`, ḢG by 0.1 s finite difference.
- **Horizon**: x north, y zenith, z east (velocity not implemented). **NTH** (for thrust): `NTH2ECI(r,v,M)` / `ECI2NTH(r,v,M)`.
- **Moon-centred / synodic**: `LCF2LCI`, `LCI2LCF`, `ECI2EarthMoonsynodic`, `EarthMoonSynodic2ECI`.
- Rotation helpers: `Rxmatrix`, `Rymatrix`, `Rzmatrix`, `coordinate_reverse`.

---

## Application programs

23 programs in `app/`; only two are argument-driven CLIs (documented in `release/readme_*.md`) — the rest are demos/batch drivers with settings hard-coded in the source.

| Program | Purpose |
|---------|---------|
| `app_orbitPrediction` | CLI: full force-model propagation (times, Keplerian elements, step, output, Cd/amrDrag/F107/Ap, order, then the 8 switches) |
| `app_observation` | Demo: one LEO + Beijing station, generates all six observation types into `data/output/demo_observation/` |
| `app_initialOrbitDetermination` | CLI: IOD from angle/range data |
| `satelliteScenario` | Multi-sat/station scenario → `obs_LEO_Beijing_Ah.txt`, `obs_LEO_Tokyo_RR.txt`, `obs_GPS_Beijing_PV.txt` |
| `observation_ranging` → `batchOD_ranging` | Paired ranging simulation (obsType 13) and batch POD from its output |
| `batch_GNSS`, `batch_position_beidou`, `batch_AzimuthElevation`, `ekf_angleData` | Batch POD from GNSS PV / SP3 (BeiDou C06) / Az-El, EKF from RA-Dec |
| `satellite`, `prediction_thruster` | `type_satellite` demos: trajectories, save/load, thruster ledger |
| `navSat_simulation` | navSat scenario: 3 MEO satellites (Walker 3/3/1, Cartesian states) + Beijing station, 6 pairwise ranging links, non-zero satellite clocks → everything under `data/navSat/` (see its `README.md`) |
| `od_RA_DEC`, `od_AZIMUTH_ELEVATION`, `od_RANGING`, `od_RANGERATE`, `od_POSITION`, `od_POSITION_VELOCITY` | Six single-observable batch POD demos sharing `app/od_demo_common.f90` (app-level module, excluded from auto-discovery and linked explicitly): one MEO + Beijing, 2-day arc, Cartesian state + amrSRP (label 3), outputs in `data/output/od_demo/` (`_obs`, `_eph`, `_residuals`, `_log`) |
| `od_MIXED` | Mixed-observation version: RA_DEC + RANGING + RANGERATE in one batch solve, same MEO/station/parameters; per-access residual files |
| `halo_orbit`, `cartesian`, `time` | Halo solution, coordinate transforms, time systems |

The classic scenario configuration files in `data/input/` (`object_*.txt`) are **dead** — their only reader (`app_orbitDetermination.f90`) was deleted; `setting_orbitDetermination.txt`, the reference for the `key = value` format, has since been removed too.

---

## Runtime data files

**All runtime data is in `data/AtlasData/` (tracked). The manual's `data/data_op/` path does not exist** — any `data/data_op/…` string is stale.

| File | Used by | Notes |
|------|---------|-------|
| `EGM2008.txt` | `module_orbitPredictionNumerical` | **The active Earth gravity field** (GGM03C is commented out) |
| `grail.txt` | same, Moon field; `test_lunar_gravity` | |
| `JPLEPH_DE440` | `PLEPH.for` | One file for every platform: `STATE` derives the record length (8,144 B) from the file size, so `FSIZER3` no longer picks a name by OS. Spans **1949-12-14 → 2150-01-22** (MJD 33264–106352); outside it, or if the file is missing/corrupt, `planetsPV` falls back to the analytical Sun/Moon with a once-per-body warning instead of aborting (planets still stop — no analytical model). See the PLEPH entry under Fortran/runtime traps |
| `Nutation.dat` | `module_space` | 20-term IAU1980 coefficients, unit 302 |
| `EOP.txt`, `leap.txt` | `module_EOPParam` | Polar motion, DUT1, leap seconds (`EOP20170101_20230426.txt` is referenced but not shipped) |
| `nrlmsise00_cheng/` | `module_atmosphere` | Density grids (`Ap=5/10/15`, F10.7 60–240) |
| `GL01.DAT`, `GLM1.DAT`, `input.txt` | nothing current | Legacy |

`data/input/` holds the inputs actually read (`isr1062*.sp3`, `eph_fit8h_refer.txt`, `obsData_alphaDelta1.txt`) plus the dead `object_*.txt` config set; `data/output/` is program output; `data/navSat/` is the navSat scenario simulation output (satellite ephemerides, ranging observations, satellite clock tables — see its `README.md`); `manual/` and `data/doc/` are documentation.

### I/O formats

- **Ephemeris**: `MJD x y z vx vy vz` (km, km/s) or the 8-column `MJD_day MJD_sec …`; `viewer/viewer.m` auto-detects both.
- **Observation / clock data**: `#`-prefixed self-describing header, then one line per epoch (see the `type_access` / `type_clock` sections). The seconds column is written `F13.6` (microsecond resolution) — a narrower field would round away the receiver-clock epoch shift of a RANGING observation (see below). Reading is list-directed, so field width never matters on input.
- **OD configuration**: the dead `data/input/setting_orbitDetermination.txt` defined a `key = value` format; the file has been deleted (the `object_*.txt` scenario configs that remain are equally dead).

---

## Viewer (MATLAB)

`viewer/` holds exactly two files: `viewer.m` and its `worldMap1.jpg` texture (missing texture degrades to a plain-coloured globe with a warning). `viewer.m` is a self-contained `classdef viewer < handle` — scenes `'Earth'`/`'Moon'`/`'Synodic'`, `touchSat(fn, name)`, `killSat(name)`, `zoomTo('L2', h)`, `zoomOut()`; the texture and all drawing helpers are private methods, so drive it from the MATLAB prompt (`v = viewer('Earth'); v.touchSat('../data/output/satellite_eph.txt')`). The old standalone helpers (`drawEarth.m`, `drawMoon.m`, `getLunarEph.m`, `getMeanSG.m`, `rotation.m`, `drawOrbit.m`, `plotOrbit.m`, `figure_3DEarth.m`) and the `main_satelliteOrbitVisualization.m` movie driver were deleted — all superseded by `viewer.m`; recover any of them from git history if needed.

---

## Development guidelines

**Adding a module:** create `src/<name>.f90` (module name matching the file where possible), add it to the right `LEVEL` in the [Makefile](Makefile), to `Atlas.vfproj` if the VS build needs it, and to [module_Atlas.f90](src/module_Atlas.f90) if it is public API.

**Conventions:** Fortran 2003+ derived types with type-bound procedures; normalized units internally, SI at I/O boundaries; `private` module variables with explicit `public`; errors via `module_ErrorHandler` (`RaiseError`/`RaiseWarning`), progress via `module_progressBar`. The code mixes `double precision`/`real(8)`, `logical`/`LOGICAL`, `endtype`/`end type`; recent modules use Chinese comments — match the file you touch.

**Terminal output is ASCII-only** (Chinese comments are fine, Chinese *strings* are not): the console here is GBK, so a UTF-8 Chinese string turns into mojibake, and the byte count mismatch also wrecks `\r`-based progress bars. Quantities printed for a human carry their unit — residual RMS in m / m/s / arcsec (see `resUnitCode`/`resRMSstr` in `module_ODor`), clock biases in s.

**Performance:** gravity-field `order` dominates (360×360 ≈ 65 000 Legendre evaluations per step; `EGM2008.txt` is cached after the first read, so raising `order` re-reads it). Typical adaptive steps: LEO 10–60 s, GPS altitude 300–600 s. Third-body terms are cheap; tidal terms and ERP (up to 64 800 grid cells) are the expensive optional ones. Use NSG2 for POD.

---

## Known gaps, traps and stale artifacts

**Build / tooling**
- `make` is green: all 39 targets (23 `app/` + 16 `test/`) build.
- `src/module_STM.f90` is not in the Makefile and never compiles; nothing references it.
- `test/test_lunar_gravity.f90` still reads `data/data_op/grail.txt` (directory gone).
- `Atlas.vfproj` omits `module_shadow.f90` and registers one test; `x64/Release/` holds a stale `module_perturbationAspherical.obj`.
- `bin/` is not tracked by git (in `.gitignore`).
- Only two targets actually call the SOFA path: `app/batch_position_beidou.f90` and `test/test_c2t06a.f90`.

**Fortran / runtime traps**
- **Never `print*` from a function that callers may use inside an I/O statement.** `print*, clk%get_clock_state(t)` where the getter also prints is *recursive I/O on the same unit*: gfortran's runtime deadlocks and the process hangs silently. Getter diagnostics must go to stderr (`write(0,*)`). `module_clock` follows this; `type_satellite%getStates` and `ErrorHandler`'s `RaiseWarning` still use `print*` and carry the same hazard.
- **A wrong direct-access record length returns zeros without any error** — that is how the ephemeris bug hid (`STATE` in [PLEPH.for](src/PLEPH.for) read 32,576-byte records with an 8,144-byte `RECL`, got zeros, and the third-body term produced NaN that spun the RKF78 loop forever). `STATE` now reads with `ACCESS='STREAM'` and derives the record length from the file size; the integrator aborts on non-finite values. **Any new binary-file reader should stream, or at least verify its record size.**
- `RECL` units are compiler-dependent: gfortran measures **bytes**, Intel is believed to use 4-byte words (inference, never verified here). The STREAM rewrite makes this moot; `NRECL` is unused by `STATE` (`FSIZER3` still returns it).
- **`PLEPH` no longer aborts: it reports status.** `STATE`/`PLEPH`/`CONST` carry `COMMON/PLEPHSTS/IPHINI,IPHSTS` (0 = ok, 1 = file/init failure — sticky, 2 = epoch outside `SS(1)..SS(2)`, 3 = record read error) and return **zeros** instead of stopping; every I/O in `STATE` has `IOSTAT=` so a missing file cannot abort inside the runtime either. Two entry points are the only supported interface for F90 callers (which cannot declare a COMMON): `PLEPH_STATUS(ISTAT)` and `PLEPH_FILENAME(NAM)` (delegates to `FSIZER3`, which owns the file name). The `IPHINI` guards are load-bearing: `FIRST` is cleared *before* the open in all three routines, so without them a failed init would leave every later call reading uninitialized `SS`/`IPT`/`PVST`. The shipped `JPLEPH_DE440` covers only **1949-12-14 → 2150-01-22** (MJD 33264–106352), so status 2 is reachable with a perfectly good file. Test: [test_jpl_fallback.f90](test/test_jpl_fallback.f90) — run it from the repo root (JPL) and from a directory without `data/` (fallback).

**Manual / code drift**
- `module_SP3Kit` (manual) is `module_sp3`; `module_perturbationAspherical` does not exist (Earth's field lives in `module_Gravity` + `module_orbitPredictionNumerical`).
- Observation constants are `OBS_TYPE_*` (11–16: RA_DEC, AZIMUTH_ELEVATION, RANGING, RANGERATE, POSITION, POSITION_VELOCITY); `type_access%dataTable` is normalized.
- `type_satelliteScenario` has no `generateObservations`; `u_*` factors convert *into* normalized units (appendix A.3 is inverted); NSG2 is implemented, not "reserved".
- The manual (`manual/atlas_manual.md`) still carries `data/data_op/` paths in ~18 places; the `files/` paths that `README.md` and `setting_orbitDetermination.txt` used to have are gone (both were fixed; the config file has since been deleted).

**Incomplete features**
- UT1 conversions print "unfinished part"; leap-second data lives in two places (see `type_t`).
- Horizon-frame velocity is not implemented; `maxEle` is not enforced.
- NRLMSISE-00 (`atmosphereModel = 2`) interface incomplete; `module_shadow` not wired into SRP/eclipse; `'SGP4'` is a placeholder.
- **Only RANGING carrys a clock term** — `calObsData1`/`calObsData2` apply `(clk_R − clk_T)·c` and the receiver-clock epoch tag for `OBS_TYPE_RANGING` only, so a clock group (label 11–13) attached to an object in no ranging link — or a link of any other observable type — gives identically zero Jacobian columns. `validateVarGroups` rejects that up front, because such a system is *rank-deficient but not zero-row*: `LUDCMP` substitutes `TINY` for the zero pivot and returns an inverse ~1e15 too large instead of reporting a singularity.
- **Clock groups are batch-only** — `invoke`/`numericalSTM` (EKF/sequential) build their own design matrix and never call `calB`; `validateVarGroups` rejects clock groups for methods 2 and 3.

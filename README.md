# Atlas Orbitalis (Atlas)

A comprehensive Fortran toolkit for satellite orbit computation and analysis, developed by Shawn Hu, Nanjing University.


## Overview

Atlas Orbitalis (Atlas) is a collection of tools for satellite orbit computation and analysis. Atlas provides essential functionalities for satellite mission analysis, orbit prediction, and tracking data processing. The name "Atlas Orbitalis" combines Greek mythology (Atlas, the Titan who carries the celestial sphere) with Latin (Orbitalis, meaning "orbital"), symbolizing the tool's purpose in celestial navigation and orbital mechanics.

Internally everything is computed in normalized units derived from Earth parameters (length 6378.1363 km, time 806.8109913067326 s, mass 5.9742272407074334×10²⁴ kg); file I/O converts at the boundary, so users work in km, km/s and degrees.

Atlas的实现依赖以下参考文献：
1. 《航天动力学引论》，刘林，胡松杰，王歆
2. 《卫星轨道理论与应用》，刘林，汤靖师
3. 《轨道力学基础》，刘林，侯锡云
4. 《Satellite Orbits Models, Methods and Applications》 Oliver Montenbruck, Eberhard Gill
5. 《Statistic Orbit Determination》, Byron Tapley, Bob Schutz, George H. Born
6. 《参考坐标系及航天应用》，黄珹，刘林
7. 《航天器定轨理论及应用》，刘林，胡松杰，曹建峰，汤靖师
8. 《导航卫星精密定轨技术》，王小亚，胡小工，蒋虎

## Key Features

- **Orbit Prediction**: Numerical (RKF78, adaptive step, full force model) and analytical (two-body, NSG1/NSG2 element propagation) orbit propagation, plus a fast approximate propagator for repeated whole-arc runs
- **Observation Simulation**: Generate synthetic tracking data — angles (RA/Dec, azimuth/elevation), ranging, range-rate/Doppler, and GNSS position/velocity
- **Initial Orbit Determination (IOD)**: General Laplace method from angle or range data
- **Precise Orbit Determination (POD)**: Batch least squares, extended Kalman filter and sequential filter, including estimation of drag/SRP area-to-mass ratios, ECOM solar radiation pressure coefficients and satellite clock offsets
- **Kalman Filtering**: Real-time orbit estimation and prediction
- **Multiple Perturbation Models**: Earth's aspherical gravity (EGM2008, to degree/order 100; a 360×360 GGM03C field ships alongside), third-body attraction of the Sun, Moon and planets with optional lunar and solar tide terms, atmospheric drag (tabular density grids), solar radiation pressure (cannonball, ECOM or Earth radiation pressure), post-Newtonian relativistic correction, and thrust
- **Multi-Body Dynamics**: Circular restricted three-body problem (CRTBP) with Lagrange points, 3rd-order halo orbits, Lissajous orbits and a Lambert solver
- **Coordinate and Time Systems**: IAU 1976/1980 and IAU 2006 (SOFA) Earth-fixed transformations, polar motion and DUT1 from EOP files, UTC/TT/UT1/GPST/BDT time systems, Moon-centred and Earth-Moon synodic frames
- **Analysis Utilities**: SP3 precise ephemeris parsing, JPL DE440 planetary/lunar ephemeris reader, thruster fuel-and-impulse ledger, GNSS-style clock models

## Repository Layout

```
Atlas/
├── Makefile                   # primary build (GNU Make + gfortran)
├── Atlas.sln / Atlas.vfproj   # secondary Visual Studio + Intel Fortran build
├── app/                       # application programs (17)
├── src/                       # library modules, layered L0–L6 below module_Atlas
├── test/                      # standalone test programs (21)
├── sofa/                      # vendored IAU SOFA sources, compiled into the library
├── viewer/                    # MATLAB orbit viewer (viewer.m + world map texture)
├── data/
│   ├── AtlasData/             # runtime data: gravity fields, ephemerides, EOP, leap seconds
│   ├── input/                 # sample inputs
│   ├── output/                # program output
│   ├── navSat/                # generated navSat scenario data
│   ├── manual/                # user manual (Chinese)
│   └── doc/                   # algorithm notes
├── release/                   # prebuilt command-line executables + their references
├── build/, bin/               # build outputs (created by make)
└── lib/, x64/                 # legacy binaries, unused
```

## Requirements

- **GNU Make** and a **gfortran ≥ 9** toolchain — the primary, supported build
- Optionally **Intel `ifx` ≥ 2021** with Visual Studio: `Atlas.sln` / `Atlas.vfproj` is a secondary build that lags the Makefile
- Optionally **MATLAB** for the orbit viewer under `viewer/`

The IAU SOFA library is vendored in `sofa/` and compiled from source as part of the build — no external dependency is needed.

## Building

```bash
make                     # library + every app and test (auto-discovered from app/, test/)
make build/libatlas.a    # library only
make bin/test_elements   # a single target
make clean               # remove build/ and bin/
```

The first build takes ~50 s because SOFA (249 F77 files) is compiled too; afterwards it is cached.

> **After changing a derived type that other modules `use`** (`type_access`, `type_satellite`, `type_ODor`, …), run `make clean` before `make`. Library objects carry no module dependency edges, so a stale `.mod` file otherwise causes either a confusing `Mismatch in components of derived type` error or objects that silently disagree about the type layout.

## Running

**Run every executable from the repository root.** All runtime data paths are relative (`data/AtlasData/…`), so running from `bin/` silently fails to find the ephemeris, gravity or EOP files.

```bash
./bin/test_elements.exe                       # a self-contained test
./bin/app_orbitPrediction.exe <args...>       # a command-line program
```

On Windows/MinGW the binaries carry the `.exe` suffix, as above; on Linux they do not.

## Application Programs

16 programs live in `app/`. Three are argument-driven command-line tools with dedicated references; the rest are demos and batch drivers whose settings are edited in the source.

| Program | Purpose |
|---------|---------|
| `app_orbitPrediction` | Full force-model orbit propagation — [reference](release/readme_orbitPrediction.md) |
| `app_observation` | Ground-station observation simulation — [reference](release/readme_observation.md) |
| `app_initialOrbitDetermination` | IOD from angle or range data — [reference](release/readme_initialOrbitDetermination.md) |
| `batch_GNSS` | Batch POD from GNSS position/velocity data |
| `batch_position_beidou` | Batch POD from BeiDou SP3 precise ephemeris |
| `batch_AzimuthElevation` | Batch POD from azimuth/elevation data |
| `ekf_angleData` | Extended Kalman filter from RA/Dec data |
| `observation_ranging` | Paired ranging observation simulation |
| `navSat_OD` | Batch POD of the navSat scenario (satellite orbits + clocks) |
| `satelliteScenario` | Multi-satellite, multi-station access generation and observation simulation |
| `navSat_simulation` | navSat scenario simulation: 3 MEO satellites + Beijing station, 6 ranging links → see [data/navSat/README.md](data/navSat/README.md) |
| `satellite` | `type_satellite` demonstration: trajectories, save/load |
| `prediction_thruster` | Low-thrust propagation with the thruster ledger |
| `halo_orbit` | Halo orbit solution in the Earth-Moon system |
| `cartesian`, `time` | Coordinate transformation and time system demonstrations |

## Tests

`test/` holds 16 standalone programs — no test framework, each prints results for human inspection (some also write into `data/output/`). Build them all with `make`, then run from the repository root:

```bash
./bin/test_fastPrediction.exe
```

| Group | Tests |
|-------|-------|
| Dynamics and multi-body | `test_CRTBP`, `test_halo`, `test_LissajousOrbit2`, `test_elements` |
| Propagators | `test_fastPrediction`, `BDS_orbitPrediction`, `orbit_prediction_ECOM` |
| Force models | `test_thruster`, `test_lunar_gravity` |
| Coordinates and time | `test_c2t06a`, `test_moon_coordinate_transformation`, `test_moon_coordinate_transformation2`, `test_OS` |
| Ephemerides and data I/O | `test_sp3`, `test_jpl_fallback` |
| State-transition matrices | `test_DXDR_nsg2` |

`test_jpl_fallback` is the one test that also runs *without* the ephemeris data — start it from a directory other than the repository root to watch the analytical Sun/Moon fallback in action.

## Runtime Data

All runtime data is tracked under `data/AtlasData/`: the `EGM2008.txt` Earth gravity field, the `grail.txt` lunar field, the JPL DE440 ephemeris (`JPLEPH_DE440`), `EOP.txt` and `leap.txt` for Earth orientation and leap seconds, and the `nrlmsise00_cheng/` atmospheric density grids.

Sample inputs live in `data/input/`, program output is written to `data/output/`, and `data/navSat/` holds a complete simulated scenario (ephemerides, ranging observations, clock tables) with [its own description](data/navSat/README.md).

## Orbit Viewer (MATLAB)

`viewer/viewer.m` post-processes ephemeris text files — either `MJD x y z vx vy vz` or the 8-column `MJD_day MJD_sec …`, in km and km/s. With the repository root as the working directory:

```matlab
addpath('viewer');

v = viewer('Earth');                                  % 3-D Earth, J2000 frame
v.touchSat('data/output/ISS_trajectory.txt');         % add an orbit
v.touchSat('data/output/BDS_MEO_ephemeris.txt', 'BDS');
v.killSat('BDS');  v.zoomOut();

s = viewer('Synodic');                                % Earth-Moon rotating frame
s.touchSat('data/output/Lissajous_ephemeris.txt', 'L2');
s.zoomTo('L2');  s.zoomOut();                         % or zoomTo([x y z], halfRange)
```

The three scenes are `'Earth'`, `'Moon'` (Moon-centred inertial frame) and `'Synodic'` (Earth-Moon rotating frame, km input, drawn in normalized units).

## Documentation

- **User manual**: `manual/atlas_manual.md` (Chinese, PDF alongside). It lags the code in places — when in doubt, trust the source.
- **Command-line references**: `release/readme_*.md` for the three CLI programs
- **Scenario data notes**: `data/navSat/README.md`

## Credits

Atlas builds on third-party code and data — see [THIRD_PARTY.md](THIRD_PARTY.md) for the
full notices:

- **IAU SOFA** — the Earth-rotation and time-scale routines in `sofa/` (SOFA Software License, reproduced in every file there)
- **JPL DE440** — the planetary and lunar ephemeris in `data/AtlasData/`, read by `src/PLEPH.for` (Caltech/JPL/NASA)
- **NRLMSISE-00** — the atmospheric model (Picone et al., 2002)
- **Numerical Recipes** LU routines (F90 translation by J-P Moreau) and the **Montenbruck & Gill** orbital-elements routine, both marked in place
- **Gravity fields and Earth-orientation data** — EGM2008 (NGA), GGM03C (UT CSR), GRAIL (NASA/JPL), IERS EOP and leap seconds
- **NASA Blue Marble** — the Earth texture used by the MATLAB viewer

## Author

**Shawn Hu**, Nanjing University — lxiahu@hotmail.com

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
Third-party components and data keep their own terms; see [THIRD_PARTY.md](THIRD_PARTY.md).

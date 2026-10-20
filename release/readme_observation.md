# Observation Demo Program

A program based on Atlas, developed by Shawn Hu (<lxiahu@hotmail.com>) in Nanjing University.

## Overview

`app_observation` is a self-contained observation simulation demo. It builds
one LEO satellite and one ground station, then generates **all six observation
types** supported by `module_access` over a one-day arc:

| Constant | Value | Observable | Components |
|---|---:|---|---:|
| `OBS_TYPE_RA_DEC` | 11 | Right ascension, declination | 2 |
| `OBS_TYPE_AZIMUTH_ELEVATION` | 12 | Azimuth, elevation | 2 |
| `OBS_TYPE_RANGING` | 13 | Range | 1 |
| `OBS_TYPE_RANGERATE` | 14 | Range-rate | 1 |
| `OBS_TYPE_POSITION` | 15 | x, y, z | 3 |
| `OBS_TYPE_POSITION_VELOCITY` | 16 | Position and velocity | 6 |

All scenario settings are hard-coded in `app/app_observation.f90` (time span,
satellite orbit, station coordinates, minimum elevation), following the
batch-program convention. Modify the source and recompile to change them.

## Build and Run

From the repository root:

```bash
make
./bin/app_observation.exe
```

Run from the repository root: all runtime paths are relative, so the
ephemeris and gravity-field files are not found from `bin/`.

## Scenario

- Epoch: 2025-01-01 00:00:00 UTC, one-day arc
- Step: 30 s
- Satellite: ~600 km altitude, 60 deg inclination, near-circular LEO
- Force model: Earth gravity 8x8, drag, Sun and Moon third-body
- Station: Beijing (116.4 deg E, 39.9 deg N, 50 m)
- Minimum elevation: 10 deg
- Constraint: optical types (`RA_DEC`, `AZIMUTH_ELEVATION`) skip epochs in the
  Earth shadow; radio types do not

`POSITION` and `POSITION_VELOCITY` are state observations and are therefore
not filtered by elevation or shadow (the existing `calObsData1` convention),
so they contain one record per simulation epoch.

## Output

Everything is written to `data/output/demo_observation/`:

| File | Contents |
|---|---|
| `obs_RA_DEC.txt` | MJD day, seconds, RA/deg, Dec/deg |
| `obs_AZIMUTH_ELEVATION.txt` | MJD day, seconds, Az/deg, El/deg |
| `obs_RANGING.txt` | MJD day, seconds, range/km |
| `obs_RANGERATE.txt` | MJD day, seconds, range-rate/(km/s) |
| `obs_POSITION.txt` | MJD day, seconds, x,y,z/km |
| `obs_POSITION_VELOCITY.txt` | MJD day, seconds, x,y,z/km, vx,vy,vz/(km/s) |
| `eph_LEO-600km.txt` | Satellite ephemeris used for the simulation |

Each observation file starts with a self-describing `#` header (generator,
access, transmitter, receiver, data type and its integer code, epoch count,
time span, column layout). Data lines are whitespace-separated and use SI
units at the file boundary (km, km/s, degrees); internally `type_access`
stores normalized units.

`readAccessData` reads these files back, adopting the header's observation
type and device names. The header type name is checked against the integer,
so files written with the pre-redesign type codes are reported and ignored
instead of being silently reinterpreted.

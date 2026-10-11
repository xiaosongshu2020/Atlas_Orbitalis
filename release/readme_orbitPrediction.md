# Orbit Prediction Tool

Satellite Orbit Prediction Program Based on Atlas Orbit Algorithm Library

## Program Overview

`orbitPrediction.exe` is a Fortran console application for satellite orbit prediction. This program is implemented based on the Atlas orbit algorithm library developed by Dr. Shawn Hu from Nanjing University, capable of accurately calculating and predicting satellite orbital evolution under complex force models.

## Key Features

- **High-precision orbit prediction**: Supports numerical integration with multiple perturbation force models
- **Multiple output formats**: Can output Keplerian orbital elements or position/velocity vectors
- **Comprehensive force models**: Includes Earth gravity field, lunisolar perturbations, atmospheric drag, solar radiation pressure, etc.
- **Flexible parameter configuration**: Supports fine control of various physical parameters and model switches

## System Requirements

- Windows/Linux operating system
- Necessary Fortran runtime libraries (if required)

## Usage

### Command Line Syntax

```bash
orbit_prediction.exe <start_time> <kepler_elements> <end_time> <step> <filename> <ephType> <Cd> <amrDrag> <F107> <Ap> <reflectivity> <amrSRP> <degree> <order> <K_lunar> <K_lunarTide> <K_solor> <K_solarTide> <K_SRP> <K_drag> <K_PN>
```

### Parameter Description

| Parameter | Description | Format/Units | Typical Value |
|-----------|-------------|--------------|---------------|
| `start_time` | Start time | "YYYY MM DD HH MM SS.SSS" | "2023 01 01 12 00 00.000" |
| `kepler_elements` | Keplerian elements | 6 space-separated values (km, deg) | "7000 0.01 45 30 60 90" |
| `end_time` | End time | "YYYY MM DD HH MM SS.SSS" | "2023 01 02 12 00 00.000" |
| `step` | Time step | seconds | 60 |
| `filename` | Output ephemeris filename | string | "output.txt" |
| `ephType` | Ephemeris type | 1: Keplerian, 2: Position/Velocity | 2 |
| `Cd` | Drag coefficient | dimensionless | 2.2 |
| `amrDrag` | Area-to-mass ratio for drag | m²/kg | 0.02 |
| `F107` | Solar flux index | SFU | 130 |
| `Ap` | Geomagnetic index | dimensionless | 15 |
| `reflectivity` | Surface reflectivity | dimensionless | 0 |
| `amrSRP` | Area-to-mass ratio for SRP | m²/kg | 0.02 |
| `degree` | Gravity field degree n (阶数) | integer | 2 |
| `order` | Gravity field order m (次数); 0 = zonal only | integer | = `degree` |
| `K_lunar` | Lunar perturbation flag | 0/1 | 0 |
| `K_lunarTide` | Lunar tide flag | 0/1 | 0 |
| `K_solor` | Solar perturbation flag | 0/1 | 0 |
| `K_solarTide` | Solar tide flag | 0/1 | 0 |
| `K_SRP` | Solar radiation pressure flag | 0/1 | 0 |
| `K_drag` | Atmospheric drag flag | 0/1 | 0 |
| `K_PN` | Post-Newtonian flag | 0/1 | 0 |

### Keplerian Elements Format

Six parameters in order:
1. Semi-major axis (km)
2. Eccentricity
3. Inclination (deg)
4. Right Ascension of Ascending Node (deg)
5. Argument of Perigee (deg)
6. Mean Anomaly (deg)

### Example

```bash
orbit_prediction.exe "2023 01 01 12 00 00.000" "7000 0.01 45 30 60 90" "2023 01 02 12 00 00.000" 60 output.txt 2 2.2 0.02 130 15 0 0.02 2 0 0 0 0 0 0 0 0
```

This example will calculate:
- From 2023-01-01 12:00:00 to 2023-01-02 12:00:00
- Time step of 60 seconds
- Initial orbit: semi-major axis 7000km, eccentricity 0.01, inclination 45°, etc.
- Output position/velocity format ephemeris to output.txt file

## Output File

The generated ephemeris file contains satellite states at specified time epochs, with format determined by `ephType` parameter:
- `ephType=1`: Keplerian orbital elements
- `ephType=2`: Position and velocity vectors (recommended)

## Algorithm Description

This program is based on the Atlas orbit algorithm library, which employs advanced numerical integration methods and precise force models, including:
- Earth non-spherical gravity field (with tesseral harmonics)
- Third-body perturbations (lunisolar gravity)
- Atmospheric drag models
- Solar radiation pressure effects
- Tidal effects
- Relativistic effects (post-Newtonian corrections)

## Developer Information

**Algorithm Library Author**: Shawn Hu  
**Affiliation**: Nanjing University  
**Contact**: hushnju@163.com

## Notes

1. Ensure all 21 parameters are provided in correct order
2. Time format must strictly follow "YYYY MM DD HH MM SS.SSS"
3. Orbital element units: distance in km, angles in degrees
4. For most applications, position/velocity output format (ephType=2) is recommended
5. Perturbation force switches (K_* parameters) are typically set to 0, enable when needed

## Error Handling

If parameter count is incorrect or format is invalid, the program will display detailed usage instructions and terminate execution. Please carefully check the error messages for guidance.
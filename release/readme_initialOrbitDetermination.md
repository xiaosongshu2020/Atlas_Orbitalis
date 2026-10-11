# Initial Orbit Determination Tool

Satellite Initial Orbit Determination Program Based on Atlas Orbit Algorithm Library

## Program Overview

`initialOrbitDetermination.exe` is a Fortran console application for satellite initial orbit determination (IOD). This program implements classical orbit determination methods using ground-based observation data to estimate a satellite's orbital elements from tracking measurements.

## Key Features

- **Multiple Observation Types**: Supports azimuth/elevation, right ascension/declination, and range/range-rate measurements
- **Classical IOD Methods**: Implements General Laplace method for orbit determination
- **Ground Station Support**: Configurable station coordinates (longitude, latitude, altitude)
- **Flexible Data Input**: Reads observation data from text files
- **Comprehensive Output**: Displays determined orbital elements and saves results to file

## System Requirements

- Windows/Linux operating system
- Necessary Fortran runtime libraries (if required)

## Usage

### Command Line Syntax

```bash
initialOrbitDetermination.exe <longitude> <latitude> <altitude> <obs_type> <data_file>
```

### Parameter Description

| Parameter | Description | Format/Units | Example |
|-----------|-------------|--------------|---------|
| `longitude` | Station longitude | Degrees (-180 to 180) | 116.385 |
| `latitude` | Station latitude | Degrees (-90 to 90) | 39.905 |
| `altitude` | Station altitude | Meters | 50.0 |
| `obs_type` | Observation type | String (AH, AD, RRD) | "AH" |
| `data_file` | Observation data file path | String | "observations.dat" |

### Observation Types

- **AH**: Azimuth and Elevation measurements
- **AD**: Right Ascension and Declination measurements  
- **RRD**: Range and Range-Rate measurements (Note: Limited implementation)

### Example

```bash
initialOrbitDetermination.exe 116.385 39.905 50.0 AH observations.dat
```

This example will:
- Use a station at longitude 116.385°, latitude 39.905°, altitude 50m
- Process azimuth/elevation observation data from "observations.dat"
- Perform initial orbit determination using the General Laplace method
- Output results to "iod_results.txt"

## Observation Data Format

The observation data file should contain time-stamped measurement data. The specific format depends on the observation type:

### AH Format (Azimuth/Elevation)
```
YYYY MM DD HH MM SS.SSS AZIMUTH ELEVATION
```
- Time: Year, Month, Day, Hour, Minute, Second
- Azimuth: Degrees
- Elevation: Degrees

### AD Format (Right Ascension/Declination)
```
YYYY MM DD HH MM SS.SSS RA DEC
```
- Time: Year, Month, Day, Hour, Minute, Second  
- RA: Right Ascension (degrees)
- DEC: Declination (degrees)

### RRD Format (Range/Range-Rate)
```
YYYY MM DD HH MM SS.SSS RANGE RANGE_RATE
```
- Time: Year, Month, Day, Hour, Minute, Second
- Range: Kilometers
- Range-Rate: km/s

## Output

The program generates two types of output:

1. **Console Display**: Shows the determined orbital elements including:
   - Semi-major axis
   - Eccentricity  
   - Inclination
   - Right ascension of ascending node
   - Argument of perigee
   - Mean anomaly

2. **Results File**: Saves detailed results to "iod_results.txt" including:
   - Input parameters
   - Determined orbital elements
   - Residual information
   - Covariance data

## Algorithm Description

This program uses the General Laplace method for initial orbit determination, which:

- Processes angular measurements (azimuth/elevation or RA/dec) from a single station
- Solves for the orbit that best fits the observed angular data
- Provides an initial estimate of the six Keplerian orbital elements
- Can be used as a starting point for more precise orbit determination methods

## Developer Information

**Algorithm Library Author**: Shawn Hu  
**Affiliation**: Nanjing University  
**Contact**: hushnju@163.com

## Notes

1. Ensure the observation data file exists and is properly formatted
2. Station coordinates should be in the WGS84 reference frame
3. For RRD observation type, functionality may be limited in current implementation
4. The minimum elevation angle is set to 5° by default
5. Results provide an initial orbit estimate; refinement may be needed for precise applications

## Error Handling

The program will display error messages and usage instructions if:
- Incorrect number of parameters provided
- Invalid coordinate values
- Unknown observation type specified  
- Observation file not found
- No valid observation data in file
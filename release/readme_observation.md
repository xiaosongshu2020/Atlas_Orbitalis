# Satellite Observation Program

A program based on Atlas, developed by Shawn Hu (<lxiahu@hotmail.com>) in Nanjing University.

## Overview

The Satellite Observation Program is a Fortran-based tool designed for simulating satellite observation scenarios. It calculates visibility and observation data based on a ground station's position and various observation parameters.

## Features

- Time-based simulation with configurable start time, end time, and step size
- Customizable observer location (longitude, latitude, altitude)
- Support for different observation types
- Constraints for minimum elevation angle
- Solar and lunar avoidance angles for observation planning
- Ephemeris file reading for satellite position data
- Data output to user-specified files


## Usage

Run the program with the following command:

```
./satellite_observation ts step te longitude latitude altitude obsType minEle solarDistance lunarDistance fnEph [fnData]
```

### Parameters

| Parameter | Description | Units |
|-----------|-------------|-------|
| ts | Initial time | Modified Julian Date (MJD) |
| step | Time step | Seconds |
| te | End time | Modified Julian Date (MJD) |
| longitude | Observer longitude | Degrees |
| latitude | Observer latitude | Degrees |
| altitude | Observer altitude | Meters |
| obsType | Observation type | Integer code |
| minEle | Minimum elevation angle | Degrees |
| solarDistance | Minimum angular distance from Sun | Degrees |
| lunarDistance | Minimum angular distance from Moon | Degrees |
| fnEph | Satellite ephemeris filename | String |
| fnData | Observation data output filename (optional) | String |

### Example

```
./satellite_observation 59000.0 60.0 59001.0 139.0 35.0 100.0 1 10.0 30.0 15.0 satellite.eph results.dat
```

This command runs a 1-day simulation (MJD 59000.0 to 59001.0) with a 60-second step size for an observer in Tokyo (139°E, 35°N) at 100m altitude, using observation type 1, minimum elevation of 10°, minimum solar distance of 30°, and minimum lunar distance of 15°.

## Observation Types

The program supports different observation types through the `obsType` parameter:

| obsType | Description |
|---------|-------------|
| 11 | isAlphaDelta - Angles |
| 12 | isRadar - A, h, range and velocity |
| 13 | isRanging - range |
| 14 | isGNSS - position and velocity |
| 15 | isPosition - position |
| 16 | isAh - A, h |
| 17 | isRRD - range and velocity |
| 19 | isOther - Other observation types |

## File Formats

### Ephemeris File Format

The ephemeris file should contain the satellite's position and velocity data in Cartesian coordinates, with each line formatted as:

```
MJD(day) MJD(sec in day) X Y Z VX VY VZ
```

Where:
- MJD(day) is the integer part of Modified Julian Date
- MJD(sec in day) is the seconds part of the day
- X, Y, Z are position coordinates in kilometers
- VX, VY, VZ are velocity components in kilometers per second

## Error Handling

The program will display error messages and exit in the following cases:
- Insufficient command-line arguments
- Invalid parameter formats
- Non-existent ephemeris file

## Notes

- All angles are in degrees
- Times are in Modified Julian Date (MJD)
- Distances are in kilometers or meters as specified



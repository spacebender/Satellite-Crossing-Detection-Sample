# Satellite Crossing and Visibility Simulation

This project simulates geometric crossing and visibility events between:

- A target satellite propagated using SGP4 (TLE)
- A tracker satellite defined via Keplerian elements

## Features
- TEME frame consistency
- 30° FOV crossing detection
- 1000 km proximity constraint
- Cylindrical Earth shadow eclipse model
- Above-horizon Earth occultation test
- 24-hour propagation with 5-second step

## Requirements
- numpy
- astropy
- sgp4
- hapsira
- matplotlib

Run with:
python Satellite_Visibility_Main.py


## Sample Output
Simulating 24 hours...

    ------------------------------
    Detected 1 crossing moments.
    First Crossing: 2025-09-01 00:29:35+00:00
    Last Crossing: 2025-09-01 00:29:35+00:00
    ------------------------------
    Detected 1 visible moments.
    Visible at 2025-09-01 00:29:35+00:00 | Distance: 116.43 km

    ========================================
    CROSSING INTERVALS
    ========================================
    Interval: 00:29:35 to 00:29:35

    ========================================
    VISIBLE INTERVALS (<1000km & SUNLIT & Above Horizon)
    ========================================
    Detected: 00:29:35 to 00:29:35

### Installation
pip install -r Requirements.txt
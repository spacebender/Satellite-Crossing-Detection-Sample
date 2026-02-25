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

pip install -r Requirements.txt
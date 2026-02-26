# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 21:05:31 2026

@author: Gilari Ramachandran Karthi

Use Python Version 3.11.4 or less


    
Sample Output:
    
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
 
**************************************************************************************************
Note:
        If only "poliastro" is available,  then replace "hapsira" by "poliastro". In which case
        use Python Version 3.10 or less
        
        
        Example : Replace "from hapsira.twobody import Orbit" to "from poliastro.twobody import Orbit"
        
        Runtime is slow possibly because of the GCRS to TEME conversion for every "step_sec = 5" seconds. 
****************************************************************************************************   
    
    
"""
import numpy as np
from sgp4.api import WGS84, Satrec, jday
from sgp4.conveniences import sat_epoch_datetime
from datetime import datetime, timezone
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun, TEME, GCRS, CartesianRepresentation, CartesianDifferential
from hapsira.twobody import Orbit 
from hapsira.bodies import Earth
from datetime import timedelta


# Required Functions 
def group_events(timestamps, step_sec):
    if not timestamps: return []
    events = []
    start_time = timestamps[0]
    for i in range(1, len(timestamps)):
        if (timestamps[i] - timestamps[i-1]).total_seconds() > step_sec:
            events.append((start_time, timestamps[i-1]))
            start_time = timestamps[i]
    events.append((start_time, timestamps[-1]))
    return events

def is_sunlit(rt, rs, R):
    rt = np.asarray(rt, dtype=float)  # Earth->target
    rs = np.asarray(rs, dtype=float)  # Earth->sun

    dot = np.dot(rt, rs)
    rs2 = np.dot(rs, rs)
    rt2 = np.dot(rt, rt)

    fe1 = (dot > 0)
    fe2 = ( (dot*dot)/rs2 - rt2 + R*R < 0 )

    return fe1 and fe2

def is_above_horizon(r_obs_km, r_tgt_km, R_earth):

    r_o = np.asarray(r_obs_km, dtype=float)
    r_t = np.asarray(r_tgt_km, dtype=float)

    delta_r = r_t - r_o

    ro_norm = np.linalg.norm(r_o)
    delta_norm = np.linalg.norm(delta_r)

   
    r_o_e = np.sqrt(ro_norm**2 - R_earth**2)

    
    f_o = (delta_norm * r_o_e + np.dot(delta_r, r_o)) > 0

    return f_o

def state_to_teme(state_trk, obstime):
    """
    Convert propagated state (inertial ~GCRS) to TEME at the given obstime.
    Returns r_teme_km, v_teme_km_s as numpy arrays.
    """
    r = state_trk.r.to(u.km)
    v = state_trk.v.to(u.km / u.s)

    rep = CartesianRepresentation(r).with_differentials(CartesianDifferential(v))
    sc_gcrs = GCRS(rep, obstime=obstime)

    sc_teme = sc_gcrs.transform_to(TEME(obstime=obstime))

    r_teme = sc_teme.cartesian.xyz.to(u.km).value
    v_teme = sc_teme.velocity.d_xyz.to(u.km / u.s).value
    return r_teme, v_teme

# ------------------------------- MAIN Section -------------------------------------------------------
R_EARTH_KM = 6371.0
# ---- ------------------------- Object Details from TLE --------------------------------------------- 
tle1 = "1 63223U 25052P   25244.59601767  .00010814  00000-0  51235-3 0  9991"
tle2 = "2 63223  97.4217 137.0451 0006365  74.2830 285.9107 15.19475170 25990"
obj = Satrec.twoline2rv(tle1,tle2,WGS84)


epoch_time = sat_epoch_datetime(obj)
year, month, day, hour, minute, seconds = [epoch_time.year,epoch_time.month,
                                            epoch_time.day,epoch_time.hour,epoch_time.minute,epoch_time.second]
jd_obj_at_epoch, fr_obj_at_epoch   = jday(year, month, day, hour, minute, seconds)

e2, r_obj_at_epoch, v_obj_at_epoch = obj.sgp4(jd_obj_at_epoch,fr_obj_at_epoch)

# start epoch 2025-09-01 00:00:00 UTC (based on tracker)
epoch_start = datetime(2025, 9, 1, 00, 00, 0, tzinfo=timezone.utc)
jd_obj_tgt_epoch, fr_obj_tgt_epoch = jday(epoch_start.year, epoch_start.month, epoch_start.day, 
                                          epoch_start.hour, epoch_start.minute, epoch_start.second)


# Position and Velocity of the target at the tracker start epoch
e2, robj_tgt_epoch, vobj_tgt_epoch = obj.sgp4(jd_obj_tgt_epoch, fr_obj_tgt_epoch)

# ----------------------------Tracker Details from the Orbital Elements ----------------------------- 
SMA_track_tgt_epoch, Incl_track_tgt_epoch,Eccen_track_tgt_epoch,RAAN_track_tgt_epoch, Argp__track_tgt_epoch,MeanAn__track_tgt_epoch =[6878,97.4,0,72.628,331.7425,0]

ss_tracker = Orbit.from_classical(
    Earth,
    SMA_track_tgt_epoch * u.km,
    Eccen_track_tgt_epoch * u.one,
    Incl_track_tgt_epoch * u.deg,
    RAAN_track_tgt_epoch * u.deg,
    Argp__track_tgt_epoch * u.deg,
    MeanAn__track_tgt_epoch * u.deg,
    epoch=Time(epoch_start)
)

# --- ---------------------------------- Propagation  Loop ------------------------------------------
duration_hrs = 24
step_sec = 5 # increase if faster result is needed 
steps = int((duration_hrs * 3600) / step_sec)

crossing_events = []
visible_events = []
r_trk_list = []
r_tgt_list = []
rel_dist_list = []

print(f"\nSimulating {duration_hrs} hours...")
epoch_start_dt = epoch_start  
epoch_start_ast = Time(epoch_start_dt) #for astropy mode
  
for i in range(steps):
    print(f"Step {i+1}/{steps}  |  Progress: {100*(i+1)/steps:.2f}%")
    current_time = epoch_start_dt + timedelta(seconds=i * step_sec)
    astropy_time = Time(current_time)
    
    # Propagate Tracker (hapsira)
    
    dt = (astropy_time - epoch_start_ast).to(u.s)
    state_trk = ss_tracker.propagate(dt)

    r_trk, v_trk = state_to_teme(state_trk, astropy_time) # tracker states
    
    # Propagate Object (SGP4)
    
    jd, fr = jday(current_time.year, current_time.month, current_time.day, 
                  current_time.hour, current_time.minute, current_time.second)
    error, r_tgt, v_tgt = obj.sgp4(jd, fr) # object states 
    
   
    #-------Relative Position--------- ---
    rel_pos = np.array(r_tgt) - np.array(r_trk)
    
    dist = np.linalg.norm(rel_pos)
    
    
    
    # Angle between Tracker Velocity and Relative Position
    unit_v_trk = v_trk / np.linalg.norm(v_trk)
    unit_rel_pos = rel_pos / dist
    angle_rad = np.arccos(np.clip(np.dot(unit_v_trk, unit_rel_pos), -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    
    # Crossing  (30 deg FOV ==> 15 deg half-angle)
    if angle_deg <= 15.0:
        crossing_events.append(current_time)
        
        # ---  Visibility/Detection  ---
        """for the visibility/detection all A, B, C muct hold """
        if dist < 1000: # A) Proximity
            sun_pos = get_sun(astropy_time).transform_to(TEME(obstime=astropy_time)).cartesian.xyz.to(u.km).value
            sunlit = is_sunlit(r_tgt, sun_pos, R_EARTH_KM) # B) Sun-lit
            above_horizon = is_above_horizon(r_trk, r_tgt, R_EARTH_KM) # C) Above Horizon
            if sunlit and above_horizon:
                visible_events.append((current_time, dist))
                
    r_trk_list.append(np.array(r_trk, dtype=float)) # Accumulate the tracker position
    r_tgt_list.append(np.array(r_tgt, dtype=float)) # Accumulate the object position
    rel_dist_list.append(np.array(r_tgt, dtype=float)) # Accumulate the relative distance
                
# ---------------- Group Intervals --------------------------------------------
crossing_intervals = group_events(crossing_events, step_sec)
visible_timestamps = [item[0] for item in visible_events]
visible_intervals = group_events(visible_timestamps, step_sec)



# --- -------------------------- Output Results ------------------------------
print("-" * 30)
if crossing_events: # Display total crossing
    print(f"Detected {len(crossing_events)} crossing moments.")
    print(f"First Crossing: {crossing_events[0]}")
    print(f"Last Crossing: {crossing_events[-1]}")
else:
    print("No crossings found.")

print("-" * 30)
if visible_events:# Display total detectionn
    print(f"Detected {len(visible_events)} visible moments.")
    for time, d in visible_events[:5]: # Print first 5 for brevity
        print(f"Visible at {time} | Distance: {d:.2f} km")
else:
    print("No visible events detected.")

# ----------------------- Display Intervals -----------------
print("\n" + "="*40)
print("CROSSING INTERVALS")
print("="*40)
for start, end in crossing_intervals: 
    print(f"Interval: {start.strftime('%H:%M:%S')} to {end.strftime('%H:%M:%S')}")

print("\n" + "="*40)
print("VISIBLE INTERVALS (<1000km & SUNLIT & Above Horizon)")
print("="*40)
if not visible_intervals:
    print("No events met detection criteria.")
else:
    for start, end in visible_intervals:
        print(f"Detected: {start.strftime('%H:%M:%S')} to {end.strftime('%H:%M:%S')}")









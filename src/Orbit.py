import numpy as np
import pandas as pd
# https://doi.org/10.1016/j.solener.2003.12.003
# Based on VSOP87; disable the jit to avoid the bug in this package
import sunposition
sunposition.disable_jit()

# detectors positions
# SK: https://www-sk.icrr.u-tokyo.ac.jp/~berns_s/SUPERK/GPS/
detectors = pd.DataFrame({ 
    'name': ['SK', 'Jinping', 'JUNO'],
    'latitude': [36 + 25/60 + 32.282/3600, 0, 0],
    'longitude': [137 + 18/60 + 37.224/3600, 0, 0],
    'altitude': [374, 0, 0],
    })

class Orbit():
    def __init__(self, det=detectors.set_index('name').loc['SK']):
        self.latitude, self.longitude, self.altitude = det['latitude'], det['longitude'], det['altitude']

    def rotationAngle(self):
        # The longitude=0 relative to the 24:00 posistion
        pass

    def revolutionAngle(self):
        pass

    def getZenith(self, timestamp):
        # timestamp may be 'now', which returns the current time using time.time(),
        # or an ISO-8601 formatted date & time, e.g. '2024-04-08T11:09:34-07:00'
        # return the Zenith angle of sun in the coodinate of detector (longitude, latitude)
        # delta_longitude = longitude - self.rotationAngle(self)
        # t, latitude, longitude, elevation
        # The sunpos in sunposition package consider the refraction of atmosphere, which is useless
        ## start time: Y=-4712,M=1,D=1,12:00. Equ 4 in 10.1016/j.solener.2003.12.003
        ## The magic number explanation: https://www.zhihu.com/question/21698547
        jd = sunposition.julian_day(timestamp)
        ## 
        latitude, longitude = sunposition._norm_lat_lon(self.latitude, self.longitude)
        elevation, Delta_t = 0, 0
        alpha_prime, delta_prime, H_prime = sunposition._sun_topo_ra_decl_hour(latitude, longitude, elevation, jd, Delta_t)
        # _sun_topo_azimuth_zenith
        phi = np.deg2rad(latitude)
        delta_prime_rad, H_prime_rad = np.deg2rad(delta_prime), np.deg2rad(H_prime)

        # eq 41, topocentric elevation angle, uncorrected
        e0 = np.rad2deg(np.arcsin(np.sin(phi)*np.sin(delta_prime_rad) + np.cos(phi)*np.cos(delta_prime_rad)*np.cos(H_prime_rad)))
        Gamma = np.rad2deg(np.arctan2(np.sin(H_prime_rad), np.cos(H_prime_rad)*np.sin(phi) - np.tan(delta_prime_rad)*np.cos(phi))) % 360
        Phi = (Gamma + 180) % 360 #azimuth from north

        return Phi, 90 - e0

    def getZenithNow(self):
        now = sunposition.time_to_datetime64('now')
        az, zen = self.getZenith(now)
        return zen

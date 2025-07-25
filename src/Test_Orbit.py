import Orbit
import time

orbit = Orbit.Orbit()
print(orbit.getZenithNow())
t = time.time()
for i in range(48):
    print(orbit.getZenith(t + i * 3600))

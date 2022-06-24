import folium
import numpy as np
import time
from astropy import units as u
from astropy.coordinates import CartesianRepresentation, EarthLocation, SphericalRepresentation
from .sampleearth import sample_to_gdf, sample_to_earthloc

def greateCircle(lat1, lon1, lat2, lon2):

    return 6371 * u.km * (
        np.arccos(np.sin(lat1.si.value) * np.sin(lat2.si.value) + np.cos(lat1.si.value) * np.cos(lat2.si.value) * np.cos(lon1.si.value - lon2.si.value)))

nList = (1, 10, 25)
#nList = (1, 10, 25, 50, 100, 250)

# sample distance
if True:
    for n in nList:
        t0 = time.perf_counter()
        points = sample_to_earthloc(n)
        t1 = time.perf_counter()

        p0 = CartesianRepresentation(points[0].geocentric)
        p0 = p0 / p0.norm()
        p1 = None
        angle = 2*np.pi
        for p in points[1:]:
            v = CartesianRepresentation(p.geocentric)
            a = np.arccos(p0.dot(v / v.norm()))
            if p1 is None or a < angle:
                p1 = p
                angle = a

        sampleDistance = greateCircle(points[0].lat, points[0].lon, p1.lat, p1.lon)
        print(f"n:{n} sample distance:{sampleDistance:.3f} numPoints={len(points)} time={t1 - t0:.3f} seconds")

# print maps
if True:
    for n in nList:
        gdf = sample_to_gdf(n)
        
        map = folium.Map()
        folium.GeoJson(gdf.to_json()).add_to(map)
        map.save(f"sampled-n{n}.html")
        
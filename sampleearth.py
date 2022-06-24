from icosphere import icosphere
from astropy.coordinates import CartesianRepresentation, SphericalRepresentation
from astropy import units as u
from geopandas import GeoDataFrame
from shapely.geometry import Point
import numpy as np


def sample(nu:int=1):
    (verts, edges) = icosphere(nu)

    p0=None
    p1=None
    sp0=None
    sp1=None
    dist=0
    points=[]
    for v in verts:
        p = CartesianRepresentation(v)
        sp = SphericalRepresentation.from_cartesian(p)

        if p0 is None:
            p0 = p
            sp0 = sp
        else:
            angleBetween = np.arccos(p0.dot(p))
            if p1 is None or angleBetween < dist:
                p1 = p
                sp1 = sp
                dist = angleBetween


        points.append(Point(sp.lon.wrap_at('180d').degree, sp.lat.degree))

    sampleDistance = greateCircle(sp0, sp1)
    
    return (GeoDataFrame(geometry=points, crs="+proj=longlat +datum=WGS84 +no_defs"), sampleDistance)

def greateCircle(s1:SphericalRepresentation, s2:SphericalRepresentation):
    lon1 = s1.lon.si.value
    lat1 = s1.lat.si.value
    lon2 = s2.lon.si.value
    lat2 = s2.lat.si.value

    return 6371 * u.km * (
        np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)))

if __name__ == "__main__":
    import folium
    import time
    
    for n in [1, 10, 25, 50, 100, 250]:
        t0 = time.perf_counter()
        (gdf, sampleDistance) = sample(n)
        t1 = time.perf_counter()
        print(f"n:{n} sample distance:{sampleDistance:.3f} numPoints={len(gdf.geometry)} time={t1 - t0:.3f} seconds")
    
        map = folium.Map()
        folium.GeoJson(gdf.to_json()).add_to(map)
        map.save(f"sampled-n{n}.html")
'''
Use an icosphere to create a point mesh of approximately evenly spaced points on the
surface of the globe.

Author: Joe Reed
'''
from icosphere import icosphere
from astropy.coordinates import CartesianRepresentation, EarthLocation, SphericalRepresentation
from astropy import units as u
from geopandas import GeoDataFrame
from shapely.geometry import Point
import numpy as np

def sample_to_gdf(nu:int=1) -> GeoDataFrame:
    '''
    Sample the globe into a mesh of approximately evently spaced points.

    This mesh is created by creating a geodesic icosahedron with the specified subdivision frequency and 
    projecting those points to the surface of the ellipsoid.

    Parameters
    ----------
    nu (int): The icosophere subdivision frequency; must be greater than zero. Defaults to 1

    Returns
    -------
    gdf (geopandas.GeoDataFrame): A data frame holding the point mesh
    '''
    (verts, edges) = icosphere(nu)

    points=[]
    for v in verts:
        p = CartesianRepresentation(v)
        sp = SphericalRepresentation.from_cartesian(p)

        points.append(Point(sp.lon.wrap_at('180d').degree, sp.lat.degree))
    
    return GeoDataFrame(geometry=points, crs="+proj=longlat +datum=WGS84 +no_defs")

def sample_to_earthloc(nu:int=1) -> list[EarthLocation]:
    '''
    Sample the globe into a mesh of approximately evently spaced points.

    This mesh is created by creating a geodesic icosahedron with the specified subdivision frequency and 
    projecting those points to the surface of the ellipsoid.

    Parameters
    ----------
    nu (int): The icosophere subdivision frequency; must be greater than zero. Defaults to 1

    Returns
    -------
    points (list[EarthLocation]): A list of mesh points
    '''
    (verts, edges) = icosphere(nu)
    points=[]
    for v in verts:
        p = CartesianRepresentation(v)
        sp = SphericalRepresentation.from_cartesian(p)

        points.append(EarthLocation.from_geodetic(lon=sp.lon, lat=sp.lat))
    
    return points

def sample(nu:int=1, unit:u.Quantity=u.deg) -> list[tuple[float, float]]:
    '''
    Sample the globe into a mesh of approximately evently spaced points.

    This mesh is created by creating a geodesic icosahedron with the specified subdivision frequency and 
    projecting those points to the surface of the ellipsoid.

    Parameters
    ----------
    nu (int): The icosophere subdivision frequency; must be greater than zero. Defaults to 1
    unit (astropy.units.Quantity, Optional): The unit in which to represent the longitude and latitude. Defaults to degrees.

    Returns
    -------
    points (list[tuple[float,float]]): A list of mesh points, represented as tuples of (lon,lat)
    '''
    (verts, edges) = icosphere(nu)
    points=[]
    for v in verts:
        p = CartesianRepresentation(v)
        sp = SphericalRepresentation.from_cartesian(p)

        points.append((sp.lon.to_value(unit), sp.lat.to_value(unit)))
    
    return points

if __name__ == "__main__":
    import folium
    import time

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

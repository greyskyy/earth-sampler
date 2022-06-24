import pytest

from .. import sampleearth

class TestSampleEarth_n1:
    nu = 1
    expectedPoints=12

    def test_sample(self):
        points = sampleearth.sample(self.nu)

        assert len(points) == self.expectedPoints
    
    def test_sampleGdf(self):
        gdf = sampleearth.sample_to_gdf(self.nu)

        assert len(gdf.geometry) == self.expectedPoints

    def test_sampleLocs(self):
        points = sampleearth.sample_to_earthloc(self.nu)

        assert len(points) == self.expectedPoints

class TestSampleEarth_n10:
    nu = 10
    expectedPoints=1002

    def test_sample(self):
        points = sampleearth.sample(self.nu)

        assert len(points) == self.expectedPoints
    
    def test_sampleGdf(self):
        gdf = sampleearth.sample_to_gdf(self.nu)

        assert len(gdf.geometry) == self.expectedPoints

    def test_sampleLocs(self):
        points = sampleearth.sample_to_earthloc(self.nu)

        assert len(points) == self.expectedPoints
# Sample-earth

Sample the earth using an [icosphere](https://github.com/vedranaa/icosphere).

## tl;dr

```
from sampleearth import sample

(gdf, sampleDistance) = sample(n=25)
```

## output

```
$ conda env create -f environment.yaml -n sample-earth
$ conda activate sample-earth
$ python sampleearth.py
n:1 sample distance:7053.644 km numPoints=12 time=0.019 seconds
n:10 sample distance:601.390 km numPoints=1002 time=0.375 seconds
n:25 sample distance:232.986 km numPoints=6252 time=2.303 seconds
n:50 sample distance:115.229 km numPoints=25002 time=9.651 seconds
n:100 sample distance:57.299 km numPoints=100002 time=37.222 seconds
n:250 sample distance:22.844 km numPoints=625002 time=234.010 seconds
```

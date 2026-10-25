# Geodesics on an Ellipsoid of Revolution

C++20 implementations of direct and inverse geodesic solutions on an oblate ellipsoid of revolution, developed directly on the ellipsoid without introducing an auxiliary sphere.

This repository accompanies the manuscript **New formulations and solutions for geodesics on an oblate ellipsoid of revolution without an auxiliary sphere**, by Jia-Chun Guo.

## Method

The formulation takes the northern geodesic vertex as the origin and uses the nonnegative Clairaut constant to describe the geodesic. A change of variable expresses the distance and longitude relations through standard Legendre elliptic integrals with a real modulus.

- **Direct problem:** a non-iterative solution using distance-to-amplitude series reversion, explicit branch selection, and periodic distance reduction. Long geodesics spanning multiple complete latitude oscillations are supported.
- **Inverse problem:** determination of the Clairaut constant through five geometric cases: meridional, equatorial-endpoint, non-equatorial symmetric, general, and near-antipodal configurations.
- **Sensitive inverse configurations:** safeguarded Newton iteration with analytical derivatives, bracket preservation, and midpoint fallback.
- **Numerical evaluation:** Boost.Math elliptic integrals, finite series evaluated with Horner's method and Clenshaw recurrence, and nominal 100-decimal-digit arithmetic using Boost.Multiprecision.

The direct and inverse programs share the same definitions and numerical kernels.

## Source files

| File | Purpose |
| --- | --- |
| `GeodesicDefinitions.hpp` | Ellipsoid parameters, numerical type, angular conventions, and geometric utilities |
| `GeodesicCore.hpp` | Shared elliptic-integral and series kernels for both problems |
| `NearAntipodalClairaut.hpp` | Safeguarded solution of branch-dependent longitude equations |
| `GeodesicDirect.cpp` | Single-point direct calculation |
| `GeodesicInverse.cpp` | Single-point shortest-path inverse calculation |
| `GeodesicDirect_Input_Karney_test_data.cpp` | Batch direct tests using Karney-format records |
| `GeodesicInverse_Input_Karney_test_data.cpp` | Batch inverse tests using Karney-format records |

## Build and run

Requirements:

- A C++20 compiler, such as GCC or Microsoft Visual C++.
- Boost.Math and Boost.Multiprecision headers. The project was developed using Boost 1.89.

Place the shared headers alongside the source files. Build each `.cpp` as a separate executable; do not link all programs into one executable.

Example GCC commands (adjust the Boost include path):

```bash
g++ -std=c++20 -O2 -I /path/to/boost GeodesicDirect.cpp -o GeodesicDirect
g++ -std=c++20 -O2 -I /path/to/boost GeodesicInverse.cpp -o GeodesicInverse
g++ -std=c++20 -O2 -I /path/to/boost GeodesicDirect_Input_Karney_test_data.cpp -o GeodesicDirectTest
g++ -std=c++20 -O2 -I /path/to/boost GeodesicInverse_Input_Karney_test_data.cpp -o GeodesicInverseTest
```

Omit `-I /path/to/boost` if Boost is already on the compiler's include path. In Visual Studio 2022, select C++20, add the Boost include directory, and use a Release configuration for timing tests.

Run a single-point program without arguments for interactive input, or supply decimal-degree coordinates:

```bash
# lat1 lon1 azi1 distance_in_metres
./GeodesicDirect 30 0 45 1000000

# lat1 lon1 lat2 lon2
./GeodesicInverse 30 0 -29.999999 179.5
```

Batch programs accept an input file and an optional CSV output path:

```bash
./GeodesicDirectTest GeodTest-short.dat direct_results.csv
./GeodesicInverseTest GeodTest-short.dat inverse_results.csv
```

Without arguments, they provide a file-selection dialog on Windows or a terminal prompt on other platforms. Decompress the reference `.dat.gz` file before running a test.

## Units and conventions

- Program input and output angles are in degrees; internal calculations use radians.
- Distances are in metres.
- The default ellipsoid is WGS84: `a = 6378137 m`, `f = 1/298.257223563`.
- The initial azimuth is forward. The inverse solver's terminal azimuth is backward; batch comparisons convert it to Karney's continuing forward convention.
- Final longitudes use the principal longitude range. For long direct geodesics, accumulated longitude contraction is retained before normalization.
- At exact poles, the inverse programs assign meridional azimuths of 0 or 180 degrees, independently of the supplied polar longitude.

## Reference data and tests

The standard batch tests read the original ten-column Karney format:

```text
lat1 lon1 azi1 lat2 lon2 azi2 s12 a12 m12 S12
```

Reference data:

- [GeodTest.dat.gz — 500,000 records](https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest.dat.gz)
- [GeodTest-short.dat.gz — 10,000-record subset](https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest-short.dat.gz)
- [Dataset archive](https://doi.org/10.5281/zenodo.32156)

For the tested 500,000 WGS84 records:

| Comparison with the tabulated reference values | Result |
| --- | --- |
| Direct terminal latitude, longitude, and forward azimuth | Every absolute difference below `1e-18` degrees |
| Inverse distance | Every absolute difference below `0.1 pm` (`1e-13 m`) |
| Both inverse azimuths simultaneously below `1e-12` degrees | 90.0280% of records |
| Both inverse azimuths simultaneously below `1e-18` degrees | 52.9528% of records |

These results measure agreement for the tested configurations. Nominal working precision and iteration tolerances do not directly specify final accuracy. Global convergence over arbitrary ellipsoid parameters and comparative computational efficiency have not been established.

## References

- Karney, C. F. F. (2013). *Algorithms for geodesics*. Journal of Geodesy, 87, 43–55. [Paper](https://arxiv.org/abs/1109.4448).
- Karney's reference dataset, linked above. See the accompanying manuscript for the derivations and further references.

## Author

Jia-Chun Guo, Anhui Jianzhu University, China.

Contact: [jchguo@ahjzu.edu.cn](mailto:jchguo@ahjzu.edu.cn).

## License

The source code is released under the MIT License. See `LICENSE` for details.

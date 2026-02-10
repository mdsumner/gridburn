
## gridburn 0.1.0

First working version. Sparse polygon rasterization with exact coverage
fractions, powered by vendored C++ code from 
[exactextract](https://github.com/isciences/exactextract) algorithm.

### Features

* `burn_sparse()` computes exact cell coverage fractions for
  polygon-grid intersections and returns a compact two-table
  representation: run-length encoded interior cells and individually
  weighted boundary cells.

* Accepts geometry input from sf (`sfc`), geos (`geos_geometry`), or
  raw WKB lists. Geometry is passed via WKB serialization, decoupling
  from any single geometry backend.

* Automatic tiling for large grids (`tile_size` argument, default 4096).
  Bounds memory usage to ~64 MB per tile regardless of grid dimensions.
  Tested to 256,000 × 128,000 grids (~32 billion cells).

* `materialise_chunk()` expands sparse results to dense matrix for
  visualisation. Supports filtering by polygon `id` and direct use
  with `ximage::ximage()` via the stored `extent` attribute.

* Output uses `c(xmin, xmax, ymin, ymax)` extent convention, matching
  `ximage` and `rasterImage()`.

### Internals

* Vendors 27 files (3159 lines) from exactextract C++ core. Zero
  dependency on GDAL, R stats framework, or exactextract's processing
  pipeline. Only the cell intersection algorithm and its geometry
  helpers.

* Links to GEOS via the
 [libgeos](https://github.com/paleolimbot/libgeos) R package (no
  system GEOS dependency, no configure script). Five `#include`
  patches (`geos_c.h` → `libgeos.h`), otherwise upstream code is
  unmodified.

* Uses [cpp11](https://cpp11.r-lib.org/) for the R–C++ bridge.

* All upstream GEOS version guards preserved — fast paths activate
  for GEOS ≥ 3.7/3.8/3.10, with fallbacks for older versions.

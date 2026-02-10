
# gridburn

Sparse polygon rasterization with exact coverage fractions.

`gridburn` computes exact fractional coverage for polygon-grid intersections
and returns results in a compact sparse format — not a dense matrix. The
algorithm is vendored from Daniel Baston's
[exactextract](https://github.com/isciences/exactextract) C++ library (the
same code behind [exactextractr](https://github.com/isciences/exactextractr)).

## Installation

```r
# install.packages("pak")
pak::pak("hypertidy/gridburn")
```

Requires [libgeos](https://github.com/paleolimbot/libgeos) and
[cpp11](https://cpp11.r-lib.org/) (installed automatically as dependencies).

## Usage

```r
library(gridburn)
library(geos)

# A polygon that partially covers a 10x10 grid
poly <- as_geos_geometry(
  "POLYGON ((1.5 1.5, 7.3 1.5, 7.3 8.2, 1.5 8.2, 1.5 1.5))"
)

result <- burn_sparse(poly, extent = c(0, 10, 0, 10), dimension = c(10, 10))
result
#> <gridburn> 10 x 10 grid, 1 geometry
#>   runs:  6 (30 interior cells)
#>   edges: 20 boundary cells
#>   sparsity: 50.0% empty

# Interior cells (weight = 1) as run-length encoded rows
result$runs

# Boundary cells with exact fractional coverage
result$edges

# Expand to dense matrix for visualisation
image(materialise_chunk(result), useRaster = TRUE)
```

## Output Format

`burn_sparse()` returns a list with two data.frames:

- **`runs`**: `(row, col_start, col_end, id)` — contiguous interior cells
  where coverage fraction ≈ 1.0, compressed by row
- **`edges`**: `(row, col, weight, id)` — boundary cells with exact
  partial coverage (0 < weight < 1)

This is far more compact than a dense matrix for typical polygon-on-grid
operations, especially when polygons are large relative to cell size.

## Design Notes

The sparse two-table format is designed for downstream use in zonal
statistics, polygon overlay analysis, and similar operations where you
need to know *which cells* a polygon covers and *by how much*, but don't
need to materialise the full grid.

`gridburn` uses [libgeos](https://github.com/paleolimbot/libgeos) for
GEOS access (no system GEOS dependency, no configure script). The
exactextract algorithm handles all geometry types (polygons,
multipolygons, lines) and correctly accounts for holes.

## License

Apache License 2.0. The vendored exactextract code is Copyright
ISciences, LLC, also Apache 2.0.

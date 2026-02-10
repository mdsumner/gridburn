#' Sparse polygon rasterization with exact coverage fractions
#'
#' Computes exact coverage fractions for polygon-grid intersections and returns
#' results in a sparse two-table format: run-length encoded interior cells and
#' individually weighted boundary cells.
#'
#' @param x geometry input, one of:
#'   - an `sfc` geometry column (from sf)
#'   - a `geos_geometry` vector (from geos)
#'   - a list of raw vectors containing WKB
#' @param extent numeric vector `c(xmin, xmax, ymin, ymax)` defining the raster extent
#' @param dimension integer vector `c(ncol, nrow)` defining the grid dimensions
#'
#' @return A list with class `"gridburn"` containing:
#'   \describe{
#'     \item{`runs`}{data.frame with columns `row`, `col_start`, `col_end`, `id` —
#'       run-length encoded interior cells (coverage fraction ≈ 1.0)}
#'     \item{`edges`}{data.frame with columns `row`, `col`, `weight`, `id` —
#'       boundary cells with partial coverage (0 < weight < 1)}
#'     \item{`extent`}{the raster extent as supplied}
#'     \item{`dimension`}{the grid dimensions as supplied}
#'   }
#'
#'   Row and column indices are 1-based. Row 1 is the top (ymax) row.
#'   The `id` column is a 1-based index into the input geometry vector.
#'
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   library(geos)
#'   poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
#'   result <- burn_sparse(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#'   result$runs
#'   result$edges
#' }
burn_sparse <- function(x, extent, dimension) {
  wkb <- as_wkb_list(x)
  extent <- as.double(extent)
  dimension <- as.integer(dimension)

  stopifnot(
    length(extent) == 4,
    length(dimension) == 2,
    dimension[1] > 0,
    dimension[2] > 0,
    extent[2] > extent[1],  # xmax > xmin
    extent[4] > extent[3]   # ymax > ymin
  )

  result <- cpp_burn_sparse(
    wkb,
    extent[1], extent[3], extent[2], extent[4],  # xmin, ymin, xmax, ymax to C++
    dimension[1], dimension[2]
  )

  result$extent <- extent
  result$dimension <- dimension
  class(result) <- "gridburn"
  result
}

#' Materialise a gridburn result to a dense matrix or vector
#'
#' Expands the sparse two-table representation into a dense coverage fraction
#' matrix. Primarily for visualisation and testing.
#'
#' @param x a `"gridburn"` object from [burn_sparse()]
#' @param id integer polygon id to materialise, or `NULL` (default) for all
#' @param type character, one of `"matrix"` (default) or `"vector"`
#'
#' @return A numeric matrix (nrow × ncol) or vector (length nrow*ncol) of
#'   coverage fractions. Values range from 0 (outside) to 1 (fully inside).
#'
#' @export
materialise_chunk <- function(x, id = NULL, type = c("matrix", "vector")) {
  stopifnot(inherits(x, "gridburn"))
  type <- match.arg(type)

  ncol <- x$dimension[1]
  nrow <- x$dimension[2]

  mat <- matrix(0, nrow = nrow, ncol = ncol)

  runs <- x$runs
  edges <- x$edges

  if (!is.null(id)) {
    runs <- runs[runs$id %in% id, , drop = FALSE]
    edges <- edges[edges$id %in% id, , drop = FALSE]
  }

  # Fill interior runs
  if (nrow(runs) > 0) {
    for (i in seq_len(nrow(runs))) {
      r <- runs$row[i]
      cs <- runs$col_start[i]
      ce <- runs$col_end[i]
      mat[r, cs:ce] <- 1
    }
  }

  # Fill edge cells
  if (nrow(edges) > 0) {
    for (i in seq_len(nrow(edges))) {
      r <- edges$row[i]
      c <- edges$col[i]
      w <- edges$weight[i]
      # For overlapping polygons, sum weights (or take max depending on use case)
      mat[r, c] <- mat[r, c] + w
    }
  }

  if (type == "vector") {
    as.vector(t(mat))  # row-major order
  } else {
    mat
  }
}

#' @export
print.gridburn <- function(x, ...) {
  ncol <- x$dimension[1]
  nrow <- x$dimension[2]
  n_runs <- nrow(x$runs)
  n_edges <- nrow(x$edges)
  n_ids <- length(unique(c(x$runs$id, x$edges$id)))

  # Compute total cells represented
  total_interior <- if (n_runs > 0) sum(x$runs$col_end - x$runs$col_start + 1L) else 0L
  total_cells <- total_interior + n_edges
  grid_cells <- as.numeric(ncol) * as.numeric(nrow)
  sparsity <- 1 - total_cells / grid_cells

  cat(sprintf("<gridburn> %d x %d grid, %d geometr%s\n",
              ncol, nrow, n_ids, if (n_ids == 1) "y" else "ies"))
  cat(sprintf("  runs:  %d (%d interior cells)\n", n_runs, total_interior))
  cat(sprintf("  edges: %d boundary cells\n", n_edges))
  cat(sprintf("  sparsity: %.1f%% empty\n", sparsity * 100))
  invisible(x)
}

# ---- internal: geometry coercion to WKB list ----

as_wkb_list <- function(x) {
  UseMethod("as_wkb_list")
}

#' @export
as_wkb_list.list <- function(x) {
  # Assume already a list of raw vectors (WKB)
  stopifnot(all(vapply(x, is.raw, logical(1))))
  x
}

#' @export
as_wkb_list.sfc <- function(x) {
  # sf geometry column → WKB via sf::st_as_binary
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required to convert sfc objects", call. = FALSE)
  }
  sf::st_as_binary(x)
}

#' @export
as_wkb_list.geos_geometry <- function(x) {
  # geos geometry → WKB via geos::geos_write_wkb
  if (!requireNamespace("geos", quietly = TRUE)) {
    stop("Package 'geos' is required to convert geos_geometry objects", call. = FALSE)
  }
  geos::geos_write_wkb(x)
}

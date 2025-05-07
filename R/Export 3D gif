#' Export Rotating GIF of GRF Vectors
#'
#' Generates and animates a rotating GIF from two datasets with optional annotations.
#'
#' @param df1 First dataset.
#' @param df2 Second dataset (optional).
#' @param max_radius Sphere size multiplier.
#' @param showLines1 Logical, draw butterfly for df1.
#' @param showLines2 Logical, draw butterfly for df2.
#' @param showPeak Logical, draw peak vector.
#' @param fps Frames per second for the GIF.
#' @param out_file Path to save the GIF (e.g. "www/GRF.gif").
#' @export
export_grf_gif <- function(df1, df2 = NULL,
                           max_radius = 0.015,
                           showLines1 = FALSE,
                           showLines2 = FALSE,
                           showPeak = FALSE,
                           fps = 10,
                           out_file = "www/3D_animation.gif") {
  # create temporary frames directory
  dir.create("frames", showWarnings = FALSE)

  # headless rgl plot setup
  png_pattern <- file.path("frames", "frame%03d.png")
  plot_dynamic_3d(df1, df2, max_radius, showLines1, showLines2, showPeak)
  view3d(theta = 0, phi = 15, zoom = 0.75)
  n_frames <- 100
  for (i in seq_len(n_frames) - 1) {
    view3d(theta = i * 360 / n_frames, phi = 15, zoom = 0.75)
    rgl.snapshot(sprintf(png_pattern, i))
  }
  rgl.close()

  # stitch images
  imgs <- image_read(list.files("frames", pattern = "png$", full.names = TRUE))
  gif  <- image_animate(imgs, fps = fps)
  image_write(gif, path = out_file)
  unlink("frames", recursive = TRUE)
}

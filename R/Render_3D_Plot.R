#' Render Interactive 3D Plot of GRF Vectors
#'
#' Creates an rgl widget showing GRF vectors as spheres with optional butterfly lines and peak vector.
#'
#' @param df1 First dataset (data.frame with Fx, Fy, Fz).
#' @param df2 Optional second dataset, same format as df1.
#' @param max_radius Multiplier for sphere radii (e.g. from slider).
#' @param showLines1 Logical, draw butterfly lines for dataset 1.
#' @param showLines2 Logical, draw butterfly lines for dataset 2.
#' @param showPeak Logical, draw peak resultant vector for both datasets.
#' @param angle_phi Numeric, elevation angle in degrees (default 15).
#' @return An rglwidget object.
#' @export
plot_dynamic_3d <- function(df1, df2 = NULL,
                            max_radius = 0.015,
                            showLines1 = FALSE,
                            showLines2 = FALSE,
                            showPeak = FALSE,
                            angle_phi = 15) {
  res1 <- compute_resultant(df1)
  res2 <- if (!is.null(df2)) compute_resultant(df2) else numeric(0)
  global_max <- max(c(res1, res2), na.rm = TRUE)
  peak_z <- max(c(abs(df1$Fz), if (!is.null(df2)) abs(df2$Fz)), na.rm = TRUE)
  axis_len <- peak_z * 1.1
  half_ax  <- axis_len / 2

  # sphere radii
  radius_A <- (max(res1, na.rm = TRUE) / global_max) * (axis_len * max_radius)
  radius_B <- (if (length(res2)) max(res2, na.rm = TRUE) else 0) / global_max * (axis_len * max_radius)

  open3d(useNULL = TRUE)
  plot3d(NA,
         xlim = c(-half_ax, half_ax),
         ylim = c(-half_ax, half_ax),
         zlim = c(0, axis_len),
         xlab = "Fx", ylab = "Fy", zlab = "Fz",
         type = "n")

  # dataset 1
  cols1 <- colorRampPalette(c("orange","red","purple"))(nrow(df1))
  spheres3d(df1$Fx, df1$Fy, df1$Fz, col = cols1, radius = radius_A)
  if (showLines1) {
    for (i in seq_len(nrow(df1))) {
      lines3d(c(0, df1$Fx[i]), c(0, df1$Fy[i]), c(0, df1$Fz[i]), col = cols1[i], lwd = 1)
    }
  }
  if (showPeak && length(res1)) {
    i1 <- which.max(res1)
    lines3d(c(0, df1$Fx[i1]), c(0, df1$Fy[i1]), c(0, df1$Fz[i1]), col = "black", lwd = 4)
  }

  # dataset 2
  if (!is.null(df2)) {
    cols2 <- colorRampPalette(c("cadetblue1","cyan3","blue"))(nrow(df2))
    spheres3d(df2$Fx, df2$Fy, df2$Fz, col = cols2, radius = radius_B)
    if (showLines2) {
      for (i in seq_len(nrow(df2))) {
        lines3d(c(0, df2$Fx[i]), c(0, df2$Fy[i]), c(0, df2$Fz[i]), col = cols2[i], lwd = 1)
      }
    }
    if (showPeak && length(res2)) {
      i2 <- which.max(res2)
      lines3d(c(0, df2$Fx[i2]), c(0, df2$Fy[i2]), c(0, df2$Fz[i2]), col = "black", lwd = 4)
    }
  }

  rglwidget()
}

#' Render 2D Time-Series Plot of GRF Components
#'
#' Produces a ggplot of Fx, Fy, Fz, and resultant over stance percentage with peak markers.
#'
#' @param data A data.frame containing columns Index, Fx, Fy, Fz, Resultant, and Dataset.
#' @param displayComponents Character vector of components to plot, e.g. c("Fx","Fy","Fz","Resultant").
#' @param showPeaks Character vector of components to mark peaks for.
#' @return A ggplot2 object.
#' @export
plot_2d_grf <- function(data, displayComponents = c("Fx","Fy","Fz","Resultant"), showPeaks = NULL) {
  require(ggplot2)
  pd <- subset(data, Component %in% displayComponents)
  pd$Component <- factor(pd$Component, levels = c("Fx","Fy","Fz","Resultant"))
  upper <- max(abs(data$Resultant), na.rm = TRUE) * 1.1
  lower <- -max(abs(data$Fy), na.rm = TRUE) * 1.1

  p <- ggplot(pd, aes(Index, Value, color = Component, linetype = Dataset)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c(Fx="blue", Fy="purple", Fz="red", Resultant="black"), drop = FALSE) +
    scale_linetype_manual(values = c("A"="solid","B"="dashed"), drop = FALSE) +
    geom_hline(yintercept = 0, size = 0.5) +
    labs(title = "Fx, Fy, Fz & Resultant", x = "Stance (%)", y = "GRF (N)") +
    ylim(lower, upper) +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_blank())

  if (length(showPeaks)) {
    peaks <- do.call(rbind, lapply(showPeaks, function(comp) {
      subset(data, Component == comp)[which.max(abs(subset(data, Component == comp)$Value)), ]
    }))
    p <- p + geom_segment(data = peaks,
                          aes(x = Index, xend = Index, y = 0, yend = Value,
                              color = Component, linetype = Dataset),
                          size = 0.8, linetype = "dotted", show.legend = FALSE)
  }

  p
}

###########
# Manhatan Plots
###########

png(
  filename = file.path(c.result.dir, "manhattan_plot.png"),
  width = 1024,
  height = 720,
  units = "px"
)

qqman::manhattan(
  results1.man,
  main = "Manhattan Plot",
  col = c("blue4", "orange3"),
  ymax = 10,
  genomewideline = y.bar.intercept,
  annotatePval = y.bar.intercept,
  suggestiveline = FALSE,
  annotateTop = TRUE
)
options(error=function() traceback(3))

dev.off()
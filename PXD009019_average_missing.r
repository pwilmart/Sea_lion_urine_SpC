
# load the ggplot libraries
library(tidyverse)

# read in the data
temp <- read_tsv("average_fractionMissing.txt")

# make a basic plot
ggplot(temp, aes(x = AverageSpC, y = FracMissing)) +
  geom_point() +
  geom_smooth(method = "loess", span = 0.05) + 
  ggtitle("Fraction Missing Data versus Average SpC") + 
  labs(x = "Average SpC", y = "Fraction Missing")

# expanded x-axis plot
ggplot(temp, aes(x = AverageSpC, y = FracMissing)) +
  coord_cartesian(xlim = c(0, 50)) +
  geom_line() + 
#  geom_smooth(method = "loess", span = 0.05) +
  ggtitle("Missing versus Average SpC") + 
  labs( x = "Average SpC", y = "Fraction Missing") + 
  geom_vline(xintercept = 2.5, linetype = "dotted") +
  geom_vline(xintercept = 5.0, linetype = "dashed")

sessionInfo()



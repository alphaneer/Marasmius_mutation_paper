## Copyright (c) 2019, Johannesson lab
## Licensed under the MIT license. See LICENSE file.


## Collect mutation rates.
mutation_rates <- read.csv2("~/Ring_data/mutation_analysis/mutation_rates.csv")

# Remove 08NW because of low coverage. It's in row 37.
mutation_rates <- mutation_rates[-c(37), ]


## Reformat
library(reshape2)

summed_cell_divisions <- c(sum(mutation_rates$Cell.divisions.per.sampling.point))
rates <- mutation_rates$Mutation.rate..all.callable.sites...adjusted.
rates_max <- mutation_rates$Mutation.rate..genomewide..maximum.based.on.hyphal.cell.size..adjusted.
rates_min <- mutation_rates$Mutation.rate..genomewide..minimum.based.on.hyphal.cell.size.
rates_noncoding <- mutation_rates$Mutation.rate..noncoding.regions..adjusted.

repeats <- mutation_rates$Cell.divisions.per.sampling.point
sites <- mutation_rates$Callable.sites.post.filtering..per.nucleus.

library(spatstat)
weighted.median(rates,repeats)
sd(rates)
median(rates)

weighted.median(rates_max,repeats)
weighted.median(rates_min,repeats)

weighted.median(rates_noncoding,repeats)
sd(rates_noncoding)
median(rates_noncoding)

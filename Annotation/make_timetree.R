library(ape)
tre <- read.tree("alignment_tree32.nwk")
nodes <- c(getMRCA(tre, c("dreissena", "ruditapes")),
getMRCA(tre, c("crassostrea", "pinctada")),
39,
38,
37,
26,
30,
33,
29,
23,
25,
22
)

age = c(303,
326,
427,
428,
440,
536,
208,
130,
236,
287,
4,
543
)

calibration <- makeChronosCalib(tre, node =nodes, age.min = age, age.max=age)

pl.tree <- chronos(tre, calibration=calibration)
write.tree(pl.tree, digits = 3, file = "mollusca_timetree.nw")

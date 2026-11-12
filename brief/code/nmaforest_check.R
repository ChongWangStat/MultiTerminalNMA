suppressPackageStartupMessages({library(NMAforest); library(netmeta)})
cat("R", R.version.string, "| netmeta", as.character(packageVersion("netmeta")), "| igraph", as.character(packageVersion("igraph")), "| NMAforest", as.character(packageVersion("NMAforest")), "\n")
data(example_data)
pdf(NULL)
res <- NMAforest(data = example_data, sm = "OR", reference = "x", model = "random",
                 comparison = c("x", "y"), study = "study", treat = "t", event = "r", N = "n",
                 study_id = "id", study_path = TRUE)
out <- res$output
print(out, digits = 6)
write.csv(out, "nmaforest_example_output.csv", row.names = FALSE)

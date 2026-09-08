# README
# This file is sourced by workflow/DGE.R,
# in order to create the metadata columns required for DEGs

################################# Metadata: condition & cluster

# Condition & cluster (e.g. 10dpci_veh_02 vs 10dpci_ccpe_02)
# Typo from obj2 to obj but it doesn't make a difference
obj2@meta.data[["condition_cluster"]] <- paste(obj@meta.data[[grping]],
                                               obj@meta.data[[clustering]],
                                               sep = "_")

# Check that the expected samples are present in this object
print("The following are the names of condition_cluster levels present in this object:")
obj2$condition_cluster %>% unique() %>% print()

######################### Metadata about whether there was a surgery
map <- c("Early_active" = "no_surgery",
         "PM_naive" = "no_surgery",
         "14dpsham" = "surgery",
         "10dpci_veh" = "surgery",
         "14dpci" = "surgery")

obj2@meta.data$surgery <- factor(
  map[as.character(obj2@meta.data[["orig.ident"]])],
  levels = unique(unname(map)))

############################# Metadata: surgery & cluster

obj2@meta.data[["surgery_cluster"]] <- paste(obj2@meta.data[["surgery"]],
                                             obj2@meta.data[[clustering]],
                                             sep = "_")

print(unique(obj2$surgery_cluster))

print(table(obj2$orig.ident, obj2$surgery))
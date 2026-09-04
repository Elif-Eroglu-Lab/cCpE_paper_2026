# README
# This file is sourced by workflow/DGE.R,
# in order to create the metadata columns required for DEGs

################################# Metadata to create

# Check that the expected samples are present in this object
print("The following are the names of conditions present in this object:")
obj$condition %>% unique() %>% print()

# Condition & cluster (e.g. vehicle_13, c-CpE_13)
obj2@meta.data[["condition_cluster"]] <- paste(obj2@meta.data[["condition"]],
                                               obj2@meta.data[[clustering]],
                                               sep = "_")

# Condition batch
obj2@meta.data[["condition_batch"]] <- paste(obj2@meta.data[["condition"]],
                                             obj2@meta.data[["batch"]],
                                             sep = "_")

# Condition batch cluster
obj2@meta.data[["condition_batch_cluster"]] <- paste(obj2@meta.data[["condition"]],
                                                     obj2@meta.data[["batch"]],
                                                     obj2@meta.data[[clustering]],
                                                     sep = "_")

# functional enrichment analysis and namespace conversion with gprofiler2

# INSTALLATION
install.packages("gprofiler2")
libarary(gprofiler2)

# GENE LIST FUNCTIONAL

gostres <- gost(query = e_data("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

# VISUALIZATION

gostplot(gostres, capped = TRUE, interactive = TRUE)

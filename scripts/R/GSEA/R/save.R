#' Save gene set summary
#'
#' Save gene set summary tables as an XLSX file
#'
#' @param gene_set_summary data.frame with gene set summary information from
#' `summarise_gene_set_results()`
#' @param path Path to output file
#'
#' @return `path` invisibly
save_gene_set_summary <- function(gene_set_summary, path) {

    summary_tables <- purrr::map(
        c("TotalSets", "UpSets", "DownSets"),
        function(.type) {
            gene_set_summary %>%
                dplyr::mutate(Sets = gene_set_summary[[.type]]) %>%
                dplyr::select(-UpSets, -DownSets, -TotalSets) %>%
                tidyr::pivot_wider(names_from = Coefficient, values_from = Sets)
        }
    )
    names(summary_tables) <- c("Total", "Up", "Down")

    writexl::write_xlsx(summary_tables, path)

    invisible(path)
}

#' Save GO cluster heat map
#'
#' Save a clustering heat map produced by `simplifyEnrichment::ht_clusters()`
#'
#' @param simplified_go_results data.frame containing simplified GO results
#' @param go_similarity matrix containing the similarity between GO IDs
#' calculated using `simplifyEnrichment::GO_similarity()`
#' @param path Path to the output file
#' @param title Title string for the plot, if `NULL` the number of terms and
#' clusters is used (this is also added to the given title)
#'
#' @return `path` invisibly
save_go_cluster_heatmap <- function(simplified_go_results, go_similarity, path,
                                    title = NULL) {

    n_clusters <- length(unique(simplified_go_results$Cluster))
    subtitle <- glue::glue(
        "{sum(!is.na(simplified_go_results$GOID))} terms, ",
        "{n_clusters} clusters"
    )
    if (is.null(title)) {
        title <- subtitle
    } else {
        title <- glue::glue("{title} ({subtitle})")
    }

    message(title)

    if (nrow(simplified_go_results) == 0) {
        message("No results, saving empty file")
        cairo_pdf(path)
        dev.off()
    } else if (n_clusters == 1) {
        message("Only one cluster, saving empty file")
        cairo_pdf(path)
        dev.off()
    } else {
        cairo_pdf(path, width = 12, height = 8)
        tryCatch(
            {
                go_ids <- simplified_go_results$GOID
                is_id_na <- is.na(go_ids)
                go_ids <- go_ids[!is_id_na]
                simplifyEnrichment::ht_clusters(
                    go_similarity[go_ids, go_ids],
                    simplified_go_results$Cluster[!is_id_na],
                    column_title = title
                )
            },
            error = function(err) {
                message("Saving heatmap failed")
            }
        )
        dev.off()
    }

    invisible(path)
}

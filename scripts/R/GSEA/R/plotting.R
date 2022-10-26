#' Plot gene set summary heat map
#'
#' Plot a heat map summarising results for a gene set collection
#'
#' @param gene_set_summary data.frame with gene set summary information from
#' `summarise_gene_set_results()`
#' @param type Which kind of sets to plot
#'
#' @return
plot_gene_set_summary_heatmap <- function(gene_set_summary,
                                          type = c("Total", "Up", "Down")) {

    type <- match.arg(type)

    gene_set_summary <- gene_set_summary %>%
        dplyr::mutate(Sets = gene_set_summary[[paste0(type, "Sets")]]) %>%
        dplyr::select(-UpSets, -DownSets, -TotalSets) %>%
        tidyr::complete(CellType, Coefficient, fill = list(Sets = NA))

    ggplot2::ggplot(
        gene_set_summary,
        ggplot2::aes(x = Coefficient, y = CellType, fill = Sets)
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(option = "magma") +
        ggplot2::labs(subtitle = glue::glue("Direction: {type}")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title  = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            panel.grid  = ggplot2::element_blank()
        )
}

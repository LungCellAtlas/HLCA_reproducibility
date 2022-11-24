#' Load results
#'
#' Loads DE modelling results for a single cell type
#'
#' @param base_dir_input Base input path for DE modelling
#' @param base_dir_output Base output path for DE modelling
#' @param cell_type The cell type to load
#'
#' @return Named list where names are coefficients and values are a `tibble`
#' with results for that coefficient
load_de_results <- function(base_dir_input, base_dir_output, cell_type) {

    input_dir  <- fs::path(base_dir_input, cell_type)
    output_dir <- fs::path(base_dir_output, cell_type)

    cli::cli_process_start("Reading modelling results")
    model_results <- readr::read_tsv(
        fs::path(output_dir, "mm_output.tsv"),
        col_types = readr::cols(.default = readr::col_double())
    )
    cli::cli_process_done()

    cli::cli_process_start("Reading gene name")
    gene_names <- readr::read_csv(
        fs::path(input_dir, "sample_gene_sums.csv"),
        n_max = 5,
        show_col_types = FALSE
    ) %>%
        colnames()

    if (gene_names[1] == "sample") {
        gene_names <- gene_names[-1]
    }

    if (length(gene_names) != nrow(model_results)) {
        cli::cli_abort(
            "Number of gene names does not match number of rows in results"
        )
    }
    cli::cli_process_done()

    cli::cli_process_start("Getting coefficients")
    coefs <- colnames(model_results)
    coefs <- coefs[stringr::str_detect(coefs, "Coef")]
    coefs <- stringr::str_remove(coefs, "Coef\\.")
    coefs <- coefs[coefs != "(Intercept)"]
    cli::cli_process_done()

    cli::cli_process_start("Tidying results")
    model_results <- model_results %>%
        dplyr::mutate(Gene = gene_names)

    results_list <- purrr::map(coefs, function(.coef) {
        res <- model_results %>%
            dplyr::select(Gene, AveExpr, tidyselect::matches(.coef))
        colnames(res) <- stringr::str_remove(
            colnames(res),
            paste0("\\.", .coef)
        )
        res
    })
    names(results_list) <- glue::glue("{cell_type}|{coefs}")
    cli::cli_process_done()

    results_list
}

#' Load uncert results
#'
#' Loads DE modelling results for uncert analysis
#'
#' @param file Path to CSV results file
#'
#' @return `tibble` with DE results
load_uncert_results <- function(file) {

    cli::cli_process_start("Reading modelling results")
    model_results <- readr::read_csv(
        file,
        show_col_types = FALSE
    )
    cli::cli_process_done()

    cli::cli_process_start("Tidying results")
    model_results <- model_results %>%
        dplyr::select(-1)
    cli::cli_process_done()

    model_results
}

#' Build GO ID map
#'
#' Build a map between GO terms and their GO IDs
#'
#' @return Named vector of GO IDs
build_goid_map <- function() {

    go_db <- AnnotationDbi::select(
        GO.db::GO.db,
        keys = AnnotationDbi::keys(GO.db::GO.db),
        columns = c("GOID", "ONTOLOGY", "TERM")
    ) %>%
        dplyr::filter(ONTOLOGY == "BP")

    goid_map <- go_db$GOID

    names(goid_map) <- go_db$TERM %>%
        stringr::str_to_upper() %>%
        stringr::str_replace_all(" ", "_") %>%
        stringr::str_replace_all("-", "_") %>%
        stringr::str_replace_all("/", "_") %>%
        stringr::str_remove_all("\\(") %>%
        stringr::str_remove_all("\\)") %>%
        stringr::str_replace_all("\\.", "_") %>%
        stringr::str_remove_all(",") %>%
        stringr::str_remove_all("'") %>%
        stringr::str_replace_all("POLYA_", "POLY_A_") %>%
        paste0("GO_", .)

    goid_map
}

#' Run CAMERA
#'
#' Run the **{limma}** CAMERA gene set test. Uses `cameraPR`.
#'
#' @param de_results data.frame with DE results
#' @param gene_sets List of gene sets
#' @param gene_map Vector where names are gene symbols and values are ENTREZ
#' gene IDs
#'
#' @return tibble of gene set results
run_camera <- function(de_results, gene_sets, gene_map) {

    gene_ids <- gene_map[de_results$Gene]

    set_indices <- limma::ids2indices(gene_sets, gene_ids)

    gene_set_names <- purrr::map(set_indices, ~ de_results$Gene[.x])

    de_up <- de_results %>%
        dplyr::filter(P.value <= 0.05) %>%
        dplyr::filter(Coef > 0)

    gene_set_names_up <- purrr::map(
        gene_set_names,
        ~ .x[.x %in% de_up$Gene]
    )

    de_down <- de_results %>%
        dplyr::filter(P.value <= 0.05) %>%
        dplyr::filter(Coef < 0)

    gene_set_names_down <- purrr::map(
        gene_set_names,
        ~ .x[.x %in% de_down$Gene]
    )

    limma::cameraPR(de_results$t, set_indices) %>%
        tibble::rownames_to_column("GeneSet") %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
            Genes  = gene_set_names[GeneSet],
            DEUp   = gene_set_names_up[GeneSet],
            DEDown = gene_set_names_down[GeneSet]
        ) %>%
        dplyr::mutate(
            Genes  = purrr::map_chr(Genes,  ~ paste(.x, collapse = ",")),
            DEUp   = purrr::map_chr(DEUp,   ~ paste(.x, collapse = ",")),
            DEDown = purrr::map_chr(DEDown, ~ paste(.x, collapse = ","))
        )
}

#' Run CAMERA diffxypy
#'
#' Run the **{limma}** CAMERA gene set test on results from **diffxpy**. Uses
#' `cameraPR`.
#'
#' @param de_results data.frame with DE results
#' @param gene_sets List of gene sets
#' @param gene_map Vector where names are gene symbols and values are ENTREZ
#' gene IDs
#'
#' @return tibble of gene set results
run_camera_diffxpy <- function(de_results, gene_sets, gene_map) {

    gene_ids <- gene_map[de_results$gene]

    set_indices <- limma::ids2indices(gene_sets, gene_ids)

    gene_set_names <- purrr::map(set_indices, ~ de_results$gene[.x])

    de_up <- de_results %>%
        dplyr::filter(qval <= 0.05) %>%
        dplyr::filter(log2fc > 0)

    gene_set_names_up <- purrr::map(
        gene_set_names,
        ~ .x[.x %in% de_up$gene]
    )

    de_down <- de_results %>%
        dplyr::filter(qval <= 0.05) %>%
        dplyr::filter(log2fc < 0)

    gene_set_names_down <- purrr::map(
        gene_set_names,
        ~ .x[.x %in% de_down$gene]
    )

    limma::cameraPR(de_results$coef_mle, set_indices) %>%
        tibble::rownames_to_column("GeneSet") %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
            Genes  = gene_set_names[GeneSet],
            DEUp   = gene_set_names_up[GeneSet],
            DEDown = gene_set_names_down[GeneSet]
        ) %>%
        dplyr::mutate(
            Genes  = purrr::map_chr(Genes,  ~ paste(.x, collapse = ",")),
            DEUp   = purrr::map_chr(DEUp,   ~ paste(.x, collapse = ",")),
            DEDown = purrr::map_chr(DEDown, ~ paste(.x, collapse = ","))
        )
}

#' Summarise gene set results
#'
#' Produce a summary table for gene set results for a given gene set collection
#'
#' @param gene_set_results List of gene set results data.frames
#' @param collection Name of the collection to summarise
#'
#' @return tibble with summary table
summarise_gene_set_results <- function(gene_set_results, collection) {

    is_collection <- stringr::str_detect(names(gene_set_results), collection)
    gene_set_results <- gene_set_results[is_collection]

    summ <- purrr::map_dfr(names(gene_set_results), function(.result_name) {
        setting <- stringr::str_split(.result_name, "\\|")[[1]]
        sig_results <- gene_set_results[[.result_name]] %>%
            dplyr::filter(FDR <= 0.05)

        tibble::tibble(
            CellType    = setting[1],
            Coefficient = setting[2],
            UpSets      = sum(sig_results$Direction == "Up"),
            DownSets    = sum(sig_results$Direction == "Down")
        )
    }) %>%
        dplyr::mutate(TotalSets = UpSets + DownSets)

    attr(summ, "Collection") <- collection

    return(summ)
}

#' Summarise gene set results by group
#'
#' Produce a summary table for gene set results for a given gene set collection
#' and a given set of groups of cell types
#'
#' @param gene_set_results List of gene set results data.frames
#' @param collection Name of the collection to summarise
#' @param groups Named list of vectors containing cell type groups
#'
#' @return tibble with summary table
summarise_gene_set_results_group <- function(gene_set_results, collection,
                                             groups) {

    is_collection <- stringr::str_detect(names(gene_set_results), collection)
    gene_set_results <- gene_set_results[is_collection]

    summ <- purrr::map_dfr(names(groups), function(.group) {
        cell_types <- stringr::str_split(names(gene_set_results), "\\|") %>%
            purrr::map_chr(~ .x[1])
        in_group <- cell_types %in% groups[[.group]]
        group_results <- gene_set_results[in_group]
        cell_types <- cell_types[in_group]

        celltype_sets <- purrr::map(unique(cell_types), function(.cell_type) {

            is_celltype <- stringr::str_detect(
                names(group_results),
                .cell_type
            )
            celltype_results <- group_results[is_celltype]

            coefficients <- stringr::str_split(
                names(celltype_results),
                "\\|"
            ) %>%
                purrr::map_chr(~ .x[2]) %>%
                unique()

            sig_sets <- purrr::map(coefficients, function(.coef) {

                result_name <- glue::glue("{.cell_type}|{.coef}|{collection}")

                up_sets <- celltype_results[[result_name]] %>%
                    dplyr::filter(
                        FDR <= 0.05,
                        Direction == "Up"
                    ) %>%
                    dplyr::pull(GeneSet)

                down_sets <- celltype_results[[result_name]] %>%
                    dplyr::filter(
                        FDR <= 0.05,
                        Direction == "Down"
                    ) %>%
                    dplyr::pull(GeneSet)

                list(Up = up_sets, Down = down_sets)
            })
            names(sig_sets) <- coefficients

            sig_sets
        })
        names(celltype_sets) <- unique(cell_types)

        coefficients <- stringr::str_split(
            names(group_results),
            "\\|"
        ) %>%
            purrr::map_chr(~ .x[2]) %>%
            unique()

        purrr::map_dfr(coefficients, function(.coef) {
            up_sets <- purrr::map(unique(cell_types), function(.cell_type) {
                celltype_sets[[.cell_type]][[.coef]][["Up"]]
            }) %>%
                purrr::reduce(intersect)

            down_sets <- purrr::map(unique(cell_types), function(.cell_type) {
                celltype_sets[[.cell_type]][[.coef]][["Down"]]
            }) %>%
                purrr::reduce(intersect)

            tibble::tibble(
                Coefficient = .coef,
                Up = length(up_sets),
                Down = length(down_sets),
                UpSets = paste(up_sets, collapse = ","),
                DownSets = paste(down_sets, collapse = ",")
            )
        }) %>%
            dplyr::mutate(CellType = .group)
    }) %>%
        dplyr::mutate(Total = Up + Down) %>%
        dplyr::relocate(CellType)

    attr(summ, "Collection") <- collection

    return(summ)
}

#' Simplify GO results
#'
#' Simplify GO gene set results by clustering them using the method from
#' **{simplifyEnrichment}**
#'
#' @param go_results data.frame of GO results
#' @param direction Direction to simplify ("Up", "Down" or "Total" (all))
#' @param goid_map Named vector mapping GO terms to GO IDs
#' @param go_similarity matrix containing the similarity between GO IDs
#' calculated using `simplifyEnrichment::GO_similarity()`
#'
#' @return data.frame of clustered GO results
simplify_go_results <- function(go_results, direction, goid_map,
                                go_similarity) {

    `%>%` <- magrittr::`%>%`

    go_terms <- get_sig_sets(go_results, direction)

    message("Mapping terms to IDs...")
    go_ids <- goid_map[go_terms]
    if (any(is.na(go_ids))) {
        message(sum(is.na(go_ids)), " terms failed to match!")
        go_ids <- go_ids[!is.na(go_ids)]
    }

    if (length(go_ids) < 2) {
        message("Less than two GO IDs found, returning empty results")
        return(
            tibble::tibble(
                id      = character(),
                term    = character(),
                cluster = numeric()
            )
        )
    }

    message("Simplifying GO terms...")
    sim_results <- simplifyEnrichment::simplifyGO(
        go_similarity[go_ids, go_ids],
        plot = FALSE
    ) %>%
        dplyr::select(GOID = id, Cluster = cluster)

    go_results %>%
        dplyr::filter(GeneSet %in% go_terms) %>%
        dplyr::mutate(GOID = goid_map[GeneSet]) %>%
        dplyr::left_join(sim_results, by = "GOID") %>%
        dplyr::relocate(GeneSet, GOID, NGenes, Direction, PValue, FDR, Cluster)
}

#' Get significant sets
#'
#' Extract significant gene sets from a data.frame containing gene set results
#'
#' @param gene_set_results data.frame containing gene set results
#' @param direction Whether to select gene sets that are "Up", "Down" or "Total"
#' (all)
#'
#' @return vector of significant gene sets
get_sig_sets <- function(gene_set_results,
                         direction = c("Total", "Up", "Down")) {

    direction <- match.arg(direction)

    if (direction %in% c("Up", "Down")) {
        gene_set_results <- gene_set_results %>%
            dplyr::filter(Direction == direction)
    }

    gene_set_results %>%
        dplyr::filter(FDR <= 0.05) %>%
        dplyr::pull(GeneSet)
}

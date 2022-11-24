library(targets)
library(here)

source(here("R", "load.R"))
source(here("R", "gene_sets.R"))
source(here("R", "plotting.R"))
source(here("R", "save.R"))

tar_option_set(packages = "magrittr")

list(
    # --- Set paths ---
    tar_target(
        base_dir,
        fs::path(
            "/",
            "..",
            "..",
            "..",
            "results",
            "covariate_modeling",
        )
    ),
    tar_target(
        base_dir_input,
        fs::path(
            base_dir,
            "input",
        )
    ),
    tar_target(
        base_dir_output,
        fs::path(
            base_dir,
            "output",
        )
    ),
    # --- Load DE results ---
    tar_target(
        cell_types,
        fs::dir_ls(base_dir_output) %>%
            fs::path_file()
    ),
    tar_target(
        celltype_de,
        load_de_results(base_dir_input, base_dir_output, cell_types),
        pattern = map(cell_types)
    ),
    tar_target(
        de_results,
        {
            res_list <- celltype_de
            names(res_list) <- stringr::str_remove(
                names(res_list),
                "celltype_de_[0-9a-f]+_"
            )
            res_list
        }
    ),
    # --- Map genes to ENTREZ IDs ---
    tar_target(
        all_genes,
        de_results %>%
            purrr::map(~ dplyr::pull(.x, Gene)) %>%
            purrr::reduce(~ c(.x, .y)) %>%
            unique()
    ),
    tar_target(
        gene_map,
        AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys    = all_genes,
            keytype = "SYMBOL",
            column  = "ENTREZID"
        )
    ),
    # --- Get gene sets ---
    tar_target(
        kegg_sets_url,
        "https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c2.cp.kegg.v7.1.entrez.rds",
        format = "url"
    ),
    tar_target(
        reactome_sets_url,
        "https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c2.cp.reactome.v7.1.entrez.rds",
        format = "url"
    ),
    tar_target(
        go_bp_sets_url,
        "https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c5.bp.v7.1.entrez.rds",
        format = "url"
    ),
    tar_target(
        hallmark_sets_url,
        "https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds",
        format = "url"
    ),
    tar_target(
        gene_set_collections,
        list(
            KEGG     = readRDS(url(kegg_sets_url)),
            REACTOME = readRDS(url(reactome_sets_url)),
            GO_BP    = readRDS(url(go_bp_sets_url)),
            HALLMARK = readRDS(url(hallmark_sets_url))
        )
    ),
    # --- Run gene set analysis ---
    tar_target(
        all_gene_set_results,
        {
            res <- list(
                run_camera(
                    de_results[[1]],
                    gene_set_collections[[1]],
                    gene_map
                )
            )
            names(res) <- paste(
                names(de_results),
                names(gene_set_collections),
                sep = "|"
            )
            res
        },
        pattern = cross(de_results, gene_set_collections),
        iteration = "list"
    ),
    tar_target(
        gene_set_results,
        {
            purrr::flatten(all_gene_set_results)
        }
    ),
    tar_target(
        gene_set_results_tsv,
        {
            out_path <- fs::path(
                "output",
                "camera_results",
                stringr::str_replace_all(names(gene_set_results), "\\|", "-"),
                ext = "tsv"
            )
            readr::write_tsv(gene_set_results[[1]], out_path)
            out_path
        },
        format = "file",
        pattern = map(gene_set_results)
    ),
    # --- Summarise gene set results ---
    tar_target(
        gene_set_summary,
        summarise_gene_set_results(
            gene_set_results,
            names(gene_set_collections)
        ),
        pattern = map(gene_set_collections),
        iteration = "list"
    ),
    tar_target(
        directions,
        c("Total", "Up", "Down")
    ),
    tar_target(
        gene_set_summary_heatmap,
        {
            collection <- attr(gene_set_summary, "Collection")
            plot <- plot_gene_set_summary_heatmap(
                gene_set_summary,
                type = directions
            ) +
                ggplot2::ggtitle(collection)
            ggplot2::ggsave(
                here::here(
                    "output",
                    "summary_heatmaps",
                    glue::glue("{collection}-{directions}.png")
                )
            )
            plot
        },
        pattern = cross(gene_set_summary, directions),
        iteration = "list"
    ),
    tar_target(
        gene_set_summary_xlsx,
        {
            collection <- attr(gene_set_summary, "Collection")
            save_gene_set_summary(
                gene_set_summary,
                here::here(
                    "output",
                    "summary_tables",
                    glue::glue("{collection}.xlsx")
                )
            )
        },
        pattern = map(gene_set_summary)
    ),
    # --- Check for common results ---
    tar_target(
        celltype_groups,
        list(
            Epithelial = c("Basal", "Multiciliated", "Secretory",
                           "Transitional_Club-AT2", "Ionocyte", "Tuft",
                           "Neuroendocrine", "SMG_serous", "SMG_mucous",
                           "SMG_duct", "AT1", "AT2"),
            Endothelial = c("EC_arterial", "EC_capillary", "EC_venous",
                            "Lymphatic_EC"),
            Stromal = c("Peribronchial_fibroblasts", "Adventitial_fibroblasts",
                        "Alveolar_fibroblasts", "Pericytes",
                        "Subpleural_fibroblasts", "Myofibroblasts",
                        "Smooth_muscle", "Fibromyocytes", "Mesothelium"),
            Immune = c("B_cell", "Plasma_cells", "T_cell_lineage",
                       "Innate_lymphoid_cell_NK", "DC1", "DC2",
                       "Migratory_DCs", "Plasmacytoid_DCs",
                       "Alveolar_macrophages", "Monocyte-derived_Mφ",
                       "Interstitial_Mφ_perivascular", "Monocytes",
                       "Mast_cells")
        )
    ),
    tar_target(
        gene_set_summary_groups,
        summarise_gene_set_results_group(
            gene_set_results,
            names(gene_set_collections),
            celltype_groups
        ),
        pattern = map(gene_set_collections),
        iteration = "list"
    ),
    # --- Simplify GO results ---
    tar_target(
        goid_map,
        build_goid_map()
    ),
    tar_target(
        go_similarity,
        simplifyEnrichment::GO_similarity(
            goid_map[names(goid_map) %in% names(gene_set_collections$GO_BP)],
            ont = "BP"
        )
    ),
    tar_target(
        gene_set_results_go,
        {
            is_go <- stringr::str_detect(names(gene_set_results), "GO_BP")
            gene_set_results[is_go]
        }
    ),
    tar_target(
        all_simplified_go_results,
        list(simplify_go_results(
            go_results    = gene_set_results_go[[1]],
            direction     = directions,
            goid_map      = goid_map,
            go_similarity = go_similarity
        )) %>%
            setNames(paste(names(gene_set_results_go), directions, sep = "|")),
        pattern = cross(gene_set_results_go, directions),
        iteration = "list"
    ),
    tar_target(
        simplified_go_results,
        {
            purrr::flatten(all_simplified_go_results) %>%
                purrr::keep(~ nrow(.x) > 0)
        }
    ),
    tar_target(
        simplified_go_tsv,
        {
            out_path <- fs::path(
                "output",
                "simplified_go",
                stringr::str_replace_all(
                    names(simplified_go_results),
                    "\\|",
                    "-"
                ),
                ext = "tsv"
            )
            readr::write_tsv(simplified_go_results[[1]], out_path)
            out_path
        },
        format = "file",
        pattern = map(simplified_go_results)
    ),
    tar_target(
        simplified_go_heatmap_pdf,
        {
            out_path <- fs::path(
                "output",
                "simplified_go",
                "heatmaps",
                stringr::str_replace_all(
                    names(simplified_go_results),
                    "\\|",
                    "-"
                ),
                ext = "pdf"
            )
            save_go_cluster_heatmap(
                simplified_go_results[[1]],
                go_similarity,
                out_path,
                title = names(simplified_go_results)
            )
            out_path
        },
        format = "file",
        pattern = map(simplified_go_results)
    ),
    # --- Load uncert DE results ---
    tar_target(
        uncert_dir,
        fs::path(
            "..",
            "..",
            "..",
            "results",
            "fibrosis_to_HLCA",
            "DEAs",
        )
    ),
    tar_target(
        all_uncert_files,
        fs::dir_ls(uncert_dir)
    ),
    tar_target(
        uncert_files,
        all_uncert_files,
        pattern = map(all_uncert_files),
        format = "file"
    ),
    tar_target(
        uncert_de,
        list(load_uncert_results(uncert_files)) %>%
            setNames(fs::path_ext_remove(fs::path_file(uncert_files))),
        pattern = map(uncert_files)
    ),
    tar_target(
        uncert_de_results,
        {
            res_list <- uncert_de
            names(res_list) <- stringr::str_remove(
                names(res_list),
                "uncert_de_[0-9a-f]+_"
            )
            res_list
        }
    ),
    tar_target(
        all_uncert_gene_set_results,
        {
            res <- list(
                run_camera_diffxpy(
                    uncert_de_results[[1]],
                    gene_set_collections[[1]],
                    gene_map
                )
            )
            names(res) <- paste(
                names(uncert_de_results),
                names(gene_set_collections),
                sep = "|"
            )
            res
        },
        pattern = cross(uncert_de_results, gene_set_collections),
        iteration = "list"
    ),
    tar_target(
        uncert_gene_set_results,
        {
            purrr::flatten(all_uncert_gene_set_results)
        }
    ),
    tar_target(
        uncert_gene_set_results_tsv,
        {
            out_path <- fs::path(
                "output",
                "uncert_camera_results",
                stringr::str_replace_all(
                    names(uncert_gene_set_results),
                    "\\|",
                    "-"
                ),
                ext = "tsv"
            )
            readr::write_tsv(uncert_gene_set_results[[1]], out_path)
            out_path
        },
        format = "file",
        pattern = map(uncert_gene_set_results)
    )
)

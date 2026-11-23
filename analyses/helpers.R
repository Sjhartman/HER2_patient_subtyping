# Shared helpers (sourced by step scripts).

# Current HGNC symbols for TCGA gene names (needed to match MSigDB / marker gene lists). A retired name is renamed only
# if its current symbol is not already a gene in the data and no other name maps to it; otherwise it keeps its TCGA
# name. Without this, alias collisions overwrite real genes (e.g. C2orf48 -> RRM2, CCRL1 -> CX3CR1).
current_symbols <- function(genes) {
  sym <- dplyr::coalesce(limma::alias2SymbolTable(genes, "Hs"), genes)
  ok <- sym == genes | !(sym %in% genes | duplicated(sym) | duplicated(sym, fromLast = TRUE))
  out <- ifelse(ok, sym, genes)
  stopifnot(!anyDuplicated(out))
  out
}

# MCP-counter marker genes (Becht 2016; analyses/ref). The source file lists CAVIN2 (Entrez 8436) as "C".
mcp_markers <- function(path = "analyses/ref/mcpcounter_genes.txt") {
  readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::transmute(gene = dplyr::if_else(`HUGO symbols` == "C" & ENTREZID == 8436, "CAVIN2", `HUGO symbols`),
                     population = `Cell population`)
}

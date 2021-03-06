# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

grab <- function(i, j) {
    .Call(`_OASW_grab`, i, j)
}

update_dys_1_and_dys_2 <- function(n, rep_ind, dys_1, dys_2, dys, max_dys) {
    invisible(.Call(`_OASW_update_dys_1_and_dys_2`, n, rep_ind, dys_1, dys_2, dys, max_dys))
}

bswap <- function(K, n, rep_ind, dys_1, dys_2, effect, dys, max_dys) {
    invisible(.Call(`_OASW_bswap`, K, n, rep_ind, dys_1, dys_2, effect, dys, max_dys))
}

hpsort <- function(n, vect) {
    invisible(.Call(`_OASW_hpsort`, n, vect))
}

silhouettes <- function(K, n, medoids, clust, clus_size, dys, dys_j, avg_dys_clus, sil) {
    invisible(.Call(`_OASW_silhouettes`, K, n, medoids, clust, clus_size, dys, dys_j, avg_dys_clus, sil))
}

hpsortint <- function(n, vect) {
    invisible(.Call(`_OASW_hpsortint`, n, vect))
}

silswap <- function(K, n, rep_ind, medoids, altmeds, clust, clus_size, dys, dys_1, dys_i, avg_dys_clus, silh, altsilh, iter) {
    invisible(.Call(`_OASW_silswap`, K, n, rep_ind, medoids, altmeds, clust, clus_size, dys, dys_1, dys_i, avg_dys_clus, silh, altsilh, iter))
}

clusanal <- function(K, n, rep_ind, dys_1, dys, max_dys, medoid, clus_size, clus_diam, clus_sep, clus_dys_1_avg, clus_dys_1_max, index_by_cluster, clus_vect, avg_dys_clust, a, b, silh, avg_clus_silh, avg_silh) {
    invisible(.Call(`_OASW_clusanal`, K, n, rep_ind, dys_1, dys, max_dys, medoid, clus_size, clus_diam, clus_sep, clus_dys_1_avg, clus_dys_1_max, index_by_cluster, clus_vect, avg_dys_clust, a, b, silh, avg_clus_silh, avg_silh))
}

sil_lab <- function(K, n, clus_lab, clus_size, disty, dys_j, avg_dys_clus, sil) {
    invisible(.Call(`_OASW_sil_lab`, K, n, clus_lab, clus_size, disty, dys_j, avg_dys_clus, sil))
}

sil_lab_swap <- function(K, n, clus_lab, alt_clus_lab, clus_size, disty, iter, dys_i, avg_dys_clus, silh, altsilh) {
    invisible(.Call(`_OASW_sil_lab_swap`, K, n, clus_lab, alt_clus_lab, clus_size, disty, iter, dys_i, avg_dys_clus, silh, altsilh))
}

clustyanlys <- function(K, n, clus_lab, clus_size, silh, avg_clus_silh, avg_clus_silhtwo, avg_silh) {
    invisible(.Call(`_OASW_clustyanlys`, K, n, clus_lab, clus_size, silh, avg_clus_silh, avg_clus_silhtwo, avg_silh))
}

sil_lab_link <- function(K, n, clus_lab, clus_size, disty, dys_j, avg_dys_clus, sil) {
    invisible(.Call(`_OASW_sil_lab_link`, K, n, clus_lab, clus_size, disty, dys_j, avg_dys_clus, sil))
}

hosil_lab_swap <- function(n, disty, clus_lab, clus_size, dys_j, avg_dys_clus, sil, best_clus_lab, all_best_clus_lab, all_best_avg_silh, copy_best_clus_lab) {
    invisible(.Call(`_OASW_hosil_lab_swap`, n, disty, clus_lab, clus_size, dys_j, avg_dys_clus, sil, best_clus_lab, all_best_clus_lab, all_best_avg_silh, copy_best_clus_lab))
}


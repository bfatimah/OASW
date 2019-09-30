#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;
#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) > (B) ? (B) : (A))
#define EPSILON 1.e-8

int grab(int i, int j);

void update_dys_1_and_dys_2(int n, IntegerVector rep_ind, NumericVector dys_1, NumericVector dys_2, NumericVector dys, double max_dys);

void bswap(int K, int n, IntegerVector rep_ind, NumericVector dys_1, NumericVector dys_2,
           NumericVector effect, NumericVector dys, double max_dys);

void hpsort(int n, NumericVector vect);

void silhouettes(int K, int n, IntegerVector medoids, IntegerVector clust, IntegerVector clus_size,
                 NumericVector dys, NumericVector dys_j, NumericVector avg_dys_clus,
                 NumericVector sil);

void hpsortint(int n, IntegerVector vect);


void silswap(int K, int n, IntegerVector rep_ind, IntegerVector medoids, IntegerVector altmeds,
             IntegerVector clust, IntegerVector clus_size,  NumericVector dys, NumericVector dys_1,
             NumericVector dys_i, NumericVector avg_dys_clus, NumericVector silh,
             NumericVector altsilh,  IntegerVector iter);

void clusanal(int K, int n, IntegerVector rep_ind, NumericVector dys_1, NumericVector dys,
              double max_dys, IntegerVector medoid, IntegerVector clus_size,
              NumericVector clus_diam, NumericVector clus_sep,
              NumericVector clus_dys_1_avg, NumericVector clus_dys_1_max,
              IntegerVector index_by_cluster,
              IntegerVector clus_vect, NumericVector avg_dys_clust, NumericVector a,
              NumericVector b, NumericVector silh, NumericVector avg_clus_silh,
              NumericVector avg_silh);

void sil_lab(int K, int n, IntegerVector  clus_lab,  IntegerVector  clus_size,  NumericVector disty,  NumericVector dys_j,  NumericVector avg_dys_clus,  NumericVector sil);

void sil_lab_swap(int K, int n, IntegerVector  clus_lab, IntegerVector  alt_clus_lab, IntegerVector  clus_size,  NumericVector disty, IntegerVector  iter,  NumericVector dys_i,  NumericVector avg_dys_clus,  NumericVector silh,  NumericVector altsilh,   NumericVector avg_clus_silh, NumericVector avg_clus_silhtwo, NumericVector avg_silh);


void clustyanlys(int K, int n, IntegerVector  clus_lab, IntegerVector  clus_size, NumericVector disty, NumericVector silh, NumericVector avg_clus_silh, NumericVector avg_clus_silhtwo, NumericVector avg_silh);

void sil_lab_link(int K, int n, IntegerVector clus_lab,  IntegerVector clus_size, NumericVector disty, NumericVector dys_j, NumericVector avg_dys_clus, NumericVector sil);

void hosil_lab_swap(int n, NumericVector disty, IntegerVector clus_lab, IntegerVector clus_size,  NumericVector dys_j, NumericVector avg_dys_clus, NumericVector sil, IntegerVector best_clus_lab, IntegerVector all_best_clus_lab, NumericVector all_best_avg_silh, IntegerVector copy_best_clus_lab);

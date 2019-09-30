//
//  trial.cpp
//  
//
//  Created by Fatima Batool on 18/09/2016.
//
//

#include "header.h"

#include <iostream>

#include <Rcpp.h>
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

// [[Rcpp::export]]
int grab(int i, int j){
    if (i == j){
        return 0;
    } else {
        if(i > j){
            int temp;
            temp = i;
            i = j;
            j = temp;
        }
        return((j*(j-1))/2)+1+i;
    }
}


// [[Rcpp::export]]

void update_dys_1_and_dys_2(int n, IntegerVector rep_ind, NumericVector dys_1, NumericVector dys_2, NumericVector dys, double max_dys)
{
    int i, j;
    
    for (i = 0; i < n; i++) {
        dys_1[i] = (1.1 * max_dys) + 1;
        dys_2[i] = (1.1 * max_dys) + 1;
        for (j = 0; j < n; j++) {
            if (rep_ind[j] != 1) continue;
            if (dys[grab(i, j)] < dys_1[i]) {
                dys_2[i] = dys_1[i];
                dys_1[i] = dys[grab(i,j)];
            } else {
                if (dys[grab(i, j)] < dys_2[i]) {
                    dys_2[i] = dys[grab(i, j)];
                }
            }
        }
    }
    
}

// [[Rcpp::export]]
void bswap(int K, int n, IntegerVector rep_ind, NumericVector dys_1, NumericVector dys_2,
           NumericVector effect, NumericVector dys, double max_dys){
    int i, j, h;
    int curr_k = 0.;
    double best_effect;
    int best_object1, best_object2;
    int keep_swapping = 1;
    double obj_fun;
    obj_fun = 0;

    
// Getting k initial Medoids (built phase of PAM)
    while(curr_k < K) {
        best_effect = 0.;
        for (i = 0; i < n; i++) {
            if(rep_ind[i] == 1)
                continue;
            effect[i] = 0.;
            for(j = 0; j < n; j++) {
                effect[i] += max(dys_1[j] - dys[grab(i,j)], 0.);
            }
            if (effect[i] > best_effect) {
                best_effect = effect[i];
                best_object1 = i;
            } // if closed
        } // for closed
    rep_ind[best_object1] = 1;
    for(i = 0; i < n; i++){
        dys_1[i] = min(dys_1[i], dys[grab(i, best_object1)]);
    }
    curr_k++;
    } //while closed
    for (i = 0; i < n; i++)
        obj_fun += dys_1[i];

    if(K == 1) return;
    
    while (keep_swapping == 1) {
        update_dys_1_and_dys_2(n, rep_ind, dys_1, dys_2, dys, max_dys);
        best_effect = 0.;
        for (h = 0; h < n; h++) {
            if (rep_ind[h] == 1)	/* skip representative objects */
                continue;
            for (i = 0; i < n; i++) {
                if (rep_ind[i] == 0)	/* skip non-representative objects */
                    continue;
                effect[h] = 0.;
                for (j = 0; j < n; j++) {
                    if (fabs(dys[grab(i,j)] - dys_1[j]) > EPSILON) {
                        if (dys[grab(j,h)] < dys_1[j]) {
                            effect[h] += dys[grab(j,h)] - dys_1[j];
                        }
                    } else {
                        effect[h] += min(dys_2[j], dys[grab(j,h)]) - dys_1[j];
                    }
                }
                if (effect[h] > best_effect) continue;
                best_effect = effect[h];
                best_object1 = h;
                best_object2 = i;
            }
        }
        if (fabs(best_effect) < EPSILON) {
            keep_swapping = 0;
        } else {
            rep_ind[best_object1] = 1;
            rep_ind[best_object2] = 0;
            obj_fun += best_effect;
        }
    }
    
    
    
return;


} // bswap close

// [[Rcpp::export]]
void hpsort(int n, NumericVector vect)
{
    /* sorts the elements of a vector in place from */
    /* smallest to largest numerical value. */
    /* n is the length of the vector to be sorted. */
    /* vect is n doubles to be sorted. */
    
    int i,ir,j,l;
    double rra;
    
    if(n<2) return;
    l=(n>>1);
    ir=(n-1);
    for(;;){
        if(l>0) rra=vect[--l];
        else{
            rra=vect[ir];
            vect[ir]=vect[0];
            if(--ir==0){
                vect[0]=rra;
                break;
            }
        }
        i=l;
        j=(l+l+1);
        while(j<=ir){
            if(j<ir && vect[j]<vect[j+1]) j++;
            if(rra<vect[j]){
                vect[i]=vect[j];
                i=j;
                j<<=1;
            }
            else break;
        }
        vect[i]=rra;
    }
}


// [[Rcpp::export]]
void silhouettes(int K, int n, IntegerVector medoids, IntegerVector clust, IntegerVector clus_size,
                 NumericVector dys, NumericVector dys_j, NumericVector avg_dys_clus,
                 NumericVector sil)
{
    /* Calculates the silhouettes of each element */
    /*  based on a set of medoids and distance */
    /*  matrix dys. Also: clusters and avg_dys_clus*/
    
    int i,j;
    double temp1,temp2;
    
    if(K==1) return;
    
    for(j=0;j<n;j++){
        for(i=0;i<K;i++){
            dys_j[i]=dys[grab(j,medoids[i])];
            avg_dys_clus[(j*K)+i]=0.;          /* clear avg_dys_clus */
            if(j==0) clus_size[i]=0;           /* clear clus_size */
        }
        hpsort(K,dys_j);
        for(i=0;i<K;i++){
            if(dys[grab(j,medoids[i])]==dys_j[0]){
                clust[j]=i;
                clus_size[i]+=1;
            }
        }
    }
    
    
    for(j=0;j<n;j++){
        for(i=(j+1);i<n;i++){
            temp1=dys[grab(i,j)];
            temp2=0.;
            if(clust[i]==clust[j]) temp2=-1;
            avg_dys_clus[(j*K)+clust[i]]+=(temp1/(clus_size[clust[i]]+temp2));
            avg_dys_clus[(i*K)+clust[j]]+=(temp1/(clus_size[clust[j]]+temp2));
        }
    }
    for(j=0;j<n;j++){
        if(clus_size[clust[j]]==1){
            sil[j]=0.;
            continue;
        }
        for(i=0;i<K;i++){
            if(clust[j]==i) dys_j[i]=0.;
            else dys_j[i]=avg_dys_clus[(j*K)+i];
        }
        hpsort(K,dys_j);
        if(dys_j[0]<0) dys_j[1]=dys_j[0];
        sil[j]=(dys_j[1]-avg_dys_clus[(j*K)+clust[j]])/max(dys_j[1],avg_dys_clus[(j*K)+clust[j]]);
    }
}

// [[Rcpp::export]]
void hpsortint(int n, IntegerVector vect)
{
    /* sorts the elements of a vector in place from */
    /* smallest to largest numerical value. */
    /* n is the length of the vector to be sorted. */
    /* vect is n doubles to be sorted. */
    
    int i,ir,j,l;
    int rra;
    
    if(n<2) return;
    l=(n>>1);
    ir=(n-1);
    for(;;){
        if(l>0) rra=vect[--l];
        else{
            rra=vect[ir];
            vect[ir]=vect[0];
            if(--ir==0){
                vect[0]=rra;
                break;
            }
        }
        i=l;
        j=(l+l+1);
        while(j<=ir){
            if(j<ir && vect[j]<vect[j+1]) j++;
            if(rra<vect[j]){
                vect[i]=vect[j];
                i=j;
                j<<=1;
            }
            else break;
        }
        vect[i]=rra;
    }
}


//the following two functions are for PAMSIL. silswap and clusanal

// [[Rcpp::export]]
void silswap(int K, int n, IntegerVector rep_ind, IntegerVector medoids, IntegerVector altmeds,
             IntegerVector clust, IntegerVector clus_size,  NumericVector dys, NumericVector dys_1,
             NumericVector dys_i, NumericVector avg_dys_clus, NumericVector silh,
             NumericVector altsilh,  IntegerVector iter)

{
    int h,i,j,k;
    int besth=0,besti=0;
    int keep_swapping=1;
    double effect,best_effect;
    double obj_fun;
    
    
    if(K==1) return;
    
    for(i=0,j=0;j<n;j++){             /* define medoids,altmeds */
        if(rep_ind[j]==1){
            medoids[i]=j;
            altmeds[i]=j;
            i++;
        }
    }
    
    silhouettes(K,n,medoids,clust,clus_size,dys,dys_i,avg_dys_clus,silh);
    iter[0] = 1;
    obj_fun=0.;                       /* define obj_fun */
    for(j=0;j<n;j++){
        obj_fun += silh[j];
    }
    obj_fun = obj_fun/n;
    
    best_effect=EPSILON/10;
    while(keep_swapping==1){
        silhouettes(K,n,medoids,clust,clus_size,dys,dys_i,avg_dys_clus,silh);
        for(i=0;i<K;i++){
            for(h=0;h<n;h++){
                if(rep_ind[h]==1) continue; /* skip medoids */
                altmeds[i]=h;
                silhouettes(K,n,altmeds,clust,clus_size,dys,dys_i,avg_dys_clus,altsilh);
                effect=0.;
                for(j=0;j<n;j++){
                    effect+=altsilh[j];
                    effect-=silh[j];
                }
                effect=effect/n;
                if(effect>best_effect){
                    best_effect=effect;
                    besth=h;
                    besti=i;
                }
            }
            altmeds[i]=medoids[i];
        }
        iter[0] += 1;
        if(best_effect<EPSILON){
            keep_swapping=0;
        }
        else{                          /* swap h for i */
            rep_ind[medoids[besti]]=0;
            rep_ind[besth]=1;
            medoids[besti]=besth;
            altmeds[besti]=besth;
            obj_fun += best_effect;
            best_effect=EPSILON/10;
        }
    }
    hpsortint(K,medoids);            /* put medoids in ascending order */
    for(i=0,j=0;j<n;j++){            /* update dys_1 */
        if(j==medoids[i]){
            dys_1[j]=0.;
            if(i<(K-1)) i++;
        }
        else{
            dys_1[j]=dys[grab(j,medoids[0])];
            for(k=1;k<K;k++){
                dys_i[k]=dys[grab(j,medoids[k])];
                if(dys_i[k]<dys_1[j]) dys_1[j]=dys_i[k];
            }
        }
    }
    return;
}
// [[Rcpp::export]]
void clusanal(int K, int n, IntegerVector rep_ind, NumericVector dys_1, NumericVector dys,
              double max_dys, IntegerVector medoid, IntegerVector clus_size,
              NumericVector clus_diam, NumericVector clus_sep,
              NumericVector clus_dys_1_avg, NumericVector clus_dys_1_max,
              IntegerVector index_by_cluster,
              IntegerVector clus_vect, NumericVector avg_dys_clust, NumericVector a,
              NumericVector b, NumericVector silh, NumericVector avg_clus_silh,
              NumericVector avg_silh)
{
    
    static int i,j,h;		/* counters */
    static int current_medoid;	/* holds index of a medoid */
    static double temp1, temp2;
    current_medoid = 0;
    h = 0;
    //fillme(0., avg_dys_clust, n * K);
    
    for (i = 0; i < K; i++) {	/* basic info on the clustering */
        
        while ( rep_ind[0 + current_medoid] != 1) {
            current_medoid++;
        }
        medoid[i] = current_medoid;
        
        clus_vect[current_medoid] = i;
        
        clus_size[i] = 1;
        
        clus_diam[i] = 0.;
        
        clus_sep[i] = max_dys * 1.1;
        
        clus_dys_1_avg[i] = 0.;
        
        clus_dys_1_max[i] = 0.;
        
        index_by_cluster[0 + h] = current_medoid;
        h++;
        for (j = 0; j < n; j++) {
            if (rep_ind[j] > 0) continue;
            if (fabs(dys_1[j] - dys[grab(current_medoid,j)]) > EPSILON)
                continue;
            clus_size[i]++;
            clus_dys_1_avg[i] += dys_1[j];
            if (dys_1[j] > clus_dys_1_max[i])
                clus_dys_1_max[i] = dys_1[j];
            rep_ind[j] = 2;
            clus_vect[j] = i;
            index_by_cluster[0+h] = j;
            h++;
        }
        clus_dys_1_avg[i] /= clus_size[i];
        current_medoid++;
    }
    
   
    
    for (j = 0; j < n; j++) {	/* advanced info on clustering */
        for (i = (j + 1); i < n; i++) {
            temp1 = dys[grab(i,j)];
            temp2 = 0.;
            if ( clus_vect[i] == clus_vect[j] ) {
                temp2 = -1;
                if (temp1 > clus_diam[clus_vect[i]]) {
                    clus_diam[clus_vect[i]] = temp1;
                }
            } else {
                if (temp1 < clus_sep[clus_vect[i]]) {
                    clus_sep[clus_vect[i]] = temp1;
                }
                if (temp1 < clus_sep[clus_vect[j]]) {
                    clus_sep[clus_vect[j]] = temp1;
                }
            }
            
            avg_dys_clust[ (j * K) + clus_vect[i]] +=\
            temp1/(clus_size[clus_vect[i]] + temp2);
            avg_dys_clust[ (i * K) + clus_vect[j]] +=\
            temp1/(clus_size[clus_vect[j]] + temp2);
            
        }
    }
    
    
    for (j = 0; j < n; j++) {	/* silhouettes */
        if (rep_ind[j] == 1 && clus_size[clus_vect[j]] == 1) {
            silh[j] = 0.;
            continue;
        }
        if (clus_vect[j] == 0) {
            b[j] = avg_dys_clust[(j * K) + 1];
        } else {
            b[j] = avg_dys_clust[(j * K) + 0];
        }
        for (i = 0; i < K; i++) {
            if (clus_vect[j] == i) {
                a[j] = avg_dys_clust[(j * K) + i];
            } else {
                if (avg_dys_clust[(j * K) + i] < b[j]) {
                    b[j] = avg_dys_clust[(j * K) + i];
                }
            }
        }
        silh[j] = (b[j] - a[j]) / max(b[j], a[j]);
    }
    for (i=0; i<n; i++)
    h = 0;
    for (i = 0; i < K; i++) {
        avg_clus_silh[i] = 0.;
        for (j = 0; j < clus_size[i]; j++) {
            avg_silh[0] += silh[index_by_cluster[h]] / n;
            avg_clus_silh[i] += silh[index_by_cluster[h]] / clus_size[i];
            h++;
        }
    }
    // Adding 1 to all clustering labels because R ploting functions does not take 0 as a label
    for (i = 0; i < n; i++){
        clus_vect[i] += 1;
    }
   
}



//OSil functions. The following 3 functions are for OSil. sil_lab, sil_lab_swap and clustyanlys

// [[Rcpp::export]]
void sil_lab(int K, int n, IntegerVector  clus_lab,  IntegerVector  clus_size,  NumericVector disty,  NumericVector dys_j,  NumericVector avg_dys_clus,  NumericVector sil){  //we will calculate silhouette's without medoids
    int i, j;
    double temp1,temp2;
    for(j = 0; j < n; j++){
        for(i = 0; i < K; i++){
            avg_dys_clus[(j*K)+i] = 0.;
            sil[j] = 0.;
            if (j == 0)  clus_size[i] = 0;
        }
    }
    
    for(i = 0; i < K; i++){
        for(j = 0; j < n; j++){
            if(clus_lab[j] == i){
                clus_size[i] += 1;
            }
        }
    }
    
    for(j=0;j<n;j++){
        for(i=(j+1);i<n;i++){
            temp1=disty[grab(i,j)];
            temp2=0.;
            if(clus_lab[i]==clus_lab[j]) temp2=-1;
            avg_dys_clus[(j*K)+clus_lab[i]]+=(temp1/(clus_size[clus_lab[i]]+temp2));
            avg_dys_clus[(i*K)+clus_lab[j]]+=(temp1/(clus_size[clus_lab[j]]+temp2));
        }
    }
    for(j=0;j<n;j++){
        if(clus_size[clus_lab[j]]==1){
            sil[j]=0.;
            continue;
        }
        for(i=0;i<K;i++){
            if(clus_lab[j]==i) dys_j[i]=0.;
            else dys_j[i]=avg_dys_clus[(j*K)+i];
        }
        hpsort(K,dys_j);
        if(dys_j[0]<0) dys_j[1]=dys_j[0];
        sil[j]=(dys_j[1]-avg_dys_clus[(j*K)+clus_lab[j]])/max(dys_j[1],avg_dys_clus[(j*K)+clus_lab[j]]);
    }
    
}

// [[Rcpp::export]]
void sil_lab_swap(int K, int n, IntegerVector  clus_lab, IntegerVector  alt_clus_lab, IntegerVector  clus_size,  NumericVector disty, IntegerVector  iter,  NumericVector dys_i,  NumericVector avg_dys_clus,  NumericVector silh,  NumericVector altsilh){
    int i, j, h;
    int bestlab = 0, bestobj = 0;
    int keep_swapping = 1;
    double obj_fun, best_obj_fun, alt_obj_fun;
    if (K==1) return;
    sil_lab(K, n, clus_lab,  clus_size, disty, dys_i, avg_dys_clus, silh); /*define silh instead of sil */
    obj_fun = 0.;
    for(i = 0; i < n; i++){
        obj_fun += silh[i];
    }
    obj_fun = obj_fun/n;
    best_obj_fun = obj_fun;
    alt_obj_fun = 0.;
    for(i = 0; i < n; i++){
        alt_clus_lab[i] = clus_lab[i];
    }
    iter[0] = 0;
    while(keep_swapping == 1){
        
        sil_lab(K,  n, clus_lab,  clus_size, disty, dys_i, avg_dys_clus, silh);
        obj_fun = 0;
        for (i = 0; i < n; i++){
            obj_fun += silh[i];
        }
        obj_fun = obj_fun/n;
        best_obj_fun = obj_fun;
        for(i = 0; i < n; i++){
            for(j = 0; j < K; j++){
                if(alt_clus_lab[i] == j) continue; /* skipping an object already present in present cluster*/
                alt_clus_lab[i] = j;
                sil_lab(K, n, alt_clus_lab, clus_size, disty, dys_i, avg_dys_clus, altsilh);
                alt_obj_fun = 0;
                for(h = 0; h < n; h++){
                    alt_obj_fun += altsilh[h];
                }
                alt_obj_fun = alt_obj_fun/n;
                if(alt_obj_fun > best_obj_fun){
                    best_obj_fun = alt_obj_fun;
                    bestlab = j;
                    bestobj = i;
                } /* end of if */
                alt_clus_lab[i] = clus_lab[i];
            } /* end of j */
            
        } /* end of i */
        
        if(best_obj_fun <= obj_fun){
            keep_swapping = 0;
            /*Rcpp::Rcout << "stop swapping " << std::endl;*/
        }
        else{
            clus_lab[bestobj] = bestlab;
            alt_clus_lab[bestobj] = bestlab;
            /*Rcpp::Rcout << "move point " << bestobj << std::endl;*/
        }
        
        /*Rcpp::Rcout << "clus_labs =" << clus_lab << std::endl;*/
        
        iter[0] += 1;
    } /*end of while */
    
    return;
    
} /*end of function sil_lab_swap */

// [[Rcpp::export]]
void clustyanlys(int K, int n, IntegerVector  clus_lab, IntegerVector  clus_size, NumericVector silh, NumericVector avg_clus_silh, NumericVector avg_clus_silhtwo, NumericVector avg_silh){
    int i, j;
    
    for(j = 0; j < K; j++){
        clus_size[j] = 0;
    }
    
    for(i = 0; i < K; i++){
        for(j = 0; j < n; j++){
            if(clus_lab[j] == i){
                clus_size[i] += 1;
            }
        }
    }
    
    
    for(j = 0; j < K; j++){
        for(i = 0; i < n; i++){
            if(clus_lab[i] == j){
                avg_clus_silhtwo[j] += silh[i];
            }
            avg_clus_silh[j] =  avg_clus_silhtwo[j]/clus_size[j];
        }
    }
    
    for(i = 0; i < n; i++){
        avg_silh[0] += silh[i];
    }
    avg_silh[0] = avg_silh[0]/n;
    
    
}


//The following two functions are for HOSil. sil_lab_link and hosil_lab_swap


// [[Rcpp::export]]
void sil_lab_link(int K, int n, IntegerVector clus_lab,  IntegerVector clus_size, NumericVector disty, NumericVector dys_j, NumericVector avg_dys_clus, NumericVector sil){  //we will calculate silhouette's without medoids
    int i, j;
    double temp1,temp2;
    for(j = 0; j < n; j++){
        for(i = 0; i < K; i++){
            avg_dys_clus[(j*K)+i] = 0.;
            sil[j] = 0.;
            if (j == 0)  clus_size[i] = 0;
        }
    }
    
    for(i = 0; i < K; i++){
        for(j = 0; j < n; j++){
            if(clus_lab[j] == i){
                clus_size[i] += 1;
            }
        }
    }
    
    
    
    for(j=0;j<n;j++){
        for(i=(j+1);i<n;i++){
            temp1=disty[grab(i,j)];
            temp2=0.;
            if(clus_lab[i]==clus_lab[j]) temp2=-1;
            avg_dys_clus[(j*K)+clus_lab[i]]+=(temp1/(clus_size[clus_lab[i]]+temp2));
            avg_dys_clus[(i*K)+clus_lab[j]]+=(temp1/(clus_size[clus_lab[j]]+temp2));
        }
    }
    
    for(j=0;j<n;j++){
        if(clus_size[clus_lab[j]]==1){
            sil[j]=0.;
            continue;
        }
        for(i=0;i<K;i++){
            if(clus_lab[j]==i) dys_j[i]=0.;
            else dys_j[i]=avg_dys_clus[(j*K)+i];
        }
        hpsort(K,dys_j);
        if(dys_j[0]<0) dys_j[1]=dys_j[0];
        sil[j]=(dys_j[1]-avg_dys_clus[(j*K)+clus_lab[j]])/max(dys_j[1],avg_dys_clus[(j*K)+clus_lab[j]]);
    }
    
    
}

// [[Rcpp::export]]
void hosil_lab_swap(int n, NumericVector disty, IntegerVector clus_lab, IntegerVector clus_size,  NumericVector dys_j, NumericVector avg_dys_clus, NumericVector sil, IntegerVector best_clus_lab, IntegerVector all_best_clus_lab, NumericVector all_best_avg_silh, IntegerVector copy_best_clus_lab){
    
    int K = n-1;
    double avg_silh, best_avg_silh;
    avg_silh = 0.;
    best_avg_silh = -2.;
    if (K==1) return;
    
    double smdist = DBL_MAX;
    
    int n_disty = n*(n-1)/2+1;
    int smloc = 0;
    int i, j;
    
    for(i = 1; i < n_disty; i++){ //ignoring first element which is always zero we purpossly added.Note it will not effect any other distance which is zero naturally
        if(disty[i] < smdist){
            smdist = disty[i];
            smloc = i;
        }
    }
    
    int get = 0;
    int obs1 = 0;
    int obs2 = 0;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            get = grab(i, j);
            if(get == smloc){
                obs1 = i;
                obs2 = j;
            }
        }
    }
    
    
    
    for( i = 0; i < n; i++){
        clus_lab[i] = i;
    }
    
    
    
    
    j = 0;
    
    clus_lab[obs1] = j;
    clus_lab[obs2] = j;
    for( i = 0; i < n; i++){
        if ( i == obs1 || i == obs2) continue;
        clus_lab[i] = j+1;
        j++;
    }
    
    
    
    
    // Call sil_lab_link here. This will be all_best_avg_silh[0] for first hiererachy
    sil_lab_link(K,  n, clus_lab,  clus_size, disty, dys_j, avg_dys_clus, sil);
    for(int gg = 0; gg<n; gg++){
        avg_silh += sil[gg];
    }
    avg_silh = avg_silh/n;  //here we dont have to compare it with best. as we comb two nearest and it the best we can get
    
    all_best_avg_silh[0] =  avg_silh;
    for(int gg = 0; gg<n; gg++){
        all_best_clus_lab[gg] = clus_lab[gg]; //next time you have to preserve first n values.
    }
    
    
    for(int i = 0; i < n; i++){
        best_clus_lab[i] = clus_lab[i];
    }
    
    for(int i = 0; i < n; i++){
        copy_best_clus_lab[i] = clus_lab[i];
    }
    
    
    
    //Now code for the remaing k = n-2 hierarchies.
    
    int k, r, l;
    r = 1;
    k = n-r;
    
    
    for(l = 0; l<n-3; l++){ //generating labels for all hierarchies
        best_avg_silh = -2.;
        for(i = 0; i < k+1; i++) //Cluster 1 to combine
        {
            for(j = i+1; j < k; j++){ //Cluster 2 to combine
                
                // Now producing complete set of labels for each pair of cluster
                for(int g = 0; g < n; g++){
                    if ((clus_lab[g] == i) || (clus_lab[g] == j))
                    {
                        copy_best_clus_lab[g] = 0;
                    }
                    
                }
                //This is final piece of code
                int w = 1;
                int s;
                for(s = 0; s < k; s++){
                    if(s == i || s == j) continue;
                    for(int jj = 0; jj < n; jj++){
                        if(clus_lab[jj] == s)
                        copy_best_clus_lab[jj] = w;
                    }
                    w += 1;
                }
                
                
                sil_lab_link(k-1,  n, copy_best_clus_lab,  clus_size, disty, dys_j, avg_dys_clus, sil);
                avg_silh = 0.;
                for(int gg = 0; gg<n; gg++){
                    avg_silh += sil[gg];
                }
                avg_silh = avg_silh/n;
                
                if(best_avg_silh < avg_silh)
                {
                    best_avg_silh = avg_silh;
                    for(int tt=0; tt<n; tt++)
                    best_clus_lab[tt] = copy_best_clus_lab[tt];
                }
                
            } // end of j
        } // end of i
        
        
        for(int yy = 0; yy<n; yy++){
            clus_lab[yy] = best_clus_lab[yy];
        }
        
        all_best_avg_silh[r] =  best_avg_silh;
        for(int gg = 0; gg<n; gg++){
            all_best_clus_lab[gg+r*n] = best_clus_lab[gg]; //next time you have to preserve first n values.
        }
        r+=1;
        k = n-r;
    } // end of l
} //end of function hosil_lab_swap


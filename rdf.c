#include <stdio.h>
#include <stdlib.h>

int numOfLines = 0;
double xMax = 0;

void rdf(double **pCoo){
    int i = 0;
    int bins = 100;
    double binWidth = (xMax/2)/bins;
    double dist = 0;
    int *histo;

    histo = malloc(bins * sizeof(int));

    for(i = 1; i < numOfLines; i++){
        dist = distance(pCoo[0], pCoo[i]);
        dist = sqrt(dist);
        histo[(int)(dist/binWidth)]++;
    }

    FILE *f = fopen("histo.txt", "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    for(i = 0; i < numOfLines; i++){
        fprintf(f, "%d\n", histo[i]);
    }
    fclose(f);

    free(histo);
}

int main(){
    int c;
    FILE *file;
    file = fopen("output.txt", "r");
    if (file) {
        while ((c = getc(file)) != EOF)
            putchar(c);
        fclose(file);
    }
}
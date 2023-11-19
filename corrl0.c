// with a little hep from chatgpt
// for the formulas
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Comparison function for qsort
int compare(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

double calcMean(int arr[], int n) {
    int i, sum = 0;
    for(i = 0; i < n; i++)
        sum += arr[i];
    return (double)sum / n;
}

double calcPearCorr(int arr1[], int arr2[], int n)
{
    int i;
    double mean1 = calcMean(arr1, n);
    double mean2 = calcMean(arr2, n);

    double ntor = 0.0, dtor1 = 0.0, dtor2 = 0.0;

    for (int i = 0; i < n; i++) {
        ntor += (arr1[i] - mean1) * (arr2[i] - mean2);
        dtor1 += pow(arr1[i] - mean1, 2);
        dtor2 += pow(arr2[i] - mean2, 2);
    }
    return ntor / sqrt(dtor1 * dtor2);
}

double calcSpearmanCorr(int arr1[], int arr2[], int n)
{
    // Create arrays to store ranks
    int *rank1 = malloc(n * sizeof(int));
    int *rank2 = malloc(n * sizeof(int));

    // Assign ranks to the arrays
    for (int i = 0; i < n; i++) {
        rank1[i] = i + 1;
        rank2[i] = i + 1;
    }

    // Sort the arrays based on the values in arr1 and arr2
    qsort(rank1, n, sizeof(int), compare);
    qsort(rank2, n, sizeof(int), compare);

    // Calculate the differences between ranks
    double sum_squared_diff = 0.0;
    for (int i = 0; i < n; i++) {
        int diff = rank1[i] - rank2[i];
        sum_squared_diff += diff * diff;
    }

    // Calculate Spearman correlation coefficient
    double correlation = 1.0 - (6.0 * sum_squared_diff) / (n * (n * n - 1));

    // Free allocated memory
    free(rank1);
    free(rank2);

    return correlation;
}

double calcKendallCorr(int arr1[], int arr2[], int n)
{
    int concordant = 0, discordant = 0;

    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if ((arr1[i] < arr1[j] && arr2[i] < arr2[j]) || (arr1[i] > arr1[j] && arr2[i] > arr2[j])) {
                concordant++;
            } else if ((arr1[i] < arr1[j] && arr2[i] > arr2[j]) || (arr1[i] > arr1[j] && arr2[i] < arr2[j])) {
                discordant++;
            }
        }
    }
    double correlation = (concordant - discordant) / sqrt((concordant + discordant) * (n * (n - 1) / 2));
    return correlation;
}

int main(int argc, char *argv[])
{
   /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
	if(argc!=2) {
		printf("Error. Pls supply argument (name of file).\n");
		exit(EXIT_FAILURE);
	}

    int seq1[] = {1, 2, 3, 4, 5};
    int seq2[] = {5, 4, 3, 2, 1};

    // Calculate correlation
    int n = sizeof(seq1) / sizeof(seq1[0]);
    double correl = calcCorr(seq1, seq2, n);

    // Display the result
    printf("Corr: %.4lf\n", correl);

    return 0;
}

   /* declarations */
   FILE 
   char
   int
   float
   double

   /* print stuff out */
   printf("Res_1: %s\n", "somedata");

   return 0;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "read_data.h"

#define UNUSED(expr) do { (void)(expr); } while (0) // To silence a -Wunused-but-set-variable warning

void read_data(const char *filename, int nb_lines, int nb_cols, double **res)
{
    FILE *fp;
    int i;

    // Allocate
    (*res) = (double*)malloc(nb_lines * nb_cols * sizeof(double));
    if ((*res) == NULL)
    {
        return;
    }

    // Open file
    fp = fopen(filename, "r");

    // Read data from file
    for (i=0 ; i < nb_lines ; i++)
    {
        char *token;
        int j;
        char * end;

        ssize_t read;
        char * line = NULL;
        size_t len = 0;
        double inter;

        read = getline(&line, &len, fp);
        UNUSED(read);

        /* get the first token */
        token = strtok(line, ",");

        j = 0;
        while(token != NULL)
        {
            inter = strtod(token, &end);
            (*res)[i*nb_cols + j] = inter;
            j++;
            token = strtok(NULL, ",\n");
        }

    }
    
    // Close file
    fclose(fp);

    printf("Read data, done !\n");

    return;
}

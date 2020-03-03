/************************************************************************************************************************************************
* LOOKAHEAD HEURISTIC FOR THE MINIMIZATION OF OPEN STACKS PROBLEM                                                                              *
/************************************************************************************************************************************************
* AUTHORS: MARCO ANTONIO MOREIRA DE CARVALHO (marco.opt@gmail.com) AND NEI YOSHIHIRO SOMA (soma@ita.br)                                         *
* AERONAUTICS TECHNOLOGICAL INSTITUTE - BRAZIL                                                                                                  *
* DATE OF CREATION: 13/02/2012                                                                                                                  *
* LAST MODIFICATION: 14/02/2012                                                                                                                 *
/************************************************************************************************************************************************
* TERMS OF USE                                                                                                                                  *
*                                                                                                                                               *
* This code available for free public use.                                                                                                      *
* The code in this archive is licensed gratis to all third parties under the terms of this paragraph.                                           *
* Copying and distribution of the files in this archive is unrestricted if and only if the files are not modified.                              *
* Modification of the files is encouraged, but the distribution of modifications of the files in this archive is unrestricted only              *
* if you meet the following conditions: modified files must carry a prominent notice stating (i) the author and date,                           *
* (ii) the new author and the date of release of the modification, (iii) that the work is licensed at no charge to all parties.                 *
*                                                                                                                                               *
* If you use the code extensively in your research, you are requested to provide appropriate attribution and thanks to the author of the code.  *
/************************************************************************************************************************************************
* DISCLAIMER OF WARRANTY                                                                                                                        *
*                                                                                                                                               *
* This source code is provided "as is" and without warranties as to performance or merchantability.                                             *
* The authors and/or distributors of this source code may have made statements about this source code.                                          *
* Any such statements do not constitute warranties and shall not be relied on by the user in deciding whether to use this source code.          *
*                                                                                                                                               *
* This source code is provided without any express or implied warranties whatsoever. Because of the diversity of conditions and hardware under  *
*  which this source code may be used, no warranty of fitness for a particular purpose is offered.                                                          *
* The user is advised to test the source code thoroughly before relying on it. The user must assume the entire risk of using the source code.   *
*                                                                                                                                               *
************************************************************************************************************************************************/

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "LinkedList.h"

typedef struct  {
    int degree;     //stores the degree of the vertex
    int id;         //stores the index of the vertex
    int removed;    //indicates whether the vertex has been removed in MDGh procedure
    int blocked;    //indicates whether the vertex is a neighbour of a removed vertex in MDGh procedure
} typeVertex;       //structure of a vertex

typePattern* pattern;   //stores informations about the patterns, defined in LinkedList.h file
typeVertex* vertex;     //stores informations about the vertex
typeList* dominator;    //stores the dominated patterns;

int** matrix;           //stores the input data
int** adjacencyMatrix;  //stores the adjancency matrix of the MOSP graph, used as an upper diagonal matrix
int* stack;             //simulates the open stacks
int* demand;            //stores the demands for pieces (original)
int* demand2;           //stores the demands for pieces
int* demand3;           //stores the demands for pieces
int* solution;          //stores the pattern sequencing
int* solution2;         //stores the alternative pattern sequencing
int* spe;               //stores the pieces sequencing
int* degrees;           //stores the vertices degrees on the complementary MOSP graph
int* Q;
int nPatterns;          //number of patterns
int nPieces;            //number of pieces
int nPatternsOriginal;  //original number of patterns
int nOpenStacks;        //number of open stacks
int lookahead;

/*
Initializes the structures, vectors and matrices used
*/

void initialization()
{
    int i;

    matrix = (int**) calloc(nPatterns, sizeof(int*));

    if (!matrix) {                                                  //tests the memory allocation
        printf("INSUFFICIENT MEMORY (A)\n");                        //not enough memory to allocate the pointer
        system("pause");
    }

    for (i = 0; i < nPatterns; i++) {
        matrix[i] = (int*) calloc(nPieces, sizeof(int));

        if (!matrix[i]) {
            printf("INSUFFICIENT MEMORY (B)\n");
            system("pause");
        }
    }

    pattern = (typePattern*) calloc(nPatterns, sizeof(typePattern));

    if (!pattern) {
        printf("INSUFFICIENT MEMORY (C)\n");
        system("pause");
    }

    stack = (int*) calloc(nPieces, sizeof(int));

    if (!stack) {
        printf("INSUFFICIENT MEMORY (D)\n");
        system("pause");
    }

    demand = (int*) calloc(nPieces, sizeof(int));

    if (!demand) {
        printf("INSUFFICIENT MEMORY (E)\n");
        system("pause");
    }

    solution = (int*) calloc(nPatterns, sizeof(int));

    if (!solution) {
        printf("INSUFFICIENT MEMORY (F)\n");
        system("pause");
    }

    solution2 = (int*) calloc(nPatterns, sizeof(int));

    if (!solution2) {
        printf("INSUFFICIENT MEMORY (F)\n");
        system("pause");
    }

    dominator = (typeList*) calloc(nPatterns, sizeof(typeList));

    if (!dominator) {
        printf("INSUFFICIENT MEMORY (G)\n");
        system("pause");
    }

    vertex = (typeVertex*) calloc(nPieces, sizeof(typeVertex));

    if (!vertex) {
        printf("INSUFFICIENT MEMORY (H)\n");
        system("pause");
    }

    adjacencyMatrix = (int**) calloc(nPieces, sizeof(int*));

    if (!adjacencyMatrix) {
        printf("INSUFFICIENT MEMORY (I)\n");
        system("pause");
    }

    for (i = 0; i < nPieces; i++) {
        adjacencyMatrix[i] = (int*) calloc(nPieces, sizeof(int));

        if (!adjacencyMatrix[i]) {
            printf("INSUFFICIENT MEMORY (J)\n");
            system("pause");
        }
    }

    spe = (int*) calloc(nPieces + 1, sizeof(int));

    if (!spe) {
        printf("INSUFFICIENT MEMORY (K)\n");
        system("pause");
    }

    degrees = (int*) calloc(nPieces, sizeof(int));

    if (!degrees) {
        printf("INSUFFICIENT MEMORY (L)\n");
        system("pause");
    }

    Q = (int*) calloc(nPieces * nPieces, sizeof(int));          //the queue for the breadth first search

    if (!Q) {
        printf("INSUFFICIENT MEMORY (M)\n");
        system("pause");
    }

    demand2 = (int*) calloc(nPieces, sizeof(int));

    if (!demand2) {
        printf("INSUFFICIENT MEMORY (N)\n");
        system("pause");
    }

    demand3 = (int*) calloc(nPieces, sizeof(int));

    if (!demand3) {
        printf("INSUFFICIENT MEMORY (O)\n");
        system("pause");
    }

    for (i = 0; i < nPatterns; i++) {                           //initializes some structures
        pattern[i].removed = 0;
        pattern[i].size = 0;
        pattern[i].blocked = 0;
        emptyList(&dominator[i]);
    }

    for (i = 0; i < nPieces; i++) {
        vertex[i].degree = 0;
        vertex[i].id = i;
        Q[i] = -1;                                              //initialises the positions of the queue
    }

    for (i = nPieces; i < nPieces * nPieces; i++)
        Q[i] = -1;
}

/*
Reads the problem from a file specified by fileName
*/
void readProblem(char* fileName)
{
    int i;
    int j;
    int k;

    char nome[256];
    FILE* fpIn = fopen(fileName, "r");                      //input file
    fscanf(fpIn, "%d %d", &nPieces, &nPatterns);
    printf("%d %d", nPatterns, nPieces);
    nPatternsOriginal = nPatterns;

    initialization();                                       //initializes all structures, vectors and matrices

    //for (j = nPieces - 1; j >= 0; j--) {
    for (i = 0; i < nPatterns; i++) {
        for (j = 0; j < nPieces; j++) {
            fscanf(fpIn, "%d", &matrix[i][j]);

            if (matrix[i][j] != 0) {
                matrix[i][j] = 1;
                pattern[i].size++;
            }

        }
    }

    fclose(fpIn);
}

/*
Builds de MOSP graph
*/

void buildMOSPGraph()
{
    int i, j, k;


    for (j = nPieces - 1; j >= 0; j--) {
        for (i = 0; i < nPatternsOriginal; i++) {
            if ((!pattern[i].removed) && (matrix[i][j] == 1)) {
                for (k = 0; k < nPieces; k++) {
                    if ((j != k) && (matrix[i][k] == 1)) {  //updates the adjacency matrix
                        if (k < j) {
                            adjacencyMatrix[k][j] = 1;
                        }
                        else {
                            adjacencyMatrix[j][k] = 1;
                        }
                    }
                }

                demand[j]++;
                demand2[j]++;
            }
        }
    }

    //Determines the degree of each vertex on the MOSP graph
    for (i = 0; i < nPieces; i++) {
        for (j = 0; j < nPieces; j++) {
            if (i < j) {
                vertex[i].degree += adjacencyMatrix[i][j];
            }
            else if (i > j) {
                vertex[i].degree += adjacencyMatrix[j][i];
            }
        }
    }

    //stores the degree of the vertices
    for (i = 0; i < nPieces; i++)
        degrees[i] = vertex[i].degree;
}

/*
Evaluates the current solution, simulating the opened and closed stacks
*/
int evaluation(int limit, int* solution)
{
    int i;
    int j;
    int index;
    int openStacks;
    int closedStacks;
    int maxOpenStacks;

    openStacks = 0;
    closedStacks = 0;
    maxOpenStacks = -1;

    for (i = 0; i < nPieces; i++) {
        demand3[i] = demand2[i];
        stack[i] = 0;
    }

    for (i = 0; i < limit; i++) {
        index = solution[i];                //simulates the sequencing of the current pattern

        for (j = 0; j < nPieces; j++) {
            if (matrix[index][j] > 0) {
                demand3[j]--;               //decreases the demand for each piece of the pattern

                if (demand3[j] == 0)
                    closedStacks++;         //simulates closed stacks

                if (stack[j] == 0) {
                    openStacks++;           //simulates new open stacks
                    stack[j] = 1;
                }
            }
        }

        if (openStacks > maxOpenStacks)
            maxOpenStacks = openStacks;     //stores the maximum number of open stacks

        openStacks -= closedStacks;
        closedStacks = 0;
    }

    return maxOpenStacks;                   //returns the maximum number of open stacks
}

/*
Gives the pattern sequencing correspondent to the piece sequencing, original version
*/
void patternSequencing()
{
    int i;
    int j;

    int index;

    index = nPatterns - 1;

    for (i = nPieces - 1; index >= 0; i--) {                    //traverses the pieces sequence in reversal order
        for (j = 0; (j < nPatternsOriginal) && (index >= 0); j++) { //finds the pattern that include the current piece
            if ((!pattern[j].blocked) && (!pattern[j].removed)) { //if the pattern wasn't already sequenced and wasn't removed by the preprocessing
                if (matrix[j][spe[i]] == 1) {
                    pattern[j].blocked = 1;
                    solution[index--] = j;                      //inserts it in the pattern sequence from the end to the beginning
                }
            }
        }
    }
}

/*
 Gives the pattern sequencing correspondent to the piece sequencing, as in Becceneri (2004)
 */
void patternSequencing2()
{
    int i;
    int j;
    int* completion;
    int index = 0;

    completion = (int*) calloc(nPatternsOriginal, sizeof(int));

    if (!completion) {
        printf("INSUFFICIENT MEMORY (XYZ)\n");
        system("pause");
    }

    for (i = 0; i < nPatternsOriginal; i++) {
        for (j = 0; j < nPieces; j++)
            if ((!pattern[i].removed) && (matrix[i][j] == 1))
                completion[i]++;
    }

    for (i = nPieces - 1; i >= 0; i--) { //trocar do início para o final
        for (j = 0; j < nPatternsOriginal; j++) {
            if ((!pattern[j].blocked) && (!pattern[j].removed)) {
                if (matrix[j][spe[i]] == 1) {
                    completion[j]--;

                    if (completion[j] == 0) {
                        pattern[j].blocked = 1;
                        solution[index++] = j;   //inserts it in the pattern sequence from the left to the right
                    }
                }
            }
        }
    }

    for (i = 0; i < nPatternsOriginal; i++)
        pattern[i].blocked = 0;

    free(completion);
}

int findMinimumDegree()
{
    int index  = -1;
    int i;
    int smaller = nPieces * nPieces + 1;

    for (i = 0; i < nPieces; i++) {
        if ((vertex[i].degree < smaller) && (!vertex[i].removed) && (!vertex[i].blocked)) {
            smaller = vertex[i].degree;                                                 //selects the vertex of smaller degree not sequenced yet
            index = i;                                                                  //stores its index
        }
    }

    return index;
}

void breadthFirstSearch()
{
    int i;
    typeList preQ;
    int preQindex = 0;
    int index;
    int index2 = 0;
    int indexQ = 0;
    int indexSearch = 0;
    int nVertices;

    emptyList(&preQ);
    nVertices = nPieces;

    do {
        if (preQindex > 0) {
            transfer(&preQ, Q, &indexQ);

            emptyList(&preQ);
            preQindex = 0;
        }

        count(&preQ);

        if (Q[indexSearch] == -1) {
            Q[indexQ] = findMinimumDegree();                                        //finds the vertex with maximum degree and enqueue it
            indexQ++;
        }

        spe[index2] = Q[indexSearch];                                               //dequeue and inserts into the pieces sequence
        index2++;
        index = Q[indexSearch];                                                     //stores the piece under analysis
        indexSearch++;

        for (i = 0; i < index && 0 < vertex[index].degree; i++) { //searches as long as the degree of the vertex allow
            if (adjacencyMatrix[i][index] == 1) {                                   //searches for its neighbours
                if ((!vertex[i].blocked) && (!vertex[i].removed)) {                 //if they are not in queue or pieces sequence
                    vertex[i].blocked = 1;                                          //blocks the neighbour

                    if (count(&preQ) > 0)
                        insereOrdenado(&preQ, i, degrees);
                    else
                        add(&preQ, i);
                    preQindex++;
                }
            }
        }

        for (i = index + 1; i < nPieces && 0 < vertex[index].degree; i++) { //same as the loop above, but the indexes of the adjacency matrix are reversed (remember it is an upper diagonal matrix)
            if (adjacencyMatrix[index][i] == 1) {
                if ((!vertex[i].blocked) && (!vertex[i].removed)) {
                    vertex[i].blocked = 1;

                    if (count(&preQ) > 0)
                        insereOrdenado(&preQ, i, degrees);
                    else
                        add(&preQ, i);
                    preQindex++;
                }
            }
        }

        vertex[index].removed = 1;                                                  //removes the dequeued piece from the problem

    }
    while (indexSearch < nPieces);                                                  //until all vertices are in the pieces sequence

    //  patternSequencing2();                                                           //obtains the pattern sequence related to the pieces sequence obtained
    finalize(&preQ);
}

/*
 Pre-processing procedure, detects and removes dominated patterns
 */
void dominationPreProcessing()
{
    int i;
    int j;
    int k;
    int flag = 1;
    int index;
    int index2;
    int counter = 0;

    for (i = 0; i < nPatterns; i++) {
        if (!pattern[i].removed) {
            for (j = i + 1; j < nPatterns; j++) {
                if (!pattern[j].removed) {
                    if (pattern[i].size < pattern[j].size) {                        //if pattern i is smaller than pattern j
                        index = i;
                        index2 = j;
                        for (k = 0; k < nPieces; k++) {
                            if ((matrix[i][k] == 1) && (matrix[j][k] == 0)) {       //and pattern i has a piece that pattern j doesn't
                                flag = 0;                                           //then i is not dominated by j
                                break;
                            }
                        }
                    }
                    else    if (pattern[i].size == pattern[j].size) {               //if patterns i and j have the same size
                        index = j;
                        index2 = i;
                        for (k = 0; k < nPieces; k++) {
                            if (matrix[i][k] != matrix[j][k]) {             //and they differ in any position
                                flag = 0;                                   //then no one is dominated
                                break;
                            }
                        }
                    }
                    else    if ((pattern[i].size > pattern[j].size)) {      //if pattern i is larger than pattern j
                        index = j;
                        index2 = i;
                        for (k = 0; k < nPieces; k++) {
                            if ((matrix[i][k] == 0) && (matrix[j][k] == 1)) { //and j has a piece that i doesn't
                                flag = 0;                                   //then j is not dominated by i
                                break;
                            }
                        }
                    }

                    if (flag == 1) {                                                //if the a dominance condition is met
                        pattern[index].removed = 1;                                 //the pattern is removed
                        counter++;

                        add(&dominator[index2], index);

                        if (count(&dominator[index]) > 0)                           //if the dominated pattern dominates other patterns
                            transferDominated(&dominator[index2], &dominator[index]);       //the list is transferred
                        /*
                        for(k=0; k<nPieces; k++)                                    //for each piece of the removed pattern
                            if(matrix[index][k] == 1)
                            {
                                demand2[k]--;                                       //the demand is decreased
                            }
                        */
                        if (index == i)
                            break;
                    }
                }

                flag = 1;
            }
        }
    }

    nPatternsOriginal = nPatterns;                                                  //stores the original number of patterns
    nPatterns -= counter;                                                           //updates the number of patterns

    //if (nPatternsOriginal == nPatterns)
      //  printf("instância %dx%dsem dominância\n", nPatterns, nPieces);
}

int improvementRule(int bestValue)
{
    int i;
    int j;
    int k;
    int l;
    int position1;
    int position2;
    int result;
    int pattern1;
    int pattern2;

    for (j = 0; j < nPieces; j++)
        demand[j] = demand2[j];

    for (i = 0; i < nPatterns - 1; i++) {
        for (k = 0; k < nPieces; k++) {
            if (matrix[solution[i]][k] > 0)
                demand[k]--;                                //simulates the decrease on the demand of each piece after the sequencing of each pattern of the solution
        }

        for (k = 0; k < nPieces; k++) {
            if ((demand[k] == lookahead) && (matrix[solution[i + 1]][k] != 1)) { //if a stack is about to close
                position1 = i + 1;                  //determines the position for reallocation of the pattern
                for (l = i + 2; l < nPatterns; l++) {       //searches for the pattern that will close that stack
                    if (matrix[solution[l]][k] == 1) {      //when found
                        position2 = l;                      //determines the original position of the pattern

                        pattern1 = solution[position1];
                        pattern2 = solution[position2];

                        for (j = position2; j > position1; j--) //reallocates the pattern
                            solution[j] = solution[j - 1];

                        solution[position1] = pattern2;

                        result = evaluation(nPatterns, solution);   //evaluates the modified solution

                        if (result <= bestValue) {                  //if there's no prejudice to the solution
                            bestValue = result;                     //stores the value found
                        }
                        else {                                          //if the solution was prejudiced
                            demand[k] = 0;

                            for (j = position1; j < position2; j++)     //restores the solution
                                solution[j] = solution[j + 1];

                            solution[position2] = pattern2;
                        }
                        position1++;
                    }
                }
            }
        }
    }

    return bestValue;
}

int improvementRule2(int bestValue)
{
    int i;
    int j;
    int k;
    int l;
    int position1;
    int position2;
    int result;
    int pattern1;
    int pattern2;

    for (j = 0; j < nPieces; j++)
        demand[j] = demand2[j];

    for (i = nPatterns - 1; i > 0; i--) {
        for (k = 0; k < nPieces; k++) {
            if (matrix[solution[i]][k] > 0)
                demand[k]--;                                //simulates the decrease on the demand of each piece after the sequencing of each pattern of the solution
        }

        for (k = 0; k < nPieces; k++) {
            if ((demand[k] == lookahead) && (matrix[solution[i - 1]][k] != 1)) { //if a stack is about to close
                position1 = i - 1;                  //determines the position for reallocation of the pattern
                for (l = i - 2; l >= 0; l--) {      //searches for the pattern that will close that stack
                    if (matrix[solution[l]][k] == 1) {      //when found
                        position2 = l;                      //determines the original position of the pattern

                        pattern1 = solution[position1];
                        pattern2 = solution[position2];

                        for (j = position2; j < position1; j++) //reallocates the pattern
                            solution[j] = solution[j + 1];

                        solution[position1] = pattern2;

                        result = evaluation(nPatterns, solution);   //evaluates the modified solution

                        if (result <= bestValue) {                  //if there's no prejudice to the solution
                            bestValue = result;                     //stores the value found
                        }
                        else {                                          //if the solution was prejudiced
                            demand[k] = 0;

                            for (j = position1; j > position2; j--)     //restores the solution
                                solution[j] = solution[j - 1];

                            solution[position2] = pattern2;

                            result = evaluation(nPatterns, solution);
                        }
                        position1--;
                    }
                }
            }
        }
    }

    return bestValue;
}

/*
Terminates all lists, vectors and matrices. Frees all pointers
*/
void termination()
{
    int i;

    for (i = 0; i < nPatternsOriginal; i++) {
        finalize(&dominator[i]);
        free(matrix[i]);
    }

    for (i = 0; i < nPieces; i++)
        free(adjacencyMatrix[i]);

    free(dominator);
    free(matrix);
    free(adjacencyMatrix);
    free(pattern);
    free(stack);
    free(demand);
    free(demand2);
    free(demand3);
    free(solution);
    free(solution2);
    free(vertex);
    free(spe);
    free(degrees);
    free(Q);
}

int max(int a, int b)
{
    return a > b ? a : b;
}

int min(int a, int b)
{
    return a < b ? a : b;
}
/*
Main procedure of the constructive heuristic
Parameters:
bfsResult   stores the solution value obtained by the heuristic
bfsTime     stores the exeution time of the heuristic
*/
void mainMethod(int* bfsResult, double* bfsTime, char* inputFileName)
{
    int i;
    int j;
    int start = 0;
    int end;
    int result;
    int result2;
    int resultb;
    int X;

    double startBFS;
    double endBFS;

    char outputFileName[256];


    FILE* fpSolution;

    sprintf(outputFileName, "Solution_%s", inputFileName);

    fpSolution = fopen(outputFileName, "w");                //file that contains the information about the solution of a problem instance

    readProblem(inputFileName);                             //reads the problem data
    startBFS = clock();                                     //time taking

    dominationPreProcessing();                               //preprocessing procedure, removes dominated patterns
    buildMOSPGraph();

    breadthFirstSearch();
    patternSequencing();
    resultb = evaluation(nPatterns, solution);

    X = 0.1 * max(nPatterns, 100);
    X = min(nPatterns, X);
    lookahead = 1;

    for (i = 0; i < X; i++) {
        lookahead = 1;

        for (j = 0; j < X; j++) {
            result2 = improvementRule(resultb);

            if (result2 < resultb) {
                resultb = result2;
                end = 0;
            }

            result2 = improvementRule2(resultb);

            if (result2 < resultb) {
                resultb = result2;
                end = 0;
            }

            lookahead++;
        }
    }

    for (i = 0; i < nPatterns; i++)
        solution2[i] = solution[i];

    //---------------------------------

    patternSequencing();
    result = evaluation(nPatterns, solution);
    lookahead = 1;

    for (i = 0; i < X; i++) {
        lookahead = 1;

        for (j = 0; j < X; j++) {
            result2 = improvementRule(result);

            if (result2 < result) {
                result = result2;
                end = 0;
            }

            result2 = improvementRule2(result);

            if (result2 < result) {
                result = result2;
                end = 0;
            }

            lookahead++;
        }
    }

    if (result > resultb) {
        result = resultb;

        for (i = 0; i < nPatterns; i++)
            solution[i] = solution2[i];
    }

    endBFS = clock();                                 //time taking

    *bfsResult = result;                              //stores the solution value
    *bfsTime = (endBFS - startBFS) / CLOCKS_PER_SEC;        //stores the execution time

    fprintf(fpSolution, "Maximum Open Stacks: %d\n", result);
    fprintf(fpSolution, "Number of Patterns: %d, Number of Pieces: %d, Patterns eliminated in the Pre-processing: %d\n", nPatternsOriginal, nPieces, nPatternsOriginal - nPatterns);
    fprintf(fpSolution, "Patterns Sequence:\n");

    for (i = 0; i < nPatterns; i++) {
        fprintf(fpSolution, "%d \n", solution[i] + 1);          //prints the indexes (from 0 to nPatterns-1) of the patterns, one per line

        if (count(&dominator[solution[i]]) > 0) {              //prints the dominated patterns right after the patterns they are dominated by
            print(&dominator[solution[i]], fpSolution);     //after the index of the pattern the word [dominated] is written to indicate it
        }
    }

    fprintf(fpSolution, "\nSolution Matrix:\n");

    for (i = 0; i < nPatterns; i++) {
        for (j = 0; j < nPieces; j++) {
            fprintf(fpSolution, "%d ", matrix[solution[i]][j]);
        }

        fprintf(fpSolution, "\n");

        if (count(&dominator[solution[i]]) > 0) {                               //prints the dominated patterns right after the patterns they are dominated by
            print2(&dominator[solution[i]], fpSolution, nPieces, matrix);       //after the index of the pattern the word [dominated] is written to indicate it
        }
    }

    fclose(fpSolution);

    termination();
}

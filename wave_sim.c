/**
 *******************************************************************************
 * @file:    wave_sim.c
 * @author:  Danny Soppit
 * @brief:   Simulates the physics of the wave equation through a terminal.
 *           Ensure your terminal window is wide enough!!
 *******************************************************************************
 */

//******************************************************************************
//  Include Files
//******************************************************************************

// STANDARD DEFINITONS
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

//******************************************************************************
//  Defines
//******************************************************************************

// Color and terminal defines
#define COLOR_BUFFER_SIZE 10
#define RESET             "\x1b[0m"
#define CURSOR            "\033[H"
#define RED               "\033[0;31m"
#define GREEN             "\033[0;32m"
#define YELLOW            "\033[0;33m"
#define BLUE              "\033[0;34m"
#define MAGNETA           "\033[0;35m"
#define CYAN              "\033[0;36m"
#define WHITE             "\033[0;37"
#define BLACK             "\033[0;30m"

// Mesh Parameters
#define Lx 10e-6                 // Length in x direction
#define Ly 10e-6                 // Length in y direction
#define dx 0.12e-6               // Grid size in x direction
#define dy 0.12e-6               // Grid size in y direction
#define Nx ((int)(Lx / dx) + 1)  // Number of nodes in x-direction
#define Ny ((int)(Ly / dy) + 1)  // Number of nodes in y-direction
#define n_stop 150               // Number of time steps

// Physical Constants
#define l 1.0e-6                 // Wavelength
#define w 18.0e-15               // Width of the pulse
#define T0 4.0e-15               // Initial time
#define c 299792458              // Speed of light

// Time Step Calculation
#define dt (1.0 / (c * sqrt((1.0 / (dx * dx)) + (1.0 / (dy * dy)))))

// Courant Numbers
#define Ox ((c * dt) / dx)       // Courant number in x-direction
#define Oy ((c * dt) / dy)       // Courant number in y-direction

// Source Location
#define xs1 50
#define ys1 50

//******************************************************************************
//  Functions
//******************************************************************************

/**
 *******************************************************************************
 * @brief:     Dynamic allocation of a 2D array
 * @parameter: rows: Number of rows
 * @parameter: cols: Number of cols
 * @return:    2D array pointer
 *******************************************************************************
 */
double** allocate2DArray(int rows, int cols)
{
    double** array = (double**) malloc(rows * sizeof(double*));

    for (int i = 0; i < rows; i++)
    {
        array[i] = (double*) malloc(cols * sizeof(double));
    }

    return array;
}

/**
 *******************************************************************************
 * @brief:     Free the memory of a 2D array
 * @parameter: array: 2D array pointer reference
 * @parameter: rows: Number of rows
 * @return:    N/A
 *******************************************************************************
 */
void free2DArray(double** array, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        free(array[i]);
    }

    free(array);
}

/**
 *******************************************************************************
 * @brief:     Initialize the 2D array to zeros
 * @parameter: array: 2D array pointer reference
 * @parameter: rows: Number of rows
 * @paramter:  cols: Number of columns
 * @return:    N/A
 *******************************************************************************
 */
void initializeArray(double** array, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            array[i][j] = 0.0;
        }
    }
}

/**
 *******************************************************************************
 * @brief:     Obtain the color value of a node
 * @parameter: value: Value of the node at position x,y,
 * @parameter: color: Pointer to the color buffer
 * @return:    N/A
 *******************************************************************************
 */
void getColor(double value, char* color)
{
    // Change threshold as desired
    const double blackThreshold = 0.05;

    // If value is below the threshold, return black
    if (fabs(value) < blackThreshold)
    {
        snprintf(color, 10, BLACK);
        return;
    }

    // Normalize and scale the values to any index as needed
    double norm = (value + 1.0) / 2.0;
    int colorIndex = (int)(norm * 7);

    switch (colorIndex)
    {
        case 0:  snprintf(color, COLOR_BUFFER_SIZE, RED);     break;
        case 1:  snprintf(color, COLOR_BUFFER_SIZE, GREEN);   break;
        case 2:  snprintf(color, COLOR_BUFFER_SIZE, YELLOW);  break;
        case 3:  snprintf(color, COLOR_BUFFER_SIZE, BLUE);    break;
        case 4:  snprintf(color, COLOR_BUFFER_SIZE, MAGNETA); break;
        case 5:  snprintf(color, COLOR_BUFFER_SIZE, CYAN);    break;
        case 6:  snprintf(color, COLOR_BUFFER_SIZE, WHITE);   break;
        default: snprintf(color, COLOR_BUFFER_SIZE, BLACK);   break;
    }

}

/**
 *******************************************************************************
 * @brief:     Print the 2D node plane via a terminal
 * @parameter: array: Pointer to the 2D array with the associated node values
 * @parameter: rows: The number of rows in the array
 * @parameter: cols: The number of cols in the array
 * @return:    N/A
 *******************************************************************************
 */
void printWave(double** array, int rows, int cols)
{
    printf(CURSOR);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            double value = array[i][j];
            char color[COLOR_BUFFER_SIZE];
            getColor(value, color);
            printf("%s* " RESET, color);
        }
        printf("\n");
    }

    // ~ 30 fps
    usleep(33000);
}

/**
 *******************************************************************************
 * @brief:     Main function of the file
 * @parameter: N/A
 * @return:    N/A
 *******************************************************************************
 */
int main(void)
{
    // Allocate memory for 2D arrays
    double** Un_p1 = allocate2DArray(Nx, Ny);
    double** Un0 = allocate2DArray(Nx, Ny);
    double** Un_m1 = allocate2DArray(Nx, Ny);

    initializeArray(Un_p1, Nx, Ny);
    initializeArray(Un0, Nx, Ny);
    initializeArray(Un_m1, Nx, Ny);

    // Time marchings starts here
    for (int n = 0; n < n_stop; n++)
    {
        // Compute the general wave equation solution
        for (int jj = 1; jj < Ny - 1; jj++)
        {
            for (int ii = 1; ii < Nx - 1; ii++)
            {
                Un_p1[ii][jj] = 2 * Un0[ii][jj]
                    + Ox * Ox * (Un0[ii + 1][jj] - 2 * Un0[ii][jj] + Un0[ii - 1][jj])
                    + Oy * Oy * (Un0[ii][jj + 1] - 2 * Un0[ii][jj] + Un0[ii][jj - 1])
                    - Un_m1[ii][jj];
            }
        }

        // Source nodes
        Un_p1[xs1][ys1] = 1 * exp(-pow((n * dt - T0) / (w / 2), 2.0))
                            * sin(((2 * M_PI * c) / l) * (n * dt));

        // Left nodes
        int ii = 0;
        for (int jj = 1; jj < Ny - 1; jj++)
        {
            Un_p1[ii][jj] = Un0[ii + 1][jj] + (((c * dt - dx) / (c * dt + dx))
                                            * (Un_p1[ii + 1][jj] - Un0[ii][jj]));
        }

        // Right nodes
        ii = Nx - 1;
        for (int jj = 1; jj < Ny - 1; jj++)
        {
            Un_p1[ii][jj] = Un0[ii - 1][jj] + (((c * dt - dx) / (c * dt + dx))
                                            * (Un_p1[ii - 1][jj] - Un0[ii][jj]));
        }

        // Top nodes
        int jj = Ny - 1;
        for (int ii = 1; ii < Nx - 1; ii++)
        {
            Un_p1[ii][jj] = Un0[ii][jj - 1] + (((c * dt - dy) / (c * dt + dy))
                                            * (Un_p1[ii][jj - 1] - Un0[ii][jj]));
        }

        // Bottom nodes
        jj = 0;
        for (int ii = 1; ii < Nx - 1; ii++)
        {
            Un_p1[ii][jj] = Un0[ii][jj + 1] + (((c * dt - dy) / (c * dt + dy))
                                            * (Un_p1[ii][jj + 1] - Un0[ii][jj]));
        }

        // Simply average the corner values
        Un_p1[0][0] = 0.5 * (Un_p1[1][0] + Un_p1[0][1]);
        Un_p1[Nx - 1][0] = 0.5 * (Un_p1[Nx - 2][0] + Un_p1[Nx - 1][1]);
        Un_p1[Nx - 1][Ny - 1] = 0.5 * (Un_p1[Nx - 2][Ny - 1] + Un_p1[Nx - 1][Ny - 2]);
        Un_p1[0][Ny - 1] = 0.5 * (Un_p1[0][Ny - 2] + Un_p1[1][Ny - 1]);

        // Console print :)
        printWave(Un_p1, Nx, Ny);

        // Swap references
        double** temp = Un_m1;
        Un_m1 = Un0;
        Un0 = Un_p1;
        Un_p1 = temp;
    }
    // Free the memory
    free2DArray(Un_p1, Nx);
    free2DArray(Un0, Nx);
    free2DArray(Un_m1, Nx);

    return 0;
}

// ************************************End of file******************************

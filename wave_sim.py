"""
Filename: wave_sim.py
Author:   Danny Soppit
Description: Simulates and plots the 2D wave equation
"""

# ~~~~~~~~~~ Python Libraries ~~~~~~~~~~~~~

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Lx = 10e-6     # Length in x direction
Ly = 10e-6     # Length in y direction
dx = 0.12e-6   # Grid size in x direction
dy = 0.12e-6   # Grid size in y direction
n_stop = 150   # Number of time steps
l = 1.0e-6     # Wavelength
w = 18.0e-15   # Width of the pulse
T0 = 4.0e-15   # Initial time
c = 299792458  # Speed of light

# ~~~~~~~~~~ Class Definitions ~~~~~~~~~~~~~

class WaveSimulation2D:
    def __init__(self, Lx, Ly, dx, dy, n_stop, l, w, T0, c):

        # Mesh and simulation parameters
        self.Lx = Lx
        self.Ly = Ly
        self.dx = dx
        self.dy = dy
        self.n_stop = n_stop
        self.l = l
        self.w = w
        self.T0 = T0
        self.c = c

        # Grid setup
        self.X = np.arange(0, Lx + dx, dx)
        self.Y = np.arange(-Ly / 2, Ly / 2 + dy, dy)
        self.Nx = len(self.X)
        self.Ny = len(self.Y)

        # Time step
        self.dt = 1 / (c * np.sqrt((1 / dx**2) + (1 / dy**2)))
        self.T = np.arange(2 * self.dt, self.dt * (n_stop + 2), self.dt)

        # Courant numbers
        self.Ox = (c * self.dt) / dx
        self.Oy = (c * self.dt) / dy

        # Initialize fields
        self.Un1 = np.zeros((self.Nx, self.Ny))   # Time level n+1
        self.Un0 = np.zeros((self.Nx, self.Ny))   # Time level n
        self.Un_1 = np.zeros((self.Nx, self.Ny))  # Time level n-1

        # To store the field at each time step
        self.U_value = np.zeros((self.Nx, self.Ny, n_stop))

        # Set up grid for plotting
        self.X_grid, self.Y_grid = np.meshgrid(self.X, self.Y)

    def apply_source(self, n):

        # Apply source condition at the center of the grid
        self.Un1[50, 50] = (np.exp(-1 * ((n * self.dt - self.T0) / (self.w / 2))**2)
                            * np.sin(((2 * np.pi * self.c) / self.l) * (n * self.dt)))

    def update_interior(self):

        # Update interior nodes using finite difference scheme.
        for jj in range(1, self.Ny - 1):
            for ii in range(1, self.Nx - 1):
                self.Un1[ii, jj] = (2 * self.Un0[ii, jj]
                                    + self.Ox**2 * (self.Un0[ii + 1, jj] - 2 * self.Un0[ii, jj] + self.Un0[ii - 1, jj])
                                    + self.Oy**2 * (self.Un0[ii, jj + 1] - 2 * self.Un0[ii, jj] + self.Un0[ii, jj - 1])
                                    - self.Un_1[ii, jj])

    def update_boundaries(self):

        # Left and right boundaries
        for jj in range(1, self.Ny - 1):
            self.Un1[0, jj] = (self.Un0[1, jj] + ((self.c * self.dt - self.dx) / (self.c * self.dt + self.dx)) * (self.Un1[1, jj] - self.Un0[0, jj]))
            self.Un1[self.Nx - 1, jj] = (self.Un0[self.Nx - 2, jj] + ((self.c * self.dt - self.dx) / (self.c * self.dt + self.dx)) * (self.Un1[self.Nx - 2, jj] - self.Un0[self.Nx - 1, jj]))

        # Top and bottom boundaries
        for ii in range(1, self.Nx - 1):
            self.Un1[ii, self.Ny - 1] = (self.Un0[ii, self.Ny - 2] + ((self.c * self.dt - self.dy) / (self.c * self.dt + self.dy)) * (self.Un1[ii, self.Ny - 2] - self.Un0[ii, self.Ny - 1]))
            self.Un1[ii, 0] = (self.Un0[ii, 1] + ((self.c * self.dt - self.dy) / (self.c * self.dt + self.dy)) * (self.Un1[ii, 1] - self.Un0[ii, 0]))

        # Corner nodes
        self.Un1[0, 0] = 0.5 * (self.Un1[1, 0] + self.Un1[0, 1])
        self.Un1[self.Nx - 1, 0] = 0.5 * (self.Un1[self.Nx - 2, 0] + self.Un1[self.Nx - 1, 1])
        self.Un1[self.Nx - 1, self.Ny - 1] = 0.5 * (self.Un1[self.Nx - 2, self.Ny - 1] + self.Un1[self.Nx - 1, self.Ny - 2])
        self.Un1[0, self.Ny - 1] = 0.5 * (self.Un1[1, self.Ny - 1] + self.Un1[0, self.Ny - 2])

    def store_fields(self, n):
        self.U_value[:, :, n] = self.Un1

    def step_time(self):
        self.Un_1 = self.Un0.copy()
        self.Un0 = self.Un1.copy()

    def plot_solution(self, ax, n):
        ax.clear()
        ax.plot_surface(self.Y_grid, self.X_grid, self.U_value[:, :, n], cmap='viridis')
        ax.set_zlim(-0.6, 0.6)
        ax.set_xlabel('Y')
        ax.set_ylabel('X')
        ax.set_zlabel('U')
        ax.set_title(f"Time Step {n + 1}")
        plt.pause(0.001)

    def run_simulation(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for n in range(self.n_stop):
            self.update_interior()
            self.apply_source(n)
            self.update_boundaries()
            self.store_fields(n)
            self.plot_solution(ax, n)
            self.step_time()

        plt.show()


if __name__ == "__main__":
    simulation = WaveSimulation2D(Lx, Ly, dx, dy, n_stop, l, w, T0, c)
    simulation.run_simulation()

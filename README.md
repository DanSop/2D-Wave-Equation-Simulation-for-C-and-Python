# 2D Wave Equation in C & Python
## Information & Math

The 2D wave equation is used to model how many waves, such as light or sound, propagate in a two dimensional space over time. This equation is fairly fundamental everywhere in physics and engineering. 
### Definition:

```math
\frac{\partial^2 u}{\partial t^2} = c^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)
```

Applying finite differences (FD) for the interior nodes $` (i, j) `$ and rearranging gives:
```math
\frac{1}{c^2} \frac{u_{i,j}^{n+1} - 2u_{i,j}^{n} + u_{i,j}^{n-1}}{\Delta t^2} =
\frac{u_{i+1,j}^{n} - 2u_{i,j}^{n} + u_{i-1,j}^{n}}{\Delta x^2} + 
\frac{u_{i,j+1}^{n} - 2u_{i,j}^{n} + u_{i,j-1}^{n}}{\Delta y^2}
```

Solving for $` ( u_{i,j}^{n+1} ) `$:
```math
\frac{u_{i,j}^{n+1} - 2u_{i,j}^{n} + u_{i,j}^{n-1}}{\Delta t^2} = \frac{u_{i+1,j}^{n} - 2u_{i,j}^{n} + u_{i-1,j}^{n}}{\Delta x^2} + \frac{u_{i,j+1}^{n} - 2u_{i,j}^{n} + u_{i,j-1}^{n}}{\Delta y^2}
```

```math
u_{i,j}^{n+1} - 2u_{i,j}^{n} + u_{i,j}^{n-1} = \frac{\Delta t^2}{\Delta x^2} \left( u_{i+1,j}^{n} - 2u_{i,j}^{n} + u_{i-1,j}^{n} \right) + \frac{\Delta t^2}{\Delta y^2} \left( u_{i,j+1}^{n} - 2u_{i,j}^{n} + u_{i,j-1}^{n} \right)
```

```math
\boxed{ u_{i,j}^{n+1} = 2u_{i,j}^{n} + \theta_x \left( u_{i+1,j}^{n} - 2u_{i,j}^{n} + u_{i-1,j}^{n} \right) 
+ \theta_y \left( u_{i,j+1}^{n} - 2u_{i,j}^{n} + u_{i,j-1}^{n} \right) - u_{i,j}^{n-1} }
```
where 
```math
\theta_x = \frac{c\Delta t^2}{\Delta x^2}, \quad \theta_y = \frac{c\Delta t^2}{\Delta y^2}
```

### Boundary Conditions

For boundary conditions, applying the 1D wave equation's radiating boundary conditions (RBCs):
```math
\left(\frac{\partial^2}{\partial x^2} - \frac{1}{c^2}\frac{\partial^2}{\partial t^2}\right)u = 0 \implies 
\left(\frac{\partial}{\partial x} - \frac{1}{c}\frac{\partial}{\partial t}\right)\left(\frac{\partial}{\partial x} + \frac{1}{c}\frac{\partial}{\partial t}\right) u = 0
```
Using the forward propagating wave:
```math
\frac{\partial u}{\partial x} + \frac{1}{c}\frac{\partial u}{\partial t} = 0
```

Discretizing at  $` (i-1/2, j) `$:
```math
\left.\frac{\partial u}{\partial x}\right|_{i-1/2,j}^{n+1/2} = \frac{u_{i,j}^{n+1} - u_{i-1,j}^{n+1}}{\Delta x}
```

```math
\left.\frac{\partial u}{\partial t}\right|_{i-1/2,j}^{n+1/2} = \frac{u_{i,j}^{n+1} - u_{i,j}^n}{\Delta t}
```

Combining these:
```math
\frac{u_{i,j}^{n+1} - u_{i-1,j}^{n+1}}{\Delta x} + \frac{1}{c} \frac{u_{i,j}^{n+1} - u_{i,j}^n}{\Delta t} = 0
```

Solving for  $` (u_{i,j}^{n+1}) `$ yields the following for the right boundary:
```math
\boxed{ u_{i,j}^{n+1} = u_{i-1,j}^n + \frac{c\Delta t - \Delta x}{c \Delta t + \Delta x} \left( u_{i-1,j}^{n+1} - u_{i,j}^n \right)}
```

Apply the same technique for top, left, and bottom boundaries.

## Code Usage

To install the packges: ```pip install -r requirements.txt```. To run the Python script - ```python wave_sim.py``` - you will notice it will be "laggy" due to the computations and zero optimizations in the code.

To run the C version build it via ```gcc -o sim wave_sim.c -lm``` and then ```./sim```. The simulation via C will be a lot faster but it only prints in the terminal (for now atleast)! The way it is printed is not ideal and was implemented fairly quickly. Feel free to change the colors or the values for display.

## Example Output:

Python output:

![simulation](https://github.com/user-attachments/assets/0d03d4ba-e9a2-4a28-aa39-787ba6962f42)

C Output:

![terminal](https://github.com/user-attachments/assets/c1c02d69-48d1-46bb-ae0a-3ae691418b12)


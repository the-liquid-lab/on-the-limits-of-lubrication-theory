

\[@Cortelezzi1981\] propose a derivation of the equations governing the
initial value problem for a disturbed layer of viscous fluid,
considering both the effects of gravity and surface tension. We here
propose to revisit this paper by following an alternative derivation
based on the ideas of \[@Prosperetti1976\]. More specifically we split
the dynamical variables into a potential and a viscous contribution. We
restrict the analysis to a 2D setting, and note that all the obtained
results can be directly transposed to 3D, see \[@Cortelezzi1981\].

Vorticity production
====================

The vorticity is directed along
![\\boldsymbol e\_y](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20e_y "\boldsymbol e_y"):

![\\begin{equation}
\\boldsymbol \\omega = \\left(
\\begin{array}{c}
\\omega\_x \\\\
\\omega\_y \\\\
\\omega\_z
\\end{array}
\\right) =
\\left.
\\begin{array}{c}
\\partial\_x \\\\
0 \\\\
\\partial\_z
\\end{array}
\\right\|
\\times
\\left\|
\\begin{array}{c}
u \\\\
0 \\\\
w
\\end{array}
\\right. =
\\left(
\\frac{\\partial u}{\\partial z} - 
\\frac{\\partial w}{\\partial x}
\\right) \\boldsymbol e\_y
\\end{equation} \\qquad (1)](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%20%5Comega%20%3D%20%5Cleft%28%0A%5Cbegin%7Barray%7D%7Bc%7D%0A%5Comega_x%20%5C%5C%0A%5Comega_y%20%5C%5C%0A%5Comega_z%0A%5Cend%7Barray%7D%0A%5Cright%29%20%3D%0A%5Cleft.%0A%5Cbegin%7Barray%7D%7Bc%7D%0A%5Cpartial_x%20%5C%5C%0A0%20%5C%5C%0A%5Cpartial_z%0A%5Cend%7Barray%7D%0A%5Cright%7C%0A%5Ctimes%0A%5Cleft%7C%0A%5Cbegin%7Barray%7D%7Bc%7D%0Au%20%5C%5C%0A0%20%5C%5C%0Aw%0A%5Cend%7Barray%7D%0A%5Cright.%20%3D%0A%5Cleft%28%0A%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20z%7D%20-%20%0A%5Cfrac%7B%5Cpartial%20w%7D%7B%5Cpartial%20x%7D%0A%5Cright%29%20%5Cboldsymbol%20e_y%0A%20%5Cqquad%20%281%29 "\begin{equation}
\boldsymbol \omega = \left(
\begin{array}{c}
\omega_x \\
\omega_y \\
\omega_z
\end{array}
\right) =
\left.
\begin{array}{c}
\partial_x \\
0 \\
\partial_z
\end{array}
\right|
\times
\left|
\begin{array}{c}
u \\
0 \\
w
\end{array}
\right. =
\left(
\frac{\partial u}{\partial z} - 
\frac{\partial w}{\partial x}
\right) \boldsymbol e_y
\end{equation} \qquad (1)")

Introducing the viscous shear stress ![\\tau\_{xz} = \\mu \\left(
\\frac{\\partial u}{\\partial z} + \\frac{\\partial w}{\\partial x} \\right)](https://latex.codecogs.com/png.latex?%5Ctau_%7Bxz%7D%20%3D%20%5Cmu%20%5Cleft%28%0A%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20z%7D%20%2B%20%5Cfrac%7B%5Cpartial%20w%7D%7B%5Cpartial%20x%7D%20%5Cright%29 "\tau_{xz} = \mu \left(
\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x} \right)")
we can express the vorticity in terms of viscous shear stress:

![\\begin{equation}
\\omega\_y = \\frac{\\tau\_{xz}}{\\mu} - 2 \\frac{\\partial w}{\\partial x}
\\end{equation} \\qquad (2)](https://latex.codecogs.com/png.latex?%0A%5Comega_y%20%3D%20%5Cfrac%7B%5Ctau_%7Bxz%7D%7D%7B%5Cmu%7D%20-%202%20%5Cfrac%7B%5Cpartial%20w%7D%7B%5Cpartial%20x%7D%0A%20%5Cqquad%20%282%29 "\begin{equation}
\omega_y = \frac{\tau_{xz}}{\mu} - 2 \frac{\partial w}{\partial x}
\end{equation} \qquad (2)")

At the free surface, the tangential shear stress vanishes so
![\\omega\_y = - 2 \\frac{\\partial w}{\\partial x}](https://latex.codecogs.com/png.latex?%5Comega_y%20%3D%20-%202%20%5Cfrac%7B%5Cpartial%20w%7D%7B%5Cpartial%20x%7D "\omega_y = - 2 \frac{\partial w}{\partial x}")
there. This can be further reexpressed in terms of local free surface
slope variation with the help of the linearised kinematic condition
![w =
\\frac{\\partial \\eta}{\\partial t}](https://latex.codecogs.com/png.latex?w%20%3D%0A%5Cfrac%7B%5Cpartial%20%5Ceta%7D%7B%5Cpartial%20t%7D "w =
\frac{\partial \eta}{\partial t}"), where
![\\eta(x,t)](https://latex.codecogs.com/png.latex?%5Ceta%28x%2Ct%29 "\eta(x,t)")
denotes the free surface position. This relation reveals that vorticity
production at the free surface is equal to the time variation of the
local slope:

![\\begin{equation}
\\omega\_y = - 2 \\frac{\\partial^2 \\eta}{\\partial x \\partial t}
\\end{equation} \\qquad (3)](https://latex.codecogs.com/png.latex?%0A%5Comega_y%20%3D%20-%202%20%5Cfrac%7B%5Cpartial%5E2%20%5Ceta%7D%7B%5Cpartial%20x%20%5Cpartial%20t%7D%0A%20%5Cqquad%20%283%29 "\begin{equation}
\omega_y = - 2 \frac{\partial^2 \eta}{\partial x \partial t}
\end{equation} \qquad (3)")

Problem splitting: the potential and the viscous
================================================

The equations governing the problem are associated with mass
conservation:

![\\begin{equation}
\\nabla \\cdot \\boldsymbol u = 0,
\\end{equation} \\qquad (4)](https://latex.codecogs.com/png.latex?%0A%5Cnabla%20%5Ccdot%20%5Cboldsymbol%20u%20%3D%200%2C%0A%20%5Cqquad%20%284%29 "\begin{equation}
\nabla \cdot \boldsymbol u = 0,
\end{equation} \qquad (4)")

and momentum conservation:

![\\begin{equation}
\\frac{\\partial \\boldsymbol u}{\\partial t} = - \\frac{1}{\\rho} \\nabla p + \\nu \\nabla^2\\boldsymbol u - g \\boldsymbol e\_z.
\\end{equation} \\qquad (5)](https://latex.codecogs.com/png.latex?%0A%5Cfrac%7B%5Cpartial%20%5Cboldsymbol%20u%7D%7B%5Cpartial%20t%7D%20%3D%20-%20%5Cfrac%7B1%7D%7B%5Crho%7D%20%5Cnabla%20p%20%2B%20%5Cnu%20%5Cnabla%5E2%5Cboldsymbol%20u%20-%20g%20%5Cboldsymbol%20e_z.%0A%20%5Cqquad%20%285%29 "\begin{equation}
\frac{\partial \boldsymbol u}{\partial t} = - \frac{1}{\rho} \nabla p + \nu \nabla^2\boldsymbol u - g \boldsymbol e_z.
\end{equation} \qquad (5)")

Kinematical description of the interface
----------------------------------------

The interface is tracked at all times with the help of the color
function
![F(x,z,t) = z - \\eta(x,t)](https://latex.codecogs.com/png.latex?F%28x%2Cz%2Ct%29%20%3D%20z%20-%20%5Ceta%28x%2Ct%29 "F(x,z,t) = z - \eta(x,t)")
which cancels at the free surface. We here look for *standing waves* for
which ![\\eta(x,t) = a(t)
\\cos(kx)](https://latex.codecogs.com/png.latex?%5Ceta%28x%2Ct%29%20%3D%20a%28t%29%0A%5Ccos%28kx%29 "\eta(x,t) = a(t)
\cos(kx)"). From this definition we can derive the expression of the
*linearised* outward normal:

![\\begin{equation}
\\boldsymbol n = \\frac{\\nabla F}{\\\|\\nabla F\\\|} \\simeq \\left(k a(t) \\sin (kx), 1 \\right)
\\end{equation} \\qquad (6)](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%20n%20%3D%20%5Cfrac%7B%5Cnabla%20F%7D%7B%5C%7C%5Cnabla%20F%5C%7C%7D%20%5Csimeq%20%5Cleft%28k%20a%28t%29%20%5Csin%20%28kx%29%2C%201%20%5Cright%29%0A%20%5Cqquad%20%286%29 "\begin{equation}
\boldsymbol n = \frac{\nabla F}{\|\nabla F\|} \simeq \left(k a(t) \sin (kx), 1 \right)
\end{equation} \qquad (6)")

This approximation is valid as long as the *amplitude* (not the layer
thickness!) is small in front of the wavelength (i.e.
![k a(t) \\ll 1](https://latex.codecogs.com/png.latex?k%20a%28t%29%20%5Cll%201 "k a(t) \ll 1")).
Finally we require the function
![F](https://latex.codecogs.com/png.latex?F "F") to satisfy the
*kinematical boundary condition*
![DF/Dt=0](https://latex.codecogs.com/png.latex?DF%2FDt%3D0 "DF/Dt=0").
Linearisation yields to flatten this condition at
![z=0](https://latex.codecogs.com/png.latex?z%3D0 "z=0"):

![\\begin{equation}
\\frac{\\partial F}{\\partial t} + w \\frac{\\partial F}{\\partial z} = 0 \\quad \\text{on z = 0}
\\end{equation} \\qquad (7)](https://latex.codecogs.com/png.latex?%0A%5Cfrac%7B%5Cpartial%20F%7D%7B%5Cpartial%20t%7D%20%2B%20w%20%5Cfrac%7B%5Cpartial%20F%7D%7B%5Cpartial%20z%7D%20%3D%200%20%5Cquad%20%5Ctext%7Bon%20z%20%3D%200%7D%0A%20%5Cqquad%20%287%29 "\begin{equation}
\frac{\partial F}{\partial t} + w \frac{\partial F}{\partial z} = 0 \quad \text{on z = 0}
\end{equation} \qquad (7)")

Dynamical boundary conditions
-----------------------------

At the free surface, the normal stress jump is balanced by surface
tension:

![\\begin{equation}
-p + 2 \\mu \\frac{\\partial w}{\\partial z} = - \\gamma \\nabla \\cdot
\\boldsymbol n \\quad \\text{on} \\quad F = 0
\\end{equation} \\qquad (8)](https://latex.codecogs.com/png.latex?%0A-p%20%2B%202%20%5Cmu%20%5Cfrac%7B%5Cpartial%20w%7D%7B%5Cpartial%20z%7D%20%3D%20-%20%5Cgamma%20%5Cnabla%20%5Ccdot%0A%5Cboldsymbol%20n%20%5Cquad%20%5Ctext%7Bon%7D%20%5Cquad%20F%20%3D%200%0A%20%5Cqquad%20%288%29 "\begin{equation}
-p + 2 \mu \frac{\partial w}{\partial z} = - \gamma \nabla \cdot
\boldsymbol n \quad \text{on} \quad F = 0
\end{equation} \qquad (8)")

and the shear stress is required to vanish (*free surface* assumption):

![\\begin{equation}
\\tau\_{xz} = 0 \\quad \\text{on} \\quad F = 0
\\end{equation} \\qquad (9)](https://latex.codecogs.com/png.latex?%0A%5Ctau_%7Bxz%7D%20%3D%200%20%5Cquad%20%5Ctext%7Bon%7D%20%5Cquad%20F%20%3D%200%0A%20%5Cqquad%20%289%29 "\begin{equation}
\tau_{xz} = 0 \quad \text{on} \quad F = 0
\end{equation} \qquad (9)")

Potential and viscous parts
---------------------------

As demonstrated in \[@Prosperetti1976\] it is convenient to split the
dynamical fields into a *potential* and a *viscous* part:

![\\begin{equation}
\\boldsymbol u = \\boldsymbol u\_p + \\boldsymbol u\_v \\quad ; \\quad p = p\_p + p\_v
\\end{equation} \\qquad (10)](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%20u%20%3D%20%5Cboldsymbol%20u_p%20%2B%20%5Cboldsymbol%20u_v%20%5Cquad%20%3B%20%5Cquad%20p%20%3D%20p_p%20%2B%20p_v%0A%20%5Cqquad%20%2810%29 "\begin{equation}
\boldsymbol u = \boldsymbol u_p + \boldsymbol u_v \quad ; \quad p = p_p + p_v
\end{equation} \qquad (10)")

Bibliography
============

Cortelezzi, L. and Prosperetti, A. - Small-amplitude waves on the surface of a layer of a viscous liquid
--------------------------------------------------------------------------------------------------------

Prosperetti, Andrea - Viscous effects on small-amplitude surface waves
----------------------------------------------------------------------

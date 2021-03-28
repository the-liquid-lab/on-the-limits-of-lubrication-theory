propose a derivation of the equations governing the initial value
problem for a disturbed layer of viscous fluid, considering both the
effects of gravity and surface tension. We here propose to revisit this
paper by following an alternative derivation based on the ideas of .
More specifically we split the dynamical variables into a potential and
a viscous contribution. We restrict the analysis to a 2D setting, and
note that all the obtained results can be directly transposed to 3D, see
.

Vorticity production
====================

The vorticity is directed along **e**<sub>*y*</sub>:
$$\\boldsymbol \\omega = \\left(
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
\\right) \\boldsymbol e\_y$$
Introducing the viscous shear stress $\\tau\_{xz} = \\mu \\left(
\\frac{\\partial u}{\\partial z} + \\frac{\\partial w}{\\partial x} \\right)$
we can express the vorticity in terms of viscous shear stress:
$$\\omega\_y = \\frac{\\tau\_{xz}}{\\mu} - 2 \\frac{\\partial w}{\\partial x}$$
At the free surface, the tangential shear stress vanishes so
$\\omega\_y = - 2 \\frac{\\partial w}{\\partial x}$ there. This can be
further reexpressed in terms of local free surface slope variation with
the help of the linearised kinematic condition $w =
\\frac{\\partial \\eta}{\\partial t}$, where *η*(*x*, *t*) denotes the
free surface position. This relation reveals that vorticity production
at the free surface is equal to the time variation of the local slope:
$$\\omega\_y = - 2 \\frac{\\partial^2 \\eta}{\\partial x \\partial t}$$

Problem splitting: the potential and the viscous
================================================

The equations governing the problem are associated with mass
conservation:
∇ ⋅ **u** = 0,
and momentum conservation:
$$\\frac{\\partial \\boldsymbol u}{\\partial t} = - \\frac{1}{\\rho} \\nabla p + \\nu \\nabla^2\\boldsymbol u - g \\boldsymbol e\_z.$$

Kinematical description of the interface
----------------------------------------

The interface is tracked at all times with the help of the color
function *F*(*x*, *z*, *t*) = *z* − *η*(*x*, *t*) which cancels at the
free surface. We here look for *standing waves* for which
*η*(*x*, *t*) = *a*(*t*)cos (*k**x*). From this definition we can derive
the expression of the *linearised* outward normal:
$$\\boldsymbol n = \\frac{\\nabla F}{\\\|\\nabla F\\\|} \\simeq \\left(k a(t) \\sin (kx), 1 \\right)$$
This approximation is valid as long as the *amplitude* (not the layer
thickness!) is small in front of the wavelength (i.e. *k**a*(*t*) ≪ 1).
Finally we require the function *F* to satisfy the *kinematical boundary
condition* *D**F*/*D**t* = 0. Linearisation yields to flatten this
condition at *z* = 0:
$$\\frac{\\partial F}{\\partial t} + w \\frac{\\partial F}{\\partial z} = 0 \\quad \\text{on z = 0}$$

Dynamical boundary conditions
-----------------------------

### At the free surface

At the free surface, the normal stress jump is balanced by surface
tension:
$$-p + 2 \\mu \\frac{\\partial w}{\\partial z} = - \\gamma \\nabla \\cdot
\\boldsymbol n \\quad \\text{on} \\quad F = 0$$
and the shear stress is required to vanish (*free surface* assumption):
*τ*<sub>*x**z*</sub> = 0  on  *F* = 0

### At the bottom

Impermeability and adherence impose:
*u*(*y* =  − *h*) = *w*(*y* =  − *h*) = 0

Potential and viscous parts
---------------------------

As demonstrated in it is convenient to split the dynamical fields into a
*potential* and a *viscous* part:
**u** = **u**<sub>*p*</sub> + **u**<sub>*v*</sub>  ;  *p* = *p*<sub>*p*</sub> + *p*<sub>*v*</sub>

#### ⊳ The potential problem.

(**u**<sub>*p*</sub>, *p*<sub>*p*</sub>) is solution to the following
problem

$$\\begin{aligned}
&\\boldsymbol u\_p = \\nabla \\phi,\\\\
&\\nabla^2 \\phi = 0,\\\\
&\\frac{\\partial \\boldsymbol u\_p}{\\partial t} = -\\frac{1}{\\rho} \\nabla p\_p - g \\boldsymbol e\_z,  \\\\
&\\frac{\\partial F}{\\partial t} + w\_p \\frac{\\partial F}{\\partial z} = 0 \\quad \\text{on} \\quad z = 0, \\label{eq:kinematic\_potential} \\\\
&\\frac{\\partial \\phi}{\\partial z} = 0  \\quad \\text{on} \\quad z = -h.\\\\\\end{aligned}$$

Note that
<a href="#eq:kinematic_potential" data-reference-type="eqref" data-reference="eq:kinematic_potential">[eq:kinematic_potential]</a>
reads
$$w\_p = \\frac{\\partial \\eta}{\\partial t} \\label{eq:kinematic\_potential\_simplified}$$

#### ⊳ The viscous correction.

(**u**<sub>*v*</sub>, *p*<sub>*v*</sub>) are then governed by the
following equations:

$$\\begin{aligned}
&\\nabla \\cdot u\_v = 0,\\\\
&\\frac{\\partial \\boldsymbol u\_v}{\\partial t} = -\\frac{1}{\\rho} \\nabla p\_v + \\nu \\nabla^2 \\boldsymbol u\_v,\\\\
&w\_v \\frac{\\partial F}{\\partial z} = 0 \\quad \\text{on} \\quad z = 0,\\\\
&w\_v = u\_p + u\_v = 0  \\quad \\text{on} \\quad z = -h.\\end{aligned}$$

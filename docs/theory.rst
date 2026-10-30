======
Theory
======

Unit system
===========

The `International System of Units (SI) <https://en.wikipedia.org/wiki/International_System_of_Units>`_ is used for the definitions and computations throughout the simulation library. For convenience some of the angular quantities are specified in degrees; the computations are always done using radians.

Coordinate frames
=================

Multiple coordinate frames are defined in marine craft theory.
At this stage, the *Stonefish* library is intended for small scale simulations where Earth curvature can be neglected. Therefore, the **North-East-Down (NED) coordinate frame** is used throughout the simulation, as the world reference frame. The NED frame is a Cartesian coordinate frame, tangent to the Earth surface at a chosen geographic location. This location is called the home and defined by latitude and longitude.

Rigid body dynamics and collision
=================================

The *Stonefish* library utilises algoritms implemented in the *Bullet Physics* library for the computation of the rigid body kinematics, dynamics and collision. The simulation engine is implementing velocity-based dynamics with semi-implicit Euler integration. The dynamics are computed using the Sequential Impulse algorithm for separate dynamical bodies. Kinematic trees are handled through an implementation of the Featherstone multi-body algorithm. Rigid collisions are computed using an impulse-based approach. Soft collisions are possible by defining collision stiffeness and damping factors. Frictional forces are computed based on static and dynamic friction coefficients defined between materials. The Stribeck function is used for the transition between sticking and sliding phases.

Hydrodynamics
=============

The *Stonefish* library delivers a novel approach to simulating hydrodynamics and aerodynamics, by performing geometry based computations. The main focus is put on hydrodynamics as the library is directed towards marine robotics. The simulated effects include: added mass, buoyancy and drag. The drag is composed of 2 elements: form drag (quadratic) and skin friction. Hydrodynamic lift is not computed for general bodies but can be found in the model of a rudder actuator, which is considered an actuated hydrofoil.

Added mass
----------

The added mass effect is encountered during the acceleration of bodies submerged in liquid. It is manifested by the dynamical response of the body as if it was heavier than it actually is, because the surrounding liquid has to be moved together with it. The added mass is formulated as a 6x6 matrix, which describes how acceleration in every direction is affected by the fluid. It is a common practice to use a diagonal representation of the matrix, which is also the case here. Moreover, due to the fact that the underlying physics library is not capable of accepting different values of mass for each axis, a mean value is used for all of the axes. In terms of the added inertia tensor, a vector of 3 inertia moments is used. The values of the added mass matrix are computed based on an automatic approximation of the body geometry using one of the 3 solids: sphere, cylinder or ellipsoid.

Buoyancy
--------

The buoyancy force is calculated based on the sum of hydrostatic forces acting on the body surface. The actual geometry is used to compute force at each face of the mesh, depending on the depth of the face centre. It allows for simulating realistic buoyancy force at the surface of the ocean, with and without geometrical waves. When the body is completely submerged the buoyancy force is based on the volume of the mesh, computed automatically during loading.

Drag
----

The drag forces are calculated as a sum of forces acting on each face of the body surface. To obtain precise values of these forces it is required to solve Navier-Stokes equations, which is not possible for a general 3D case in realtime. Therefore, the computations implemented in the *Stonefish* library have to be based on the local velocity of fluid as if there was no body. The result is not quantitively correct but it gives a good approximation and allows for effects not possible when using simple formulas, e.g., a water current acting on a part of the body.

For every face of the mesh, described by its area *A*, its outward unit normal *n* and its centroid *c*, the relative velocity of the fluid with respect to the body is evaluated as

.. code-block:: text

    v_c = v_fluid(c) - (v + omega x (c - p))

where *v* and *omega* are the linear and angular velocity of the body, and *p* is the reference point used for the computation. Only the faces which are exposed to the flow, i.e. for which *v_c* · *n* < 0, contribute to the form drag. Their contributions are accumulated and scaled by the dynamic pressure factor, which yields the total form drag force

.. code-block:: text

    F_form = 0.5 * rho * Cd * sum_faces( (-v_c . n) * A * v_c )

Each face contributes the normal velocity component *-v_c* · *n*, multiplied by the relative velocity vector *v_c* and the face area *A*, so the accumulated term is quadratic in the relative velocity. The form drag coefficient *Cd* is therefore dimensionless and independent of speed: a value calibrated at one velocity produces correct drag at every other velocity, as in the standard drag equation.

For a convex closed body surface and a uniform relative velocity *v_c*, the exposed faces satisfy

.. code-block:: text

    sum_faces( (-v_c . n) * A ) = |v_c| * A_proj

In other words, after factoring out the velocity magnitude, the sum of the face areas weighted by the normal velocity component equals the projected area *A_proj*. For a convex body *A_proj* is the area of the shadow that the body casts on a plane perpendicular to *v_c* (for a sphere *A_proj* = pi * R^2). For concave meshes and compound bodies the identity does not hold exactly, because some faces can be hidden behind others: *A_proj* has to be read as the projected area counted with multiplicity, i.e. every exposed face contributes its own projected area even if it is occluded from the outside. In still water and for a purely translational motion, *v_c* = -*v*, and the formula above reduces to the standard drag equation

.. code-block:: text

    F_form = -0.5 * rho * Cd * A_proj * |v| * v

i.e. a force of magnitude ``0.5 * rho * Cd * A_proj * |v|^2`` which always opposes the motion. When the body rotates or a water current is present, the local relative velocity is used separately for every face, which additionally produces a torque and lets a non uniform flow deflect the body instead of only dragging it.

The skin friction (viscous drag) is computed from the tangential component of the relative velocity at each face, *v_t* = *v_c* - (*v_c* · *n*) *n*. The sum of *A* times *v_t* over all faces, scaled by the liquid density *rho* and the skin friction coefficient *Cf*, gives a force which is linear in the relative velocity:

.. code-block:: text

    F_friction = rho * Cf * sum_faces( A * v_t )

Both coefficients are independent of the velocity. The form drag coefficient *Cd* is dimensionless: it scales the quadratic term above and can be calibrated at any velocity. The skin friction coefficient *Cf* is the coefficient of the linear skin friction term above; it is not claimed to be dimensionless. The coefficients are estimated automatically from the geometry of the body and can be overridden per body axis with the ``quadratic_drag`` (form drag, quadratic) and ``viscous_drag`` (skin friction, linear) attributes of the ``<hydrodynamics>`` tag, or with the ``SetHydrodynamicCoefficients`` method. The three components of each coefficient are mixed according to the direction of the resulting force expressed in the body frame, which allows different drag to be defined for the principal axes of the body.
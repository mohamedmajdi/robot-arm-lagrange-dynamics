# Design of generic algorithm to calculate the lagrange dynamics of serial robot manipulators and implementing it using Matlab

Derivation of the dynamic model of a manipulator plays an important role for simulation of motion, analysis of manipulator structures, and design of control algorithms. The dynamic model of a manipulator provides a description of the relationship between the joint actuator torques and the motion of the structure. With Lagrange formulation, the equations of motion can be derived in a systematic way independently of the reference coordinate frame. Once a set of variables qi,i=1,...,n termed generalised coordinates, are chosen which
effectively describe the link positions of ann-DOF manipulator, the Lagrangian of the mechanical system can be defined as a function of the generalised coordinates.

# Algorithm Steps

1. Data collection from user
2. Forward Kinematics calculations
3. Jacobian matrices calculation
4. Finding Manipulator Inertia Matrix
5. Velocity Vector calculation
6. Determining Gravitational Vector
7. Calculating the final dynamic equation of motion

# Code demo



https://github.com/mohamedmajdi/robot-arm-lagrange-dynamics/assets/69417860/0c979a8f-f00b-4120-9661-e5c97c504e86



# Key techniques used in programming

Matlab software provides a lot of features for working with Matrices and
symbols, thus enabling writing algebraic robotic softwares easily. A cell array is
a data type with indexed data containers called cells, where each cell can
contain any type of data. We used that data type to store matrices resulting from
different calculations through the code, then advanced array indexing techniques
were used to extract desired cells from the cell arrays. For loops were used to
represent any summation or looping over certain portions of code. Also, if
conditions were implemented as logical thresholdings. The sequence of
operations was taken into consideration and the different data types used were
chosen based on its position in the program. Some data were treated as symbols
using the symbolic toolbox and others were used for mathematical comparisons
and operations and were treated as numbers.

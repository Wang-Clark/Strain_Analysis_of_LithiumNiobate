## Description
The notebook documents used to derive and calculate the anisotropic strain distribution in piezoelectric material thin films, specifically lithium niobate (TFLN) and GaAs. This part can achieve more precise computational results through multiphysics finite element simulations (solid mechanics + electrostatics). 

The Python `.ipynb` file (`Strain_to_EnergyShift.ipynb`) allows users to calculate the effects of strain on the band structure of GaAs based on the simulated strain data, with the results saved in `Strain_induced_EnergyShift.dat`. Detailed explanations are included in documents.

## Software Requirements
- **System Requirements:** A system capable of running Python 3.12 and Jupyter Notebook.  
- **Python Dependencies:** `numpy`, `scipy`, `matplotlib`.

## Eulerian Angle Rotation and Material Strain Analysis
Details of the Eulerian angle rotation and material strain analysis. The code can be easily modified to include other crystal orientations of lithium niobate or other materials. The derived strain distribution, dependent on crystal orientation angles, can be numerically solved using the file `GaAsStrain_from_LNStresstransfer.ipynb`,, with the output provided in `GaAsStrainfromLNStress_Calculation.dat`. For more precise results, this part can be calculated using multiphysics finite element simulations (combining solid mechanics and electrostatic fields), with the output provided in `GaAsStrain_fromFEM_Simulation.dat`.

Due to GitHub's Markdown not supporting MathJax formula rendering, equations may appear with ghosting. This can be resolved by adding the plugin 'GitHub with MathJax'.
### Rotation Matrices

The ZXZ Eulerian angle rotation is used to  transform the anisotropic property tensors (or the coordinate axes system) of LN, passive rotation `A = a1 . a2 . a3`, while conventional active rotation: `B = a3 . a2 . a1`.

Rotation matrices:

$$
\pmb{a1 = \begin{pmatrix} \cos\theta & \sin\theta & 0 \\ -\sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{pmatrix};}
$$

$$
\pmb{a2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos\phi & \sin\phi \\ 0 & -\sin\phi & \cos\phi \end{pmatrix};}
$$

$$
\pmb{a3 = \begin{pmatrix} \cos\varphi & \sin\varphi & 0 \\ -\sin\varphi & \cos\varphi & 0 \\ 0 & 0 & 1 \end{pmatrix};}
$$

$$
\pmb{A = a1 . a2 . a3;}
$$

Substitute specific angles into `A`:

$$
\pmb{Ax = A \Big|_{\theta \to -\theta, \phi \to -\pi/2, \varphi \to -\pi/2};}
$$

Resulting matrix for X-cut LN:

$$
\pmb{Ax = \begin{pmatrix} 0 & -\cos\theta & \sin\theta \\ 0 & -\sin\theta & -\cos\theta \\ 1 & 0 & 0 \end{pmatrix}}
$$

## Bond Transformation Matrix for Strain and Stress

Flatten the `Ax` matrix:

$$
\pmb{\{a11, a12, a13, a21, a22, a23, a31, a32, a33\} = \text{Flatten}[Ax];}
$$

Strain transformation matrix `Nstrain`:

$$
\pmb{Nstrain = \begin{pmatrix}
a11^2 & a21^2 & a31^2 & 2a21a31 & 2a31a11 & 2a11a21 \\
a12^2 & a22^2 & a32^2 & 2a22a32 & 2a32a12 & 2a12a22 \\
a13^2 & a23^2 & a33^2 & 2a23a33 & 2a33a13 & 2a13a23 \\
a12a13 & a22a23 & a32a33 & a22a33 + a32a23 & a12a33 + a32a13 & a22a13 + a12a23 \\
a13a11 & a23a21 & a33a31 & a21a33 + a31a23 & a31a13 + a11a33 & a11a23 + a21a13 \\
a11a12 & a21a22 & a31a32 & a21a32 + a31a22 & a31a12 + a11a32 & a11a22 + a21a12
\end{pmatrix};}
$$

Stress transformation matrix:

$$
\pmb{Mstress = \text{FullSimplify}[\text{Inverse}[Nstrain]];}
$$

Simplified strain transformation matrix:

$$
\pmb{\text{FullSimplify}[Nstrain] = \begin{pmatrix}
0 & 0 & 1 & 0 & 0 & 0 \\
\cos^2\theta & \sin^2\theta & 0 & 0 & 0 & \sin(2\theta) \\
\sin^2\theta & \cos^2\theta & 0 & 0 & 0 & -2\cos\theta\sin\theta \\
-\cos\theta\sin\theta & \cos\theta\sin\theta & 0 & 0 & 0 & \cos(2\theta) \\
0 & 0 & 0 & -\cos\theta & \sin\theta & 0 \\
0 & 0 & 0 & -\sin\theta & -\cos\theta & 0
\end{pmatrix}}
$$

## Piezoelectric Coupling Tensor for X-cut Lithium Niobate (LN)

Piezoelectric tensor `ePiezo` for X-cut LN:

$$
\pmb{ePiezo = \begin{pmatrix} 0 & 0 & 0 & 0 & e15 & -e22 \\ -e22 & e22 & 0 & e15 & 0 & 0 \\ e31 & e31 & e33 & 0 & 0 & 0 \end{pmatrix};}
$$

Substitute specific values:

$$
\pmb{ecaxis = ePiezo \Big|_{e15 \to 3.7, e22 \to 2.5, e31 \to 0.2, e33 \to 1.3};}
$$

Resulting matrix for c-axis:

$$
\pmb{ecaxis = \begin{pmatrix} 0 & 0 & 0 & 0 & 3.7 & -2.5 \\ -2.5 & 2.5 & 0 & 3.7 & 0 & 0 \\ 0.2 & 0.2 & 1.3 & 0 & 0 & 0 \end{pmatrix}}
$$

Transformed piezoelectric tensor for X-cut LN:

$$
\pmb{eXcut = \text{FullSimplify}[Ax . ecaxis . \text{Transpose}[Mstress]];}
$$

## Stress Tensor at X-cut LN

Electric field vector:

$$
\pmb{Fp = \{0, E, 0\};}
$$

Stress tensor:

$$
\pmb{TXcutLN = \text{Transpose}[eXcut] . Fp;}
$$

Resulting stress tensor:

$$
\pmb{TXcutLN = \begin{pmatrix}
E \cos\theta (2.95 - 3.15 \cos(2\theta) - 2.5 \cos\theta \sin\theta) \\
E (-1.3 \cos^3\theta - 7.6 \cos\theta \sin^2\theta - 2.5 \sin^3\theta) \\
E (-0.2 \cos\theta + 2.5 \sin\theta) \\
0 \\
0 \\
E \sin\theta (0.55 - 3.15 \cos(2\theta) - 2.5 \cos\theta \sin\theta)
\end{pmatrix}}
$$

## Strain Tensor in GaAs

Compliance tensor for GaAs:

$$
\pmb{sGaAs = \begin{pmatrix}
12.64 & -4.23 & -4.23 & 0 & 0 & 0 \\
-4.23 & 12.64 & -4.23 & 0 & 0 & 0 \\
-4.23 & -3.58 & 12.64 & 0 & 0 & 0 \\
0 & 0 & 0 & 18.6 & 0 & 0 \\
0 & 0 & 0 & 0 & 18.6 & 0 \\
0 & 0 & 0 & 0 & 0 & 18.6
\end{pmatrix};}
$$

Strain tensor in GaAs:

$$
\pmb{StrainGaAs = \text{FullSimplify}[sGaAs . TXcutLN];}
$$

Resulting strain tensor:

$$
\pmb{StrainGaAs = \begin{pmatrix}
E (30.38\cos\theta - 26.47 \cos(3\theta) - 10.70  \sin\theta - 10.70\sin(3\theta)) \\
E (-41.38 \cos\theta + 26.47 \cos(3\theta) - 32.11  \sin\theta + 10.70\sin(3\theta)) \\
E (3.91 \cos\theta + 42.81 \sin\theta ) \\
0 \\
0 \\
E \sin\theta (10.23 - 58.59 \cos(2\theta) - 46.5 \cos\theta \sin\theta)
\end{pmatrix}}
$$

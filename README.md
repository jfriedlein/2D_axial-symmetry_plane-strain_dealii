# 2D_axial-symmetry_plane-strain_dealii
Proposal of a framework and functions to handle 2D computations such as axial symmetry or plane strain in deal.ii

## ToDo
Add a simple code that uses the functions in the correct order with a Doxygen documentation. E.g. expand on step-3 using the Rod ([Numerical examples in deal.ii](https://github.com/jfriedlein/Numerical_examples_in_dealii)) in 2D and 3D.

## Remarks
* The code outlined below is also compatible to 3D, so you don't have to differentiate between 2D and 3D. It does no harm to your 3D computation.
* Obviously the code design is far from optimal. We recompute the same values multiple times and even slightly slow down a standard 3D computaton for which none of the below steps are necessary. The focus is currently on usability and comprehensibility.
* Axisymmetry is "exact" in tangential direction. We utilise no discretisation in this "third dimension". Hence, keep in mind that when you try to compare the axisymmetric computation with a 3D model, the latter needs a fairly fine spatial discretisation in the tangential direction to converge to the axisymmetric model.
* Plane strain results can be verified by 3D computations. However, it is hardly possible to make the 3D model thick enough (extension in the third dimension) to represent the plane strain state. We were able to obtain the best results by applying symmetry constraints on both z-planes, which forces the model to acquire no normal strains in the thickness direction ("plane strain").

## Argument list
* FEValues<dim> fe_values_ref_u; // The FEValues element that corresponds to the displacement dofs, e.g. extracted via fe_values_ref[u_fe]
* SymmetricTensor<2,dim> eps_n1; // small strain symmetric second order strain tensor
* SymmetricTensor<4,3> Tangent_3D; // 3D fourth order tangent modulus
	
<a href="https://www.codecogs.com/eqnedit.php?latex=\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" title="\overset{4}{C} = \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon }}" /></a>
	
* SymmetricTensor<4,dim> Tangent; // fourth order tangent modulus in dimension dim
	
* SymmetricTensor<2,dim> Tangent_theta;

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{C}_\theta&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial&space;{\varepsilon_\theta&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{C}_\theta&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial&space;{\varepsilon_\theta&space;}}" title="\boldsymbol{C}_\theta = \frac{\partial\boldsymbol{\sigma}}{\partial {\varepsilon_\theta }}" /></a>

* unsigned int type_2D; // 0: plane strain; 1: axial symmetry
* Vector<double> current_solution; // Vector containing the displacements for all dofs
* unsigned int k; // number of the current quadrature point (e.g. 0..3 for 4 QPs in 2D per cell)
* unsigned int j; // number of the current dof (innermost loop k->i->j; e.g. 0...7 for 8 dofs per cell in 2D)


## Plane strain
The plane strain case is handled straightforward. We expand the 2D strain tensor by the 3D components, which are all zero. Then we input this 3D strain tensor into a 3D material model and compute the full 3D stress tensor, which is, for instance, needed to compute the true von Mises stress (Exception: For elasticity the out-of-plane stress can directly be computed from Hooke's law). We extract only the relevant 2D components for the computation of the residuum and tangent.

To implement this you only require steps 1, 2 and 3 as outlined below for the more general case of axial symmetry. Because plane strain is a special case of axial symmetry, you can kill two birds with one stone by implementing steps 1 to 6 (compatible to 2D plane strain, axisymmetry and 3D).

<img src="https://github.com/jfriedlein/2D_axial-symmetry_plane-strain_dealii/blob/master/images/plane%20strain%20-%20sketch.jpg" width="1000">


## Axial symmetry
### Background
To start with, we have to integrate over a circular domain as outlined for deal.ii here: [mail-archive-msg08250](https://www.mail-archive.com/dealii@googlegroups.com/msg08250.html) (and in the following messages). In short: We have to multiply the standard JxW-value by the factor (2 * pi * r), where r is the radial coordinate of the current quadrature point. As a consequence, the results we obtain, e.g. the global force acting on the loaded face, represents the force that acts on the full model (the 2D cross section rotated by the above factor (2 * pi) producing the 360° model). If you leave the factor of (2 * pi) out of the integration, your results correspond to a section of 1 rad of the full model.

Besides this minor addition, the following major extensions are needed in the assembly routine when mechanical problems (displacement, deformation) shall be solved that use kinematics (e.g. small strain, deformation gradient).

The 2D displacement gradient computed from the displacement field (components u and w) needs to be expanded by the normal strain eps_theta. The latter is, in contrast to 2D plane strain, non-zero as for instance explained by Petr Krysl in ["A pragmatic Introduction to the Finite Element Method"](https://www.semanticscholar.org/paper/A-Pragmatic-Introduction-to-the-Finite-Element-for-Krysl/1a1e70b1dae4e971ba21af396f54b7fbaac14ffd) in chapter 15.3 (excerpt in \images folder) or sketched in the figure below.

@todo add screenshot from Krysl with citation to the image folder

<img src="https://github.com/jfriedlein/2D_axial-symmetry_plane-strain_dealii/blob/master/images/axialsymmetry%20-%20sketch.jpg" width="1000">

To sum up, radial displacements u_r produce tangential strains eps_theta = u_r / r.

### Small strains

This strain in the third dimension must be inserted as the out-of-plane component in the expanded 2D strain tensor as shown in the following figure.

@todo add this figure

#### Steps:

1. This can be done by calling the following function after the strain tensor eps_n1 was constructed (e.g. as symmetric part of the gradient of the deformations) and before calling the material model. The latter needs to handle the 3D case (3D strain and stress tensors) and must be called with the expanded strain 'eps_n1_3D'. (Everything located inside the loop over the quadrature points)

@todo explain parameter.type_2D and add enumerator to the handling_2D.h header, add enumerator axiSym instead of type_2D in input

@todo change fe_values_ref to ..._u, if possible

		SymmetricTensor<2,3> eps_n1_3D = prepare_strain<dim> (eps_n1, type_2D, fe_values_ref, current_solution, k);

2. Call the 3D material model with the properly prepared 3D strain tensor.

		{eps_n1_3D, history_n} -> material model -> {stress_3D, Tangent_3D, history_tmp}

3. Extract the 2D tangent modulus 'Tangent' from the full 3D quantity. Moreover, the 2D stress tensor can be extracted, which is, for instance, needed for the residuum. For (dim=3) the 'extract_dim' functions simply return the input argument unchanged. Call:

		SymmetricTensor<2,dim> Tangent = extract_dim<dim> (Tangent_3D);
		SymmetricTensor<2,dim> stress = extract_dim<dim> (stress_3D);
		
4. Extract from the full 3D tangent modulus 'Tangent_3D' (fourth order tensor) the components that contain the dependency of the stress tensor on the out-of-plane strain eps_theta summarised in 'Tangent_theta'. Call:

		SymmetricTensor<2,dim> Tangent_theta = extract_theta<dim> (Tangent_3D);

5. Modify the JxW-value for the axial symmetric domain (see factor 2 * pi * r above) by calling

		const double JxW = get_JxW<dim> (type_2D, fe_values_ref, k);

6. In the loop over the dofs 'j' (to assemble the cell_matrix(i,j)) compute the theta contribution 'Tangent_axisym' and add it to the linearisation of the stress tensor (typically 'Tangent * sym_shape_gradient_wrt_ref_config_j')

		SymmetricTensor<2,dim> Tangent_axisym = get_Tangent_axisym_addOn(Tangent_theta, fe_values_ref[u_fe],k,j);
		cell_matrix(i,j) += ( sym_shape_gradient_wrt_ref_config_i * ( Tangent * sym_shape_gradient_wrt_ref_config_j + Tangent_axisym ) * JxW;

### Finite strains

repeat using the deformation gradient and the function prepare_DefoGrad



# 2D_axial-symmetry_plane-strain_dealii
Proposal of a framework and functions to handle 2D computations such as axial symmetry or plane strain in deal.ii

## ToDo
Add a simple code that uses the functions in the correct order with a Doxygen documentation. E.g. expand on step-3 using the Rod ([Numerical examples in deal.ii](https://github.com/jfriedlein/Numerical_examples_in_dealii)) in 2D and 3D.

## Remarks
* Axisymmetry is "exact" in tangential direction. We utilise no discretisation in this "third dimension". Hence, keep in mind that when you try to compare the axisymmetric computation with a 3D model, the latter needs a fairly fine spatial discretisation in the tangential direction to converge to the axisymmetric model.
* Plane strain results can be verified by 3D computations. However, it is hardly possible to make the 3D model thick enough (extension in the third dimension) to represent the plane strain state. We were able to obtain the best results by applying symmetry constraints on both z-planes, which forces the model to acquire no normal strains in the thickness direction ("plane strain").


## Plane strain
<img src="https://github.com/jfriedlein/2D_axial-symmetry_plane-strain_dealii/blob/master/images/plane%20strain%20-%20sketch.jpg" width="1000">


## Axial symmetry
To start with, we have to integrate over a circular domain as outlined for deal.ii here: [mail-archive-msg08250](https://www.mail-archive.com/dealii@googlegroups.com/msg08250.html) (and in the following messages). In short: We have to multiply the standard JxW-value with the factor (2 * pi * r), where r is the radial coordinate of the current quadrature point. As a consequence, the results we obtain, e.g. the global force acting on the loaded face, represents the force that acts on the full model (the 2D cross section rotated by the above factor (2 * pi) producing the 360° model). If you leave the factor of (2 * pi) out of the integration, your results correspond to a section of 1 rad of the full model.

Besides this minor addition, the following major extensions are needed in the assembly routine when mechanical problems shall be solved that use kinematics (e.g. small strain, deformation gradient).

The 2D displacement gradient computed from the displacement field (components u and w) needs to be expanded by the normal strain eps_theta. The latter is, in contrast to 2D plane strain, not zero as for instance explained by Petr Krysl in ["A pragmatic Introduction to the Finite Element Method"](https://www.semanticscholar.org/paper/A-Pragmatic-Introduction-to-the-Finite-Element-for-Krysl/1a1e70b1dae4e971ba21af396f54b7fbaac14ffd) in chapter 15.3 (excerpt in \images folder).

@todo add screenshot from Krysl with citation to the image folder

<img src="https://github.com/jfriedlein/2D_axial-symmetry_plane-strain_dealii/blob/master/images/axialsymmetry%20-%20sketch.jpg" width="1000">

To sum up, radial displacements u_r produce tangential strains eps_theta = u_r / r.

### Small strains

This strain in the third dimension must be inserted as the out-of-plane component in the expanded 2D strain tensor as shown in the following figure.

@todo add this figure

1. This can be done by calling the following function after the strain tensor eps_n1 was constructed (e.g. as symmetric part of the gradient of the deformations) and before calling the material model. The latter needs to handle the 3D case (3D strain and stress tensors) and must be called with the expanded strain 'eps_n1_3D'. (Everything located inside the loop over the quadrature points)

@todo explain parametert.type_2D and add enumerator to the handling_2D.h header

@todo change fe_values_ref to ..._u, if possible

		SymmetricTensor<2,3> eps_n1_3D = prepare_strain<dim> (eps_n1, type_2D, fe_values_ref, current_solution, k);

2. Call the 3D material model with the properly prepared 3D strain tensor.

		{eps_n1_3D, history_n} -> material model -> {stress_3D, Tangent_3D, history_tmp}

3. Extract from the full 3D tangent modulus 'Tangent_3D' (fourth order tensor) the components that contain the dependency of the stress tensor on the out-of-plane strain eps_theta summarised in 'Tangent_theta'. Call:

		SymmetricTensor<2,dim> Tangent_theta = extract_theta<dim> (Tangent_3D);

4. Modify the JxW-value for the axial symmetric domain by calling

		const double JxW = get_JxW<dim> (type_2D, fe_values_ref, k);

5. In the loop over the dofs 'j' (to assemble the cell_matrix(i,j)) compute the theta contribution 'Tangent_axisym' and add it to the linearisation of the stress tensor (typically 'Tangent * sym_shape_gradient_wrt_ref_config_j')

		SymmetricTensor<2,dim> Tangent_axisym = get_Tangent_axisym_addOn(Tangent_theta, fe_values_ref[u_fe],k,j);
		cell_matrix(i,j) += ( sym_shape_gradient_wrt_ref_config_i * ( Tangent * sym_shape_gradient_wrt_ref_config_j + Tangent_axisym ) * JxW;

### Finite strains

repeat using the deformation gradient and the function prepare_DefoGrad



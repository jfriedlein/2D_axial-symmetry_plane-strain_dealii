# 2D_axial-symmetry_plane-strain_dealii
Proposal of a framework and functions to handle 2D computations such as axial symmetry or plane strain in deal.ii

## Axial symmetry
To start with, we have to integrate over a circular domain as outlined for deal.ii here: https://www.mail-archive.com/dealii@googlegroups.com/msg08250.html (and in the following messages). In short: We have to multiply the standard JxW-value with the factor (2*pi*r), where r is the radial coordinate of the current quadrature point.

Besides this minor addition, the following major extensions are needed in the assembly routine when mechanical problems shall be solved that use kinematics (e.g. small strain, deformation gradient).

The 2D displacement gradient computed from the displacement field (components u and w) needs to be expanded by the normal strain eps_theta. The latter is, in contrast to 2D plane strain, not zero as for instance explained by Petr Krysl in https://www.semanticscholar.org/paper/A-Pragmatic-Introduction-to-the-Finite-Element-for-Krysl/1a1e70b1dae4e971ba21af396f54b7fbaac14ffd in chapter 15.3.

@todo add sketch and maybe explanation in own words

To sum up, radial displacements u_r produce tangential strains eps_theta = u_r / r.

### Small strains

This strain in the third dimension must be inserted as the out-of-plane component in the expanded 2D strain tensor as shown in the following figure.

@todo add this figure

1. This can be done by calling the following function after the strain tensor eps_n1 was constructed and before calling the material model with this strain (in the loop over the quadrature points).

@todo explain parametert.type_2D and add enumerator to the handling_2D.h header

		SymmetricTensor<2,3> eps_n1_3D = prepare_strain<dim> (eps_n1, parameter.type_2D, fe_values_ref, current_solution, k);

2. Call the 3D material model with the properly prepared 3D strain tensor.

3. Extract from the full 3D tangent modulus (fourth order tensor Tangent_3D) the components that contain the dependency of the stress tensor on the out-of-plane strain eps_theta summarised in Tangent_theta. Call:

		SymmetricTensor<2,2> Tangent_theta = extract_theta<dim> (Tangent_3D);

4. Modify the JxW-value for the axial symmetric domain by calling

		const double JxW = get_JxW<dim> (parameter.type_2D, fe_values_ref, k);

5. In the loop over the dofs j (to assemble the cell_matrix(i,j)) compute the value of the shape function for the u-component of the displacements

		double shape_fnc_j_u = fe_values_ref[u_fe].value(j,k)[enums::u];

6. Add the theta contribution to the tangent as

		cell_matrix(i,j) += ( sym_shape_gradient_wrt_ref_config_i * ( Tangent * sym_shape_gradient_wrt_ref_config_j + Tangent_theta/get_radial_x<dim>(fe_values_ref,k) * shape_fnc_j_u ) * JxW;

### Finite strains

repeat using the deformation gradient and the function prepare_DefoGrad



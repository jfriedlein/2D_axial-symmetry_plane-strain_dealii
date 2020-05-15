#ifndef handling_2D_h
#define handling_2D_h

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>

#include "../MA-Code/enumerator_list.h"

// @todo Check whether the following three headers are needed at all
#include <iostream>
#include <fstream>
#include <cmath>

/*!
 * Extract the dim components from a full 3D tensor
 */
template<int dim>
SymmetricTensor<2,dim> extract_dim ( SymmetricTensor<2,3> &symTensor_3D )
{
    SymmetricTensor<2,dim> symTensor_dim;
    for ( unsigned int i=0; i<dim; i++)
        for ( unsigned int j=i; j<dim; j++)
            symTensor_dim[i][j] = symTensor_3D[i][j];

    return symTensor_dim;
}
template<int dim>
SymmetricTensor<4,dim> extract_dim ( SymmetricTensor<4,3> &symTensor_3D )
{
	SymmetricTensor<4,dim> symTensor_dim;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			for ( unsigned int k=0; k<dim; ++k )
				for ( unsigned int l=k; l<dim; ++l )
				  symTensor_dim[i][j][k][l] = symTensor_3D[i][j][k][l];

    return symTensor_dim;
}


/*!
 * Expand a dim component tensor to a full 3D tensor
 */
template<int dim>
SymmetricTensor<2,3> expand_3D ( SymmetricTensor<2,dim> &symTensor_dim )
{
    SymmetricTensor<2,3> symTensor_3D;
    for ( unsigned int i=0; i<dim; i++)
        for ( unsigned int j=i; j<dim; j++)
        	symTensor_3D[i][j] = symTensor_dim[i][j];

    return symTensor_3D;
}
template<int dim>
Tensor<2,3> expand_3D ( Tensor<2,dim> &tensor_dim )
{
	Tensor<2,3> tensor_3D;
    for ( unsigned int i=0; i<dim; i++)
        for ( unsigned int j=0; j<dim; j++)
        	tensor_3D[i][j] = tensor_dim[i][j];

    return tensor_3D;
}
template<int dim>
SymmetricTensor<4,3> expand_3D ( SymmetricTensor<4,dim> &symTensor_dim )
{
	SymmetricTensor<4,dim> symTensor_3D;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			for ( unsigned int k=0; k<dim; ++k )
				for ( unsigned int l=k; l<dim; ++l )
					symTensor_3D[i][j][k][l] = symTensor_dim[i][j][k][l];

    return symTensor_3D;
}


/*!
 * Extract the tangent contribution that belongs to the theta-strain
 */
template<int dim>
SymmetricTensor<2,dim> extract_theta ( SymmetricTensor<4,3> &symTensor_3D )
{
	SymmetricTensor<2,dim> symTensor_theta;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			symTensor_theta[i][j] = symTensor_3D[i][j][enums::theta][enums::theta];

    return symTensor_theta;
}


/*
 * ################################## Small strain ###################################
 */


/*!
 * Prepare the 2D strain tensor for plane strain, where plane strain is defined
 * by the out-of plane strains being zero
 */
template <int dim>
SymmetricTensor<2,3> prepare_planeStrain( SymmetricTensor<2,dim> &strain )
{
	return expand_3D<dim> (strain);
}


template <int dim>
const double get_radial_x( const FEValues<dim> &fe_values_ref_u, const unsigned int &current_QP )
{
	return (fe_values_ref_u.quadrature_point(current_QP)[enums::r]);
}


template <int dim>
const double get_radial_u( const Vector<double> &current_solution, const FEValues<dim> &fe_values_ref_u,
						   const unsigned int &current_QP )
{
	// ToDo-optimize: Here we evaluate all QPs of the cell but just require the data at the current QP
	 std::vector< Tensor<1,dim> > cell_solution(fe_values_ref_u.get_quadrature().size());
	 // ToDo: use instead of extractor 0 the same value as for \a u_fe in MA-Code.cc
	 fe_values_ref_u.get_function_values(current_solution, cell_solution);
	 return cell_solution[current_QP][enums::u];
}

/*!
 * Prepare the 2D strain tensor for an axial symmetric computation by incorporating
 * the dependency of the ouf-of plane normal strain on the radial displacement
 */
template <int dim>
SymmetricTensor<2,3> prepare_axiSym( SymmetricTensor<2,dim> &strain, const FEValues<dim> &fe_values_ref_u,
									 const Vector<double> &current_solution, const unsigned int &current_QP )
{
	// Get the radial coordinate of the current QP and extract its radial displacement
	 const double radial_x = get_radial_x<dim>(fe_values_ref_u,current_QP);
	 const double radial_u = get_radial_u<dim>(current_solution, fe_values_ref_u, current_QP);
	 
	// Enter the out-of plane normal strain into the 3D strain tensor
	 SymmetricTensor<2,3> strain_3D = expand_3D<dim>(strain);
	 strain_3D[enums::theta][enums::theta] = radial_u/radial_x;

	return strain_3D;
}


// Handling 2D plane strain and axial symmetry
template <int dim>
SymmetricTensor<2,3> prepare_strain ( SymmetricTensor<2,dim> &strain, const unsigned int &type_2D,
									  const FEValues<dim> &fe_values_ref_u, const Vector<double> &current_solution,
									  const unsigned int &current_QP )
{
	if ( dim==3 )
		return expand_3D<dim> (strain); // does nothing for 3D
	else
		switch ( enums::enum_type_2D(type_2D) )
		{
			case enums::planeStrain:
				return prepare_planeStrain<dim> (strain);
			break;
			case enums::axiSym:
				return prepare_axiSym<dim> (strain, fe_values_ref_u, current_solution, current_QP);
			break;
		}
}


// ToDo: check whether we can declare the enumerator for the 2D type (plane strain, axisym) in here in the same namespace enums::
template<int dim>
SymmetricTensor<2,dim> get_Tangent_axisym_addOn( const SymmetricTensor<2,dim> &Tangent_theta,
		  	  	  	  	  	  	  	  	  	  	 const FEValues<dim> &fe_values_ref_u,
												 const unsigned int &current_QP, const unsigned int &j )
{
	// \f$ \frac{\partial \sigma}{\partial \varepsilon_\theta} / r \cdot N_j^u \f$
	return Tangent_theta/get_radial_x<dim>(fe_values_ref_u,k) * fe_values_ref_u.value(j,k)[enums::u];
}

/*
 * ################################## Finite strain ###################################
 */


/*!
 * Prepare the 2D strain tensor for plane strain, where plane strain is defined
 * by the out-of plane strains being zero
 */
template <int dim>
Tensor<2,3> prepare_planeStrain( Tensor<2,dim> &F )
{
	Tensor<2,3> F_3D = expand_3D<dim>(F);
	F_3D[enums::z][enums::z] = 1.;
	return F_3D;
}


/*!
 * Prepare the 2D strain tensor for an axial symmetric computation by incorporating
 * the dependency of the ouf-of plane normal strain on the radial displacement
 */
template <int dim>
Tensor<2,3> prepare_axiSym( Tensor<2,dim> &F, const FEValues<dim> &fe_values_ref_u, const Vector<double> &current_solution,
							const unsigned int &current_QP )
{
	// Get the radial coordinate of the current QP and extract its radial displacement
	 const double radial_x = get_radial_x<dim>(fe_values_ref_u,current_QP);
	 const double radial_u = get_radial_u<dim>(current_solution, fe_values_ref_u, current_QP);

	// Enter the out-of plane normal strain into the 3D strain tensor
	 Tensor<2,3> F_3D = expand_3D<dim>(F);
	 // small strain?
	 F_3D[enums::theta][enums::theta] = 1. + radial_u/radial_x;
	 // deformed configuration?
	 //F_3D[enums::theta][enums::theta] = 1. + radial_u/(radial_x+radial_u);

	return F_3D;
}


// Handling 2D plane strain and axial symmetry
template <int dim>
Tensor<2,3> prepare_DefoGrad( Tensor<2,dim> &F, const unsigned int &type_2D, const FEValues<dim> &fe_values_ref_u,
							  const Vector<double> &current_solution, const unsigned int &current_QP )
{
	if ( dim==3 )
		return expand_3D<dim>(F);  // does nothing for 3D
	else
		switch ( type_2D )
		{
			case enums::planeStrain:
				return prepare_planeStrain<dim> (F);
			break;
			case enums::axiSym:
				return prepare_axiSym<dim> (F, fe_values_ref_u, current_solution, current_QP);
			break;
		}
}


template<int dim>
const double get_JxW ( const unsigned int &type_2D, const FEValues<dim> &fe_values_ref_u, const unsigned int &current_QP )
{
	if ( dim==2 && enums::enum_type_2D(type_2D)==enums::axiSym )
	{
		const double radial_x = get_radial_x<dim>(fe_values_ref_u,current_QP);
		// 2 * pi * r: pi = 4 * arctan(1 rad)
		return (fe_values_ref_u.JxW(current_QP) * (2. * (4.*std::atan(1.)) * radial_x));
	}
	else
		return fe_values_ref_u.JxW(current_QP);
}


#endif // handling_2D_h
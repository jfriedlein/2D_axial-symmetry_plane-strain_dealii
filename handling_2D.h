#ifndef handling_2D_h
#define handling_2D_h

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>

// @todo Check whether the following three headers are needed at all
#include <iostream>
#include <fstream>
#include <cmath>

namespace enums
{
     enum enum_coord_ax
	 {
		r = 0,
		theta = 2,
		u = 0,
		w = 1
	 };

     /**
      * 2D model types:
      * - plane strain: uses the 3D material models but sets all out-of-plane strains to zero for the input
      * - axisymmetric: for axial symmetry, which requires a special out-of-plane strain or deformation gradient component
      * (see <a href="https://github.com/jfriedlein/2D_axial-symmetry_plane-strain_dealii">2D axial symmetry and plane strain</a>)
      */
      enum enum_type_2D
 	 {
 		 planeStrain = 0,//!< planeStrain
 		 axiSym = 1      //!< axiSym
 	 };
}

/*!
 * Extract the dim components from a full 3D tensor
 */
/**
 * @todo-optimize Find a  more efficient way to avoid the transformation for 3D
 * @param symTensor_3D
 * @return
 */
template<int dim>
SymmetricTensor<2,dim> extract_dim ( const SymmetricTensor<2,3> &symTensor_3D )
{
	SymmetricTensor<2,dim> symTensor_dim;
	for ( unsigned int i=0; i<dim; i++)
		for ( unsigned int j=i; j<dim; j++)
			symTensor_dim[i][j] = symTensor_3D[i][j];

	return symTensor_dim;
}
template<int dim>
Tensor<2,dim> extract_dim ( const Tensor<2,3> &Tensor_3D )
{
	Tensor<2,dim> Tensor_dim;
	for ( unsigned int i=0; i<dim; i++)
		for ( unsigned int j=0; j<dim; j++)
			Tensor_dim[i][j] = Tensor_3D[i][j];

	return Tensor_dim;
}
template<int dim>
SymmetricTensor<4,dim> extract_dim ( const SymmetricTensor<4,3> &symTensor_3D )
{
	SymmetricTensor<4,dim> symTensor_dim;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			for ( unsigned int k=0; k<dim; ++k )
				for ( unsigned int l=k; l<dim; ++l )
				  symTensor_dim[i][j][k][l] = symTensor_3D[i][j][k][l];

	return symTensor_dim;
}
template<int dim>
Tensor<4,dim> extract_dim ( const Tensor<4,3> &tensor_3D )
{
	Tensor<4,dim> tensor_dim;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=0; j<dim; ++j )
			for ( unsigned int k=0; k<dim; ++k )
				for ( unsigned int l=0; l<dim; ++l )
					tensor_dim[i][j][k][l] = tensor_3D[i][j][k][l];

	return tensor_dim;
}

/*!
 * Expand a dim component tensor to a full 3D tensor
 */
template<int dim>
SymmetricTensor<2,3> expand_3D ( const SymmetricTensor<2,dim> &symTensor_dim )
{
    SymmetricTensor<2,3> symTensor_3D;
    for ( unsigned int i=0; i<dim; i++)
        for ( unsigned int j=i; j<dim; j++)
        	symTensor_3D[i][j] = symTensor_dim[i][j];

    return symTensor_3D;
}
template<int dim>
Tensor<2,3> expand_3D ( const Tensor<2,dim> &tensor_dim )
{
	Tensor<2,3> tensor_3D;
    for ( unsigned int i=0; i<dim; i++)
        for ( unsigned int j=0; j<dim; j++)
        	tensor_3D[i][j] = tensor_dim[i][j];

    return tensor_3D;
}
template<int dim>
SymmetricTensor<4,3> expand_3D ( const SymmetricTensor<4,dim> &symTensor_dim )
{
	SymmetricTensor<4,dim> symTensor_3D;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			for ( unsigned int k=0; k<dim; ++k )
				for ( unsigned int l=k; l<dim; ++l )
					symTensor_3D[i][j][k][l] = symTensor_dim[i][j][k][l];

    return symTensor_3D;
}
template<int dim>
Tensor<1,3> expand_3D ( const Tensor<1,dim> &tensor_dim )
{
	Tensor<1,3> tensor_3D;
    for ( unsigned int i=0; i<dim; i++)
        	tensor_3D[i] = tensor_dim[i];

    return tensor_3D;
}


/*!
 * Extract the tangent contribution that belongs to the theta-strain
 */
template<int dim>
SymmetricTensor<2,dim> extract_theta ( const SymmetricTensor<4,3> &symTensor_3D )
{
	SymmetricTensor<2,dim> symTensor_theta;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			symTensor_theta[i][j] = symTensor_3D[i][j][enums::theta][enums::theta];

    return symTensor_theta;
}
template<int dim>
Tensor<2,dim> extract_theta ( const Tensor<4,3> &Tensor_3D )
{
	Tensor<2,dim> Tensor_theta;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=0; j<dim; ++j )
			Tensor_theta[i][j] = Tensor_3D[i][j][enums::theta][enums::theta];

    return Tensor_theta;
}

template<int dim>
SymmetricTensor<2,dim> extract_theta_secondPair ( const SymmetricTensor<4,3> &symTensor_3D )
{
	SymmetricTensor<2,dim> symTensor_theta;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			symTensor_theta[i][j] = symTensor_3D[enums::theta][enums::theta][i][j];

    return symTensor_theta;
}
template<int dim>
Tensor<2,dim> extract_theta_secondPair ( const Tensor<4,3> &Tensor_3D )
{
	Tensor<2,dim> Tensor_theta;
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=0; j<dim; ++j )
			Tensor_theta[i][j] = Tensor_3D[enums::theta][enums::theta][i][j];

    return Tensor_theta;
}

template<int dim>
SymmetricTensor<2,dim> extract_comps_secondPair ( const SymmetricTensor<4,dim> &symTensor_3D, const unsigned int i, const unsigned int j )
{
	SymmetricTensor<2,dim> symTensor_theta;
	for ( unsigned int k=0; k<dim; ++k )
		for ( unsigned int l=k; l<dim; ++l )
			symTensor_theta[k][l] = symTensor_3D[i][j][k][l];

    return symTensor_theta;
}
template<int dim>
Tensor<2,dim> extract_comps_secondPair ( const Tensor<4,dim> &Tensor_3D, const unsigned int i, const unsigned int j )
{
	Tensor<2,dim> Tensor_theta;
	for ( unsigned int k=0; k<dim; ++k )
		for ( unsigned int l=0; l<dim; ++l )
			Tensor_theta[k][l] = Tensor_3D[i][j][k][l];

    return Tensor_theta;
}


template <int dim>
double get_radial_x( const FEValues<dim> &fe_values_ref, const unsigned int &current_QP )
{
	return (fe_values_ref.quadrature_point(current_QP)[enums::r]);
}
/**
 * @todo try to merge FEValues and FEFaceValues cases
 * @param fe_values_ref
 * @param current_QP
 * @return
 */
template <int dim>
double get_radial_x( const FEFaceValues<dim> &fe_face_values_ref, const unsigned int &current_QP )
{
	return (fe_face_values_ref.quadrature_point(current_QP)[enums::r]);
}


// @ToDo: avoid using the FEextractor maybe use fe_values[u_fe] as input
template <int dim>
double get_radial_u( const Vector<double> &current_solution, const FEValues<dim> &fe_values_ref,
					 const unsigned int &current_QP )
{
	// ToDo-optimize: Here we evaluate all QPs of the cell but just require the data at the current QP
	 std::vector< Tensor<1,dim> > cell_solution(fe_values_ref.get_quadrature().size());
	 // ToDo: use instead of extractor 0 the same value as for \a u_fe in MA-Code.cc
	 fe_values_ref[(FEValuesExtractors::Vector) 0].get_function_values(current_solution, cell_solution);
	 return cell_solution[current_QP][enums::u];
}


/*!
 * Return the JxW value, possibly scaled by factor 2*pi*r for axisymmetry
 */
template<int dim>
double get_JxW ( const unsigned int &type_2D, const FEValues<dim> &fe_values_ref, const unsigned int &current_QP )
{
	if ( dim==2 && enums::enum_type_2D(type_2D)==enums::axiSym )
	{
		const double radial_x = get_radial_x<dim>(fe_values_ref,current_QP);
		// 2 * pi * r: pi = 4 * arctan(1 rad)
		return (fe_values_ref.JxW(current_QP) * (2. * (4.*std::atan(1.)) * radial_x));
	}
	else
		return fe_values_ref.JxW(current_QP);
}
/**
 * @todo Try to merge the FEValues and FEFaceValue
 * @param type_2D
 * @param fe_face_values_ref
 * @param current_QP
 * @return
 */
template<int dim>
double get_JxW ( const unsigned int &type_2D, const FEFaceValues<dim> &fe_face_values_ref, const unsigned int &current_QP )
{
	if ( dim==2 && enums::enum_type_2D(type_2D)==enums::axiSym )
	{
		const double radial_x = get_radial_x<dim>(fe_face_values_ref,current_QP);
		// 2 * pi * r: pi = 4 * arctan(1 rad)
		return (fe_face_values_ref.JxW(current_QP) * (2. * (4.*std::atan(1.)) * radial_x));
	}
	else
		return fe_face_values_ref.JxW(current_QP);
}


/*
 * ################################## Small strain ###################################
 */


/*!
 * Prepare the 2D strain tensor for plane strain, where plane strain is defined
 * by the out-of plane strains being zero
 */
template <int dim>
SymmetricTensor<2,3> prepare_planeStrain( const SymmetricTensor<2,dim> &strain )
{
	return expand_3D<dim> (strain);
}


/*!
 * Prepare the 2D strain tensor for an axial symmetric computation by incorporating
 * the dependency of the ouf-of plane normal strain on the radial displacement
 */
template <int dim>
SymmetricTensor<2,3> prepare_axiSym( const SymmetricTensor<2,dim> &strain, const FEValues<dim> &fe_values_ref,
									 const Vector<double> &current_solution, const unsigned int &current_QP )
{
	// Get the radial coordinate of the current QP and extract its radial displacement
	 const double radial_x = get_radial_x<dim>(fe_values_ref,current_QP);
	 const double radial_u = get_radial_u<dim>(current_solution, fe_values_ref, current_QP);

	// Enter the out-of plane normal strain into the 3D strain tensor
	 SymmetricTensor<2,3> strain_3D = expand_3D<dim>(strain);
	 strain_3D[enums::theta][enums::theta] = radial_u/radial_x;

	return strain_3D;
}


// Handling 2D plane strain and axial symmetry
template <int dim>
SymmetricTensor<2,3> prepare_strain ( const SymmetricTensor<2,dim> &strain, const unsigned int &type_2D,
									  const FEValues<dim> &fe_values_ref, const Vector<double> &current_solution,
									  const unsigned int &current_QP )
{
	if ( dim==3 )
		return expand_3D<dim> (strain); // does nothing for 3D @todo try SymTensor<2,3>(strain)
	else
		switch ( enums::enum_type_2D(type_2D) )
		{
			case enums::planeStrain:
				return prepare_planeStrain<dim> (strain);
			break;
			case enums::axiSym:
				return prepare_axiSym<dim> (strain, fe_values_ref, current_solution, current_QP);
			break;
			default:
				AssertThrow(false,ExcMessage("sd"));
				break;
		}
}


template<int dim>
SymmetricTensor<2,dim> get_dS_theta_axisym_sstrain( const SymmetricTensor<2,dim> &Tangent_theta,
													 const FEValues<dim> &fe_values_ref,
													 const unsigned int current_QP, const unsigned int j )
{
	// \f$ \frac{\partial \sigma}{\partial \varepsilon_\theta} / r \cdot N_j^u \f$
	return Tangent_theta / get_radial_x<dim>(fe_values_ref,current_QP)
			* fe_values_ref[(FEValuesExtractors::Vector) 0].value(j,current_QP)[enums::u];
}


/*
 * ################################## Finite strain ###################################
 */


/*!
 * Prepare the 2D strain tensor for plane strain, where plane strain is defined
 * by the out-of plane strains being zero
 */
template <int dim>
Tensor<2,3> prepare_planeStrain( const Tensor<2,dim> &F )
{
	Tensor<2,3> F_3D = expand_3D<dim>(F);
	F_3D[enums::theta][enums::theta] = 1.;
	return F_3D;
}


/*!
 * Prepare the 2D strain tensor for an axial symmetric computation by incorporating
 * the dependency of the ouf-of plane normal strain on the radial displacement
 */
template <int dim>
Tensor<2,3> prepare_axiSym( const Tensor<2,dim> &F, const FEValues<dim> &fe_values_ref, const Vector<double> &current_solution,
							const unsigned int &current_QP )
{
	// Get the radial coordinate of the current QP and extract its radial displacement
	 const double radial_x = get_radial_x<dim>(fe_values_ref,current_QP); // this deliveres undeformed values
	 const double radial_u = get_radial_u<dim>(current_solution, fe_values_ref, current_QP);

	// Enter the out-of plane normal strain into the 3D strain tensor
	 Tensor<2,3> F_3D = expand_3D<dim>(F);
	 F_3D[enums::theta][enums::theta] = 1. + radial_u/radial_x;
//	 const double s=5;
//	 const double w=1;
//	 F_3D[enums::theta][enums::theta] = 1. + radial_u/radial_x * (std::tanh(s*radial_x-w)+std::tanh(w)) / (1.+std::tanh(w));

	return F_3D;
}


// Handling 2D plane strain and axial symmetry
template <int dim>
Tensor<2,3> prepare_DefoGrad( const Tensor<2,dim> &F, const unsigned int &type_2D, const FEValues<dim> &fe_values_ref,
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
				return prepare_axiSym<dim> (F, fe_values_ref, current_solution, current_QP);
			break;
			default:
				AssertThrow(false,ExcMessage("sd"));
				break;
		}
}


/**
 * @note Here we desire the pure dS_dC derivative part not "2*" in Tangent_theta.
 * @param Tangent_theta
 * @param fe_values_ref
 * @param current_solution
 * @param current_QP
 * @param j
 * @return
 */
template<int dim>
SymmetricTensor<2,dim> get_dS_theta_axisym_fstrain( const SymmetricTensor<2,dim> &Tangent_theta,
		  	  	  	  	  	  	  	  	  	  	 	const FEValues<dim> &fe_values_ref, const Vector<double> &current_solution,
												 	const unsigned int current_QP, const unsigned int j )
{
	double shape_fnc_j_u = fe_values_ref[(FEValuesExtractors::Vector) 0].value(j,current_QP)[enums::u];
	const double radial_u = get_radial_u<dim>(current_solution, fe_values_ref, current_QP);
	const double radial_x = get_radial_x<dim>(fe_values_ref,current_QP);

	double dCtheta_dur = 2. * ( 1. + radial_u / radial_x ) / radial_x;

	return Tangent_theta * dCtheta_dur * shape_fnc_j_u;
}

template<int dim>
SymmetricTensor<2,dim> get_dS_theta_dF_axisym_fstrain( const Tensor<2,dim> &Tangent_theta,
		  	  	  	  	  	  	  	  	  	  	 	const FEValues<dim> &fe_values_ref, const Vector<double> &current_solution,
												 	const unsigned int current_QP, const unsigned int j )
{
	double shape_fnc_j_u = fe_values_ref[(FEValuesExtractors::Vector) 0].value(j,current_QP)[enums::u];
	const double radial_u = get_radial_u<dim>(current_solution, fe_values_ref, current_QP);
	const double radial_x = get_radial_x<dim>(fe_values_ref,current_QP);

	double dFtheta_dur = 1. / radial_x;

	return symmetrize( Tangent_theta * dFtheta_dur * shape_fnc_j_u );
}


// cell_rhs(i) -= ( R + u_r ) / std::pow( R,2 ) * shape_fnc_i_u[0] * stress_S_3D[2][2] * JxW;
// Equivalent spatial description: cell_rhs(i) -= 1. / ( R + u_r ) * shape_fnc_i_u[0] * Cauchy_stress_3D[2][2] * determinant(DefoGradient_3D) * JxW;
double get_axisym_residual_contribution( const double &R, const double &u_r, const SymmetricTensor<2,3> &stress_S_3D,
										 const double &shape_fnc_i_u, const double &JxW )
{
	return ( R + u_r ) / std::pow( R,2 ) * shape_fnc_i_u * stress_S_3D[enums::theta][enums::theta] * JxW;
}

/**
 * with deltaF
 * @param R
 * @param u_r
 * @param stress_S_3D
 * @param deltaF
 * @param d_S_thetatheta_d_F
 * @param Tangent_theta_theta
 * @param shape_fnc_i_u
 * @param shape_fnc_j_u
 * @param JxW
 * @return
 */
template<int dim>
double get_axisym_linearisation_contribution( const double &R, const double &u_r,
											  const SymmetricTensor<2,3> &stress_S_3D, const Tensor<2,dim> &deltaF,
											  const Tensor<2,dim> &d_S_thetatheta_d_F, const double &Tangent_theta_theta,
											  const double &shape_fnc_i_u, const double &shape_fnc_j_u, const double &JxW )
{
	return shape_fnc_i_u / std::pow( R,2 )
		   * (
				 ( stress_S_3D[enums::theta][enums::theta]  ) * shape_fnc_j_u
				 + (R + u_r)
				   * (
						   double_contract<0,0,1,1>( d_S_thetatheta_d_F, deltaF )
						 + Tangent_theta_theta * shape_fnc_j_u / R
					 )
		   )   * JxW;
}
/**
 * with deltaC
 * @param R
 * @param u_r
 * @param stress_S_3D
 * @param deltaRCG
 * @param d_S_thetatheta_d_C
 * @param Tangent_theta_theta
 * @param shape_fnc_i_u
 * @param shape_fnc_j_u
 * @param JxW
 * @return
 */
template<int dim>
double get_axisym_linearisation_contribution( const double &R, const double &u_r,
											  const SymmetricTensor<2,3> &stress_S_3D, const SymmetricTensor<2,dim> &deltaRCG,
											  const SymmetricTensor<2,dim> &d_S_thetatheta_d_C, const double &Tangent_theta_theta,
											  const double &shape_fnc_i_u, const double &shape_fnc_j_u, const double &JxW )
{
	return shape_fnc_i_u / std::pow( R,2 )
		   * (
				 ( stress_S_3D[enums::theta][enums::theta]  ) * shape_fnc_j_u
				 + (R + u_r)
				   * (
						 d_S_thetatheta_d_C * deltaRCG
						 + Tangent_theta_theta * 2. * ( 1. + u_r / R ) / R * shape_fnc_j_u
					 )
		     ) * JxW;
}



//template<int dim>
//SymmetricTensor<2,dim> get_Tangent_axisym_addOn_fstrain( const SymmetricTensor<2,dim> &Tangent_theta,
//														 const FEValues<dim> &fe_values_ref, const Vector<double> &current_solution,
//														 const unsigned int &current_QP, const unsigned int &j )
//{
//	const double radial_u = get_radial_u<dim>(current_solution, fe_values_ref, current_QP);
//	const double radial_x = get_radial_x<dim>(fe_values_ref,current_QP);
//	// small strain:?
//	 double d_Ctheta_d_ur = 2. * ( 1. + radial_u / radial_x ) / radial_x;
//	// deformed configuration?
////	 double d_Ctheta_d_ur = 2. * ( radial_x + 2.*radial_u ) * radial_x / std::pow(radial_u+radial_x,3);
//
//	// \f$ \frac{\partial S}{\partial C_\theta} \cdot \frac{\partial C_\theta}{\partial u_r} \cdot N_j^u \f$
//	return 0.5 * Tangent_theta * d_Ctheta_d_ur * fe_values_ref[(FEValuesExtractors::Vector) 0].value(j,current_QP)[enums::u];
//}



#endif // handling_2D_h

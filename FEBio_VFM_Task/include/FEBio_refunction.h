#pragma once
#include "common_FEBio.h"
#include "JacobianPolicy.h"

struct PhysicalDeformationGradient
{
	mat3d gradient;
	double determinant;
	double integration_jacobian;
};

struct VirtualFieldGradient
{
	mat3d gradient;
	double determinant;
};

double domain_invjac0(FESolidDomain& domain, const FESolidElement& el, double Ji[3][3], int n);

PhysicalDeformationGradient domain_physical_deformation_gradient(
	FESolidDomain& domain, FESolidElement& el, int n);

VirtualFieldGradient domain_virtual_field_gradient(
	FESolidDomain& domain, FESolidElement& el, int n);

void domain_init(FEMesh& mesh, FESolidElement& el, FESolidDomain& domain);

mat3ds FENeoHookean_Stress(FENeoHookean& fe, FEMaterialPoint& mp);

mat3ds FEIsotropicElastic_Stress(FEIsotropicElastic& fe, FEMaterialPoint& mp);

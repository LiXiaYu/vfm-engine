#include "FEBio_refunction.h"

double domain_invjac0(FESolidDomain& domain, const FESolidElement& el, double Ji[3][3], int n)
{
	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	domain.GetReferenceNodalCoordinates(el, r0);

	// calculate Jacobian
	double J[3][3] = { 0 };
	int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		const double& Gri = el.Gr(n)[i];
		const double& Gsi = el.Gs(n)[i];
		const double& Gti = el.Gt(n)[i];

		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;

		J[0][0] += Gri * x; J[0][1] += Gsi * x; J[0][2] += Gti * x;
		J[1][0] += Gri * y; J[1][1] += Gsi * y; J[1][2] += Gti * y;
		J[2][0] += Gri * z; J[2][1] += Gsi * z; J[2][2] += Gti * z;
	}

	// calculate the determinant
	double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
		+ J[0][1] * (J[1][2] * J[2][0] - J[2][2] * J[1][0])
		+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

	// make sure the determinant is positive
	if (!vfm::is_admissible_jacobian(det, vfm::JacobianRole::ReferenceMapping)) throw NegativeJacobian(el.GetID(), n + 1, det);

	// calculate the inverse jacobian
	double deti = 1.0 / det;

	Ji[0][0] = deti * (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
	Ji[1][0] = deti * (J[1][2] * J[2][0] - J[1][0] * J[2][2]);
	Ji[2][0] = deti * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

	Ji[0][1] = deti * (J[0][2] * J[2][1] - J[0][1] * J[2][2]);
	Ji[1][1] = deti * (J[0][0] * J[2][2] - J[0][2] * J[2][0]);
	Ji[2][1] = deti * (J[0][1] * J[2][0] - J[0][0] * J[2][1]);

	Ji[0][2] = deti * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
	Ji[1][2] = deti * (J[0][2] * J[1][0] - J[0][0] * J[1][2]);
	Ji[2][2] = deti * (J[0][0] * J[1][1] - J[0][1] * J[1][0]);

	return det;
}

namespace
{
struct EvaluatedDeformationGradient
{
	mat3d gradient;
	double determinant;
	double target_mapping_jacobian;
};

EvaluatedDeformationGradient evaluate_domain_deformation_gradient(FESolidDomain& domain, FESolidElement& el, int n)
{
	mat3d F;
	// nodal points
	vec3d r[FEElement::MAX_NODES];
	domain.GetCurrentNodalCoordinates(el, r);

	// calculate inverse jacobian
//    double Ji[3][3];
//    invjac0(el, Ji, n);
	mat3d& Ji = el.m_J0i[n];

	// shape function derivatives
	double* Grn = el.Gr(n);
	double* Gsn = el.Gs(n);
	double* Gtn = el.Gt(n);

	mat3d F1(0, 0, 0, 0, 0, 0, 0, 0, 0);

	// calculate deformation gradient
	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		double Gri = Grn[i];
		double Gsi = Gsn[i];
		double Gti = Gtn[i];

		double x = r[i].x;
		double y = r[i].y;
		double z = r[i].z;

		F1[0][0] += x * Gri; F1[0][1] += x * Gsi; F1[0][2] += x * Gti;
		F1[1][0] += y * Gri; F1[1][1] += y * Gsi; F1[1][2] += y * Gti;
		F1[2][0] += z * Gri; F1[2][1] += z * Gsi; F1[2][2] += z * Gti;

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		double GX = Ji[0][0] * Gri + Ji[1][0] * Gsi + Ji[2][0] * Gti;
		double GY = Ji[0][1] * Gri + Ji[1][1] * Gsi + Ji[2][1] * Gti;
		double GZ = Ji[0][2] * Gri + Ji[1][2] * Gsi + Ji[2][2] * Gti;

		// calculate deformation gradient F
		F[0][0] += GX * x; F[0][1] += GY * x; F[0][2] += GZ * x;
		F[1][0] += GX * y; F[1][1] += GY * y; F[1][2] += GZ * y;
		F[2][0] += GX * z; F[2][1] += GY * z; F[2][2] += GZ * z;
	}

	return { F, F.det(), F1.det() };
}
}

PhysicalDeformationGradient domain_physical_deformation_gradient(FESolidDomain& domain, FESolidElement& el, int n)
{
	auto evaluated = evaluate_domain_deformation_gradient(domain, el, n);
	if (!vfm::is_admissible_jacobian(evaluated.determinant, vfm::JacobianRole::PhysicalDeformation))
	{
		throw NegativeJacobian(el.GetID(), n, evaluated.determinant, &el);
	}
	return { evaluated.gradient, evaluated.determinant, evaluated.target_mapping_jacobian };
}

VirtualFieldGradient domain_virtual_field_gradient(FESolidDomain& domain, FESolidElement& el, int n)
{
	auto evaluated = evaluate_domain_deformation_gradient(domain, el, n);
	if (!vfm::is_admissible_jacobian(evaluated.determinant, vfm::JacobianRole::VirtualField))
	{
		throw NegativeJacobian(el.GetID(), n, evaluated.determinant, &el);
	}
	return { evaluated.gradient, evaluated.determinant };
}
void domain_init(FEMesh& mesh, FESolidElement& el, FESolidDomain& domain) {

	// evaluate nodal coordinates
	const int NELN = FEElement::MAX_NODES;
	vec3d r0[NELN], r[NELN], v[NELN], a[NELN];
	int neln = el.Nodes();
	for (int j = 0; j < neln; ++j)
	{
		FENode& node = mesh.Node(el.m_node[j]);
		r0[j] = node.m_r0;
	}
	//mesh.Domain(0).
	// initialize reference Jacobians
	double Ji[3][3];



	// loop over the integration points
	int nint = el.GaussPoints();
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		// initiali Jacobian

		try
		{


			//mp.m_J0 = static_cast<FESolidDomain&>(mesh.Domain(0)).invjac0(el, Ji, n);
			mp.m_J0 = domain_invjac0(domain, el, Ji, n);
		}
		catch (NegativeJacobian& e)
		{
			throw e;
		}

		el.m_J0i[n] = mat3d(Ji);

		// material point coordinates
		mp.m_r0 = el.Evaluate(r0, n);
	}
}

mat3ds FENeoHookean_Stress(FENeoHookean& fe, FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double detF = pt.m_J;
	double detFi = 1.0 / detF;
	double lndetF = log(detF);

	// get the material parameters
	double E = fe.m_E(mp);
	double v = fe.m_v(mp);

	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// lame parameters
	double lam = v * E / ((1 + v) * (1 - 2 * v));
	double mu = 0.5 * E / (1 + v);

	// Identity
	mat3dd I(1);

	// calculate stress
	mat3ds s = (b - I) * (mu * detFi) + I * (lam * lndetF * detFi);

	return s;
}

mat3ds FEIsotropicElastic_Stress(FEIsotropicElastic& fe, FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d& F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

	double E = fe.m_E(mp);
	double v = fe.m_v(mp);

	// lame parameters
	double lam = Ji * (v * E / ((1 + v) * (1 - 2 * v)));
	double mu = Ji * (0.5 * E / (1 + v));

	// calculate left Cauchy-Green tensor (ie. b-matrix)
	mat3ds b = pt.LeftCauchyGreen();

	// calculate trace of Green-Lagrance strain tensor
	double trE = 0.5 * (b.tr() - 3);

	// calculate square of b-matrix
	// (we commented out the matrix components we do not need)
	mat3ds b2 = b.sqr();

	// calculate stress
	mat3ds s = b * (lam * trE - mu) + b2 * mu;

	return s;
}

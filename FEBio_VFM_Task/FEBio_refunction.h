#pragma once
#include "common_FEBio.h"

double domain_invjac0(FESolidDomain& domain, const FESolidElement& el, double Ji[3][3], int n);

double domain_defgrad_GJ(FESolidDomain& domain, FESolidElement& el, mat3d& F, int n);

void domain_init(FEMesh& mesh, FESolidElement& el, FESolidDomain& domain);

mat3ds FENeoHookean_Stress(FENeoHookean& fe, FEMaterialPoint& mp);

mat3ds FEIsotropicElastic_Stress(FEIsotropicElastic& fe, FEMaterialPoint& mp);
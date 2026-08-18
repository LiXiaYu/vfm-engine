#include "JacobianPolicy.h"

#include <iostream>
#include <limits>

namespace
{
int failures = 0;

void expect(bool condition, const char* description)
{
	if (!condition)
	{
		std::cerr << "FAILED: " << description << '\n';
		++failures;
	}
}
}

int main()
{
	using vfm::JacobianRole;
	using vfm::is_admissible_jacobian;

	expect(is_admissible_jacobian(1.0, JacobianRole::ReferenceMapping),
		"positive reference mapping is admissible");
	expect(!is_admissible_jacobian(0.0, JacobianRole::ReferenceMapping),
		"singular reference mapping is rejected");
	expect(!is_admissible_jacobian(-1.0, JacobianRole::ReferenceMapping),
		"reversed reference mapping is rejected");

	expect(is_admissible_jacobian(1.0, JacobianRole::PhysicalDeformation),
		"positive physical deformation is admissible");
	expect(!is_admissible_jacobian(0.0, JacobianRole::PhysicalDeformation),
		"singular physical deformation is rejected");
	expect(!is_admissible_jacobian(-1.0, JacobianRole::PhysicalDeformation),
		"negative physical deformation is rejected");

	expect(is_admissible_jacobian(1.0, JacobianRole::VirtualField),
		"positive virtual-field determinant is admissible");
	expect(is_admissible_jacobian(0.0, JacobianRole::VirtualField),
		"singular virtual-field determinant is admissible when no inverse is used");
	expect(is_admissible_jacobian(-1.0, JacobianRole::VirtualField),
		"negative virtual-field determinant is admissible");

	expect(!is_admissible_jacobian(std::numeric_limits<double>::infinity(),
		JacobianRole::VirtualField), "infinite determinant is rejected");
	expect(!is_admissible_jacobian(std::numeric_limits<double>::quiet_NaN(),
		JacobianRole::VirtualField), "NaN determinant is rejected");

	if (failures == 0)
	{
		std::cout << "Jacobian policy tests passed\n";
	}
	return failures == 0 ? 0 : 1;
}

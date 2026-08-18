#pragma once

#include <cmath>

namespace vfm
{
enum class JacobianRole
{
	ReferenceMapping,
	PhysicalDeformation,
	VirtualField
};

inline bool is_admissible_jacobian(double determinant, JacobianRole role) noexcept
{
	if (!std::isfinite(determinant)) return false;

	switch (role)
	{
	case JacobianRole::ReferenceMapping:
	case JacobianRole::PhysicalDeformation:
		return determinant > 0.0;
	case JacobianRole::VirtualField:
		// The virtual gradient is used to form a test strain, not an
		// invertible physical deformation. Its orientation is unrestricted.
		return true;
	}

	return false;
}
}

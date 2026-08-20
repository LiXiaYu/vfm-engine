# Gaussian integration across solid element types

This validation checks the integration measure used by VFM-engine independently
of a constitutive optimization. Each model contains one affine solid element
with a known analytical volume and one planar pressure surface.

The VFM run injects two stress states, zero followed by a constant isotropic
stress, and uses the linear virtual displacement `u*_x = c x`. The validator
checks three independent identities:

1. the sum of the per-Gauss-point volume measures equals the analytical volume;
2. the internal virtual work equals `sigma : epsilon*` times that volume;
3. the pressure-surface virtual work equals the analytical planar-surface
   integral.

The matrix covers `hex8`, `tet4`, `penta6`, and `pyra5`, including both `quad4`
and `tri3` pressure surfaces. Generated data is written below
`out/validation/gauss_integration_solid_elements/`.

## Integration rule under test

The volume and surface integrals must use the complete FEBio quadrature
measures

```text
dV_n = det(dx/dxi)_n * solid.GaussWeights()[n]
dA_n = |dx/dr x dx/ds|_n * surface.GaussWeights()[n]
```

FEBio's `H`, `Gr`, `Gs`, and `Gt` values are shape functions and shape-function
derivatives. They do not contain the quadrature weight. This distinction is
invisible for rules whose weights happen to be one, but it is essential for
tetrahedral, wedge, pyramid, and triangular-surface rules.

The reviewed run produces the following values. The virtual field scale is
`c = 1e-4`, the imposed normal stress is `2`, and the expected internal work is
`2 * (c + c^2 / 2) * volume`.

| Element | Integrated volume | Internal virtual work | Surface virtual work |
| --- | ---: | ---: | ---: |
| `hex8` | `1` | `2.0001e-4` | `2.0e-4` |
| `tet4` | `1/6` | `3.3335e-5` | `3.333333333333e-5` |
| `penta6` | `0.5` | `1.00005e-4` | `-1.0e-4` |
| `pyra5` | `4/3` | `2.6668e-4` | `1.333333333333e-4` |

Before the surface-weight correction, the `tet4` and `pyra5` triangular-face
results were respectively `2.0e-4` and `8.0e-4`: exactly six times their
analytical values. Their three surface points each have weight `1/6`. The
validator compares signed work as well as magnitude, so a reversed pressure
normal also fails.

This case establishes the integration identities for the four linear solid
families above. Additional higher-order element formulations should be added
to this matrix before claiming coverage for those formulations. It isolates
quadrature and does not by itself validate a constitutive update, time
discretization, or parameter optimizer; those remain covered by their own
cases.

Run from the repository root:

```powershell
.\validation\cases\gauss_integration_solid_elements\run.ps1
```

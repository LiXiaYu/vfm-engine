# Tests

This directory contains automated tests for VFM-engine itself.

Tests should be small, deterministic, and suitable for local development and continuous integration. As the test suite grows, organize it into focused subdirectories such as `unit/`, `integration/`, `regression/`, and `fixtures/`.

Paper-specific experiments and large research datasets do not belong here.

The current `jacobian_policy` unit test documents and enforces the different
admissibility rules for reference mappings, physical deformations, and virtual
test fields. Configure the project with testing enabled, build it, and run CTest
from the configured build directory.

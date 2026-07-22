# Affine Coxeter Explorer

A desktop application for exploring conjugacy classes and coconjugation sets in affine Weyl groups of dimensions 2 and 3.

## Supported groups
- Ã₁ × Ã₁
- Ã₂, B̃₂, C̃₂, G̃₂
- Ã₃, B̃₃, C̃₃

## Requirements
- Java 11 or later
- SageMath 9.0 or later

## Running the app
Double-click the appropriate launcher for your OS:
- `Launch (macOS).command`
- `Launch (Windows).bat`
- `Launch (Linux).sh`

Or run manually:
```bash
java --source 11 AffineCoxeterExplorer.java
```

## How to read the pictures
Alcove colors encode the spherical direction (finite Weyl part) of each
element. In dimension 2, every element of the finite Weyl group has its own
color, fixed per type, so colors can be compared between figures. In
dimension 3, where the finite Weyl groups are too large for that:

- **Conjugacy-class figures** use distinct colors within the class, so every
  spherical direction appearing in the figure is distinguishable.
- **Coconjugation figures** use one color per conjugacy class of the finite
  Weyl group, so equal color means conjugate finite parts.

In both dimensions, the identity alcove is cyan and labeled `e`, striping
marks the identity as a member of the computed set, the chosen element(s)
are outlined in red (and blue, for coconjugation), and in dimension 3 the
uncolored tessellation is drawn as a gray wireframe.

## Authors
Amy Herron, Anne Thomas

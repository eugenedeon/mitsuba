Mitsuba 0.x with additional BSDFs and phase functions
===================================

This is my fork of the 0.x mitsuba repo:
http://mitsuba-renderer.org/

## BRDFs

A batch script in `scenes/lambert_sphere/Figure8` tests a variety of diffuse BRDFs, which can be compared using the provided Mathematica notebook.

### Lambert-sphere BRDF
<img src="http://www.eugenedeon.com/wp-content/uploads/2021/06/LS_BRDF_teaser-scaled.jpg" alt="Lambert-sphere BRDF teaser">
This models the reflectance from a half space comprised of sparsely distributed spherical Lambertian particles from the EGSR 2021 paper.
Three versions implemented: `lambert_sphere`, `lambert_sphere_fast` and `lambert_sphere_hapke81`.

### Chandrasekhar's BRDF
This BRDF is for a half space with isotropic scattering.  Albedo inversion and importance sampling implemented.
Implemented as: `chandra`.

### Diffusion-transport BDRF
This models the reflectance from a non-exponential half space with isotropic scattering where the distances between collisions are drawn from a Gamma/Erlang-2 distribution.  Implemented as: `diffusion_transport`.

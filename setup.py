from distutils.core import setup
from distutils.extension import Extension

try:
    import Cython
    USE_CYTHON = True
    print "using cython"
except ImportError:
    USE_CYTHON = False
    print "cython not present, use the default C file"

extension = ".pyx" if USE_CYTHON else ".c"

sourcefiles = [
# cython files
"src/lyapunov_exponents" + extension,
# pure C file
"src/random.c",
"src/monodromy.c",
"src/lin_alg_complex.c",
"src/geodesic.c"]

extensions = [Extension("lyapunov_exponents", sourcefiles)]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name="lyapunov_exponents",
    version = "0.alpha1",
    author = "C. Fougeron",
    author_email = "charles.fougeron@ens.fr",
    description="Computing Lyapunov of K problem",
    long_description=
"""LEKZ: Lyapunov exponents of the Kontzevich-Zorich cocycles

Each affine SL(2,R)-ergodic measure on the moduli space of quadratic
differentials on Riemman surfaces give rise to Lyapunov exponents. This package
allows to compute approximations of them for various affine measures:
- the Masur-Veech measure
- coverings construction (which contain the so called pillowcase covers and
  origamis)""",
    ext_modules = extensions)

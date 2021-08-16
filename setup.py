from setuptools import setup, find_packages
from Cython.Build import cythonize


with open("README.md") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

extras = {
    'dev': ['bump2version'],
}

setup(
    name="optbeam",
    version="2.1.1",
    description=("Simulation of reflection and refraction of polarized "
                 "opticial beams at plane and curved dielectric interfaces"),
    long_description=readme,
    long_description_content_type='text/markdown',
    url="https://github.com/DanielKotik/Optical-beams-MEEP",
    author="Daniel Kotik et al.",
    author_email="kotik@physics.org",
    license=license,
    packages=find_packages(exclude=("scripts")),
    ext_modules=cythonize(["optbeam/**/beams.py",
                           "optbeam/**/helpers.py"],
                          compiler_directives={'language_level': 3}),
    zip_safe=False,
    include_package_data=True,
    install_requires=["scipy", "cython"],
    extras_require=extras,
)

from setuptools import setup, find_packages


with open("README.md") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

extras = {
    'dev': ['bump2version'],
}

setup(
    name="optbeam",
    version="2.0.1",
    description=("Simulation of reflection and refraction of polarized "
                 "opticial beams at plane and curved dielectric interfaces"),
    long_description=readme,
    long_description_content_type='text/markdown',
    url="https://github.com/DanielKotik/Optical-beams-MEEP",
    author="Daniel Kotik et al.",
    author_email="kotik@physics.org",
    license=license,
    packages=find_packages(exclude=("scripts")),
    include_package_data=True,
    install_requires=["scipy"],
    extras_require=extras,
)
